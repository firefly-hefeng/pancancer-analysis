# cell_communication_analysis.R
# 基于功能模块和空间共定位的细胞通讯分析 - 修复基因名称匹配问题

library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(networkD3)
library(igraph)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(future)
library(parallel)
library(reshape2)
library(stringr)
library(ggraph)

# 设置并行计算 - 增加全局对象大小限制
options(future.globals.maxSize = 50 * 1024^3)  # 设置为 50 GB
future::plan(multicore, workers = 6)

# 参数设置
output_base_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/fine_annotation"
communication_output_dir <- file.path(output_base_dir, "cell_communication_analysis")
dir.create(communication_output_dir, showWarnings = FALSE, recursive = TRUE)

# 创建子目录
sub_dirs <- c("functional_modules", "spatial_colocation", "pathway_analysis", 
              "network_analysis", "plots", "tables", "interactive_plots", "qc_reports")
for (dir in sub_dirs) {
  dir.create(file.path(communication_output_dir, dir), showWarnings = FALSE, recursive = TRUE)
}

message("=== 细胞通讯分析框架 ===")
message("分析策略:")
message("1. 功能模块分层 - 按生物学功能对细胞类型分组")
message("2. 空间共定位分层 - 基于细胞空间分布概率")
message("3. 信号通路分层 - 按通路类型分类分析")
message("4. 多维度网络构建与可视化")
message(paste("Future全局对象大小限制:", format(options()$future.globals.maxSize / 1024^3, digits = 2), "GB"))

# 定义功能模块
functional_modules <- list(
  "Immune_Regulation" = c("T_cells", "B_cell", "NK_cell", "Myeloid"),
  "Tissue_Structure" = c("Epithelial_cells", "Endothelial_cells", "Fibroblasts"),
  "Microenvironment_Support" = c("CAFs", "Pericytes", "Smooth_muscle"),
  "Antigen_Presentation" = c("Dendritic_cells", "Macrophages", "B_cell"),
  "Cytotoxicity" = c("T_cells", "NK_cell", "CTL"),
  "Angiogenesis" = c("Endothelial_cells", "Pericytes", "Macrophages"),
  "ECM_Remodeling" = c("Fibroblasts", "CAFs", "Macrophages")
)

# 定义空间共定位可能性权重
spatial_colocalization_weights <- list(
  "High_Probability" = list(
    pairs = list(
      c("T_cells", "B_cell"),
      c("T_cells", "Dendritic_cells"),
      c("Endothelial_cells", "Pericytes"),
      c("Fibroblasts", "Epithelial_cells"),
      c("Macrophages", "T_cells"),
      c("CAFs", "Fibroblasts")
    ),
    weight = 1.0
  ),
  "Medium_Probability" = list(
    pairs = list(
      c("NK_cell", "T_cells"),
      c("Macrophages", "Endothelial_cells"),
      c("Fibroblasts", "Endothelial_cells"),
      c("B_cell", "Dendritic_cells")
    ),
    weight = 0.7
  ),
  "Low_Probability" = list(
    pairs = list(
      c("Epithelial_cells", "T_cells"),
      c("Epithelial_cells", "B_cell"),
      c("NK_cell", "Fibroblasts")
    ),
    weight = 0.4
  )
)

# 定义信号通路分类
pathway_categories <- list(
  "Immune_Signaling" = c("CXCL", "CCL", "TNF", "IFNG", "IL", "TGFb", "CD40"),
  "Growth_Factors" = c("VEGF", "EGF", "FGF", "PDGF", "HGF", "IGF"),
  "ECM_Interaction" = c("COLLAGEN", "LAMININ", "FN1", "VTN", "THBS"),
  "Cell_Adhesion" = c("ICAM", "VCAM", "CADM", "CDH", "ITGB"),
  "Metabolic_Signaling" = c("GAS", "ANGPTL", "SEMA", "PLAU"),
  "Developmental" = c("WNT", "NOTCH", "BMP", "SHH", "JAG")
)

# 错误收集器
error_collector <- list()

# 安全执行函数
safe_execute <- function(expr, step, cancer_type = "unknown", default_return = NULL) {
  tryCatch({
    result <- expr
    message(paste("SUCCESS [", step, "]:", cancer_type))
    return(result)
  }, error = function(e) {
    error_info <- list(
      timestamp = Sys.time(),
      step = step,
      cancer_type = cancer_type,
      error = e$message,
      traceback = capture.output(traceback())
    )
    error_collector[[paste(step, cancer_type, sep = "_")]] <<- error_info
    message(paste("ERROR [", step, "]:", cancer_type, "-", e$message))
    return(default_return)
  })
}

# 数据加载函数 - 优化版本
load_cancer_data <- function(cancer_type) {
  message(paste("加载数据:", cancer_type))
  
  # 检查quickfix文件
  quickfix_file <- file.path(output_base_dir, cancer_type, 
                            paste0(cancer_type, "_fine_annotated_quickfix.rds"))
  main_file <- file.path(output_base_dir, cancer_type, 
                        paste0(cancer_type, "_fine_annotated.rds"))
  
  if (file.exists(quickfix_file)) {
    use_file <- quickfix_file
    message("  使用quickfix版本")
  } else if (file.exists(main_file)) {
    use_file <- main_file
    message("  使用原始版本")
  } else {
    message(paste("  数据文件不存在:", cancer_type))
    return(NULL)
  }
  
  return(safe_execute({
    seurat_obj <- readRDS(use_file)
    
    # 确定细胞类型列
    if ("fine_cell_type" %in% colnames(seurat_obj@meta.data)) {
      cell_type_col <- "fine_cell_type"
    } else if ("fine_cell_type_fixed" %in% colnames(seurat_obj@meta.data)) {
      cell_type_col <- "fine_cell_type_fixed"
    } else {
      stop("未找到细胞类型注释列")
    }
    
    # 检查数据标准化状态
    has_normalized_data <- FALSE
    
    if ("data" %in% slotNames(seurat_obj@assays$RNA)) {
      data_matrix <- slot(seurat_obj@assays$RNA, "data")
      if (!is.null(data_matrix) && nrow(data_matrix) > 0 && ncol(data_matrix) > 0) {
        counts_matrix <- slot(seurat_obj@assays$RNA, "counts")
        if (!identical(data_matrix, counts_matrix)) {
          has_normalized_data <- TRUE
          message("  数据已标准化")
        }
      }
    }
    
    if (!has_normalized_data) {
      message("  数据未标准化，执行标准化...")
      seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", 
                                   scale.factor = 10000, verbose = FALSE)
    }
    
    # 简化细胞类型名称
    seurat_obj$main_cell_type <- gsub("_.*", "", seurat_obj@meta.data[[cell_type_col]])
    
    # 过滤低质量细胞和基因
    seurat_obj <- subset(seurat_obj, 
                        subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & 
                                percent.mt < 25)
    
    # 确保有足够的细胞进行分析
    cell_counts <- table(seurat_obj$main_cell_type)
    valid_types <- names(cell_counts)[cell_counts >= 50]
    
    if (length(valid_types) < 3) {
      stop(paste("有效细胞类型过少:", length(valid_types)))
    }
    
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj$main_cell_type %in% valid_types])
    
    # 如果数据集太大，进行下采样
    if (ncol(seurat_obj) > 50000) {
      message("  数据集较大，进行下采样...")
      Idents(seurat_obj) <- "main_cell_type"
      seurat_obj <- subset(seurat_obj, downsample = 5000)
      message(paste("  下采样后细胞数:", ncol(seurat_obj)))
    }
    
    # 清理不需要的数据以减少内存占用
    if ("scale.data" %in% slotNames(seurat_obj@assays$RNA)) {
      slot(seurat_obj@assays$RNA, "scale.data") <- matrix()
    }
    
    # 清理不需要的reductions
    seurat_obj@reductions <- list()
    
    # 清理graphs和neighbors
    seurat_obj@graphs <- list()
    seurat_obj@neighbors <- list()
    
    message(paste("  成功加载:", ncol(seurat_obj), "细胞,", length(valid_types), "种细胞类型"))
    message(paste("  细胞类型:", paste(valid_types, collapse = ", ")))
    
    # 估算对象大小
    obj_size <- object.size(seurat_obj)
    message(paste("  对象大小:", format(obj_size, units = "auto")))
    
    # **检查并修复基因名称格式**
    gene_names <- rownames(seurat_obj)
    message(paste("  基因名称示例:", paste(head(gene_names, 5), collapse = ", ")))
    
    # 检查基因名称格式（是否全大写、全小写或首字母大写）
    all_upper <- all(gene_names == toupper(gene_names))
    all_lower <- all(gene_names == tolower(gene_names))
    
    if (all_upper) {
      message("  检测到基因名称为全大写格式")
    } else if (all_lower) {
      message("  检测到基因名称为全小写格式")
    } else {
      message("  检测到基因名称为混合格式（可能是首字母大写）")
    }
    
    return(seurat_obj)
    
  }, "data_loading", cancer_type, NULL))
}

# 创建CellChat对象 - 修复基因名称匹配问题
create_cellchat_object <- function(seurat_obj, cancer_type) {
  message(paste("创建CellChat对象:", cancer_type))
  
  return(safe_execute({
    # 提取标准化数据
    data.input <- slot(seurat_obj@assays$RNA, "data")
    
    # 转换为普通矩阵以减少future传输大小（如果不是太大）
    if (ncol(data.input) < 20000 && nrow(data.input) < 20000) {
      message("  转换为普通矩阵以优化传输")
      data.input <- as.matrix(data.input)
    }
    
    # 检查数据
    if (is.null(data.input) || nrow(data.input) == 0 || ncol(data.input) == 0) {
      stop("提取的表达数据为空")
    }
    
    message(paste("  表达矩阵维度:", nrow(data.input), "x", ncol(data.input)))
    message(paste("  数据类型:", class(data.input)))
    
    # 检查数据范围
    non_zero_values <- data.input[data.input > 0]
    if (length(non_zero_values) > 0) {
      message(paste("  非零值范围:", round(min(non_zero_values), 3), "到", round(max(non_zero_values), 3)))
    }
    
    # **检查基因名称格式**
    gene_names <- rownames(data.input)
    message(paste("  原始基因名称示例:", paste(head(gene_names, 5), collapse = ", ")))
    
    meta <- seurat_obj@meta.data
    
    # 创建CellChat对象
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "main_cell_type")
    
    # 设置配体-受体数据库
    CellChatDB <- CellChatDB.human
    
    # **关键修复：检查并统一基因名称格式**
    db_genes <- unique(c(CellChatDB$interaction$ligand, CellChatDB$interaction$receptor))
    expr_genes <- rownames(data.input)
    
    message(paste("  数据库基因名称示例:", paste(head(db_genes, 5), collapse = ", ")))
    message(paste("  表达数据基因名称示例:", paste(head(expr_genes, 5), collapse = ", ")))
    
    # 初始检查重叠
    initial_overlap <- intersect(db_genes, expr_genes)
    message(paste("  初始基因重叠数:", length(initial_overlap)))
    
    # 如果重叠很少，尝试转换基因名称格式
    if (length(initial_overlap) < 100) {
      message("  基因重叠过少，尝试转换基因名称格式...")
      
      # 尝试不同的转换策略
      # 策略1：全部转为大写
      expr_genes_upper <- toupper(expr_genes)
      overlap_upper <- intersect(db_genes, expr_genes_upper)
      
      # 策略2：全部转为首字母大写
      expr_genes_title <- paste0(toupper(substr(expr_genes, 1, 1)), 
                                 tolower(substr(expr_genes, 2, nchar(expr_genes))))
      overlap_title <- intersect(db_genes, expr_genes_title)
      
      # 策略3：数据库基因转为小写与表达基因匹配
      db_genes_lower <- tolower(db_genes)
      expr_genes_lower <- tolower(expr_genes)
      
      message(paste("  大写转换后重叠:", length(overlap_upper)))
      message(paste("  首字母大写转换后重叠:", length(overlap_title)))
      
      # 选择最佳转换策略
      if (length(overlap_upper) > length(initial_overlap) && 
          length(overlap_upper) > length(overlap_title)) {
        message("  使用全大写转换策略")
        
        # 转换表达数据的基因名称为大写
        new_gene_names <- toupper(rownames(cellchat@data))
        
        # 检查是否有重复
        if (any(duplicated(new_gene_names))) {
          message("  警告：转换后有重复基因名，需要合并")
          # 对于重复的基因，保留第一个或者取平均值
          # 这里简单保留第一个
          keep_idx <- !duplicated(new_gene_names)
          cellchat@data <- cellchat@data[keep_idx, ]
          rownames(cellchat@data) <- new_gene_names[keep_idx]
        } else {
          rownames(cellchat@data) <- new_gene_names
        }
        
      } else if (length(overlap_title) > length(initial_overlap)) {
        message("  使用首字母大写转换策略")
        
        new_gene_names <- paste0(toupper(substr(rownames(cellchat@data), 1, 1)), 
                                tolower(substr(rownames(cellchat@data), 2, 
                                              nchar(rownames(cellchat@data)))))
        
        if (any(duplicated(new_gene_names))) {
          message("  警告：转换后有重复基因名")
          keep_idx <- !duplicated(new_gene_names)
          cellchat@data <- cellchat@data[keep_idx, ]
          rownames(cellchat@data) <- new_gene_names[keep_idx]
        } else {
          rownames(cellchat@data) <- new_gene_names
        }
        
      } else {
        message("  警告：无法找到有效的基因名称转换策略")
        message("  将继续使用原始基因名称，但可能导致检测到的通讯较少")
      }
    }
    
    # 使用完整数据库
    cellchat@DB <- CellChatDB
    message("  使用完整CellChatDB数据库")
    message(paste("  数据库包含", nrow(CellChatDB$interaction), "个配体-受体对"))
    
    # 最终检查数据库基因与表达数据的重叠
    db_genes <- unique(c(CellChatDB$interaction$ligand, CellChatDB$interaction$receptor))
    expr_genes <- rownames(cellchat@data)
    overlap_genes <- intersect(db_genes, expr_genes)
    message(paste("  最终数据库与表达数据基因重叠:", length(overlap_genes), "/", length(db_genes)))
    
    if (length(overlap_genes) < 100) {
      warning("数据库基因与表达数据重叠较少，可能影响分析结果")
      message(paste("  重叠基因示例:", paste(head(overlap_genes, 10), collapse = ", ")))
    }
    
    return(cellchat)
    
  }, "create_cellchat", cancer_type, NULL))
}

# 预处理和分析CellChat - 优化版本（禁用future避免大对象传输）
analyze_cellchat <- function(cellchat, cancer_type) {
  message(paste("分析细胞通讯:", cancer_type))
  
  return(safe_execute({
    # 临时禁用future以避免大对象传输问题
    old_plan <- future::plan()
    future::plan(sequential)
    
    on.exit({
      future::plan(old_plan)
    })
    
    # 预处理数据
    message("  步骤1: 子集化数据...")
    cellchat <- subsetData(cellchat)
    
    # 检查数据是否正确
    if (is.null(cellchat@data.signaling) || nrow(cellchat@data.signaling) == 0) {
      # 提供更详细的诊断信息
      message("  诊断信息:")
      message(paste("    - 原始数据基因数:", nrow(cellchat@data)))
      message(paste("    - 数据库配体-受体对数:", nrow(cellchat@DB$interaction)))
      
      db_genes <- unique(c(cellchat@DB$interaction$ligand, cellchat@DB$interaction$receptor))
      expr_genes <- rownames(cellchat@data)
      overlap <- intersect(db_genes, expr_genes)
      
      message(paste("    - 数据库基因数:", length(db_genes)))
      message(paste("    - 表达数据基因数:", length(expr_genes)))
      message(paste("    - 重叠基因数:", length(overlap)))
      
      if (length(overlap) > 0) {
        message(paste("    - 重叠基因示例:", paste(head(overlap, 20), collapse = ", ")))
      }
      
      stop("subsetData后数据为空，可能数据库基因与表达数据不匹配")
    }
    
    message(paste("  信号基因数:", nrow(cellchat@data.signaling)))
    message(paste("  细胞数:", ncol(cellchat@data.signaling)))
    
    # 识别过表达基因
    message("  步骤2: 识别过表达基因...")
    cellchat <- identifyOverExpressedGenes(cellchat)
    
    # 识别过表达相互作用
    message("  步骤3: 识别过表达相互作用...")
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # 计算通讯概率 - 使用单线程
    message("  步骤4: 计算通讯概率...")
    cellchat <- computeCommunProb(cellchat, 
                                 type = "triMean", 
                                 trim = 0.1,
                                 population.size = TRUE,
                                 nboot = 100)
    
    # 过滤通讯
    message("  步骤5: 过滤低置信度通讯...")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    # 检查是否有足够的通讯
    if (is.null(cellchat@net$prob) || sum(cellchat@net$prob > 0, na.rm = TRUE) == 0) {
      warning("未检测到显著的细胞通讯")
      return(cellchat)
    }
    
    message(paste("  检测到", sum(cellchat@net$prob > 0, na.rm = TRUE), "个细胞对通讯"))
    
    # 推断细胞通讯网络
    message("  步骤6: 计算通路水平的通讯...")
    cellchat <- computeCommunProbPathway(cellchat)
    
    message("  步骤7: 聚合细胞通讯网络...")
    cellchat <- aggregateNet(cellchat)
    
    # 识别信号通路角色
    message("  步骤8: 分析信号通路角色...")
    cellchat <- netAnalysis_signalingRole_network(cellchat, slot.name = "netP")
    
    message(paste("  最终检测到通讯数:", sum(cellchat@net$count > 0, na.rm = TRUE)))
    message(paste("  涉及的信号通路数:", length(unique(cellchat@netP$pathways))))
    
    # 恢复future设置
    future::plan(old_plan)
    
    return(cellchat)
    
  }, "analyze_cellchat", cancer_type, NULL))
}

# 功能模块分析
analyze_functional_modules <- function(cellchat, cancer_type) {
  message(paste("功能模块分析:", cancer_type))
  
  return(safe_execute({
    if (is.null(cellchat@net$prob) || sum(cellchat@net$prob > 0, na.rm = TRUE) == 0) {
      message("  无有效通讯网络，跳过功能模块分析")
      return(list())
    }
    
    cell_types <- levels(cellchat@idents)
    module_results <- list()
    
    for (module_name in names(functional_modules)) {
      module_cells <- intersect(functional_modules[[module_name]], cell_types)
      
      if (length(module_cells) >= 2) {
        module_net <- subsetCommunication(cellchat, 
                                        sources.use = module_cells, 
                                        targets.use = module_cells)
        
        if (nrow(module_net) > 0) {
          module_stats <- list(
            total_interactions = nrow(module_net),
            avg_prob = mean(module_net$prob, na.rm = TRUE),
            unique_pathways = length(unique(module_net$pathway_name)),
            cell_types = module_cells,
            interactions = module_net
          )
          
          module_results[[module_name]] <- module_stats
          message(paste("  ", module_name, ":", module_stats$total_interactions, "个相互作用"))
        }
      }
    }
    
    if (length(module_results) == 0) {
      message("  未检测到功能模块间的通讯")
    }
    
    return(module_results)
    
  }, "functional_modules", cancer_type, list()))
}

# 空间共定位分析
analyze_spatial_colocation <- function(cellchat, cancer_type) {
  message(paste("空间共定位分析:", cancer_type))
  
  return(safe_execute({
    if (is.null(cellchat@net$prob) || sum(cellchat@net$prob > 0, na.rm = TRUE) == 0) {
      message("  无有效通讯网络，跳过空间共定位分析")
      return(list())
    }
    
    cell_types <- levels(cellchat@idents)
    spatial_results <- list()
    
    for (prob_level in names(spatial_colocalization_weights)) {
      prob_info <- spatial_colocalization_weights[[prob_level]]
      weight <- prob_info$weight
      
      level_interactions <- data.frame()
      
      for (pair in prob_info$pairs) {
        if (all(pair %in% cell_types)) {
          pair_net <- subsetCommunication(cellchat, 
                                        sources.use = pair[1], 
                                        targets.use = pair[2])
          
          if (nrow(pair_net) > 0) {
            pair_net$spatial_weight <- weight
            pair_net$spatial_probability <- prob_level
            level_interactions <- rbind(level_interactions, pair_net)
          }
        }
      }
      
      if (nrow(level_interactions) > 0) {
        spatial_results[[prob_level]] <- list(
          interactions = level_interactions,
          total_count = nrow(level_interactions),
          avg_prob = mean(level_interactions$prob, na.rm = TRUE),
          weighted_prob = mean(level_interactions$prob * level_interactions$spatial_weight, na.rm = TRUE)
        )
        
        message(paste("  ", prob_level, ":", nrow(level_interactions), "个相互作用"))
      }
    }
    
    if (length(spatial_results) == 0) {
      message("  未检测到预定义的空间共定位通讯")
    }
    
    return(spatial_results)
    
  }, "spatial_colocation", cancer_type, list()))
}

# 信号通路分类分析
analyze_pathway_categories <- function(cellchat, cancer_type) {
  message(paste("信号通路分类分析:", cancer_type))
  
  return(safe_execute({
    if (is.null(cellchat@netP$pathways) || length(cellchat@netP$pathways) == 0) {
      message("  无有效通路信息，跳过通路分类分析")
      return(list())
    }
    
    detected_pathways <- cellchat@netP$pathways
    category_results <- list()
    
    for (category_name in names(pathway_categories)) {
      category_patterns <- pathway_categories[[category_name]]
      
      # 匹配通路
      matched_pathways <- detected_pathways[
        sapply(detected_pathways, function(p) {
          any(sapply(category_patterns, function(pattern) {
            grepl(pattern, p, ignore.case = TRUE)
          }))
        })
      ]
      
      if (length(matched_pathways) > 0) {
        # 提取这些通路的通讯
        category_net <- subsetCommunication(cellchat, 
                                           signaling = matched_pathways)
        
        if (nrow(category_net) > 0) {
          category_results[[category_name]] <- list(
            pathways = matched_pathways,
            total_interactions = nrow(category_net),
            avg_prob = mean(category_net$prob, na.rm = TRUE),
            interactions = category_net
          )
          
          message(paste("  ", category_name, ":", 
                       length(matched_pathways), "个通路,", 
                       nrow(category_net), "个相互作用"))
        }
      }
    }
    
    if (length(category_results) == 0) {
      message("  未检测到预定义类别的通路")
    }
    
    return(category_results)
    
  }, "pathway_categories", cancer_type, list()))
}

# 网络中心性分析
analyze_network_centrality <- function(cellchat, cancer_type) {
  message(paste("网络中心性分析:", cancer_type))
  
  return(safe_execute({
    if (is.null(cellchat@net$weight) || sum(cellchat@net$weight > 0, na.rm = TRUE) == 0) {
      message("  无有效网络权重，跳过中心性分析")
      return(list())
    }
    
    # 获取网络权重矩阵
    weight_matrix <- cellchat@net$weight
    
    # 创建igraph对象
    g <- graph_from_adjacency_matrix(weight_matrix, mode = "directed", weighted = TRUE)
    
    # 计算各种中心性指标
    centrality_metrics <- list(
      degree_in = degree(g, mode = "in"),
      degree_out = degree(g, mode = "out"),
      degree_total = degree(g, mode = "all"),
      betweenness = betweenness(g, directed = TRUE),
      closeness_in = closeness(g, mode = "in"),
      closeness_out = closeness(g, mode = "out"),
      eigenvector = eigen_centrality(g, directed = TRUE)$vector,
      pagerank = page_rank(g)$vector
    )
    
    # 转换为数据框
    centrality_df <- data.frame(
      cell_type = V(g)$name,
      in_degree = centrality_metrics$degree_in,
      out_degree = centrality_metrics$degree_out,
      total_degree = centrality_metrics$degree_total,
      betweenness = centrality_metrics$betweenness,
      closeness_in = centrality_metrics$closeness_in,
      closeness_out = centrality_metrics$closeness_out,
      eigenvector = centrality_metrics$eigenvector,
      pagerank = centrality_metrics$pagerank
    )
    
    # 识别hub细胞（高度中心性）
    hub_threshold <- quantile(centrality_df$total_degree, 0.75, na.rm = TRUE)
    centrality_df$is_hub <- centrality_df$total_degree >= hub_threshold
    
    message(paste("  Hub细胞类型:", sum(centrality_df$is_hub)))
    message(paste("  平均度数:", round(mean(centrality_df$total_degree), 2)))
    
    return(list(
      metrics = centrality_df,
      graph = g,
      hub_cells = centrality_df$cell_type[centrality_df$is_hub]
    ))
    
  }, "network_centrality", cancer_type, list()))
}

# 可视化功能模块网络
plot_functional_module_network <- function(module_results, cancer_type) {
  message(paste("绘制功能模块网络:", cancer_type))
  
  return(safe_execute({
    if (length(module_results) == 0) {
      message("  无功能模块结果，跳过可视化")
      return(NULL)
    }
    
    # 准备模块统计数据
    module_stats <- data.frame(
      Module = names(module_results),
      Interactions = sapply(module_results, function(x) x$total_interactions),
      AvgProb = sapply(module_results, function(x) x$avg_prob),
      Pathways = sapply(module_results, function(x) x$unique_pathways)
    )
    
    # 创建气泡图
    p1 <- ggplot(module_stats, aes(x = Module, y = Interactions)) +
      geom_point(aes(size = AvgProb, color = Pathways)) +
      scale_size_continuous(name = "平均概率") +
      scale_color_viridis_c(name = "通路数") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(cancer_type, "- 功能模块通讯概览"),
           x = "功能模块", y = "相互作用数")
    
    # 保存图片
    ggsave(file.path(communication_output_dir, "functional_modules", 
                    paste0(cancer_type, "_module_overview.pdf")),
           p1, width = 10, height = 6)
    
    # 创建详细的网络图
    all_interactions <- do.call(rbind, lapply(names(module_results), function(m) {
      df <- module_results[[m]]$interactions
      df$module <- m
      return(df)
    }))
    
    if (nrow(all_interactions) > 0) {
      # 创建边数据框
      edges <- all_interactions %>%
        group_by(source, target, module) %>%
        summarise(weight = sum(prob), .groups = "drop")
      
      # 创建节点数据框
      all_cells <- unique(c(edges$source, edges$target))
      nodes <- data.frame(
        name = all_cells,
        module = sapply(all_cells, function(cell) {
          for (m in names(functional_modules)) {
            if (cell %in% functional_modules[[m]]) return(m)
          }
          return("Other")
        })
      )
      
      # 创建igraph对象
      g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
      
      # 使用ggraph绘制网络
      p2 <- ggraph(g, layout = "fr") +
        geom_edge_link(aes(width = weight, alpha = weight), 
                      arrow = arrow(length = unit(2, 'mm')), 
                      end_cap = circle(3, 'mm')) +
        geom_node_point(aes(color = module), size = 5) +
        geom_node_text(aes(label = name), repel = TRUE, size = 3) +
        scale_edge_width(range = c(0.5, 2)) +
        scale_edge_alpha(range = c(0.3, 0.8)) +
        scale_color_brewer(palette = "Set3") +
        theme_void() +
        labs(title = paste(cancer_type, "- 功能模块通讯网络"))
      
      ggsave(file.path(communication_output_dir, "functional_modules", 
                      paste0(cancer_type, "_module_network.pdf")),
             p2, width = 12, height = 10)
      
      return(list(overview = p1, network = p2))
    }
    
    return(list(overview = p1))
    
  }, "plot_functional_modules", cancer_type, NULL))
}

# 可视化空间共定位分析
plot_spatial_colocation <- function(spatial_results, cancer_type) {
  message(paste("绘制空间共定位图:", cancer_type))
  
  return(safe_execute({
    if (length(spatial_results) == 0) {
      message("  无空间共定位结果，跳过可视化")
      return(NULL)
    }
    
    # 准备数据
    spatial_stats <- data.frame(
      Probability_Level = names(spatial_results),
      Interactions = sapply(spatial_results, function(x) x$total_count),
      AvgProb = sapply(spatial_results, function(x) x$avg_prob),
      WeightedProb = sapply(spatial_results, function(x) x$weighted_prob)
    )
    
    # 设置因子顺序
    spatial_stats$Probability_Level <- factor(
      spatial_stats$Probability_Level,
      levels = c("High_Probability", "Medium_Probability", "Low_Probability")
    )
    
    # 创建条形图
    p1 <- ggplot(spatial_stats, aes(x = Probability_Level, y = Interactions, fill = WeightedProb)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = Interactions), vjust = -0.5) +
      scale_fill_gradient(low = "lightblue", high = "darkblue", name = "加权概率") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(cancer_type, "- 空间共定位通讯分布"),
           x = "共定位概率", y = "相互作用数")
    
    ggsave(file.path(communication_output_dir, "spatial_colocation", 
                    paste0(cancer_type, "_spatial_distribution.pdf")),
           p1, width = 8, height = 6)
    
    # 创建热图展示细胞对的通讯强度
    all_interactions <- do.call(rbind, lapply(names(spatial_results), function(level) {
      df <- spatial_results[[level]]$interactions
      df$probability_level <- level
      return(df)
    }))
    
    if (nrow(all_interactions) > 0) {
      # 汇总细胞对的通讯
      pair_summary <- all_interactions %>%
        group_by(source, target) %>%
        summarise(
          total_prob = sum(prob * spatial_weight),
          count = n(),
          .groups = "drop"
        )
      
      # 创建矩阵
      all_cells <- unique(c(pair_summary$source, pair_summary$target))
      comm_matrix <- matrix(0, nrow = length(all_cells), ncol = length(all_cells),
                           dimnames = list(all_cells, all_cells))
      
      for (i in 1:nrow(pair_summary)) {
        comm_matrix[pair_summary$source[i], pair_summary$target[i]] <- pair_summary$total_prob[i]
      }
      
      # 绘制热图
      pdf(file.path(communication_output_dir, "spatial_colocation", 
                   paste0(cancer_type, "_spatial_heatmap.pdf")),
          width = 10, height = 9)
      
      col_fun <- colorRamp2(c(0, max(comm_matrix)/2, max(comm_matrix)), 
                           c("white", "yellow", "red"))
      
      ht <- Heatmap(comm_matrix,
                   name = "加权通讯强度",
                   col = col_fun,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   column_title = paste(cancer_type, "- 空间共定位通讯热图"),
                   row_names_gp = gpar(fontsize = 10),
                   column_names_gp = gpar(fontsize = 10))
      
      draw(ht)
      dev.off()
      
      return(list(distribution = p1, heatmap_data = comm_matrix))
    }
    
    return(list(distribution = p1))
    
  }, "plot_spatial_colocation", cancer_type, NULL))
}

# 可视化通路分类
plot_pathway_categories <- function(pathway_results, cancer_type) {
  message(paste("绘制通路分类图:", cancer_type))
  
  return(safe_execute({
    if (length(pathway_results) == 0) {
      message("  无通路分类结果，跳过可视化")
      return(NULL)
    }
    
    # 准备数据
    pathway_stats <- data.frame(
      Category = names(pathway_results),
      Pathways = sapply(pathway_results, function(x) length(x$pathways)),
      Interactions = sapply(pathway_results, function(x) x$total_interactions),
      AvgProb = sapply(pathway_results, function(x) x$avg_prob)
    )
    
    # 创建堆叠条形图
    p1 <- ggplot(pathway_stats, aes(x = reorder(Category, -Interactions), y = Interactions)) +
      geom_bar(aes(fill = AvgProb), stat = "identity") +
      geom_text(aes(label = Pathways), vjust = -0.5, size = 3) +
      scale_fill_viridis_c(name = "平均概率") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(cancer_type, "- 信号通路分类分析"),
           subtitle = "数字表示通路数量",
           x = "通路类别", y = "相互作用数")
    
    ggsave(file.path(communication_output_dir, "pathway_analysis", 
                    paste0(cancer_type, "_pathway_categories.pdf")),
           p1, width = 10, height = 6)
    
    # 创建详细的通路网络图
    # 选择top通路类别
    top_categories <- pathway_stats %>%
      arrange(desc(Interactions)) %>%
      head(5) %>%
      pull(Category)
    
    for (category in top_categories) {
      cat_data <- pathway_results[[category]]
      
      if (nrow(cat_data$interactions) > 0) {
        # 创建通路网络
        pathway_edges <- cat_data$interactions %>%
          group_by(source, target, pathway_name) %>%
          summarise(weight = sum(prob), .groups = "drop")
        
        # 限制显示的边数
        if (nrow(pathway_edges) > 50) {
          pathway_edges <- pathway_edges %>%
            arrange(desc(weight)) %>%
            head(50)
        }
        
        # 创建节点
        all_cells <- unique(c(pathway_edges$source, pathway_edges$target))
        nodes <- data.frame(name = all_cells)
        
        # 创建igraph
        g <- graph_from_data_frame(pathway_edges, directed = TRUE, vertices = nodes)
        
        # 绘图
        p_cat <- ggraph(g, layout = "stress") +
          geom_edge_link(aes(width = weight, color = pathway_name, alpha = weight),
                        arrow = arrow(length = unit(2, 'mm')),
                        end_cap = circle(3, 'mm')) +
          geom_node_point(size = 4, color = "steelblue") +
          geom_node_text(aes(label = name), repel = TRUE, size = 3) +
          scale_edge_width(range = c(0.5, 2)) +
          scale_edge_alpha(range = c(0.3, 0.8)) +
          theme_void() +
          labs(title = paste(cancer_type, "-", category, "通路网络"),
               subtitle = paste("Top", nrow(pathway_edges), "个相互作用"))
        
        ggsave(file.path(communication_output_dir, "pathway_analysis", 
                        paste0(cancer_type, "_", category, "_network.pdf")),
               p_cat, width = 12, height = 10)
      }
    }
    
    return(list(overview = p1))
    
  }, "plot_pathway_categories", cancer_type, NULL))
}

# 可视化网络中心性
plot_network_centrality <- function(centrality_results, cancer_type) {
  message(paste("绘制网络中心性图:", cancer_type))
  
  return(safe_execute({
    if (is.null(centrality_results) || is.null(centrality_results$metrics)) {
      message("  无中心性结果，跳过可视化")
      return(NULL)
    }
    
    centrality_df <- centrality_results$metrics
    
    # 创建雷达图数据
    radar_data <- centrality_df %>%
      select(cell_type, total_degree, betweenness, eigenvector, pagerank) %>%
      mutate(across(where(is.numeric), ~scales::rescale(.x, to = c(0, 1))))
    
    # 创建散点图矩阵
    p1 <- ggplot(centrality_df, aes(x = in_degree, y = out_degree)) +
      geom_point(aes(size = betweenness, color = is_hub)) +
      geom_text(aes(label = cell_type), size = 3, vjust = -0.5) +
      scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                        name = "Hub细胞") +
      scale_size_continuous(name = "中介中心性") +
      theme_minimal() +
      labs(title = paste(cancer_type, "- 网络中心性分析"),
           x = "入度", y = "出度")
    
    ggsave(file.path(communication_output_dir, "network_analysis", 
                    paste0(cancer_type, "_centrality_scatter.pdf")),
           p1, width = 10, height = 8)
    
    # 创建中心性指标条形图
    centrality_long <- centrality_df %>%
      select(cell_type, total_degree, betweenness, eigenvector, pagerank) %>%
      pivot_longer(cols = -cell_type, names_to = "metric", values_to = "value") %>%
      group_by(metric) %>%
      mutate(value_scaled = scales::rescale(value, to = c(0, 1)))
    
    p2 <- ggplot(centrality_long, aes(x = reorder(cell_type, -value_scaled), 
                                      y = value_scaled, fill = metric)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_brewer(palette = "Set2", name = "中心性指标",
                       labels = c("总度数", "中介", "特征向量", "PageRank")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(cancer_type, "- 中心性指标比较"),
           x = "细胞类型", y = "标准化值")
    
    ggsave(file.path(communication_output_dir, "network_analysis", 
                    paste0(cancer_type, "_centrality_comparison.pdf")),
           p2, width = 12, height = 6)
    
    # 绘制网络图，标注hub细胞
    g <- centrality_results$graph
    V(g)$size <- scales::rescale(centrality_df$total_degree, to = c(3, 15))
    V(g)$color <- ifelse(centrality_df$is_hub, "red", "lightblue")
    
    pdf(file.path(communication_output_dir, "network_analysis", 
                 paste0(cancer_type, "_hub_network.pdf")),
        width = 12, height = 10)
    
    plot(g,
         vertex.label = V(g)$name,
         vertex.label.cex = 0.8,
         edge.arrow.size = 0.3,
         edge.width = E(g)$weight * 2,
         layout = layout_with_fr(g),
         main = paste(cancer_type, "- Hub细胞网络 (红色为Hub)"))
    
    dev.off()
    
    return(list(scatter = p1, comparison = p2))
    
  }, "plot_centrality", cancer_type, NULL))
}

# 创建交互式网络图
create_interactive_network <- function(cellchat, cancer_type) {
  message(paste("创建交互式网络:", cancer_type))
  
  return(safe_execute({
    if (is.null(cellchat@net$weight) || sum(cellchat@net$weight > 0, na.rm = TRUE) == 0) {
      message("  无有效网络数据，跳过交互式可视化")
      return(NULL)
    }
    
    # 获取网络数据
    weight_matrix <- cellchat@net$weight
    cell_types <- rownames(weight_matrix)
    
    # 创建边列表
    edges <- data.frame()
    for (i in 1:nrow(weight_matrix)) {
      for (j in 1:ncol(weight_matrix)) {
        if (weight_matrix[i, j] > 0) {
          edges <- rbind(edges, data.frame(
            source = i - 1,  # networkD3使用0-based索引
            target = j - 1,
            value = weight_matrix[i, j]
          ))
        }
    }
    }
    
    # 创建节点列表
    nodes <- data.frame(
      name = cell_types,
      group = 1:length(cell_types)
    )
    
    # 创建交互式网络图
    network <- forceNetwork(
      Links = edges,
      Nodes = nodes,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      Group = "group",
      opacity = 0.9,
      fontSize = 14,
      zoom = TRUE,
      legend = TRUE,
      bounded = FALSE,
      colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);")
    )
    
    # 保存为HTML
    saveNetwork(network, 
                file.path(communication_output_dir, "interactive_plots", 
                         paste0(cancer_type, "_interactive_network.html")))
    
    message("  交互式网络图已保存")
    return(network)
    
  }, "interactive_network", cancer_type, NULL))
}

# 生成综合报告
generate_comprehensive_report <- function(cancer_type, all_results) {
  message(paste("生成综合报告:", cancer_type))
  
  return(safe_execute({
    report_file <- file.path(communication_output_dir, "qc_reports", 
                            paste0(cancer_type, "_analysis_report.txt"))
    
    sink(report_file)
    
    cat("=" , rep("=", 70), "=\n", sep = "")
    cat("细胞通讯分析综合报告\n")
    cat("癌症类型:", cancer_type, "\n")
    cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("=" , rep("=", 70), "=\n\n", sep = "")
    
    # 基础统计
    if (!is.null(all_results$cellchat)) {
      cellchat <- all_results$cellchat
      cat("1. 基础网络统计\n")
      cat(rep("-", 72), "\n", sep = "")
      cat("细胞类型数:", length(levels(cellchat@idents)), "\n")
      
      if (!is.null(cellchat@net$count)) {
        cat("总通讯对数:", sum(cellchat@net$count > 0, na.rm = TRUE), "\n")
        cat("平均通讯强度:", round(mean(cellchat@net$count[cellchat@net$count > 0], na.rm = TRUE), 2), "\n")
      }
      
      if (!is.null(cellchat@netP$pathways)) {
        cat("检测到的信号通路:", length(cellchat@netP$pathways), "\n")
        cat("主要通路:", paste(head(cellchat@netP$pathways, 10), collapse = ", "), "\n")
      }
      cat("\n")
    }
    
    # 功能模块分析
    if (!is.null(all_results$functional_modules) && length(all_results$functional_modules) > 0) {
      cat("2. 功能模块分析\n")
      cat(rep("-", 72), "\n", sep = "")
      
      for (module in names(all_results$functional_modules)) {
        mod_data <- all_results$functional_modules[[module]]
        cat(sprintf("  %s:\n", module))
        cat(sprintf("    - 相互作用数: %d\n", mod_data$total_interactions))
        cat(sprintf("    - 平均概率: %.4f\n", mod_data$avg_prob))
        cat(sprintf("    - 涉及通路: %d\n", mod_data$unique_pathways))
        cat(sprintf("    - 细胞类型: %s\n", paste(mod_data$cell_types, collapse = ", ")))
      }
      cat("\n")
    }
    
    # 空间共定位分析
    if (!is.null(all_results$spatial_colocation) && length(all_results$spatial_colocation) > 0) {
      cat("3. 空间共定位分析\n")
      cat(rep("-", 72), "\n", sep = "")
      
      for (prob_level in names(all_results$spatial_colocation)) {
        spatial_data <- all_results$spatial_colocation[[prob_level]]
        cat(sprintf("  %s:\n", prob_level))
        cat(sprintf("    - 相互作用数: %d\n", spatial_data$total_count))
        cat(sprintf("    - 平均概率: %.4f\n", spatial_data$avg_prob))
        cat(sprintf("    - 加权概率: %.4f\n", spatial_data$weighted_prob))
      }
      cat("\n")
    }
    
    # 通路分类分析
    if (!is.null(all_results$pathway_categories) && length(all_results$pathway_categories) > 0) {
      cat("4. 信号通路分类\n")
      cat(rep("-", 72), "\n", sep = "")
      
      for (category in names(all_results$pathway_categories)) {
        cat_data <- all_results$pathway_categories[[category]]
        cat(sprintf("  %s:\n", category))
        cat(sprintf("    - 通路数: %d\n", length(cat_data$pathways)))
        cat(sprintf("    - 相互作用数: %d\n", cat_data$total_interactions))
        cat(sprintf("    - 平均概率: %.4f\n", cat_data$avg_prob))
        cat(sprintf("    - 通路列表: %s\n", paste(cat_data$pathways, collapse = ", ")))
      }
      cat("\n")
    }
    
    # 网络中心性分析
    if (!is.null(all_results$centrality) && !is.null(all_results$centrality$metrics)) {
      cat("5. 网络中心性分析\n")
      cat(rep("-", 72), "\n", sep = "")
      
      cent_df <- all_results$centrality$metrics
      hub_cells <- all_results$centrality$hub_cells
      
      cat("Hub细胞类型:\n")
      for (hub in hub_cells) {
        hub_data <- cent_df[cent_df$cell_type == hub, ]
        cat(sprintf("  %s:\n", hub))
        cat(sprintf("    - 总度数: %.2f\n", hub_data$total_degree))
        cat(sprintf("    - 中介中心性: %.4f\n", hub_data$betweenness))
        cat(sprintf("    - PageRank: %.4f\n", hub_data$pagerank))
      }
      cat("\n")
    }
    
    cat("=" , rep("=", 70), "=\n", sep = "")
    cat("报告结束\n")
    
    sink()
    
    message("  综合报告已生成")
    return(report_file)
    
  }, "generate_report", cancer_type, NULL))
}

# 保存所有结果
save_analysis_results <- function(cancer_type, all_results) {
  message(paste("保存分析结果:", cancer_type))
  
  return(safe_execute({
    # 保存RDS对象
    saveRDS(all_results, 
            file.path(communication_output_dir, "tables", 
                     paste0(cancer_type, "_complete_results.rds")))
    
    # 保存功能模块结果
    if (!is.null(all_results$functional_modules) && length(all_results$functional_modules) > 0) {
      for (module in names(all_results$functional_modules)) {
        mod_data <- all_results$functional_modules[[module]]
        if (!is.null(mod_data$interactions) && nrow(mod_data$interactions) > 0) {
          write.csv(mod_data$interactions,
                   file.path(communication_output_dir, "functional_modules",
                            paste0(cancer_type, "_", module, "_interactions.csv")),
                   row.names = FALSE)
        }
      }
    }
    
    # 保存空间共定位结果
    if (!is.null(all_results$spatial_colocation) && length(all_results$spatial_colocation) > 0) {
      for (prob_level in names(all_results$spatial_colocation)) {
        spatial_data <- all_results$spatial_colocation[[prob_level]]
        if (!is.null(spatial_data$interactions) && nrow(spatial_data$interactions) > 0) {
          write.csv(spatial_data$interactions,
                   file.path(communication_output_dir, "spatial_colocation",
                            paste0(cancer_type, "_", prob_level, "_interactions.csv")),
                   row.names = FALSE)
        }
      }
    }
    
    # 保存通路分类结果
    if (!is.null(all_results$pathway_categories) && length(all_results$pathway_categories) > 0) {
      for (category in names(all_results$pathway_categories)) {
        cat_data <- all_results$pathway_categories[[category]]
        if (!is.null(cat_data$interactions) && nrow(cat_data$interactions) > 0) {
          write.csv(cat_data$interactions,
                   file.path(communication_output_dir, "pathway_analysis",
                            paste0(cancer_type, "_", category, "_interactions.csv")),
                   row.names = FALSE)
        }
      }
    }
    
    # 保存中心性分析结果
    if (!is.null(all_results$centrality) && !is.null(all_results$centrality$metrics)) {
      write.csv(all_results$centrality$metrics,
               file.path(communication_output_dir, "network_analysis",
                        paste0(cancer_type, "_centrality_metrics.csv")),
               row.names = FALSE)
    }
    
    message("  所有结果已保存")
    
  }, "save_results", cancer_type, NULL))
}

# 主分析流程
analyze_cancer_communications <- function(cancer_type) {
  message("\n")
  message(rep("=", 80))
  message(paste("开始分析:", cancer_type))
  message(rep("=", 80))
  
  # 存储所有结果
  all_results <- list()
  
  # 1. 加载数据
  seurat_obj <- load_cancer_data(cancer_type)
  if (is.null(seurat_obj)) {
    message(paste("跳过", cancer_type, "- 数据加载失败"))
    return(NULL)
  }
  
  # 2. 创建CellChat对象
  cellchat <- create_cellchat_object(seurat_obj, cancer_type)
  if (is.null(cellchat)) {
    message(paste("跳过", cancer_type, "- CellChat对象创建失败"))
    return(NULL)
  }
  
  # 清理Seurat对象以释放内存
  rm(seurat_obj)
  gc()
  
  # 3. 分析CellChat
  cellchat <- analyze_cellchat(cellchat, cancer_type)
  if (is.null(cellchat)) {
    message(paste("跳过", cancer_type, "- CellChat分析失败"))
    return(NULL)
  }
  
  all_results$cellchat <- cellchat
  
  # 4. 功能模块分析
  all_results$functional_modules <- analyze_functional_modules(cellchat, cancer_type)
  
  # 5. 空间共定位分析
  all_results$spatial_colocation <- analyze_spatial_colocation(cellchat, cancer_type)
  
  # 6. 通路分类分析
  all_results$pathway_categories <- analyze_pathway_categories(cellchat, cancer_type)
  
  # 7. 网络中心性分析
  all_results$centrality <- analyze_network_centrality(cellchat, cancer_type)
  
  # 8. 可视化
  plot_functional_module_network(all_results$functional_modules, cancer_type)
  plot_spatial_colocation(all_results$spatial_colocation, cancer_type)
  plot_pathway_categories(all_results$pathway_categories, cancer_type)
  plot_network_centrality(all_results$centrality, cancer_type)
  
  # 9. 创建交互式网络
  create_interactive_network(cellchat, cancer_type)
  
  # 10. 生成报告
  generate_comprehensive_report(cancer_type, all_results)
  
  # 11. 保存结果
  save_analysis_results(cancer_type, all_results)
  
  message(paste("完成分析:", cancer_type))
  message(rep("=", 80))
  message("\n")
  
  # 清理内存
  gc()
  
  return(all_results)
}

# ============================================================================
# 主执行部分
# ============================================================================

# 获取所有癌症类型
cancer_types <- list.dirs(output_base_dir, full.names = FALSE, recursive = FALSE)
cancer_types <- cancer_types[cancer_types != "cell_communication_analysis"]

message(paste("发现", length(cancer_types), "个癌症类型"))
message(paste("癌症类型:", paste(cancer_types, collapse = ", ")))

# 分析每个癌症类型
all_cancer_results <- list()

for (cancer_type in cancer_types) {
  tryCatch({
    result <- analyze_cancer_communications(cancer_type)
    if (!is.null(result)) {
      all_cancer_results[[cancer_type]] <- result
    }
    
    # 每个癌症分析后强制垃圾回收
    gc()
    
  }, error = function(e) {
    message(paste("ERROR: 分析", cancer_type, "时发生严重错误:", e$message))
    error_collector[[paste("main", cancer_type, sep = "_")]] <- list(
      timestamp = Sys.time(),
      step = "main_analysis",
      cancer_type = cancer_type,
      error = e$message
    )
  })
}

# 保存所有错误日志
if (length(error_collector) > 0) {
  error_log_file <- file.path(communication_output_dir, "qc_reports", "error_log.rds")
  saveRDS(error_collector, error_log_file)
  
  # 生成错误摘要
  error_summary_file <- file.path(communication_output_dir, "qc_reports", "error_summary.txt")
  sink(error_summary_file)
  
  cat("=" , rep("=", 70), "=\n", sep = "")
  cat("错误日志摘要\n")
  cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("=" , rep("=", 70), "=\n\n", sep = "")
  
  for (err_name in names(error_collector)) {
    err_info <- error_collector[[err_name]]
    cat(sprintf("错误ID: %s\n", err_name))
    cat(sprintf("  时间: %s\n", err_info$timestamp))
    cat(sprintf("  步骤: %s\n", err_info$step))
    cat(sprintf("  癌症类型: %s\n", err_info$cancer_type))
    cat(sprintf("  错误信息: %s\n\n", err_info$error))
  }
  
  cat("=" , rep("=", 70), "=\n", sep = "")
  cat(sprintf("总错误数: %d\n", length(error_collector)))
  
  sink()
  
  message(paste("发现", length(error_collector), "个错误，详见:", error_log_file))
}

# 生成跨癌症比较摘要
if (length(all_cancer_results) > 1) {
  message("生成跨癌症比较摘要...")
  
  comparison_file <- file.path(communication_output_dir, "qc_reports", "cross_cancer_summary.txt")
  sink(comparison_file)
  
  cat("=" , rep("=", 70), "=\n", sep = "")
  cat("跨癌症细胞通讯比较摘要\n")
  cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("=" , rep("=", 70), "=\n\n", sep = "")
  
  for (cancer in names(all_cancer_results)) {
    cat(sprintf("\n%s\n", cancer))
    cat(rep("-", 72), "\n", sep = "")
    
    result <- all_cancer_results[[cancer]]
    
    if (!is.null(result$cellchat) && !is.null(result$cellchat@net$count)) {
      cat(sprintf("  总通讯数: %d\n", sum(result$cellchat@net$count > 0, na.rm = TRUE)))
    }
    
    if (!is.null(result$functional_modules)) {
      cat(sprintf("  功能模块数: %d\n", length(result$functional_modules)))
    }
    
    if (!is.null(result$spatial_colocation)) {
      cat(sprintf("  空间共定位类型: %d\n", length(result$spatial_colocation)))
    }
    
    if (!is.null(result$pathway_categories)) {
      cat(sprintf("  通路类别数: %d\n", length(result$pathway_categories)))
    }
    
    if (!is.null(result$centrality) && !is.null(result$centrality$hub_cells)) {
      cat(sprintf("  Hub细胞类型: %s\n", paste(result$centrality$hub_cells, collapse = ", ")))
    }
  }
  
  cat("\n")
  cat("=" , rep("=", 70), "=\n", sep = "")
  
  sink()
  
  message("跨癌症比较摘要已生成")
}

message("\n")
message(rep("=", 80))
message("所有分析完成!")
message(paste("成功分析:", length(all_cancer_results), "个癌症类型"))
message(paste("输出目录:", communication_output_dir))
message(rep("=", 80))