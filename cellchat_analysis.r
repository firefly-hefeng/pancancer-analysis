# cell_communication_analysis.R
# 基于功能模块和空间共定位的细胞通讯分析

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

# 设置并行计算
future::plan(multicore, workers = 6)

# 参数设置
output_base_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation"
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
    # 高共定位可能性组合
    pairs = list(
      c("T_cells", "B_cell"),           # 淋巴滤泡
      c("T_cells", "Dendritic_cells"),  # 免疫突触
      c("Endothelial_cells", "Pericytes"), # 血管壁
      c("Fibroblasts", "Epithelial_cells"), # 基质-上皮界面
      c("Macrophages", "T_cells"),      # 炎症微环境
      c("CAFs", "Fibroblasts")          # 基质细胞群
    ),
    weight = 1.0
  ),
  "Medium_Probability" = list(
    pairs = list(
      c("NK_cell", "T_cells"),          # 免疫效应区
      c("Macrophages", "Endothelial_cells"), # 血管周围
      c("Fibroblasts", "Endothelial_cells"), # 血管支持
      c("B_cell", "Dendritic_cells")    # 抗原呈递
    ),
    weight = 0.7
  ),
  "Low_Probability" = list(
    pairs = list(
      c("Epithelial_cells", "T_cells"), # 上皮-免疫
      c("Epithelial_cells", "B_cell"),  # 上皮-体液免疫
      c("NK_cell", "Fibroblasts")       # NK-基质
    ),
    weight = 0.4
  )
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
      error = e$message
    )
    error_collector[[paste(step, cancer_type, sep = "_")]] <<- error_info
    message(paste("ERROR [", step, "]:", cancer_type, "-", e$message))
    return(default_return)
  })
}

# 数据加载函数
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
    
    # 简化细胞类型名称（提取主要类型）
    seurat_obj$main_cell_type <- gsub("_.*", "", seurat_obj@meta.data[[cell_type_col]])
    
    # 过滤低质量细胞和基因
    seurat_obj <- subset(seurat_obj, 
                        subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & 
                                percent.mt < 25)
    
    # 确保有足够的细胞进行分析
    cell_counts <- table(seurat_obj$main_cell_type)
    valid_types <- names(cell_counts)[cell_counts >= 50]
    
    if (length(valid_types) < 3) {
      stop("有效细胞类型过少")
    }
    
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj$main_cell_type %in% valid_types])
    
    message(paste("  成功加载:", ncol(seurat_obj), "细胞,", length(valid_types), "种细胞类型"))
    
    return(seurat_obj)
    
  }, "data_loading", cancer_type, NULL))
}

# 创建CellChat对象
create_cellchat_object <- function(seurat_obj, cancer_type) {
  message(paste("创建CellChat对象:", cancer_type))
  
  return(safe_execute({
    # 提取数据
    data.input <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
    meta <- seurat_obj@meta.data
    
    # 创建CellChat对象
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "main_cell_type")
    
    # 设置配体-受体数据库
    CellChatDB <- CellChatDB.human
    
    # 根据细胞类型选择合适的数据库子集
    available_cell_types <- unique(seurat_obj$main_cell_type)
    
    if (any(c("T_cells", "B_cell", "NK_cell", "Myeloid") %in% available_cell_types)) {
      # 如果有免疫细胞，使用完整数据库
      cellchat@DB <- CellChatDB
      message("  使用完整CellChatDB数据库")
    } else {
      # 否则使用分泌信号数据库
      cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")
      message("  使用分泌信号数据库")
    }
    
    return(cellchat)
    
  }, "create_cellchat", cancer_type, NULL))
}

# 预处理和分析CellChat
analyze_cellchat <- function(cellchat, cancer_type) {
  message(paste("分析细胞通讯:", cancer_type))
  
  return(safe_execute({
    # 预处理数据
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # 计算通讯概率
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    # 推断细胞通讯网络
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    # 识别信号通路
    cellchat <- netAnalysis_signalingRole_network(cellchat, slot.name = "netP", width = 8, height = 2.5)
    
    return(cellchat)
    
  }, "analyze_cellchat", cancer_type, NULL))
}

# 功能模块分析
analyze_functional_modules <- function(cellchat, cancer_type) {
  message(paste("功能模块分析:", cancer_type))
  
  return(safe_execute({
    # 获取细胞类型
    cell_types <- levels(cellchat@idents)
    
    # 为每个功能模块计算通讯强度
    module_results <- list()
    
    for (module_name in names(functional_modules)) {
      module_cells <- intersect(functional_modules[[module_name]], cell_types)
      
      if (length(module_cells) >= 2) {
        # 提取模块内的通讯
        module_net <- subsetCommunication(cellchat, 
                                        sources.use = module_cells, 
                                        targets.use = module_cells)
        
        if (nrow(module_net) > 0) {
          # 计算模块通讯统计
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
    
    return(module_results)
    
  }, "functional_modules", cancer_type, list()))
}

# 空间共定位分析
analyze_spatial_colocalization <- function(cellchat, cancer_type) {
  message(paste("空间共定位分析:", cancer_type))
  
  return(safe_execute({
    # 获取通讯网络
    net <- cellchat@net
    cell_types <- levels(cellchat@idents)
    
    # 创建空间权重矩阵
    spatial_weights <- matrix(0.1, nrow = length(cell_types), ncol = length(cell_types))
    rownames(spatial_weights) <- colnames(spatial_weights) <- cell_types
    
    # 应用空间共定位权重
    for (prob_level in names(spatial_colocalization_weights)) {
      pairs <- spatial_colocalization_weights[[prob_level]]$pairs
      weight <- spatial_colocalization_weights[[prob_level]]$weight
      
      for (pair in pairs) {
        if (all(pair %in% cell_types)) {
          spatial_weights[pair[1], pair[2]] <- weight
          spatial_weights[pair[2], pair[1]] <- weight
        }
      }
    }
    
    # 计算空间加权的通讯强度
    weighted_prob <- net$prob * spatial_weights[rownames(net$prob), colnames(net$prob)]
    weighted_count <- net$count * spatial_weights[rownames(net$count), colnames(net$count)]
    
    # 分层分析结果
    colocalization_results <- list(
      "High_Colocalization" = list(),
      "Medium_Colocalization" = list(),
      "Low_Colocalization" = list()
    )
    
    # 按共定位概率分组
    for (prob_level in names(spatial_colocalization_weights)) {
      pairs <- spatial_colocalization_weights[[prob_level]]$pairs
      level_name <- gsub("_Probability", "_Colocalization", prob_level)
      
      level_interactions <- c()
      level_strength <- c()
      
      for (pair in pairs) {
        if (all(pair %in% cell_types)) {
          # 双向通讯强度
          strength1 <- weighted_prob[pair[1], pair[2]]
          strength2 <- weighted_prob[pair[2], pair[1]]
          
          level_interactions <- c(level_interactions, 
                                paste(pair[1], "->", pair[2]),
                                paste(pair[2], "->", pair[1]))
          level_strength <- c(level_strength, strength1, strength2)
        }
      }
      
      if (length(level_interactions) > 0) {
        colocalization_results[[level_name]] <- list(
          interactions = level_interactions,
          strength = level_strength,
          avg_strength = mean(level_strength, na.rm = TRUE),
          n_interactions = length(level_interactions)
        )
      }
    }
    
    # 添加加权网络到结果
    colocalization_results$weighted_networks <- list(
      prob = weighted_prob,
      count = weighted_count,
      spatial_weights = spatial_weights
    )
    
    return(colocalization_results)
    
  }, "spatial_colocalization", cancer_type, list()))
}

# 信号通路分层分析
analyze_pathway_layers <- function(cellchat, cancer_type) {
  message(paste("信号通路分层分析:", cancer_type))
  
  return(safe_execute({
    # 获取所有通讯
    df.net <- subsetCommunication(cellchat)
    
    if (nrow(df.net) == 0) {
      return(list())
    }
    
    # 定义通路类别
    pathway_categories <- list(
      "Cytokine" = c("IL", "TNF", "IFN", "CSF", "CXCL", "CCL"),
      "Growth_Factor" = c("VEGF", "PDGF", "FGF", "TGF", "EGF", "IGF"),
      "Cell_Adhesion" = c("CDH", "ITGA", "ITGB", "NCAM", "PECAM"),
      "ECM_Interaction" = c("COL", "FN", "VTN", "THBS", "SPP1"),
      "Immune_Checkpoint" = c("CD28", "CD80", "CD86", "PDCD1", "CD274"),
      "Metabolism" = c("LDLR", "APOE", "LIPID")
    )
    
    # 为每个通路分类统计
    pathway_results <- list()
    
    for (category in names(pathway_categories)) {
      keywords <- pathway_categories[[category]]
      
      # 查找匹配的通路
      matching_pathways <- c()
      for (keyword in keywords) {
        matches <- grep(keyword, df.net$pathway_name, ignore.case = TRUE, value = TRUE)
        matching_pathways <- c(matching_pathways, matches)
      }
      matching_pathways <- unique(matching_pathways)
      
      if (length(matching_pathways) > 0) {
        category_df <- df.net[df.net$pathway_name %in% matching_pathways, ]
        
        pathway_results[[category]] <- list(
          pathways = matching_pathways,
          n_pathways = length(matching_pathways),
          n_interactions = nrow(category_df),
          avg_prob = mean(category_df$prob, na.rm = TRUE),
          top_interactions = category_df[order(category_df$prob, decreasing = TRUE)[1:min(10, nrow(category_df))], ],
          cell_type_involvement = table(c(category_df$source, category_df$target))
        )
        
        message(paste("  ", category, ":", length(matching_pathways), "通路,", nrow(category_df), "相互作用"))
      }
    }
    
    return(pathway_results)
    
  }, "pathway_analysis", cancer_type, list()))
}

# 网络分析和hub识别
analyze_communication_network <- function(cellchat, cancer_type) {
  message(paste("网络分析:", cancer_type))
  
  return(safe_execute({
    # 获取网络
    net <- cellchat@net
    
    # 计算网络指标
    cell_types <- levels(cellchat@idents)
    n_types <- length(cell_types)
    
    # 出度和入度
    outgoing <- rowSums(net$prob, na.rm = TRUE)
    incoming <- colSums(net$prob, na.rm = TRUE)
    
    # 网络中心性指标
    # 创建igraph对象
    adj_matrix <- net$prob
    adj_matrix[is.na(adj_matrix)] <- 0
    
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE)
    
    # 计算中心性指标
    centrality_metrics <- list(
      betweenness = betweenness(g, weights = E(g)$weight),
      closeness = closeness(g, weights = E(g)$weight),
      eigenvector = eigen_centrality(g, weights = E(g)$weight)$vector,
      page_rank = page_rank(g, weights = E(g)$weight)$vector
    )
    
    # 识别hub细胞类型
    hub_threshold <- 0.7  # 前30%
    
    hub_analysis <- list()
    for (metric in names(centrality_metrics)) {
      values <- centrality_metrics[[metric]]
      threshold <- quantile(values, hub_threshold, na.rm = TRUE)
      hubs <- names(values)[values >= threshold]
      
      hub_analysis[[metric]] <- list(
        hubs = hubs,
        values = values[hubs],
        threshold = threshold
      )
    }
    
    # 综合网络统计
    network_stats <- list(
      n_cell_types = n_types,
      total_interactions = sum(net$count > 0, na.rm = TRUE),
      avg_prob = mean(net$prob[net$prob > 0], na.rm = TRUE),
      network_density = sum(net$prob > 0, na.rm = TRUE) / (n_types * (n_types - 1)),
      outgoing_strength = outgoing,
      incoming_strength = incoming,
      centrality_metrics = centrality_metrics,
      hub_analysis = hub_analysis,
      igraph_object = g
    )
    
    return(network_stats)
    
  }, "network_analysis", cancer_type, list()))
}

# 生成功能模块可视化
plot_functional_modules <- function(module_results, cancer_type) {
  message(paste("绘制功能模块图:", cancer_type))
  
  return(safe_execute({
    if (length(module_results) == 0) {
      return(NULL)
    }
    
    # 准备数据
    module_data <- data.frame(
      Module = names(module_results),
      Total_Interactions = sapply(module_results, function(x) x$total_interactions),
      Avg_Probability = sapply(module_results, function(x) x$avg_prob),
      Unique_Pathways = sapply(module_results, function(x) x$unique_pathways),
      stringsAsFactors = FALSE
    )
    
    # 1. 相互作用数量条形图
    p1 <- ggplot(module_data, aes(x = reorder(Module, Total_Interactions), y = Total_Interactions)) +
      geom_col(aes(fill = Module), alpha = 0.8) +
      coord_flip() +
      scale_fill_viridis_d(option = "plasma") +
      labs(title = paste(cancer_type, "- Functional Module Interactions"),
           x = "Functional Module", y = "Total Interactions") +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    # 2. 通讯概率散点图
    p2 <- ggplot(module_data, aes(x = Unique_Pathways, y = Avg_Probability)) +
      geom_point(aes(size = Total_Interactions, color = Module), alpha = 0.7) +
      scale_color_viridis_d(option = "plasma") +
      scale_size_continuous(range = c(3, 10), name = "Total\nInteractions") +
      labs(title = "Module Communication Quality",
           x = "Unique Pathways", y = "Average Probability") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    
    # 3. 模块网络图
    # 创建模块间连接数据
    module_connections <- data.frame()
    
    for (i in 1:(length(module_results)-1)) {
      for (j in (i+1):length(module_results)) {
        module1 <- names(module_results)[i]
        module2 <- names(module_results)[j]
        
        cells1 <- module_results[[i]]$cell_types
        cells2 <- module_results[[j]]$cell_types
        
        # 检查是否有共同细胞类型
        shared_cells <- intersect(cells1, cells2)
        connection_strength <- length(shared_cells)
        
        if (connection_strength > 0) {
          module_connections <- rbind(module_connections, 
                                    data.frame(from = module1, to = module2, 
                                             strength = connection_strength))
        }
      }
    }
    
    if (nrow(module_connections) > 0) {
      # 创建网络图
      g_modules <- graph_from_data_frame(module_connections, directed = FALSE)
      
      # 设置节点属性
      V(g_modules)$size <- sapply(names(module_results), 
                                function(x) sqrt(module_results[[x]]$total_interactions) * 3)
      V(g_modules)$color <- viridis(length(module_results))
      
      # 设置边属性
      E(g_modules)$width <- module_connections$strength * 2
      
      # 网络布局
      layout_coords <- layout_with_fr(g_modules)
      
      p3 <- ggraph(g_modules, layout = layout_coords) +
        geom_edge_link(aes(width = strength), alpha = 0.6, color = "gray50") +
        geom_node_point(aes(size = size, color = name), alpha = 0.8) +
        geom_node_text(aes(label = name), size = 3, vjust = -1.2) +
        scale_color_viridis_d(option = "plasma", name = "Module") +
        scale_size_identity() +
        scale_edge_width_identity() +
        labs(title = "Functional Module Network") +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
              legend.position = "none")
    } else {
      p3 <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "No module connections found") +
        theme_void()
    }
    
    # 组合图形
    combined_plot <- (p1 | p2) / p3
    combined_plot <- combined_plot + 
      plot_annotation(title = paste(cancer_type, "- Functional Module Analysis"),
                     theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
    
    # 保存图形
    plot_file <- file.path(communication_output_dir, "functional_modules",
                          paste0(cancer_type, "_functional_modules.png"))
    ggsave(plot_file, combined_plot, width = 16, height = 12, dpi = 300, bg = "white")
    
    plot_file_pdf <- file.path(communication_output_dir, "functional_modules",
                              paste0(cancer_type, "_functional_modules.pdf"))
    ggsave(plot_file_pdf, combined_plot, width = 16, height = 12, bg = "white")
    
    return(combined_plot)
    
  }, "plot_functional_modules", cancer_type, NULL))
}

# 生成空间共定位可视化
plot_spatial_colocalization <- function(colocalization_results, cancer_type) {
  message(paste("绘制空间共定位图:", cancer_type))
  
  return(safe_execute({
    if (length(colocalization_results) == 0 || 
        is.null(colocalization_results$weighted_networks)) {
      return(NULL)
    }
    
    # 获取权重矩阵和概率矩阵
    spatial_weights <- colocalization_results$weighted_networks$spatial_weights
    weighted_prob <- colocalization_results$weighted_networks$prob
    
    # 1. 空间权重热图
    weight_df <- melt(spatial_weights)
    colnames(weight_df) <- c("Source", "Target", "Weight")
    
    p1 <- ggplot(weight_df, aes(x = Source, y = Target, fill = Weight)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma", name = "Spatial\nWeight") +
      labs(title = "Spatial Colocalization Weights", x = "Source", y = "Target") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    
    # 2. 加权通讯概率热图
    prob_df <- melt(weighted_prob)
    colnames(prob_df) <- c("Source", "Target", "Probability")
    prob_df$Probability[is.na(prob_df$Probability)] <- 0
    
    p2 <- ggplot(prob_df, aes(x = Source, y = Target, fill = Probability)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = median(prob_df$Probability[prob_df$Probability > 0]),
                          name = "Weighted\nProbability") +
      labs(title = "Spatially Weighted Communication", x = "Source", y = "Target") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    
    # 3. 分层统计柱状图
    layer_stats <- data.frame(
      Layer = character(),
      Avg_Strength = numeric(),
      N_Interactions = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (layer in c("High_Colocalization", "Medium_Colocalization", "Low_Colocalization")) {
      if (layer %in% names(colocalization_results) && 
          !is.null(colocalization_results[[layer]]$avg_strength)) {
        layer_stats <- rbind(layer_stats, 
                           data.frame(
                             Layer = layer,
                             Avg_Strength = colocalization_results[[layer]]$avg_strength,
                             N_Interactions = colocalization_results[[layer]]$n_interactions
                           ))
      }
    }
    
    if (nrow(layer_stats) > 0) {
      layer_stats$Layer <- factor(layer_stats$Layer, 
                                levels = c("High_Colocalization", "Medium_Colocalization", "Low_Colocalization"))
      
      p3 <- ggplot(layer_stats, aes(x = Layer, y = Avg_Strength)) +
        geom_col(aes(fill = Layer), alpha = 0.8) +
        geom_text(aes(label = paste0("n=", N_Interactions)), vjust = -0.5) +
        scale_fill_manual(values = c("High_Colocalization" = "#440154", 
                                   "Medium_Colocalization" = "#31688e", 
                                   "Low_Colocalization" = "#35b779")) +
        labs(title = "Communication Strength by Spatial Layer", 
             x = "Colocalization Layer", y = "Average Strength") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
              legend.position = "none")
    } else {
      p3 <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "No layer data available") +
        theme_void()
    }
    
    # 组合图形
    combined_plot <- (p1 | p2) / p3
    combined_plot <- combined_plot + 
      plot_annotation(title = paste(cancer_type, "- Spatial Colocalization Analysis"),
                     theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
    
    # 保存图形
    plot_file <- file.path(communication_output_dir, "spatial_colocation",
                          paste0(cancer_type, "_spatial_colocalization.png"))
    ggsave(plot_file, combined_plot, width = 16, height = 12, dpi = 300, bg = "white")
    
    plot_file_pdf <- file.path(communication_output_dir, "spatial_colocation",
                              paste0(cancer_type, "_spatial_colocalization.pdf"))
    ggsave(plot_file_pdf, combined_plot, width = 16, height = 12, bg = "white")
    
    return(combined_plot)
    
  }, "plot_spatial_colocalization", cancer_type, NULL))
}

# 生成网络分析可视化
plot_network_analysis <- function(network_stats, cancer_type) {
  message(paste("绘制网络分析图:", cancer_type))
  
  return(safe_execute({
    if (length(network_stats) == 0) {
      return(NULL)
    }
    
    # 1. 细胞类型强度图
    strength_data <- data.frame(
      Cell_Type = names(network_stats$outgoing_strength),
      Outgoing = network_stats$outgoing_strength,
      Incoming = network_stats$incoming_strength,
      stringsAsFactors = FALSE
    )
    
    strength_long <- melt(strength_data, id.vars = "Cell_Type", 
                         variable.name = "Direction", value.name = "Strength")
    
    p1 <- ggplot(strength_long, aes(x = reorder(Cell_Type, Strength), y = Strength, fill = Direction)) +
      geom_col(position = "dodge", alpha = 0.8) +
      coord_flip() +
      scale_fill_manual(values = c("Outgoing" = "#e31a1c", "Incoming" = "#1f78b4")) +
      labs(title = "Cell Type Communication Strength", 
           x = "Cell Type", y = "Communication Strength") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    
    # 2. 中心性指标图
    centrality_data <- data.frame(
      Cell_Type = names(network_stats$centrality_metrics$betweenness),
      Betweenness = network_stats$centrality_metrics$betweenness,
      Closeness = network_stats$centrality_metrics$closeness,
      Eigenvector = network_stats$centrality_metrics$eigenvector,
      PageRank = network_stats$centrality_metrics$page_rank,
      stringsAsFactors = FALSE
    )
    
    centrality_long <- melt(centrality_data, id.vars = "Cell_Type",
                           variable.name = "Metric", value.name = "Value")
    
    p2 <- ggplot(centrality_long, aes(x = Cell_Type, y = Value, fill = Metric)) +
      geom_col(position = "dodge", alpha = 0.8) +
      scale_fill_viridis_d(option = "plasma") +
      labs(title = "Network Centrality Metrics", x = "Cell Type", y = "Centrality Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    
    # 3. 网络图
    g <- network_stats$igraph_object
    
    # 设置节点大小和颜色
    V(g)$size <- (network_stats$outgoing_strength + network_stats$incoming_strength) * 10
    V(g)$color <- viridis(length(V(g)))
    
    # 设置边宽度
    edge_weights <- E(g)$weight
    E(g)$width <- (edge_weights / max(edge_weights, na.rm = TRUE)) * 5
    
    # 网络布局
    layout_coords <- layout_with_fr(g)
    
    p3 <- ggraph(g, layout = layout_coords) +
      geom_edge_link(aes(width = width), alpha = 0.6, color = "gray50") +
      geom_node_point(aes(size = size, color = name), alpha = 0.8) +
      geom_node_text(aes(label = name), size = 3, vjust = -1.2) +
      scale_color_viridis_d(option = "plasma", name = "Cell Type") +
      scale_size_identity() +
      scale_edge_width_identity() +
      labs(title = "Cell Communication Network") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "none")
    
    # 组合图形
    combined_plot <- (p1 | p2) / p3
    combined_plot <- combined_plot + 
      plot_annotation(title = paste(cancer_type, "- Network Analysis"),
                     theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
    
    # 保存图形
    plot_file <- file.path(communication_output_dir, "network_analysis",
                          paste0(cancer_type, "_network_analysis.png"))
    ggsave(plot_file, combined_plot, width = 16, height = 12, dpi = 300, bg = "white")
    
    plot_file_pdf <- file.path(communication_output_dir, "network_analysis",
                              paste0(cancer_type, "_network_analysis.pdf"))
    ggsave(plot_file_pdf, combined_plot, width = 16, height = 12, bg = "white")
    
    return(combined_plot)
    
  }, "plot_network_analysis", cancer_type, NULL))
}

# 生成综合报告
generate_communication_report <- function(cancer_type, cellchat, module_results, 
                                        colocalization_results, pathway_results, 
                                        network_stats) {
  message(paste("生成通讯分析报告:", cancer_type))
  
  return(safe_execute({
    report_file <- file.path(communication_output_dir, "qc_reports",
                           paste0(cancer_type, "_communication_report.txt"))
    
    sink(report_file)
    
    cat("=== 细胞通讯分析报告 ===\n")
    cat("癌种:", cancer_type, "\n")
    cat("分析时间:", as.character(Sys.time()), "\n\n")
    
    # 基本统计
    cat("=== 基本统计 ===\n")
    if (!is.null(cellchat)) {
      cat("分析细胞类型数:", length(levels(cellchat@idents)), "\n")
      cat("检测到的配体-受体对数:", nrow(cellchat@LR$LRsig), "\n")
      cat("总通讯相互作用数:", sum(cellchat@net$count > 0, na.rm = TRUE), "\n")
      cat("平均通讯概率:", round(mean(cellchat@net$prob[cellchat@net$prob > 0], na.rm = TRUE), 4), "\n")
    }
    
    # 功能模块分析结果
    cat("\n=== 功能模块分析 ===\n")
    if (length(module_results) > 0) {
      for (module in names(module_results)) {
        stats <- module_results[[module]]
        cat(paste0(module, ":\n"))
        cat(paste0("  - 相互作用数: ", stats$total_interactions, "\n"))
        cat(paste0("  - 平均概率: ", round(stats$avg_prob, 4), "\n"))
        cat(paste0("  - 独特通路数: ", stats$unique_pathways, "\n"))
        cat(paste0("  - 涉及细胞类型: ", paste(stats$cell_types, collapse = ", "), "\n\n"))
      }
    } else {
      cat("未检测到功能模块相互作用\n")
    }
    
    # 空间共定位分析结果
    cat("=== 空间共定位分析 ===\n")
    if (length(colocalization_results) > 0) {
      for (layer in c("High_Colocalization", "Medium_Colocalization", "Low_Colocalization")) {
        if (layer %in% names(colocalization_results) && 
            !is.null(colocalization_results[[layer]]$avg_strength)) {
          cat(paste0(layer, ":\n"))
          cat(paste0("  - 相互作用数: ", colocalization_results[[layer]]$n_interactions, "\n"))
          cat(paste0("  - 平均强度: ", round(colocalization_results[[layer]]$avg_strength, 4), "\n\n"))
        }
      }
    } else {
      cat("未检测到空间共定位模式\n")
    }
    
    # 信号通路分析结果
    cat("=== 信号通路分析 ===\n")
    if (length(pathway_results) > 0) {
      for (category in names(pathway_results)) {
        stats <- pathway_results[[category]]
        cat(paste0(category, ":\n"))
        cat(paste0("  - 通路数: ", stats$n_pathways, "\n"))
        cat(paste0("  - 相互作用数: ", stats$n_interactions, "\n"))
        cat(paste0("  - 平均概率: ", round(stats$avg_prob, 4), "\n\n"))
      }
    } else {
      cat("未检测到特定信号通路\n")
    }
    
    # 网络分析结果
    cat("=== 网络分析 ===\n")
    if (length(network_stats) > 0) {
      cat(paste0("网络密度: ", round(network_stats$network_density, 4), "\n"))
      cat(paste0("总相互作用数: ", network_stats$total_interactions, "\n"))
      
      # Hub细胞类型
      cat("\nHub细胞类型 (基于PageRank):\n")
      if ("page_rank" %in% names(network_stats$hub_analysis)) {
        hubs <- network_stats$hub_analysis$page_rank$hubs
        if (length(hubs) > 0) {
          for (hub in hubs) {
            value <- network_stats$hub_analysis$page_rank$values[hub]
            cat(paste0("  - ", hub, ": ", round(value, 4), "\n"))
          }
        } else {
          cat("  未识别到明显的hub细胞类型\n")
        }
      }
    }
    
    cat("\n=== 报告结束 ===\n")
    sink()
    
    message(paste("  报告已保存:", report_file))
    
  }, "generate_report", cancer_type, NULL))
}

# 保存分析结果表格
save_analysis_tables <- function(cancer_type, module_results, colocalization_results, 
                                pathway_results, network_stats) {
  message(paste("保存分析表格:", cancer_type))
  
  return(safe_execute({
    # 1. 功能模块结果表
    if (length(module_results) > 0) {
      module_table <- data.frame(
        Module = names(module_results),
        Total_Interactions = sapply(module_results, function(x) x$total_interactions),
        Avg_Probability = sapply(module_results, function(x) x$avg_prob),
        Unique_Pathways = sapply(module_results, function(x) x$unique_pathways),
        Cell_Types = sapply(module_results, function(x) paste(x$cell_types, collapse = "; ")),
        stringsAsFactors = FALSE
      )
      
      write.csv(module_table, 
                file.path(communication_output_dir, "tables", 
                         paste0(cancer_type, "_functional_modules.csv")),
                row.names = FALSE)
    }
    
    # 2. 空间共定位结果表
    if (length(colocalization_results) > 0) {
      coloc_table <- data.frame()
      for (layer in names(colocalization_results)) {
        if (layer != "weighted_networks" && !is.null(colocalization_results[[layer]]$interactions)) {
          layer_data <- data.frame(
            Layer = layer,
            Interaction = colocalization_results[[layer]]$interactions,
            Strength = colocalization_results[[layer]]$strength,
            stringsAsFactors = FALSE
          )
          coloc_table <- rbind(coloc_table, layer_data)
        }
      }
      
      if (nrow(coloc_table) > 0) {
        write.csv(coloc_table,
                  file.path(communication_output_dir, "tables",
                           paste0(cancer_type, "_spatial_colocalization.csv")),
                  row.names = FALSE)
      }
    }
    
    # 3. 信号通路结果表
    if (length(pathway_results) > 0) {
      pathway_table <- data.frame(
        Category = names(pathway_results),
        N_Pathways = sapply(pathway_results, function(x) x$n_pathways),
        N_Interactions = sapply(pathway_results, function(x) x$n_interactions),
        Avg_Probability = sapply(pathway_results, function(x) x$avg_prob),
        Top_Pathways = sapply(pathway_results, function(x) paste(x$pathways[1:min(3, length(x$pathways))], collapse = "; ")),
        stringsAsFactors = FALSE
      )
      
      write.csv(pathway_table,
                file.path(communication_output_dir, "tables",
                         paste0(cancer_type, "_pathway_analysis.csv")),
                row.names = FALSE)
    }
    
    # 4. 网络分析结果表
    if (length(network_stats) > 0) {
      network_table <- data.frame(
        Cell_Type = names(network_stats$outgoing_strength),
        Outgoing_Strength = network_stats$outgoing_strength,
        Incoming_Strength = network_stats$incoming_strength,
        Betweenness = network_stats$centrality_metrics$betweenness,
        Closeness = network_stats$centrality_metrics$closeness,
        Eigenvector = network_stats$centrality_metrics$eigenvector,
        PageRank = network_stats$centrality_metrics$page_rank,
        stringsAsFactors = FALSE
      )
      
      write.csv(network_table,
                file.path(communication_output_dir, "tables",
                         paste0(cancer_type, "_network_analysis.csv")),
                row.names = FALSE)
    }
    
  }, "save_tables", cancer_type, NULL))
}

# 主分析函数（单个癌种）
analyze_single_cancer_communication <- function(cancer_type) {
  message(paste("\n", paste(rep("=", 80), collapse = "")))
  message(paste("开始细胞通讯分析:", cancer_type))
  message(paste(rep("=", 80), collapse = ""))
  
  # 1. 数据加载
  message("步骤1: 数据加载")
  seurat_obj <- load_cancer_data(cancer_type)
  
  if (is.null(seurat_obj)) {
    message(paste("数据加载失败，跳过:", cancer_type))
    return(NULL)
  }
  
  # 2. 创建CellChat对象
  message("步骤2: 创建CellChat对象")
  cellchat <- create_cellchat_object(seurat_obj, cancer_type)
  
  if (is.null(cellchat)) {
    message(paste("CellChat对象创建失败，跳过:", cancer_type))
    return(NULL)
  }
  
  # 3. CellChat分析
  message("步骤3: CellChat分析")
  cellchat <- analyze_cellchat(cellchat, cancer_type)
  
  if (is.null(cellchat)) {
    message(paste("CellChat分析失败，跳过:", cancer_type))
    return(NULL)
  }
  
  # 4. 功能模块分析
  message("步骤4: 功能模块分析")
  module_results <- analyze_functional_modules(cellchat, cancer_type)
  
  # 5. 空间共定位分析
  message("步骤5: 空间共定位分析")
  colocalization_results <- analyze_spatial_colocalization(cellchat, cancer_type)
  
  # 6. 信号通路分层分析
  message("步骤6: 信号通路分析")
  pathway_results <- analyze_pathway_layers(cellchat, cancer_type)
  
  # 7. 网络分析
  message("步骤7: 网络分析")
  network_stats <- analyze_communication_network(cellchat, cancer_type)
  
  # 8. 生成可视化
  message("步骤8: 生成可视化")
  plot_functional_modules(module_results, cancer_type)
  plot_spatial_colocalization(colocalization_results, cancer_type)
  plot_network_analysis(network_stats, cancer_type)
  
  # 9. 保存结果表格
  message("步骤9: 保存结果表格")
  save_analysis_tables(cancer_type, module_results, colocalization_results, 
                      pathway_results, network_stats)
  
  # 10. 生成报告
  message("步骤10: 生成报告")
  generate_communication_report(cancer_type, cellchat, module_results, 
                              colocalization_results, pathway_results, network_stats)
  
  # 保存完整结果
  results <- list(
    cancer_type = cancer_type,
    cellchat = cellchat,
    module_results = module_results,
    colocalization_results = colocalization_results,
    pathway_results = pathway_results,
    network_stats = network_stats,
    analysis_time = Sys.time()
  )
  
  results_file <- file.path(communication_output_dir, 
                           paste0(cancer_type, "_communication_results.rds"))
  saveRDS(results, results_file, compress = "xz")
  
  message(paste("分析完成:", cancer_type))
  return(results)
}

# 主分析函数
main_communication_analysis <- function() {
  message("开始细胞通讯分析...")
  message(paste("输出目录:", communication_output_dir))
  
  # 读取分析汇总
  analysis_summary_file <- file.path(output_base_dir, "analysis_summary.csv")
  if (!file.exists(analysis_summary_file)) {
    stop("分析汇总文件不存在:", analysis_summary_file)
  }
  
  analysis_summary <- read.csv(analysis_summary_file, stringsAsFactors = FALSE)
  cancer_types <- analysis_summary$Cancer_Type
  
  message(paste("发现癌种:", length(cancer_types)))
  
  all_results <- list()
  success_count <- 0
  
  # 主分析循环
  for (cancer_type in cancer_types) {
    result <- analyze_single_cancer_communication(cancer_type)
    
    if (!is.null(result)) {
      all_results[[cancer_type]] <- result
      success_count <- success_count + 1
    }
    
    # 强制垃圾回收释放内存
    gc()
  }
  
  message(paste("\n分析完成! 成功:", success_count, "/ 总计:", length(cancer_types)))
  
  # 保存错误汇总
  if (length(error_collector) > 0) {
    error_summary_file <- file.path(communication_output_dir, "error_summary.rds")
    saveRDS(error_collector, error_summary_file)
    message(paste("错误详情保存在:", error_summary_file))
  }
  
  # 保存所有结果
  final_results_file <- file.path(communication_output_dir, "all_communication_results.rds")
  saveRDS(all_results, final_results_file, compress = "xz")
  message(paste("所有结果保存在:", final_results_file))
  
  return(list(
    results = all_results,
    errors = error_collector,
    success_rate = success_count / length(cancer_types)
  ))
}

# 执行分析
if (!interactive()) {
  message("开始执行细胞通讯分析...")
  communication_output <- main_communication_analysis()
  
  message("\n=== 最终统计 ===")
  message(paste("成功分析:", length(communication_output$results)))
  message(paste("错误数量:", length(communication_output$errors)))
  message(paste("成功率:", round(communication_output$success_rate * 100, 1), "%"))
  
} else {
  message("脚本在交互模式下运行")
  message("请运行: communication_output <- main_communication_analysis()")
}

message("\n细胞通讯分析脚本加载完成!")
message("分析框架:")
message("1. 功能模块分层 - 免疫调节、组织结构、微环境支持等")
message("2. 空间共定位分层 - 高/中/低共定位可能性")
message("3. 信号通路分层 - 细胞因子、生长因子、细胞粘附等")
message("4. 网络中心性分析 - Hub细胞识别")
message("5. 多维度可视化和报告生成")