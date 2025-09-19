library(Matrix)
library(SummarizedExperiment) 
library(BiocParallel)
register(MulticoreParam(workers = 16)) 
library(SingleR)
library(Seurat)
library(dplyr)
library(stringr)
library(celldex)
library(org.Hs.eg.db)
library(harmony)

# 设置路径
input_dir <- "/mnt/public7/pancancercol/zhaolingyu/panCNV/cnv_output"
output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation"
sample_info_file <- "/mnt/public7/pancancercol/zhaolingyu/panCNV/inferCNV_ref/sample_info.csv"

# 创建输出目录
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 检查目录创建是否成功
if (!dir.exists(output_dir)) {
  stop("无法创建输出目录: ", output_dir)
}
message("✓ 输出目录创建成功: ", output_dir)

# 读取样本信息
sample_info <- read.csv(sample_info_file, stringsAsFactors = FALSE)
message("读取到 ", nrow(sample_info), " 个样本信息")

# 准备参考数据 - 按照您原代码的方式处理
message("准备参考数据...")
prepare_reference_data <- function() {
  # Monaco参考数据
  ref_monaco <- celldex::MonacoImmuneData()
  gene_symbols_monaco <- rownames(ref_monaco)
  
  ensembl_ids_monaco <- mapIds(
    org.Hs.eg.db,
    keys = gene_symbols_monaco,
    keytype = "SYMBOL",
    column = "ENSEMBL",
    multiVals = "first"
  )
  
  keep_genes_monaco <- !is.na(ensembl_ids_monaco)
  ensembl_ids_monaco <- ensembl_ids_monaco[keep_genes_monaco]
  ref_monaco <- ref_monaco[keep_genes_monaco, ]
  rownames(ref_monaco) <- ensembl_ids_monaco
  
  # Blueprint参考数据
  ref_blueprint <- celldex::BlueprintEncodeData()
  gene_symbols_blueprint <- rownames(ref_blueprint)
  
  ensembl_ids_blueprint <- mapIds(
    org.Hs.eg.db,
    keys = gene_symbols_blueprint,
    keytype = "SYMBOL", 
    column = "ENSEMBL",
    multiVals = "first"
  )
  
  keep_genes_blueprint <- !is.na(ensembl_ids_blueprint)
  ensembl_ids_blueprint <- ensembl_ids_blueprint[keep_genes_blueprint]
  ref_blueprint <- ref_blueprint[keep_genes_blueprint, ]
  rownames(ref_blueprint) <- ensembl_ids_blueprint
  
  message("Monaco参考数据: ", nrow(ref_monaco), " 个基因")
  message("Blueprint参考数据: ", nrow(ref_blueprint), " 个基因")
  
  return(list(monaco = ref_monaco, blueprint = ref_blueprint))
}

ref_data <- prepare_reference_data()

# 定义细胞类型和marker基因 - 转换为ENSEMBL ID
convert_markers_to_ensembl <- function(marker_list) {
  ensembl_markers <- list()
  
  for (cell_type in names(marker_list)) {
    ensembl_markers[[cell_type]] <- list()
    
    for (subtype in names(marker_list[[cell_type]])) {
      symbols <- marker_list[[cell_type]][[subtype]]
      
      ensembl_ids <- mapIds(
        org.Hs.eg.db,
        keys = symbols,
        keytype = "SYMBOL",
        column = "ENSEMBL", 
        multiVals = "first"
      )
      
      # 移除NA值
      ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]
      
      if (length(ensembl_ids) > 0) {
        ensembl_markers[[cell_type]][[subtype]] <- ensembl_ids
      }
    }
  }
  
  return(ensembl_markers)
}

# 定义marker基因（SYMBOL格式）
cell_markers_symbol <- list(
  T_cells = list(
    "CD4_T_cells" = c("CD4", "IL7R", "CCR7"),
    "CD8_T_cells" = c("CD8A", "CD8B", "GZMK"),
    "Regulatory_T_cells" = c("FOXP3", "IL2RA", "CTLA4"),
    "Th1_cells" = c("TBX21", "IFNG", "IL2"),
    "Th2_cells" = c("GATA3", "IL4", "IL13"),
    "Th17_cells" = c("RORC", "IL17A", "IL17F"),
    "Exhausted_T_cells" = c("PDCD1", "TIGIT", "LAG3", "HAVCR2"),
    "Memory_T_cells" = c("IL7R", "CCR7", "SELL", "TCF7"),
    "Effector_T_cells" = c("GZMA", "GZMB", "PRF1", "IFNG"),
    "gamma_delta_T_cells" = c("TRGC1", "TRGC2", "TRDV1"),
    "MAIT_cells" = c("SLC4A10", "NCR3", "KLRB1")
  ),
  
  Myeloid = list(
    "Classical_Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9"),
    "Non_classical_Monocytes" = c("FCGR3A", "MS4A7", "CDKN1C"),
    "Intermediate_Monocytes" = c("CD14", "FCGR3A", "CX3CR1"),
    "M1_Macrophages" = c("CD68", "CD86", "IL1B", "TNF"),
    "M2_Macrophages" = c("CD68", "MRC1", "ARG1", "IL10"),
    "Dendritic_cells" = c("CD1C", "FCER1A", "CLEC10A"),
    "Plasmacytoid_DC" = c("LILRA4", "IRF7", "CLEC4C"),
    "Neutrophils" = c("FCGR3B", "CEACAM8", "CSF3R")
  ),
  
  B_cell = list(
    "Naive_B_cells" = c("MS4A1", "TCL1A", "FCER2"),
    "Memory_B_cells" = c("MS4A1", "CD27", "IGHD"),
    "Plasma_cells" = c("MZB1", "SDC1", "BLIMP1", "XBP1"),
    "Germinal_center_B_cells" = c("BCL6", "AICDA", "MEF2B")
  ),
  
  NK_cell = list(
    "NK_cells" = c("KLRD1", "NCAM1", "NKG7"),
    "CD56bright_NK" = c("NCAM1", "CD160", "KLRC1"),
    "CD56dim_NK" = c("FCGR3A", "CX3CR1", "FGFBP2")
  ),
  
  Epithelial_cells = list(
    "Basal_epithelial" = c("KRT14", "KRT5", "TP63"),
    "Luminal_epithelial" = c("KRT8", "KRT18", "EPCAM"),
    "Secretory_epithelial" = c("MUC1", "SCGB1A1", "SCGB3A1")
  ),
  
  Endothelial_cells = list(
    "Vascular_endothelial" = c("PECAM1", "VWF", "ENG"),
    "Lymphatic_endothelial" = c("PROX1", "LYVE1", "FLT4"),
    "Tip_cells" = c("KDR", "FLT1", "ANGPT2")
  ),
  
  Fibroblasts = list(
    "CAFs" = c("ACTA2", "FAP", "PDGFRB"),
    "Normal_fibroblasts" = c("COL1A1", "COL3A1", "FN1"),
    "Myofibroblasts" = c("ACTA2", "TAGLN", "MYH11")
  )
)

# 转换marker基因到ENSEMBL ID
message("转换marker基因到ENSEMBL ID...")
cell_markers <- convert_markers_to_ensembl(cell_markers_symbol)

# 输出转换结果统计
for (cell_type in names(cell_markers)) {
  message("  ", cell_type, ": ", length(cell_markers[[cell_type]]), " 个亚型")
}

# 定义细胞类型映射关系
cell_type_mapping <- list(
  "T_cells" = "T_cells",
  "Monocyte" = "Myeloid", 
  "Macrophage" = "Myeloid",
  "DC" = "Myeloid",
  "Neutrophils" = "Myeloid",
  "GMP" = "Myeloid",
  "Pro-Myelocyte" = "Myeloid",
  "CMP" = "Myeloid",
  "B_cell" = "B_cell",
  "Pre-B_cell_CD34-" = "B_cell",
  "Pro-B_cell_CD34+" = "B_cell",
  "NK_cell" = "NK_cell",
  "Epithelial_cells" = "Epithelial_cells",
  "Endothelial_cells" = "Endothelial_cells",
  "Fibroblasts" = "Fibroblasts",
  "Smooth_muscle_cells" = "Fibroblasts",
  "Chondrocytes" = "Other",
  "HSC_-G-CSF" = "Other",
  "HSC_CD34+" = "Other",
  "Tissue_stem_cells" = "Other"
)

# 检查基因重叠的函数
check_gene_overlap <- function(test_genes, ref_genes, data_name = "") {
  overlap <- intersect(test_genes, ref_genes)
  overlap_pct <- length(overlap) / length(test_genes) * 100
  
  message(data_name, " 基因重叠情况:")
  message("  测试数据基因数: ", length(test_genes))
  message("  参考数据基因数: ", length(ref_genes))
  message("  重叠基因数: ", length(overlap))
  message("  重叠比例: ", round(overlap_pct, 2), "%")
  
  return(list(overlap = overlap, overlap_pct = overlap_pct))
}

# 计算marker基因得分的函数
calculate_marker_scores <- function(seurat_obj, markers_list) {
  scores_list <- list()
  
  for (cell_type in names(markers_list)) {
    genes <- markers_list[[cell_type]]
    available_genes <- intersect(genes, rownames(seurat_obj))
    
    if (length(available_genes) > 0) {
      message("  计算 ", cell_type, " marker得分，可用基因: ", length(available_genes), "/", length(genes))
      seurat_obj <- AddModuleScore(
        seurat_obj, 
        features = list(available_genes), 
        name = paste0(cell_type, "_score"),
        assay = "RNA"
      )
      scores_list[[cell_type]] <- paste0(cell_type, "_score1")
    } else {
      message("  跳过 ", cell_type, "：无可用marker基因")
    }
  }
  
  return(list(seurat_obj = seurat_obj, score_names = scores_list))
}

# SingleR注释函数 - 使用ENSEMBL ID
perform_singleR_annotation <- function(seurat_obj, ref_data, method = "cluster") {
  # 获取表达矩阵 - 保持ENSEMBL ID格式
  counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  
  results <- list()
  
  # Monaco注释
  tryCatch({
    # 检查基因重叠
    monaco_overlap <- check_gene_overlap(rownames(counts_mat), rownames(ref_data$monaco), "Monaco")
    
    if (monaco_overlap$overlap_pct > 10) {  # 至少10%重叠
      message("开始Monaco注释...")
      
      if (method == "cluster" && "seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
        monaco_results <- SingleR(
          test = counts_mat,
          ref = ref_data$monaco,
          labels = ref_data$monaco$label.fine,
          method = "cluster",
          clusters = seurat_obj$seurat_clusters,
          BPPARAM = bpparam()
        )
        seurat_obj$monaco_fine <- monaco_results$labels[seurat_obj$seurat_clusters]
      } else {
        # 对于大量细胞，使用采样
        if (ncol(counts_mat) > 10000) {
          message("细胞数量较多，使用采样进行注释...")
          sample_cells <- sample(colnames(counts_mat), 5000)
          monaco_results <- SingleR(
            test = counts_mat[, sample_cells],
            ref = ref_data$monaco,
            labels = ref_data$monaco$label.fine,
            method = "single",
            BPPARAM = bpparam()
          )
          # 使用聚类水平注释扩展到所有细胞
          temp_obj <- subset(seurat_obj, cells = sample_cells)
          temp_obj$monaco_labels <- monaco_results$labels
          
          # 基于聚类分配标签到所有细胞
          cluster_labels <- sapply(levels(seurat_obj$seurat_clusters), function(cluster) {
            cluster_cells <- names(seurat_obj$seurat_clusters)[seurat_obj$seurat_clusters == cluster]
            sample_cluster_cells <- intersect(cluster_cells, sample_cells)
            if (length(sample_cluster_cells) > 0) {
              most_common <- names(sort(table(temp_obj$monaco_labels[sample_cluster_cells]), decreasing = TRUE))[1]
              return(most_common)
            } else {
              return("Unknown")
            }
          })
          
          seurat_obj$monaco_fine <- cluster_labels[seurat_obj$seurat_clusters]
        } else {
          monaco_results <- SingleR(
            test = counts_mat,
            ref = ref_data$monaco,
            labels = ref_data$monaco$label.fine,
            method = "single",
            BPPARAM = bpparam()
          )
          seurat_obj$monaco_fine <- monaco_results$labels
        }
      }
      results$monaco <- TRUE
      message("Monaco注释完成")
    } else {
      message("Monaco基因重叠不足，跳过注释")
      seurat_obj$monaco_fine <- "Unknown"
      results$monaco <- FALSE
    }
  }, error = function(e) {
    message("Monaco注释失败: ", e$message)
    seurat_obj$monaco_fine <- "Unknown"
    results$monaco <- FALSE
  })
  
  # Blueprint注释
  tryCatch({
    # 检查基因重叠
    blueprint_overlap <- check_gene_overlap(rownames(counts_mat), rownames(ref_data$blueprint), "Blueprint")
    
    if (blueprint_overlap$overlap_pct > 10) {  # 至少10%重叠
      message("开始Blueprint注释...")
      
      if (method == "cluster" && "seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
        blueprint_results <- SingleR(
          test = counts_mat,
          ref = ref_data$blueprint,
          labels = ref_data$blueprint$label.fine,
          method = "cluster",
          clusters = seurat_obj$seurat_clusters,
          BPPARAM = bpparam()
        )
        seurat_obj$blueprint_fine <- blueprint_results$labels[seurat_obj$seurat_clusters]
      } else {
        # 对于大量细胞，使用采样
        if (ncol(counts_mat) > 10000) {
          message("细胞数量较多，使用采样进行注释...")
          sample_cells <- sample(colnames(counts_mat), 5000)
          blueprint_results <- SingleR(
            test = counts_mat[, sample_cells],
            ref = ref_data$blueprint,
            labels = ref_data$blueprint$label.fine,
            method = "single",
            BPPARAM = bpparam()
          )
          # 使用聚类水平注释扩展到所有细胞
          temp_obj <- subset(seurat_obj, cells = sample_cells)
          temp_obj$blueprint_labels <- blueprint_results$labels
          
          # 基于聚类分配标签到所有细胞
          cluster_labels <- sapply(levels(seurat_obj$seurat_clusters), function(cluster) {
            cluster_cells <- names(seurat_obj$seurat_clusters)[seurat_obj$seurat_clusters == cluster]
            sample_cluster_cells <- intersect(cluster_cells, sample_cells)
            if (length(sample_cluster_cells) > 0) {
              most_common <- names(sort(table(temp_obj$blueprint_labels[sample_cluster_cells]), decreasing = TRUE))[1]
              return(most_common)
            } else {
              return("Unknown")
            }
          })
          
          seurat_obj$blueprint_fine <- cluster_labels[seurat_obj$seurat_clusters]
        } else {
          blueprint_results <- SingleR(
            test = counts_mat,
            ref = ref_data$blueprint,
            labels = ref_data$blueprint$label.fine,
            method = "single",
            BPPARAM = bpparam()
          )
          seurat_obj$blueprint_fine <- blueprint_results$labels
        }
      }
      results$blueprint <- TRUE
      message("Blueprint注释完成")
    } else {
      message("Blueprint基因重叠不足，跳过注释")
      seurat_obj$blueprint_fine <- "Unknown"
      results$blueprint <- FALSE
    }
  }, error = function(e) {
    message("Blueprint注释失败: ", e$message)
    seurat_obj$blueprint_fine <- "Unknown"
    results$blueprint <- FALSE
  })
  
  return(list(seurat_obj = seurat_obj, results = results))
}

# 精细化注释函数（修正版）
perform_fine_annotation <- function(major_cell_type, cell_subset) {  # 修正：移除第一个参数
  message("对 ", major_cell_type, " 进行精细化注释，细胞数: ", ncol(cell_subset))
  
  # 数据预处理
  cell_subset <- NormalizeData(cell_subset)
  cell_subset <- FindVariableFeatures(cell_subset, nfeatures = 2000)
  cell_subset <- ScaleData(cell_subset)
  cell_subset <- RunPCA(cell_subset, npcs = 30, verbose = FALSE)
  
  # 根据细胞数调整聚类参数
  ncells <- ncol(cell_subset)
  if (ncells < 100) {
    resolution <- 0.2
    k_param <- min(10, ncells - 1)
  } else if (ncells < 500) {
    resolution <- 0.3
    k_param <- 20
  } else {
    resolution <- 0.5
    k_param <- 30
  }
  
  cell_subset <- FindNeighbors(cell_subset, dims = 1:20, k.param = k_param)
  cell_subset <- FindClusters(cell_subset, resolution = resolution)
  cell_subset <- RunUMAP(cell_subset, dims = 1:20, verbose = FALSE)
  
  # SingleR注释
  annotation_results <- perform_singleR_annotation(cell_subset, ref_data, method = "cluster")
  cell_subset <- annotation_results$seurat_obj
  
  # 计算marker基因得分
  if (major_cell_type %in% names(cell_markers)) {
    message("计算marker基因得分...")
    marker_results <- calculate_marker_scores(cell_subset, cell_markers[[major_cell_type]])
    cell_subset <- marker_results$seurat_obj
    score_names <- marker_results$score_names
    
    # 基于marker得分进行注释
    cell_subset <- assign_cell_types_by_markers(cell_subset, major_cell_type, score_names)
  }
  
  # 综合注释结果
  cell_subset <- integrate_annotations(cell_subset, major_cell_type)
  
  return(cell_subset)
}

# 基于marker基因得分分配细胞类型
assign_cell_types_by_markers <- function(seurat_obj, major_cell_type, score_names) {
  if (length(score_names) == 0) {
    seurat_obj$marker_based_annotation <- "Unknown"
    return(seurat_obj)
  }
  
  # 获取所有得分列
  score_cols <- unlist(score_names)
  score_data <- seurat_obj@meta.data[, score_cols, drop = FALSE]
  
  # 为每个细胞分配最高得分的细胞类型
  max_scores <- apply(score_data, 1, function(x) {
    max_idx <- which.max(x)
    if (length(max_idx) > 0 && !is.na(max_idx) && x[max_idx] > 0) {
      return(names(score_names)[max_idx])
    } else {
      return("Unassigned")
    }
  })
  
  seurat_obj$marker_based_annotation <- max_scores
  
  return(seurat_obj)
}

# 整合多种注释结果
integrate_annotations <- function(seurat_obj, major_cell_type) {
  meta_data <- seurat_obj@meta.data
  
  # 获取聚类、SingleR和marker注释结果
  cluster_info <- paste0("Cluster_", meta_data$seurat_clusters)
  monaco_info <- if("monaco_fine" %in% colnames(meta_data)) meta_data$monaco_fine else rep("Unknown", nrow(meta_data))
  blueprint_info <- if("blueprint_fine" %in% colnames(meta_data)) meta_data$blueprint_fine else rep("Unknown", nrow(meta_data))
  marker_info <- if("marker_based_annotation" %in% colnames(meta_data)) meta_data$marker_based_annotation else rep("Unknown", nrow(meta_data))
  
  # 综合注释逻辑
  final_annotation <- sapply(1:nrow(meta_data), function(i) {
    cluster <- cluster_info[i]
    monaco <- monaco_info[i]
    blueprint <- blueprint_info[i]
    marker <- marker_info[i]
    
    # 优先使用marker基因结果，然后是SingleR结果
    if (!is.na(marker) && marker != "Unassigned" && marker != "Unknown") {
      return(paste0(major_cell_type, "_", marker))
    } else if (!is.na(monaco) && monaco != "Unknown") {
      return(paste0(major_cell_type, "_", gsub("[^A-Za-z0-9_]", "_", monaco)))
    } else if (!is.na(blueprint) && blueprint != "Unknown") {
      return(paste0(major_cell_type, "_", gsub("[^A-Za-z0-9_]", "_", blueprint)))  
    } else {
      return(paste0(major_cell_type, "_", cluster))
    }
  })
  
  seurat_obj$fine_cell_type <- final_annotation
  
  return(seurat_obj)
}

# 主分析函数（修正版）
process_cancer_type <- function(cancer_type, sample_info_subset) {
  message("\n==== 开始处理癌种: ", cancer_type, " ====")
  
  # 创建癌种输出目录
  cancer_output_dir <- file.path(output_dir, cancer_type)
  dir.create(cancer_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 检查目录是否创建成功
  if (!dir.exists(cancer_output_dir)) {
    stop("无法创建输出目录: ", cancer_output_dir)
  }
  message("✓ 输出目录创建成功: ", cancer_output_dir)
  
  # 收集所有样本数据
  all_seurat_objects <- list()
  
  for (i in 1:nrow(sample_info_subset)) {
    sample_row <- sample_info_subset[i, ]
    sample_id <- sample_row$ID
    file_prefix <- sample_row$file_name
    
    # 构建文件路径
    seurat_file <- file.path(input_dir, cancer_type, "seurat_objects", 
                           paste0(sample_id, "_", file_prefix, "_annotated.rds"))
    
    message("检查文件: ", seurat_file)
    
    if (file.exists(seurat_file)) {
      message("✓ 读取样本: ", sample_id, "_", file_prefix)
      
      tryCatch({
        seurat_obj <- readRDS(seurat_file)
        
        # 添加样本标识
        seurat_obj$sample_id <- sample_id
        seurat_obj$file_prefix <- file_prefix
        seurat_obj$batch <- paste0(sample_id, "_", file_prefix)
        
        all_seurat_objects[[paste0(sample_id, "_", file_prefix)]] <- seurat_obj
        
      }, error = function(e) {
        warning("读取文件失败 ", seurat_file, ": ", e$message)
      })
      
    } else {
      warning("✗ 文件不存在: ", seurat_file)
    }
  }
  
  if (length(all_seurat_objects) == 0) {
    warning("癌种 ", cancer_type, " 没有找到有效样本")
    return(NULL)
  }
  
  message("成功读取 ", length(all_seurat_objects), " 个样本")
  
  # 合并所有样本
  message("合并样本...")
  combined_seurat <- merge(all_seurat_objects[[1]], 
                          y = if(length(all_seurat_objects) > 1) all_seurat_objects[-1] else NULL, 
                          add.cell.ids = names(all_seurat_objects))
  
  message("合并后总细胞数: ", ncol(combined_seurat))
  
  # 添加主要细胞类型分组
  combined_seurat$major_cell_type <- sapply(combined_seurat$singleR_labels, function(x) {
    if (x %in% names(cell_type_mapping)) {
      return(cell_type_mapping[[x]])
    } else {
      return("Other")
    }
  })
  
  # 统计主要细胞类型
  major_type_counts <- table(combined_seurat$major_cell_type)
  message("主要细胞类型分布:")
  for (mt in names(major_type_counts)) {
    message("  ", mt, ": ", major_type_counts[mt], " cells")
  }
  
  # 按主要细胞类型分组进行精细化注释
  major_types <- names(major_type_counts)
  major_types <- major_types[major_types != "Other" & major_type_counts >= 50]  # 至少50个细胞
  
  message("将要处理的主要细胞类型: ", paste(major_types, collapse = ", "))
  
  annotated_objects <- list()
  
  for (major_type in major_types) {
    message("\n--- 处理主要细胞类型: ", major_type, " ---")
    
    # 提取该细胞类型的所有细胞
    cells_of_type <- colnames(combined_seurat)[combined_seurat$major_cell_type == major_type]
    
    message("提取 ", length(cells_of_type), " 个 ", major_type, " 细胞")
    
    # 提取子集
    cell_subset <- subset(combined_seurat, cells = cells_of_type)
    message("子集细胞数验证: ", ncol(cell_subset))
    
    # 批次效应校正（如果有多个样本）
    if (length(unique(cell_subset$batch)) > 1) {
      message("进行批次效应校正...")
      cell_subset <- NormalizeData(cell_subset)
      cell_subset <- FindVariableFeatures(cell_subset, nfeatures = 2000)
      cell_subset <- ScaleData(cell_subset)
      cell_subset <- RunPCA(cell_subset, npcs = 30, verbose = FALSE)
      
      # 使用Harmony进行批次校正
      tryCatch({
        cell_subset <- RunHarmony(cell_subset, group.by.vars = "batch", 
                                 dims.use = 1:20, verbose = FALSE)
        cell_subset <- RunUMAP(cell_subset, reduction = "harmony", 
                              dims = 1:20, verbose = FALSE)
        message("✓ Harmony批次校正完成")
      }, error = function(e) {
        message("⚠️ Harmony失败，使用标准PCA: ", e$message)
        cell_subset <- RunUMAP(cell_subset, dims = 1:20, verbose = FALSE)
      })
    } else {
      message("只有一个批次，跳过批次校正")
    }
    
    # 进行精细化注释 - 修正：传入正确的参数
    tryCatch({
      cell_subset <- perform_fine_annotation(major_type, cell_subset)  # 修正这里
      message("✓ ", major_type, " 精细化注释完成")
      
      # 验证注释结果
      if ("fine_cell_type" %in% colnames(cell_subset@meta.data)) {
        fine_types <- table(cell_subset$fine_cell_type)
        message("  ", major_type, " 精细类型: ", paste(names(fine_types), collapse = ", "))
        annotated_objects[[major_type]] <- cell_subset
      } else {
        warning("⚠️ ", major_type, " 缺少fine_cell_type字段")
      }
    }, error = function(e) {
      warning("⚠️ ", major_type, " 精细化注释失败: ", e$message)
    })
  }
  
  # 合并所有注释结果
  if (length(annotated_objects) > 0) {
    message("\n合并注释结果...")
    
    # 检查annotated_objects
    message("annotated_objects包含: ", paste(names(annotated_objects), collapse = ", "))
    
    # 将精细注释结果添加回原始对象
    fine_annotations <- rep("Other", ncol(combined_seurat))
    names(fine_annotations) <- colnames(combined_seurat)
    
    for (major_type in names(annotated_objects)) {
      obj <- annotated_objects[[major_type]]
      if (!is.null(obj) && "fine_cell_type" %in% colnames(obj@meta.data)) {
        fine_annotations[colnames(obj)] <- obj$fine_cell_type
        message("✓ 已合并 ", major_type, " 的注释结果 (", ncol(obj), " 个细胞)")
      } else {
        warning("⚠️ ", major_type, " 的注释结果为空或缺少fine_cell_type")
      }
    }
    
    combined_seurat$fine_cell_type <- fine_annotations
    
    # 检查最终结果
    message("最终细胞类型统计:")
    final_types <- table(combined_seurat$fine_cell_type)
    for (ft in names(final_types)) {
      message("  ", ft, ": ", final_types[ft])
    }
    
    # 保存结果
    message("保存结果...")
    tryCatch({
      # 保存主要结果
      main_file <- file.path(cancer_output_dir, paste0(cancer_type, "_fine_annotated.rds"))
      saveRDS(combined_seurat, main_file)
      message("✓ 主要结果已保存")
      
      # 检查文件是否真的保存成功
      if (file.exists(main_file)) {
        file_size <- file.size(main_file)
        message("✓ 文件保存确认: ", main_file, " (", round(file_size/1024/1024, 2), " MB)")
      } else {
        warning("⚠️ 文件保存失败: ", main_file)
      }
      
      # 保存细胞类型统计
      cell_type_stats <- table(combined_seurat$fine_cell_type, combined_seurat$sample_id)
      stats_file <- file.path(cancer_output_dir, paste0(cancer_type, "_cell_type_stats.csv"))
      write.csv(cell_type_stats, stats_file)
      message("✓ 细胞类型统计已保存: ", stats_file)
      
      # 保存每个主要细胞类型的详细结果
      for (major_type in names(annotated_objects)) {
        detail_file <- file.path(cancer_output_dir, paste0(cancer_type, "_", major_type, "_detailed.rds"))
        saveRDS(annotated_objects[[major_type]], detail_file)
        if (file.exists(detail_file)) {
          detail_size <- file.size(detail_file)
          message("✓ ", major_type, " 详细结果已保存: ", detail_file, " (", round(detail_size/1024/1024, 2), " MB)")
        } else {
          warning("⚠️ ", major_type, " 详细结果保存失败")
        }
      }
      
    }, error = function(e) {
      warning("保存结果时出错: ", e$message)
      return(NULL)
    })
    
    return(combined_seurat)
  } else {
    warning("没有annotated_objects，跳过保存")
    return(NULL)
  }
}

# 主执行流程
main <- function() {
  message("开始癌种内精细化注释分析...")
  
  # 检查输入数据
  message("检查输入目录结构...")
  if (!dir.exists(input_dir)) {
    stop("输入目录不存在: ", input_dir)
  }
  
  cancer_dirs <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
  message("发现癌种目录: ", paste(cancer_dirs, collapse = ", "))
  
  # 按癌种分组处理
  cancer_types <- unique(sample_info$Cancer_general)
  message("sample_info中的癌种: ", paste(cancer_types, collapse = ", "))
  
  results_summary <- list()
  
  for (cancer_type in cancer_types) {
    tryCatch({
      sample_subset <- sample_info[sample_info$Cancer_general == cancer_type, ]
      
      message("\n", paste(rep("=", 60), collapse = ""))
      message("处理癌种: ", cancer_type, " (", nrow(sample_subset), " 个样本)")
      message(paste(rep("=", 60), collapse = ""))
      
      # 检查癌种目录是否存在
      cancer_dir <- file.path(input_dir, cancer_type)
      if (!dir.exists(cancer_dir)) {
        warning("癌种目录不存在: ", cancer_dir)
        next
      }
      
      seurat_dir <- file.path(cancer_dir, "seurat_objects")
      if (!dir.exists(seurat_dir)) {
        warning("seurat_objects目录不存在: ", seurat_dir)
        next
      }
      
      result <- process_cancer_type(cancer_type, sample_subset)
      
      if (!is.null(result)) {
        results_summary[[cancer_type]] <- list(
          n_samples = length(unique(result$sample_id)),
          n_cells = ncol(result),
          n_cell_types = length(unique(result$fine_cell_type))
        )
        message("✓ 癌种 ", cancer_type, " 处理完成")
      } else {
        message("✗ 癌种 ", cancer_type, " 处理失败")
      }
      
    }, error = function(e) {
      warning("处理癌种 ", cancer_type, " 时出错: ", e$message)
      print(traceback())
    })
  }
  
  # 保存分析汇总
  if (length(results_summary) > 0) {
    summary_df <- do.call(rbind, lapply(names(results_summary), function(ct) {
      data.frame(
        Cancer_Type = ct,
        N_Samples = results_summary[[ct]]$n_samples,
        N_Cells = results_summary[[ct]]$n_cells,
        N_CellTypes = results_summary[[ct]]$n_cell_types
      )
    }))
    
    summary_file <- file.path(output_dir, "analysis_summary.csv")
    write.csv(summary_df, summary_file, row.names = FALSE)
    message("✓ 分析汇总已保存: ", summary_file)
    
    message("\n==== 分析汇总 ====")
    print(summary_df)
  } else {
    message("⚠️ 没有成功处理的癌种")
  }
  
  message("\n==== 所有癌种处理完成 ====")
  message("结果保存在: ", output_dir)
  
  return(if(exists("summary_df")) summary_df else NULL)
}

# 运行主函数
if (!interactive()) {
  summary_result <- main()
  if (!is.null(summary_result)) {
    print(summary_result)
  }
}