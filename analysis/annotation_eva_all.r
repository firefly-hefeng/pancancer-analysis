# ============================================================================
# 高级单细胞注释质量评估系统 - 完整版 v2.4
# Advanced Quality Assessment System for Single-cell Annotation
# Version: 2.4 - Enhanced with Publication-Quality Figures
# ============================================================================

# ============================================================================
# 0. 包管理和环境初始化
# ============================================================================

# 必需包列表
required_packages <- c(
  # 核心分析
  "Seurat", "Matrix", "SingleR", "celldex",
  # 数据处理
  "dplyr", "tidyr", "data.table", "reshape2",
  # 统计分析
  "fossil", "cluster", "mclust", "irr",
  # 功能富集
  "clusterProfiler", "org.Hs.eg.db", "GSVA", "msigdbr",
  # 可视化
  "ggplot2", "ComplexHeatmap", "circlize", "ggrepel", 
  "patchwork", "viridis", "ggsci", "ggridges", "ggalluvial",
  # 性能优化
  "future", "future.apply", "parallel", "RANN",
  # 其他
  "scales", "RColorBrewer", "cowplot", "pheatmap", "gridExtra"
)

# 检查并安装缺失的包
install_if_missing <- function(packages) {
  missing <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing) > 0) {
    message("正在安装缺失的包: ", paste(missing, collapse = ", "))
    
    # 分别处理 CRAN 和 Bioconductor 包
    cran_packages <- c("fossil", "cluster", "mclust", "irr", "ggridges", 
                       "ggalluvial", "ggsci", "RANN", "future", "future.apply",
                       "dplyr", "tidyr", "data.table", "reshape2", "ggplot2",
                       "ggrepel", "patchwork", "viridis", "scales", 
                       "RColorBrewer", "cowplot", "pheatmap", "gridExtra")
    bioc_packages <- c("clusterProfiler", "org.Hs.eg.db", "GSVA", 
                       "msigdbr", "ComplexHeatmap", "celldex", "SingleR",
                       "Seurat", "circlize")
    
    missing_cran <- intersect(missing, cran_packages)
    missing_bioc <- intersect(missing, bioc_packages)
    
    if (length(missing_cran) > 0) {
      install.packages(missing_cran, dependencies = TRUE)
    }
    
    if (length(missing_bioc) > 0) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(missing_bioc, update = FALSE)
    }
  }
  message("✓ 所有必需的包已就绪")
}

# 加载所有包
load_packages <- function(packages) {
  invisible(lapply(packages, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
  message("✓ 所有包加载成功")
}

# 初始化环境
initialize_environment <- function(n_cores = 16, memory_limit = 500000) {
  # 性能优化
  library(future)
  plan(multicore, workers = n_cores)
  options(future.globals.maxSize = memory_limit * 1024^2)
  
  library(data.table)
  setDTthreads(n_cores)
  
  # 设置随机种子
  set.seed(42)
  
  message("✓ 环境初始化完成:")
  message("  CPU核心数: ", n_cores)
  message("  内存限制: ", memory_limit, " MB")
}

# ============================================================================
# 1. 数据加载模块
# ============================================================================

# 自动检测所有癌种
detect_cancer_types <- function(annotation_dir) {
  message("正在检测癌种...")
  
  # 查找所有RDS文件
  rds_files <- list.files(annotation_dir, pattern = "\\.rds$", 
                         full.names = FALSE, ignore.case = TRUE)
  
  if (length(rds_files) == 0) {
    stop("未找到任何RDS文件: ", annotation_dir)
  }
  
  # 提取癌种名称（去除后缀）
  cancer_types <- gsub("_with_umap\\.rds$", "", rds_files, ignore.case = TRUE)
  cancer_types <- gsub("\\.rds$", "", cancer_types, ignore.case = TRUE)
  
  message("✓ 检测到 ", length(cancer_types), " 个癌种:")
  for (ct in cancer_types) {
    message("  - ", ct)
  }
  
  return(cancer_types)
}

# 数据加载函数 - 优化版
load_annotation_data_optimized <- function(cancer_type, 
                                          annotation_dir,
                                          max_cells = 100000,
                                          sample_fraction = 1.0) {
  
  message("正在加载 ", cancer_type, " 的数据...")
  
  # 查找RDS文件
  rds_pattern <- paste0(cancer_type, ".*\\.rds$")
  rds_files <- list.files(annotation_dir, pattern = rds_pattern, 
                         full.names = TRUE, ignore.case = TRUE)
  
  if (length(rds_files) == 0) {
    message("错误: 未找到 ", cancer_type, " 的RDS文件")
    message("查找目录: ", annotation_dir)
    message("匹配模式: ", rds_pattern)
    return(NULL)
  }
  
  # 如果有多个文件，选择第一个
  if (length(rds_files) > 1) {
    message("找到多个RDS文件，使用: ", basename(rds_files[1]))
  }
  
  rds_file <- rds_files[1]
  message("加载文件: ", basename(rds_file))
  
  # 加载数据
  seurat_obj <- tryCatch({
    readRDS(rds_file)
  }, error = function(e) {
    message("加载RDS文件时出错: ", e$message)
    return(NULL)
  })
  
  if (is.null(seurat_obj)) {
    return(NULL)
  }
  
  # 检查对象类型
  if (!inherits(seurat_obj, "Seurat")) {
    message("错误: 对象不是Seurat对象")
    return(NULL)
  }
  
  original_n_cells <- ncol(seurat_obj)
  message("原始细胞数: ", original_n_cells)
  
  # 检查必要的列
  required_cols <- c("fine_cell_type")
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  
  if (length(missing_cols) > 0) {
    message("警告: 缺少必需列: ", paste(missing_cols, collapse = ", "))
    message("可用列: ", paste(colnames(seurat_obj@meta.data), collapse = ", "))
  }
  
  # 采样策略
  is_sampled <- FALSE
  
  # 1. 如果细胞数超过max_cells，进行采样
  if (original_n_cells > max_cells) {
    message("从 ", original_n_cells, " 个细胞中采样 ", max_cells, " 个")
    
    # 分层采样（如果有细胞类型注释）
    if ("fine_cell_type" %in% colnames(seurat_obj@meta.data)) {
      
      set.seed(42)
      sample_cells <- c()
      
      for (ct in unique(seurat_obj$fine_cell_type)) {
        ct_cells <- colnames(seurat_obj)[seurat_obj$fine_cell_type == ct]
        n_sample <- min(length(ct_cells), 
                       ceiling(max_cells * length(ct_cells) / original_n_cells))
        n_sample <- max(10, n_sample)  # 每个类型至少10个细胞
        
        if (length(ct_cells) > n_sample) {
          sample_cells <- c(sample_cells, sample(ct_cells, n_sample))
        } else {
          sample_cells <- c(sample_cells, ct_cells)
        }
      }
      
      # 确保不超过max_cells
      if (length(sample_cells) > max_cells) {
        set.seed(42)
        sample_cells <- sample(sample_cells, max_cells)
      }
      
      seurat_obj <- subset(seurat_obj, cells = sample_cells)
      
    } else {
      # 随机采样
      set.seed(42)
      sample_cells <- sample(colnames(seurat_obj), max_cells)
      seurat_obj <- subset(seurat_obj, cells = sample_cells)
    }
    
    is_sampled <- TRUE
    message("采样后细胞数: ", ncol(seurat_obj))
  }
  
  # 2. 按比例采样
  if (sample_fraction < 1.0 && !is_sampled) {
    n_sample <- floor(ncol(seurat_obj) * sample_fraction)
    message("按 ", sample_fraction * 100, "% 采样 (", n_sample, " 个细胞)")
    
    set.seed(42)
    sample_cells <- sample(colnames(seurat_obj), n_sample)
    seurat_obj <- subset(seurat_obj, cells = sample_cells)
    
    is_sampled <- TRUE
  }
  
  # 确保必要的降维已完成
  if (!"umap" %in% names(seurat_obj@reductions)) {
    message("警告: 未找到UMAP，将跳过空间分析")
  }
  
  if (!"pca" %in% names(seurat_obj@reductions)) {
    message("警告: 未找到PCA")
  }
  
  # 打印数据摘要
  message("✓ 数据加载成功:")
  message("  最终细胞数: ", ncol(seurat_obj))
  message("  基因数: ", nrow(seurat_obj))
  
  if ("fine_cell_type" %in% colnames(seurat_obj@meta.data)) {
    message("  细胞类型数: ", length(unique(seurat_obj$fine_cell_type)))
  }
  
  return(list(
    main = seurat_obj,
    original_n_cells = original_n_cells,
    is_sampled = is_sampled
  ))
}

# ============================================================================
# 2. 辅助函数
# ============================================================================

# 安全采样函数
safe_sample <- function(x, size, replace = FALSE, prob = NULL) {
  if (length(x) == 0) {
    return(character(0))
  }
  
  actual_size <- min(size, length(x))
  
  if (actual_size == 0) {
    return(character(0))
  }
  
  if (actual_size == length(x)) {
    return(x)
  }
  
  return(sample(x, actual_size, replace = replace, prob = prob))
}

# 获取全面的marker基因列表
get_comprehensive_markers <- function() {
  list(
    T_cells = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "IL7R", "TRAC", "TRBC1", "TRBC2"),
    CD4_T = c("CD4", "IL7R", "CCR7", "SELL", "TCF7", "LEF1"),
    CD8_T = c("CD8A", "CD8B", "GZMK", "GZMH", "CCL5", "NKG7"),
    Treg = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "TNFRSF18"),
    NK_cell = c("NCAM1", "NKG7", "KLRD1", "GNLY", "GZMB", "PRF1", "KLRB1", "KLRF1"),
    B_cell = c("CD19", "MS4A1", "CD79A", "CD79B", "IGHM", "IGHD", "IGKC"),
    Plasma = c("JCHAIN", "MZB1", "DERL3", "SDC1", "TNFRSF17", "XBP1"),
    Myeloid = c("CD14", "CD68", "LYZ", "FCGR3A", "CST3", "S100A8", "S100A9"),
    Macrophage = c("CD68", "CD163", "MSR1", "MRC1", "C1QA", "C1QB", "C1QC"),
    Monocyte = c("CD14", "VCAN", "S100A8", "S100A9", "FCN1", "CD16"),
    DC = c("FCER1A", "CD1C", "CLEC9A", "CLEC10A", "IRF8", "BATF3"),
    pDC = c("LILRA4", "IL3RA", "CLEC4C", "IRF7", "TCF4"),
    Neutrophil = c("FCGR3B", "CXCR2", "CSF3R", "FPR1", "CXCR1"),
    Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"),
    Epithelial_cells = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"),
    Endothelial_cells = c("PECAM1", "VWF", "CDH5", "ENG", "CD34", "CLDN5"),
    Fibroblasts = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "PDGFRB"),
    Cancer = c("MKI67", "TOP2A", "PCNA", "CDK1", "AURKA"),
    Cycling = c("MKI67", "TOP2A", "PCNA", "STMN1", "HMGB2", "TUBB")
  )
}

# 提取主要细胞类型
extract_major_type <- function(cell_type_label) {
  
  cell_type_label <- tolower(as.character(cell_type_label))
  
  # 精确匹配规则
  if (grepl("^t[- _]?cell|^cd[48][+ ]?t|tcell", cell_type_label)) {
    if (grepl("cd4|helper|th[0-9]", cell_type_label)) return("CD4_T")
    if (grepl("cd8|cytotoxic|ctl", cell_type_label)) return("CD8_T")
    if (grepl("reg|treg", cell_type_label)) return("Treg")
    return("T_cells")
  }
  
  if (grepl("^nk[- _]?cell|natural.*killer", cell_type_label)) return("NK_cell")
  
  if (grepl("^b[- _]?cell|bcell", cell_type_label)) return("B_cell")
  
  if (grepl("plasma", cell_type_label)) return("Plasma")
  
  if (grepl("macro", cell_type_label)) return("Macrophage")
  
  if (grepl("mono", cell_type_label)) return("Monocyte")
  
  if (grepl("dendritic|^dc[^a-z]|pdc", cell_type_label)) {
    if (grepl("pdc|plasmacytoid", cell_type_label)) return("pDC")
    return("DC")
  }
  
  if (grepl("neutro|granulocyte", cell_type_label)) return("Neutrophil")
  
  if (grepl("mast", cell_type_label)) return("Mast")
  
  if (grepl("myeloid|monocyte|macrophage", cell_type_label)) return("Myeloid")
  
  if (grepl("epithe|tumor|cancer|malig|carcinoma", cell_type_label)) return("Epithelial_cells")
  
  if (grepl("endothe|^ec[^a-z]", cell_type_label)) return("Endothelial_cells")
  
  if (grepl("fibro|caf|mesench", cell_type_label)) return("Fibroblasts")
  
  if (grepl("prolif|cycling|^s[- ]?phase|^g2", cell_type_label)) return("Cycling")
  
  return("Other")
}

# 计算归一化互信息(NMI)
calculate_nmi <- function(labels1, labels2) {
  
  if (length(labels1) != length(labels2)) {
    stop("标签长度必须相同")
  }
  
  dt <- data.table(l1 = as.character(labels1), l2 = as.character(labels2))
  
  joint_counts <- dt[, .N, by = .(l1, l2)]
  total <- nrow(dt)
  
  margin1 <- dt[, .N, by = l1]
  margin2 <- dt[, .N, by = l2]
  
  setkey(joint_counts, l1, l2)
  setkey(margin1, l1)
  setkey(margin2, l2)
  
  joint_counts <- merge(joint_counts, margin1, by = "l1")
  setnames(joint_counts, c("l1", "l2", "N_joint", "N_l1"))
  joint_counts <- merge(joint_counts, margin2, by = "l2")
  setnames(joint_counts, c("l2", "l1", "N_joint", "N_l1", "N_l2"))
  
  joint_counts[, mi_term := (N_joint/total) * log2((N_joint * total) / (N_l1 * N_l2))]
  mi <- sum(joint_counts$mi_term, na.rm = TRUE)
  
  margin1[, entropy_term := -(N/total) * log2(N/total)]
  margin2[, entropy_term := -(N/total) * log2(N/total)]
  
  entropy1 <- sum(margin1$entropy_term)
  entropy2 <- sum(margin2$entropy_term)
  
  if (entropy1 + entropy2 == 0) return(0)
  
  nmi <- 2 * mi / (entropy1 + entropy2)
  
  return(nmi)
}

# 计算基尼系数
calculate_gini <- function(x) {
  x <- sort(x[x >= 0])
  n <- length(x)
  if (n == 0 || sum(x) == 0) return(0)
  G <- 2 * sum((1:n) * x) / (n * sum(x)) - (n + 1) / n
  return(G)
}

# 计算空间统计量
calculate_spatial_statistics <- function(coords, labels, k = 50) {
  
  if (!requireNamespace("RANN", quietly = TRUE)) {
    message("RANN包不可用，跳过空间统计")
    return(list(
      within_type_dist = NA,
      between_type_dist = NA,
      separation_index = NA,
      morans_i = NA
    ))
  }
  
  library(RANN)
  
  k <- min(k, nrow(coords) - 1)
  
  if (k < 2) {
    return(list(
      within_type_dist = NA,
      between_type_dist = NA,
      separation_index = NA,
      morans_i = NA
    ))
  }
  
  nn_result <- nn2(coords, coords, k = k + 1)
  
  # 类内vs类间距离
  within_dists <- c()
  between_dists <- c()
  
  for (i in 1:min(1000, nrow(coords))) {
    neighbors <- nn_result$nn.idx[i, -1]
    neighbor_labels <- labels[neighbors]
    neighbor_dists <- nn_result$nn.dists[i, -1]
    
    same_type <- neighbor_labels == labels[i]
    if (sum(same_type) > 0) {
      within_dists <- c(within_dists, mean(neighbor_dists[same_type]))
    }
    if (sum(!same_type) > 0) {
      between_dists <- c(between_dists, mean(neighbor_dists[!same_type]))
    }
  }
  
  # Moran's I (简化版)
  label_numeric <- as.numeric(factor(labels))
  
  morans_i <- tryCatch({
    mean_label <- mean(label_numeric)
    
    numerator <- 0
    denominator <- sum((label_numeric - mean_label)^2)
    
    n_sample <- min(1000, nrow(coords))
    
    for (i in sample(1:nrow(coords), n_sample)) {
      neighbors <- nn_result$nn.idx[i, 2:min(11, k+1)]
      w <- 1 / (nn_result$nn.dists[i, 2:min(11, k+1)] + 1e-6)
      w <- w / sum(w)
      
      numerator <- numerator + sum(w * (label_numeric[i] - mean_label) * 
                                    (label_numeric[neighbors] - mean_label))
    }
    
    morans_i <- (n_sample / denominator) * numerator
    morans_i
  }, error = function(e) NA)
  
  return(list(
    within_type_dist = mean(within_dists, na.rm = TRUE),
    between_type_dist = mean(between_dists, na.rm = TRUE),
    separation_index = mean(between_dists, na.rm = TRUE) / mean(within_dists, na.rm = TRUE),
    morans_i = morans_i
  ))
}

# 计算批次混合指数
calculate_mixing_metric <- function(coords, batch_labels, k = 50) {
  
  if (!requireNamespace("RANN", quietly = TRUE)) {
    return(list(
      mixing_score = NA,
      batch_entropy = NA,
      low_mixing_pct = NA
    ))
  }
  
  library(RANN)
  
  k <- min(k, nrow(coords) - 1)
  
  if (k < 2) {
    return(list(
      mixing_score = NA,
      batch_entropy = NA,
      low_mixing_pct = NA
    ))
  }
  
  nn_result <- nn2(coords, coords, k = k + 1)
  
  # 计算每个细胞邻域的批次熵
  batch_entropy <- sapply(1:nrow(coords), function(i) {
    neighbors <- nn_result$nn.idx[i, -1]
    neighbor_batches <- batch_labels[neighbors]
    batch_props <- prop.table(table(neighbor_batches))
    -sum(batch_props * log2(batch_props + 1e-10))
  })
  
  # 标准化熵
  n_batches <- length(unique(batch_labels))
  max_entropy <- log2(n_batches)
  normalized_entropy <- batch_entropy / max_entropy
  
  mixing_score <- mean(normalized_entropy)
  
  return(list(
    mixing_score = mixing_score,
    batch_entropy = batch_entropy,
    low_mixing_pct = sum(normalized_entropy < 0.5) / length(normalized_entropy)
  ))
}

# ============================================================================
# 3. 评估模块（修复版）
# ============================================================================

# 模块 1: 注释一致性评估
assess_annotation_consistency <- function(seurat_obj) {
  message("\n=== 模块 1: 注释一致性评估 ===")
  
  results <- list()
  
  # 1.1 多方法一致性矩阵
  annotation_fields <- c("monaco_fine", "blueprint_fine", 
                         "marker_based_annotation", "fine_cell_type")
  available_fields <- annotation_fields[annotation_fields %in% 
                                        colnames(seurat_obj@meta.data)]
  
  if (length(available_fields) >= 2) {
    
    consistency_matrix <- matrix(NA, 
                                 nrow = length(available_fields),
                                 ncol = length(available_fields),
                                 dimnames = list(available_fields, available_fields))
    
    ari_matrix <- consistency_matrix
    nmi_matrix <- consistency_matrix
    kappa_matrix <- consistency_matrix
    
    for (i in 1:length(available_fields)) {
      for (j in 1:length(available_fields)) {
        if (i != j) {
          labels1 <- as.character(seurat_obj@meta.data[[available_fields[i]]])
          labels2 <- as.character(seurat_obj@meta.data[[available_fields[j]]])
          
          valid_idx <- !grepl("Unknown|Unassigned|Other", labels1) & 
                       !grepl("Unknown|Unassigned|Other", labels2)
          
          if (sum(valid_idx) > 100) {
            labels1 <- labels1[valid_idx]
            labels2 <- labels2[valid_idx]
            
            consistency_matrix[i, j] <- sum(labels1 == labels2) / length(labels1)
            
            ari_matrix[i, j] <- tryCatch({
              fossil::adj.rand.index(labels1, labels2)
            }, error = function(e) NA)
            
            nmi_matrix[i, j] <- calculate_nmi(labels1, labels2)
            
            kappa_matrix[i, j] <- tryCatch({
              irr::kappa2(data.frame(labels1, labels2))$value
            }, error = function(e) NA)
          }
        } else {
          consistency_matrix[i, j] <- 1
          ari_matrix[i, j] <- 1
          nmi_matrix[i, j] <- 1
          kappa_matrix[i, j] <- 1
        }
      }
    }
    
    results$consistency_matrices <- list(
      raw_consistency = consistency_matrix,
      ari = ari_matrix,
      nmi = nmi_matrix,
      kappa = kappa_matrix
    )
    
    results$overall_consistency <- mean(ari_matrix[upper.tri(ari_matrix)], na.rm = TRUE)
    
    message("  总体一致性 (ARI): ", round(results$overall_consistency, 3))
  }
  
  # 1.2 层级一致性检查
  if (all(c("major_cell_type", "fine_cell_type") %in% colnames(seurat_obj@meta.data))) {
    
    hierarchy_table <- table(
      Major = seurat_obj$major_cell_type,
      Fine = seurat_obj$fine_cell_type
    )
    
    hierarchy_purity <- apply(hierarchy_table, 2, function(x) max(x) / sum(x))
    
    results$hierarchy_consistency <- list(
      mean_purity = mean(hierarchy_purity),
      median_purity = median(hierarchy_purity),
      hierarchy_table = hierarchy_table,
      purity_per_type = hierarchy_purity
    )
    
    message("  层级一致性: ", round(mean(hierarchy_purity), 3))
  }
  
  # 1.3 细胞类型映射质量
  if ("fine_cell_type" %in% colnames(seurat_obj@meta.data)) {
    
    type_counts <- table(seurat_obj$fine_cell_type)
    
    singleton_types <- names(type_counts)[type_counts == 1]
    tiny_types <- names(type_counts)[type_counts < 10]
    unknown_types <- grep("Unknown|Unassigned|Other", names(type_counts), value = TRUE)
    
    results$mapping_quality <- list(
      total_types = length(type_counts),
      singleton_types = length(singleton_types),
      tiny_types = length(tiny_types),
      unknown_types = length(unknown_types),
      unknown_proportion = sum(type_counts[unknown_types]) / sum(type_counts)
    )
    
    message("  映射质量:")
    message("    细胞类型总数: ", length(type_counts))
    message("    未知比例: ", 
            round(results$mapping_quality$unknown_proportion * 100, 2), "%")
  }
  
  return(results)
}

# 模块 2: 生物学合理性评估（增强版）
assess_biological_validity <- function(seurat_obj, 
                                       max_cells = 20000,
                                       n_top_markers = 10) {
  message("\n=== 模块 2: 生物学合理性评估 ===")
  
  results <- list()
  
  if (ncol(seurat_obj) > max_cells) {
    set.seed(42)
    sample_cells <- safe_sample(colnames(seurat_obj), max_cells)
    calc_obj <- subset(seurat_obj, cells = sample_cells)
  } else {
    calc_obj <- seurat_obj
  }
  
  # 2.1 Marker基因表达验证
  message("  验证marker基因表达...")
  
  canonical_markers <- get_comprehensive_markers()
  
  if ("fine_cell_type" %in% colnames(calc_obj@meta.data)) {
    
    cell_types <- unique(calc_obj$fine_cell_type)
    cell_types <- cell_types[!grepl("Unknown|Unassigned", cell_types)]
    
    marker_validation <- list()
    
    for (ct in cell_types) {
      
      major_type <- extract_major_type(ct)
      
      if (major_type %in% names(canonical_markers)) {
        
        ct_cells <- colnames(calc_obj)[calc_obj$fine_cell_type == ct]
        other_cells <- setdiff(colnames(calc_obj), ct_cells)
        
        ct_markers <- canonical_markers[[major_type]]
        available_markers <- intersect(ct_markers, rownames(calc_obj))
        
        if (length(available_markers) >= 2 && length(ct_cells) >= 10) {
          
          expr_mat <- GetAssayData(calc_obj, assay = "RNA", layer = "data")
          
          ct_expr <- Matrix::rowMeans(expr_mat[available_markers, ct_cells, drop = FALSE])
          bg_expr <- Matrix::rowMeans(expr_mat[available_markers, other_cells, drop = FALSE])
          
          log2fc <- log2((ct_expr + 0.1) / (bg_expr + 0.1))
          pct_in <- Matrix::rowSums(expr_mat[available_markers, ct_cells, drop = FALSE] > 0) / 
                    length(ct_cells)
          pct_out <- Matrix::rowSums(expr_mat[available_markers, other_cells, drop = FALSE] > 0) / 
                     length(other_cells)
          
          p_values <- sapply(available_markers, function(gene) {
            tryCatch({
              wilcox.test(
                expr_mat[gene, ct_cells],
                expr_mat[gene, other_cells],
                alternative = "greater"
              )$p.value
            }, error = function(e) 1)
          })
          
          marker_validation[[ct]] <- list(
            cell_type = ct,
            n_cells = length(ct_cells),
            markers_tested = available_markers,
            mean_log2fc = mean(log2fc),
            median_log2fc = median(log2fc),
            n_significant = sum(p_values < 0.05),
            mean_pct_in = mean(pct_in),
            mean_pct_out = mean(pct_out),
            marker_scores = data.frame(
              gene = available_markers,
              log2fc = log2fc,
              pct_in = pct_in,
              pct_out = pct_out,
              p_value = p_values,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
    
    results$marker_validation <- marker_validation
    
    if (length(marker_validation) > 0) {
      avg_log2fc <- mean(sapply(marker_validation, function(x) x$mean_log2fc))
      avg_pct_diff <- mean(sapply(marker_validation, function(x) 
        x$mean_pct_in - x$mean_pct_out))
      
      results$marker_validation_summary <- list(
        n_types_validated = length(marker_validation),
        avg_log2fc = avg_log2fc,
        avg_pct_difference = avg_pct_diff
      )
      
      message("    验证了 ", length(marker_validation), " 个细胞类型")
      message("    平均 Log2FC: ", round(avg_log2fc, 3))
    }
  }
  
  # 2.2 聚类-注释匹配度
  message("  评估聚类-注释对齐度...")
  
  if (all(c("seurat_clusters", "fine_cell_type") %in% colnames(calc_obj@meta.data))) {
    
    dt <- data.table(
      cluster = as.character(calc_obj$seurat_clusters),
      celltype = as.character(calc_obj$fine_cell_type)
    )
    
    cluster_purity <- dt[, .(
      purity = max(.N) / .N[1],
      entropy = -sum((.N / .N[1]) * log2(.N / .N[1] + 1e-10)),
      dominant_type = celltype[which.max(.N)],
      n_types = uniqueN(celltype),
      n_cells = .N[1]
    ), by = .(cluster, total = .N)]
    
    celltype_spread <- dt[, .(
      n_clusters = uniqueN(cluster),
      entropy = -sum((.N / .N[1]) * log2(.N / .N[1] + 1e-10)),
      dominant_cluster = cluster[which.max(.N)]
    ), by = celltype]
    
    # 轮廓系数
    if ("umap" %in% names(calc_obj@reductions)) {
      
      umap_coords <- Embeddings(calc_obj, "umap")
      
      if (nrow(umap_coords) > 5000) {
        sample_idx <- safe_sample(1:nrow(umap_coords), 5000)
        umap_coords <- umap_coords[sample_idx, ]
        celltype_labels <- calc_obj$fine_cell_type[sample_idx]
      } else {
        celltype_labels <- calc_obj$fine_cell_type
      }
      
      sil <- tryCatch({
        celltype_numeric <- as.numeric(factor(celltype_labels))
        cluster::silhouette(celltype_numeric, dist(umap_coords))
      }, error = function(e) NULL)
      
      if (!is.null(sil)) {
        sil_scores <- summary(sil)$clus.avg.widths
        
        results$silhouette_analysis <- list(
          mean_silhouette = mean(sil[, "sil_width"]),
          per_type_silhouette = sil_scores
        )
        
        message("    平均轮廓系数: ", round(mean(sil[, "sil_width"]), 3))
      }
    }
    
    results$cluster_annotation_alignment <- list(
      cluster_purity = cluster_purity,
      celltype_spread = celltype_spread,
      mean_cluster_purity = mean(cluster_purity$purity),
      mean_cluster_entropy = mean(cluster_purity$entropy),
      mean_celltype_clusters = mean(celltype_spread$n_clusters)
    )
    
    message("    平均聚类纯度: ", round(mean(cluster_purity$purity), 3))
  }
  
  # 2.3 空间连续性分析
  message("  分析空间连贯性...")
  
  if ("umap" %in% names(calc_obj@reductions) && 
      "fine_cell_type" %in% colnames(calc_obj@meta.data)) {
    
    umap_coords <- Embeddings(calc_obj, "umap")
    cell_types <- calc_obj$fine_cell_type
    
    if (nrow(umap_coords) > 5000) {
      set.seed(42)
      sample_idx <- safe_sample(1:nrow(umap_coords), 5000)
      umap_coords <- umap_coords[sample_idx, ]
      cell_types <- cell_types[sample_idx]
    }
    
    spatial_stats <- calculate_spatial_statistics(umap_coords, cell_types)
    
    results$spatial_coherence <- spatial_stats
    
    if (!is.na(spatial_stats$separation_index)) {
      message("    分离指数: ", round(spatial_stats$separation_index, 3))
    }
  }
  
  # 2.4 细胞类型比例分析
  message("  分析细胞类型比例...")
  
  if ("fine_cell_type" %in% colnames(calc_obj@meta.data)) {
    
    type_counts <- table(calc_obj$fine_cell_type)
    type_props <- prop.table(type_counts)
    
    gini <- calculate_gini(as.numeric(type_props))
    
    shannon <- -sum(type_props * log2(type_props + 1e-10))
    max_shannon <- log2(length(type_props))
    normalized_shannon <- shannon / max_shannon
    
    simpson <- 1 - sum(type_props^2)
    
    results$proportion_analysis <- list(
      n_types = length(type_props),
      gini_coefficient = gini,
      shannon_diversity = shannon,
      normalized_shannon = normalized_shannon,
      simpson_diversity = simpson,
      type_proportions = type_props
    )
    
    message("    细胞类型多样性 (Shannon): ", round(normalized_shannon, 3))
    message("    基尼系数: ", round(gini, 3))
  }
  
  rm(calc_obj)
  gc()
  
  return(results)
}

# 模块 3: 统计显著性评估（修复GSVA）
assess_statistical_significance <- function(seurat_obj,
                                           test_method = "wilcox",
                                           min_pct = 0.1,
                                           logfc_threshold = 0.25,
                                           max_cells_per_type = 2000) {
  message("\n=== 模块 3: 统计显著性评估 ===")
  
  results <- list()
  
  if (!"fine_cell_type" %in% colnames(seurat_obj@meta.data)) {
    message("  未找到细胞类型注释")
    return(NULL)
  }
  
  message("  查找差异基因...")
  
  cell_types <- unique(seurat_obj$fine_cell_type)
  cell_types <- cell_types[!grepl("Unknown|Unassigned", cell_types)]
  
  type_counts <- table(seurat_obj$fine_cell_type)
  cell_types <- cell_types[cell_types %in% names(type_counts)[type_counts >= 30]]
  
  if (length(cell_types) == 0) {
    message("  没有有效的细胞类型用于DEG分析")
    return(NULL)
  }
  
  Idents(seurat_obj) <- "fine_cell_type"
  
  deg_results <- list()
  
  for (ct in cell_types[1:min(10, length(cell_types))]) {
    
    message("    测试 ", ct, "...")
    
    ct_cells <- colnames(seurat_obj)[seurat_obj$fine_cell_type == ct]
    
    if (length(ct_cells) > max_cells_per_type) {
      set.seed(42)
      ct_cells <- safe_sample(ct_cells, max_cells_per_type)
    }
    
    markers <- tryCatch({
      FindMarkers(
        seurat_obj,
        ident.1 = ct,
        min.pct = min_pct,
        logfc.threshold = logfc_threshold,
        test.use = test_method,
        verbose = FALSE
      )
    }, error = function(e) {
      message("      失败: ", e$message)
      return(NULL)
    })
    
    if (!is.null(markers) && nrow(markers) > 0) {
      
      markers$gene <- rownames(markers)
      markers <- markers[order(markers$p_val_adj, -markers$avg_log2FC), ]
      
      deg_results[[ct]] <- list(
        cell_type = ct,
        n_significant = sum(markers$p_val_adj < 0.05),
        n_upregulated = sum(markers$avg_log2FC > 0 & markers$p_val_adj < 0.05),
        n_downregulated = sum(markers$avg_log2FC < 0 & markers$p_val_adj < 0.05),
        top_genes = head(markers, 50),
        all_markers = markers
      )
      
      message("      发现 ", sum(markers$p_val_adj < 0.05), " 个显著基因")
    }
  }
  
  results$deg_analysis <- deg_results
  
  # GSVA分析（修复版）
  message("  计算基因集得分...")
  
  if (requireNamespace("GSVA", quietly = TRUE) && 
      requireNamespace("msigdbr", quietly = TRUE)) {
    
    hallmark_sets <- tryCatch({
      msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    }, error = function(e) {
      message("    无法加载Hallmark基因集，尝试使用预定义基因集")
      
      # 使用预定义的关键通路基因集
      predefined_sets <- list(
        "HALLMARK_INFLAMMATORY_RESPONSE" = c("IL1B", "IL6", "TNF", "CXCL8", "CCL2", "PTGS2"),
        "HALLMARK_INTERFERON_GAMMA_RESPONSE" = c("IRF1", "STAT1", "GBP1", "CXCL9", "CXCL10", "IDO1"),
        "HALLMARK_HYPOXIA" = c("VEGFA", "HIF1A", "LDHA", "PDK1", "SLC2A1", "ENO1"),
        "HALLMARK_APOPTOSIS" = c("BAX", "BCL2", "CASP3", "CASP8", "FAS", "FASLG"),
        "HALLMARK_P53_PATHWAY" = c("TP53", "MDM2", "CDKN1A", "BAX", "PUMA", "GADD45A"),
        "HALLMARK_MYC_TARGETS" = c("MYC", "NCL", "NPM1", "CDK4", "LDHA", "PKM"),
        "HALLMARK_E2F_TARGETS" = c("E2F1", "PCNA", "MCM2", "MCM6", "CDC6", "CDK1"),
        "HALLMARK_G2M_CHECKPOINT" = c("CDK1", "CCNB1", "CCNB2", "PLK1", "AURKA", "AURKB"),
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = c("VIM", "CDH2", "SNAI1", "TWIST1", "ZEB1", "FN1"),
        "HALLMARK_ANGIOGENESIS" = c("VEGFA", "KDR", "FLT1", "ANGPT2", "TEK", "PECAM1")
      )
      
      gene_set_df <- data.frame(
        gs_name = rep(names(predefined_sets), sapply(predefined_sets, length)),
        gene_symbol = unlist(predefined_sets),
        stringsAsFactors = FALSE
      )
      
      return(gene_set_df)
    })
    
    if (!is.null(hallmark_sets) && nrow(hallmark_sets) > 0) {
      
      hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
      
      if (ncol(seurat_obj) > 5000) {
        set.seed(42)
        sample_cells <- safe_sample(colnames(seurat_obj), 5000)
        expr_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[, sample_cells]
        cell_types_sampled <- seurat_obj$fine_cell_type[sample_cells]
      } else {
        expr_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
        cell_types_sampled <- seurat_obj$fine_cell_type
      }
      
      # 过滤基因集，只保留存在的基因
      hallmark_list <- lapply(hallmark_list, function(genes) {
        intersect(genes, rownames(expr_mat))
      })
      hallmark_list <- hallmark_list[sapply(hallmark_list, length) >= 5]
      
      if (length(hallmark_list) > 0) {
        gsva_scores <- tryCatch({
          GSVA::gsva(
            as.matrix(expr_mat),
            hallmark_list,
            method = "gsva",
            verbose = FALSE,
            parallel.sz = 1
          )
        }, error = function(e) {
          message("    GSVA计算失败: ", e$message)
          return(NULL)
        })
        
        if (!is.null(gsva_scores)) {
          
          gsva_by_celltype <- lapply(unique(cell_types_sampled), function(ct) {
            ct_cells <- which(cell_types_sampled == ct)
            if (length(ct_cells) > 0) {
              rowMeans(gsva_scores[, ct_cells, drop = FALSE])
            } else {
              NULL
            }
          })
          
          names(gsva_by_celltype) <- unique(cell_types_sampled)
          gsva_by_celltype <- gsva_by_celltype[!sapply(gsva_by_celltype, is.null)]
          
          if (length(gsva_by_celltype) > 0) {
            gsva_matrix <- do.call(cbind, gsva_by_celltype)
            
            results$gsva_analysis <- list(
              scores_matrix = gsva_matrix,
              top_pathways_per_type = apply(gsva_matrix, 2, function(x) {
                names(sort(x, decreasing = TRUE)[1:min(5, length(x))])
              })
            )
            
            message("    GSVA完成，分析了 ", ncol(gsva_matrix), " 个细胞类型")
          }
        }
      } else {
        message("    没有足够的基因集可用于GSVA分析")
      }
    }
  }
  
  return(results)
}

# 模块 4: 预测置信度评估
assess_prediction_confidence <- function(seurat_obj) {
  message("\n=== 模块 4: 预测置信度评估 ===")
  
  results <- list()
  
  score_fields <- c("monaco_scores", "blueprint_scores")
  available_scores <- score_fields[score_fields %in% colnames(seurat_obj@meta.data)]
  
  if (length(available_scores) > 0) {
    
    for (score_field in available_scores) {
      
      scores <- seurat_obj@meta.data[[score_field]]
      
      if (is.matrix(scores) || is.data.frame(scores)) {
        
        max_scores <- apply(scores, 1, max, na.rm = TRUE)
        
        score_diffs <- apply(scores, 1, function(x) {
          sorted <- sort(x, decreasing = TRUE)
          if (length(sorted) >= 2) {
            sorted[1] - sorted[2]
          } else {
            0
          }
        })
        
        score_entropy <- apply(scores, 1, function(x) {
          x <- x / sum(x)
          -sum(x * log2(x + 1e-10))
        })
        
        results[[score_field]] <- list(
          mean_max_score = mean(max_scores, na.rm = TRUE),
          median_max_score = median(max_scores, na.rm = TRUE),
          mean_score_diff = mean(score_diffs, na.rm = TRUE),
          mean_entropy = mean(score_entropy, na.rm = TRUE),
          low_confidence_pct = sum(max_scores < 0.5, na.rm = TRUE) / length(max_scores),
          high_confidence_pct = sum(max_scores > 0.8, na.rm = TRUE) / length(max_scores),
          score_distribution = max_scores
        )
        
        message("  ", score_field, ":")
        message("    平均最大得分: ", round(mean(max_scores, na.rm = TRUE), 3))
        message("    高置信度 %: ", 
                round(sum(max_scores > 0.8, na.rm = TRUE) / length(max_scores) * 100, 2))
      }
    }
  }
  
  if (all(c("monaco_fine", "blueprint_fine") %in% colnames(seurat_obj@meta.data))) {
    
    monaco <- as.character(seurat_obj$monaco_fine)
    blueprint <- as.character(seurat_obj$blueprint_fine)
    
    consensus <- monaco == blueprint
    
    valid_idx <- !grepl("Unknown|Unassigned", monaco) & 
                 !grepl("Unknown|Unassigned", blueprint)
    
    consensus_valid <- consensus[valid_idx]
    
    results$method_consensus <- list(
      overall_consensus_rate = mean(consensus),
      valid_consensus_rate = mean(consensus_valid),
      n_consensus = sum(consensus),
      n_conflict = sum(!consensus),
      consensus_by_type = table(
        Consensus = consensus,
        CellType = seurat_obj$fine_cell_type
      )
    )
    
    message("  方法一致性:")
    message("    总体一致率: ", round(mean(consensus) * 100, 2), "%")
  }
  
  if ("umap" %in% names(seurat_obj@reductions) && 
      "fine_cell_type" %in% colnames(seurat_obj@meta.data)) {
    
    message("  计算邻域纯度...")
    
    if (ncol(seurat_obj) > 5000) {
      set.seed(42)
      sample_cells <- safe_sample(colnames(seurat_obj), 5000)
      umap_coords <- Embeddings(seurat_obj, "umap")[sample_cells, ]
      cell_types <- seurat_obj$fine_cell_type[sample_cells]
    } else {
      umap_coords <- Embeddings(seurat_obj, "umap")
      cell_types <- seurat_obj$fine_cell_type
    }
    
    if (requireNamespace("RANN", quietly = TRUE)) {
      library(RANN)
      k <- min(30, nrow(umap_coords) - 1)
      nn_result <- nn2(umap_coords, umap_coords, k = k + 1)
      
      neighborhood_purity <- sapply(1:nrow(umap_coords), function(i) {
        neighbors <- nn_result$nn.idx[i, -1]
        neighbor_types <- cell_types[neighbors]
        sum(neighbor_types == cell_types[i]) / length(neighbor_types)
      })
      
      results$neighborhood_confidence <- list(
        mean_purity = mean(neighborhood_purity),
        median_purity = median(neighborhood_purity),
        low_purity_pct = sum(neighborhood_purity < 0.5) / length(neighborhood_purity),
        purity_distribution = neighborhood_purity
      )
      
      message("    平均邻域纯度: ", round(mean(neighborhood_purity), 3))
    }
  }
  
  return(results)
}

# 模块 5: 批次效应评估
assess_batch_effects <- function(seurat_obj) {
  message("\n=== 模块 5: 批次效应评估 ===")
  
  results <- list()
  
  if (!"batch" %in% colnames(seurat_obj@meta.data)) {
    message("  未找到批次信息")
    return(NULL)
  }
  
  if ("fine_cell_type" %in% colnames(seurat_obj@meta.data)) {
    
    dt <- data.table(
      batch = as.character(seurat_obj$batch),
      celltype = as.character(seurat_obj$fine_cell_type)
    )
    
    batch_celltype_table <- table(dt$batch, dt$celltype)
    batch_celltype_prop <- prop.table(batch_celltype_table, margin = 1)
    
    if (nrow(batch_celltype_prop) > 1) {
      batch_cor <- cor(t(batch_celltype_prop))
      batch_dist <- as.dist(1 - batch_cor)
      
      results$batch_similarity <- list(
        correlation_matrix = batch_cor,
        distance_matrix = as.matrix(batch_dist),
        mean_correlation = mean(batch_cor[upper.tri(batch_cor)])
      )
      
      message("  平均批次相关性: ", round(mean(batch_cor[upper.tri(batch_cor)]), 3))
    }
    
    chi_test <- suppressWarnings(chisq.test(batch_celltype_table))
    
    n <- sum(batch_celltype_table)
    min_dim <- min(dim(batch_celltype_table)) - 1
    cramers_v <- sqrt(chi_test$statistic / (n * min_dim))
    
    results$batch_independence <- list(
      chi_square = as.numeric(chi_test$statistic),
      p_value = chi_test$p.value,
      cramers_v = as.numeric(cramers_v),
      contingency_table = batch_celltype_table
    )
    
    message("  批次独立性检验:")
    message("    卡方检验 p值: ", format(chi_test$p.value, scientific = TRUE))
    message("    Cramér's V: ", round(as.numeric(cramers_v), 3))
  }
  
  if ("umap" %in% names(seurat_obj@reductions)) {
    
    if (ncol(seurat_obj) > 5000) {
      set.seed(42)
      sample_cells <- safe_sample(colnames(seurat_obj), 5000)
      umap_coords <- Embeddings(seurat_obj, "umap")[sample_cells, ]
      batch_labels <- seurat_obj$batch[sample_cells]
    } else {
      umap_coords <- Embeddings(seurat_obj, "umap")
      batch_labels <- seurat_obj$batch
    }
    
    mixing_metric <- calculate_mixing_metric(umap_coords, batch_labels)
    
    results$batch_mixing <- mixing_metric
    
    if (!is.na(mixing_metric$mixing_score)) {
      message("  批次混合得分: ", round(mixing_metric$mixing_score, 3))
    }
  }
  
  return(results)
}

# 模块 6: 稳健性评估（修复版）
assess_robustness <- function(seurat_obj,
                              n_iterations = 5,
                              sampling_rates = c(0.7, 0.8, 0.9),
                              max_cells = 10000) {
  message("\n=== 模块 6: 稳健性评估 ===")
  
  results <- list()
  
  if (!"fine_cell_type" %in% colnames(seurat_obj@meta.data)) {
    message("  未找到细胞类型注释")
    return(NULL)
  }
  
  if (ncol(seurat_obj) > max_cells) {
    set.seed(42)
    sample_cells <- safe_sample(colnames(seurat_obj), max_cells)
    calc_obj <- subset(seurat_obj, cells = sample_cells)
  } else {
    calc_obj <- seurat_obj
  }
  
  original_labels <- calc_obj$fine_cell_type
  n_cells <- ncol(calc_obj)
  
  message("  测试子采样稳定性...")
  
  sampling_stability <- list()
  
  for (rate in sampling_rates) {
    
    message("    采样率: ", rate * 100, "%")
    
    ari_scores <- sapply(1:n_iterations, function(iter) {
      
      set.seed(42 + iter)
      sample_size <- floor(n_cells * rate)
      sample_cells <- safe_sample(colnames(calc_obj), sample_size)
      
      if (length(sample_cells) == 0) return(NA)
      
      subset_obj <- subset(calc_obj, cells = sample_cells)
      
      subset_obj <- tryCatch({
        subset_obj <- FindNeighbors(subset_obj, dims = 1:20, 
                                    k.param = min(20, ncol(subset_obj)-1), 
                                    verbose = FALSE)
        subset_obj <- FindClusters(subset_obj, resolution = 0.5, 
                                  algorithm = 1, verbose = FALSE)
        subset_obj
      }, error = function(e) {
        message("      迭代 ", iter, " 失败: ", e$message)
        return(NULL)
      })
      
      if (is.null(subset_obj)) return(NA)
      
      common_cells <- intersect(sample_cells, colnames(calc_obj))
      
      if (length(common_cells) < 10) return(NA)
      
      labels1 <- as.character(original_labels[common_cells])
      labels2 <- as.character(subset_obj$seurat_clusters[common_cells])
      
      ari <- tryCatch({
        fossil::adj.rand.index(labels1, labels2)
      }, error = function(e) NA)
      
      return(ari)
    })
    
    ari_scores <- ari_scores[!is.na(ari_scores)]
    
    if (length(ari_scores) > 0) {
      sampling_stability[[as.character(rate)]] <- list(
        mean_ari = mean(ari_scores),
        sd_ari = sd(ari_scores),
        min_ari = min(ari_scores),
        max_ari = max(ari_scores),
        cv = sd(ari_scores) / mean(ari_scores),
        scores = ari_scores
      )
      
      message("      平均 ARI: ", round(mean(ari_scores), 3), 
              " (SD: ", round(sd(ari_scores), 3), ")")
    }
  }
  
  results$sampling_stability <- sampling_stability
  
  message("  测试参数敏感性...")
  
  resolution_values <- c(0.3, 0.5, 0.8, 1.0)
  resolution_stability <- list()
  
  for (res in resolution_values) {
    
    test_obj <- tryCatch({
      test_obj <- FindClusters(calc_obj, resolution = res, 
                              algorithm = 1, verbose = FALSE)
      test_obj
    }, error = function(e) NULL)
    
    if (!is.null(test_obj)) {
      
      ari <- tryCatch({
        fossil::adj.rand.index(
          as.character(original_labels),
          as.character(test_obj$seurat_clusters)
        )
      }, error = function(e) NA)
      
      if (!is.na(ari)) {
        resolution_stability[[as.character(res)]] <- list(
          resolution = res,
          n_clusters = length(unique(test_obj$seurat_clusters)),
          ari_with_annotation = ari
        )
        
        message("    分辨率 ", res, ": ARI = ", round(ari, 3), 
                " (", length(unique(test_obj$seurat_clusters)), " clusters)")
      }
    }
  }
  
  results$parameter_sensitivity <- resolution_stability
  
  rm(calc_obj)
  gc()
  
  return(results)
}

# ============================================================================
# 4. 可视化模块（增强版，专业英文标注）
# ============================================================================

generate_comprehensive_visualizations <- function(seurat_obj, 
                                                  assessment_results,
                                                  output_dir,
                                                  cancer_type,
                                                  max_cells_plot = 10000) {
  message("\n=== 生成可视化图表 ===")
  
  plot_dir <- file.path(output_dir, "figures", cancer_type)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 设置高质量图形参数
  publication_theme <- theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black", size = 1)
    )
  
  # 安全采样用于绘图
  n_cells <- ncol(seurat_obj)
  
  if (n_cells > max_cells_plot) {
    message("为绘图采样细胞 (", n_cells, " -> ", max_cells_plot, ")...")
    
    if ("fine_cell_type" %in% colnames(seurat_obj@meta.data)) {
      set.seed(42)
      sample_cells <- c()
      
      for (ct in unique(seurat_obj$fine_cell_type)) {
        ct_cells <- colnames(seurat_obj)[seurat_obj$fine_cell_type == ct]
        n_ct <- length(ct_cells)
        
        n_sample <- ceiling(max_cells_plot * n_ct / n_cells)
        n_sample <- max(10, min(n_sample, n_ct))
        
        sampled <- safe_sample(ct_cells, n_sample)
        sample_cells <- c(sample_cells, sampled)
      }
      
      if (length(sample_cells) > max_cells_plot) {
        set.seed(42)
        sample_cells <- safe_sample(sample_cells, max_cells_plot)
      }
      
      plot_obj <- subset(seurat_obj, cells = sample_cells)
      message("分层采样: ", ncol(plot_obj), " 个细胞")
      
    } else {
      set.seed(42)
      sample_cells <- safe_sample(colnames(seurat_obj), max_cells_plot)
      plot_obj <- subset(seurat_obj, cells = sample_cells)
      message("随机采样: ", ncol(plot_obj), " 个细胞")
    }
  } else {
    plot_obj <- seurat_obj
    message("使用全部 ", n_cells, " 个细胞进行绘图")
  }
  
  # 图1: 一致性矩阵热图
  if (!is.null(assessment_results$consistency$consistency_matrices)) {
    
    matrices <- assessment_results$consistency$consistency_matrices
    
    png(file.path(plot_dir, "01_Annotation_Consistency_Matrices.png"), 
        width = 3500, height = 2500, res = 300)
    
    par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
    
    for (metric_name in names(matrices)) {
      mat <- matrices[[metric_name]]
      
      if (!is.null(mat) && !all(is.na(mat))) {
        
        col_fun <- colorRamp2(
          c(min(mat, na.rm = TRUE), 
            mean(range(mat, na.rm = TRUE)), 
            max(mat, na.rm = TRUE)),
          c("#2166AC", "white", "#B2182B")
        )
        
        ht <- Heatmap(
          mat,
          name = metric_name,
          col = col_fun,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if (!is.na(mat[i, j])) {
              grid.text(sprintf("%.2f", mat[i, j]), x, y, 
                       gp = gpar(fontsize = 11, fontface = "bold"))
            }
          },
          column_title = paste("Annotation Consistency -", metric_name),
          column_title_gp = gpar(fontsize = 14, fontface = "bold"),
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 11),
          column_names_gp = gpar(fontsize = 11),
          column_names_rot = 45,
          border = TRUE,
          heatmap_legend_param = list(
            title = metric_name,
            legend_direction = "horizontal",
            legend_width = unit(6, "cm"),
            title_gp = gpar(fontsize = 12, fontface = "bold"),
            labels_gp = gpar(fontsize = 10)
          )
        )
        
        draw(ht, heatmap_legend_side = "bottom")
      }
    }
    
    dev.off()
    
    message("✓ 一致性矩阵已保存")
  }
  
  # 图2: 细胞类型分布
  if ("fine_cell_type" %in% colnames(plot_obj@meta.data)) {
    
    dt <- data.table(CellType = as.character(plot_obj$fine_cell_type))
    celltype_counts <- dt[, .N, by = CellType][order(-N)]
    
    if (nrow(celltype_counts) > 40) {
      celltype_counts <- celltype_counts[1:40]
    }
    
    celltype_counts[, Percentage := N / sum(N) * 100]
    
    p_dist <- ggplot(celltype_counts, aes(x = reorder(CellType, N), y = N)) +
      geom_bar(stat = "identity", aes(fill = Percentage), width = 0.8) +
      geom_text(aes(label = N), hjust = -0.2, size = 4, fontface = "bold") +
      scale_fill_viridis_c(option = "plasma", begin = 0.2, end = 0.9) +
      coord_flip() +
      labs(
        title = paste(cancer_type, "Cell Type Distribution"),
        subtitle = paste("Total cells:", format(sum(celltype_counts$N), big.mark = ",")),
        x = "Cell Type",
        y = "Number of Cells",
        fill = "Percentage (%)"
      ) +
      publication_theme +
      theme(
        axis.text.y = element_text(size = 10),
        legend.position = "right"
      )
    
    ggsave(file.path(plot_dir, "02_Cell_Type_Distribution.png"), p_dist,
           width = 14, height = max(10, nrow(celltype_counts) * 0.3), 
           dpi = 300, bg = "white")
    
    message("✓ 细胞类型分布已保存")
  }
  
  # 图3: Marker验证热图
  if (!is.null(assessment_results$biological_validity$marker_validation)) {
    
    marker_data <- assessment_results$biological_validity$marker_validation
    
    if (length(marker_data) > 0) {
      
      all_genes <- unique(unlist(lapply(marker_data, function(x) x$markers_tested)))
      all_celltypes <- names(marker_data)
      
      if (length(all_celltypes) > 30) {
        all_celltypes <- all_celltypes[1:30]
        marker_data <- marker_data[all_celltypes]
      }
      
      log2fc_matrix <- matrix(0, nrow = length(all_genes), 
                             ncol = length(all_celltypes))
      rownames(log2fc_matrix) <- all_genes
      colnames(log2fc_matrix) <- all_celltypes
      
      for (ct in all_celltypes) {
        scores <- marker_data[[ct]]$marker_scores
        genes <- scores$gene
        log2fc_matrix[genes, ct] <- scores$log2fc
      }
      
      keep_genes <- rowSums(log2fc_matrix > 0.5) > 0
      log2fc_matrix <- log2fc_matrix[keep_genes, , drop = FALSE]
      
      if (nrow(log2fc_matrix) > 0) {
        
        png(file.path(plot_dir, "03_Marker_Gene_Validation_Heatmap.png"), 
            width = max(3000, ncol(log2fc_matrix) * 100), 
            height = max(2500, nrow(log2fc_matrix) * 60), 
            res = 300)
        
        col_fun <- colorRamp2(
          c(-1, 0, 1, 2, 3),
          c("#2166AC", "white", "#FEE090", "#FC8D59", "#B2182B")
        )
        
        ht <- Heatmap(
          log2fc_matrix,
          name = "Log2FC",
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          show_row_names = TRUE,
          show_column_names = TRUE,
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 11, fontface = "bold"),
          column_names_rot = 45,
          column_title = "Marker Gene Expression Validation",
          column_title_gp = gpar(fontsize = 16, fontface = "bold"),
          border = TRUE,
          heatmap_legend_param = list(
            title = "Log2 Fold Change",
            legend_direction = "vertical",
            title_gp = gpar(fontsize = 12, fontface = "bold"),
            labels_gp = gpar(fontsize = 10)
          )
        )
        
        draw(ht)
        
        dev.off()
        
        message("✓ Marker验证热图已保存")
      }
    }
  }
  
  # 图4: 聚类质量
  if (!is.null(assessment_results$biological_validity$cluster_annotation_alignment)) {
    
    cluster_purity <- assessment_results$biological_validity$cluster_annotation_alignment$cluster_purity
    
    if (!is.null(cluster_purity) && nrow(cluster_purity) > 0) {
      
      cluster_purity_df <- as.data.frame(cluster_purity)
      
      p_purity <- ggplot(cluster_purity_df, aes(x = reorder(cluster, purity), y = purity)) +
        geom_bar(stat = "identity", fill = "#4575B4", width = 0.7) +
        geom_hline(yintercept = 0.7, linetype = "dashed", color = "#D73027", size = 1.2) +
        geom_text(aes(label = round(purity, 2)), hjust = -0.2, size = 3.5, fontface = "bold") +
        coord_flip() +
        ylim(0, 1.1) +
        labs(
          title = paste(cancer_type, "Cluster Purity"),
          subtitle = "Red dashed line indicates 0.7 threshold",
          x = "Cluster",
          y = "Purity"
        ) +
        publication_theme
      
      p_entropy <- ggplot(cluster_purity_df, aes(x = reorder(cluster, entropy), y = entropy)) +
        geom_bar(stat = "identity", fill = "#F46D43", width = 0.7) +
        geom_text(aes(label = round(entropy, 2)), hjust = -0.2, size = 3.5, fontface = "bold") +
        coord_flip() +
        labs(
          title = paste(cancer_type, "Cluster Entropy"),
          subtitle = "Lower entropy indicates better separation",
          x = "Cluster",
          y = "Entropy"
        ) +
        publication_theme
      
      p_combined <- p_purity / p_entropy
      
      ggsave(file.path(plot_dir, "04_Cluster_Quality_Metrics.png"), p_combined,
             width = 12, height = 14, dpi = 300, bg = "white")
      
      message("✓ 聚类质量图已保存")
    }
  }
  
  # 图5: 置信度分布（增强版）
  if (!is.null(assessment_results$confidence)) {
    
    confidence_plots <- list()
    
    for (score_field in names(assessment_results$confidence)) {
      
      if (grepl("_scores$", score_field)) {
        
        score_data <- assessment_results$confidence[[score_field]]
        
        if (!is.null(score_data$score_distribution)) {
          
          score_df <- data.frame(
            MaxScore = score_data$score_distribution
          )
          
          p <- ggplot(score_df, aes(x = MaxScore)) +
            geom_histogram(aes(y = ..density..), bins = 50, 
                          fill = "#74ADD1", color = "#313695", alpha = 0.7) +
            geom_density(color = "#D73027", size = 1.5) +
            geom_vline(xintercept = 0.5, linetype = "dashed", 
                      color = "#F46D43", size = 1) +
            geom_vline(xintercept = 0.8, linetype = "dashed", 
                      color = "#1A9850", size = 1) +
            annotate("text", x = 0.5, y = Inf, label = "Low Confidence", 
                    vjust = 2, size = 4, color = "#F46D43", fontface = "bold") +
            annotate("text", x = 0.8, y = Inf, label = "High Confidence", 
                    vjust = 2, size = 4, color = "#1A9850", fontface = "bold") +
            labs(
              title = paste("Prediction Confidence -", score_field),
              subtitle = paste0(
                "Mean: ", round(score_data$mean_max_score, 3),
                " | High Confidence: ", round(score_data$high_confidence_pct * 100, 1), "%"
              ),
              x = "Maximum Prediction Score",
              y = "Density"
            ) +
            publication_theme
          
          confidence_plots[[score_field]] <- p
        }
      }
    }
    
    if (!is.null(assessment_results$confidence$neighborhood_confidence)) {
      
      nb_data <- assessment_results$confidence$neighborhood_confidence
      
      if (!is.null(nb_data$purity_distribution)) {
        
        nb_df <- data.frame(
          NeighborhoodPurity = nb_data$purity_distribution
        )
        
        p_nb <- ggplot(nb_df, aes(x = NeighborhoodPurity)) +
          geom_histogram(aes(y = ..density..), bins = 50,
                        fill = "#91CF60", color = "#1A9850", alpha = 0.7) +
          geom_density(color = "#D73027", size = 1.5) +
          geom_vline(xintercept = nb_data$mean_purity, 
                    linetype = "dashed", color = "#313695", size = 1.2) +
          annotate("text", x = nb_data$mean_purity, y = Inf, 
                  label = paste("Mean:", round(nb_data$mean_purity, 3)), 
                  vjust = 2, size = 4, color = "#313695", fontface = "bold") +
          labs(
            title = "Neighborhood Purity Distribution",
            subtitle = paste0(
              "Mean: ", round(nb_data$mean_purity, 3),
              " | Low Purity: ", round(nb_data$low_purity_pct * 100, 1), "%"
            ),
            x = "Neighborhood Purity",
            y = "Density"
          ) +
          publication_theme
        
        confidence_plots[["neighborhood"]] <- p_nb
      }
    }
    
    if (length(confidence_plots) > 0) {
      
      combined_confidence <- wrap_plots(confidence_plots, ncol = 2)
      
      ggsave(file.path(plot_dir, "05_Confidence_Distributions.png"), 
             combined_confidence,
             width = 16, height = 7 * ceiling(length(confidence_plots) / 2),
             dpi = 300, bg = "white")
      
      message("✓ 置信度分布图已保存")
    }
  }
  
  # 图6: 质量汇总雷达图和柱状图
  quality_metrics <- data.frame(
    Metric = character(),
    Score = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(assessment_results$consistency$overall_consistency)) {
    quality_metrics <- rbind(quality_metrics, data.frame(
      Metric = "Annotation\nConsistency",
      Score = assessment_results$consistency$overall_consistency
    ))
  }
  
  if (!is.null(assessment_results$biological_validity$cluster_annotation_alignment)) {
    quality_metrics <- rbind(quality_metrics, data.frame(
      Metric = "Cluster\nPurity",
      Score = assessment_results$biological_validity$cluster_annotation_alignment$mean_cluster_purity
    ))
  }
  
  if (!is.null(assessment_results$biological_validity$spatial_coherence)) {
    sep_index <- assessment_results$biological_validity$spatial_coherence$separation_index
    if (!is.na(sep_index)) {
      sep_score <- min(1, sep_index / 3)
      quality_metrics <- rbind(quality_metrics, data.frame(
        Metric = "Spatial\nCoherence",
        Score = sep_score
      ))
    }
  }
  
  if (!is.null(assessment_results$biological_validity$marker_validation_summary)) {
    avg_fc <- assessment_results$biological_validity$marker_validation_summary$avg_log2fc
    fc_score <- min(1, max(0, avg_fc / 2))
    quality_metrics <- rbind(quality_metrics, data.frame(
      Metric = "Marker\nValidation",
      Score = fc_score
    ))
  }
  
  if (!is.null(assessment_results$robustness$sampling_stability)) {
    stab_scores <- sapply(assessment_results$robustness$sampling_stability, 
                         function(x) ifelse(is.na(x$mean_ari), 0, x$mean_ari))
    if (length(stab_scores) > 0) {
      quality_metrics <- rbind(quality_metrics, data.frame(
        Metric = "Sampling\nStability",
        Score = mean(stab_scores)
      ))
    }
  }
  
  if (!is.null(assessment_results$confidence)) {
    conf_scores <- c()
    for (field in names(assessment_results$confidence)) {
      if (grepl("_scores$", field)) {
        conf_scores <- c(conf_scores, 
                        assessment_results$confidence[[field]]$mean_max_score)
      }
    }
    if (length(conf_scores) > 0) {
      quality_metrics <- rbind(quality_metrics, data.frame(
        Metric = "Prediction\nConfidence",
        Score = mean(conf_scores)
      ))
    }
  }
  
  if (nrow(quality_metrics) > 0) {
    
    quality_metrics$Score <- pmax(0, pmin(1, quality_metrics$Score))
    quality_metrics$Metric <- factor(quality_metrics$Metric, 
                                     levels = quality_metrics$Metric)
    
    # 雷达图
    p_radar <- ggplot(quality_metrics, aes(x = Metric, y = Score, group = 1)) +
      geom_polygon(fill = "#74ADD1", alpha = 0.4, color = "#313695", size = 1.5) +
      geom_point(size = 5, color = "#D73027", shape = 16) +
      geom_hline(yintercept = c(0.5, 0.7, 0.9), 
                linetype = "dashed", color = "gray60", alpha = 0.6, size = 0.8) +
      coord_polar() +
      ylim(0, 1) +
      labs(
        title = paste(cancer_type, "Overall Quality Assessment"),
        subtitle = paste("Mean Score:", round(mean(quality_metrics$Score), 3))
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray40"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white")
      )
    
    # 柱状图
    p_bar <- ggplot(quality_metrics, aes(x = reorder(Metric, Score), y = Score)) +
      geom_bar(stat = "identity", aes(fill = Score), width = 0.7) +
      geom_text(aes(label = round(Score, 3)), hjust = -0.2, size = 5, fontface = "bold") +
      geom_hline(yintercept = 0.7, linetype = "dashed", color = "#1A9850", size = 1.2) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "#F46D43", size = 1.2) +
      scale_fill_gradient2(
        low = "#D73027", mid = "#FEE090", high = "#1A9850",
        midpoint = 0.7,
        limits = c(0, 1)
      ) +
      coord_flip() +
      ylim(0, 1.15) +
      labs(
        title = paste(cancer_type, "Quality Metrics Summary"),
        x = "",
        y = "Score (0-1)"
      ) +
      publication_theme +
      theme(
        legend.position = "none",
        axis.text.y = element_text(size = 13, face = "bold")
      )
    
    p_combined_metrics <- p_radar / p_bar
    
    ggsave(file.path(plot_dir, "06_Quality_Summary.png"), p_combined_metrics,
           width = 14, height = 16, dpi = 300, bg = "white")
    
    message("✓ 质量汇总图已保存")
  }
  
  # 图7: 批次效应
  if (!is.null(assessment_results$batch_effects)) {
    
    if (!is.null(assessment_results$batch_effects$batch_independence$contingency_table)) {
      
      cont_table <- assessment_results$batch_effects$batch_independence$contingency_table
      cont_df <- as.data.frame(cont_table)
      names(cont_df) <- c("Batch", "CellType", "Count")
      
      top_types <- cont_df %>%
        group_by(CellType) %>%
        summarise(Total = sum(Count)) %>%
        arrange(desc(Total)) %>%
        head(20) %>%
        pull(CellType)
      
      cont_df <- cont_df %>% filter(CellType %in% top_types)
      
      n_types <- length(unique(cont_df$CellType))
      color_palette <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", 
                                          "#984EA3", "#FF7F00", "#FFFF33", 
                                          "#A65628", "#F781BF", "#999999"))(n_types)
      
      p_batch_dist <- ggplot(cont_df, aes(x = Batch, y = Count, fill = CellType)) +
        geom_bar(stat = "identity", position = "fill", width = 0.8) +
        scale_fill_manual(values = color_palette) +
        scale_y_continuous(labels = scales::percent) +
        labs(
          title = "Cell Type Distribution Across Batches",
          subtitle = "Showing top 20 most abundant cell types",
          y = "Proportion",
          x = "Batch",
          fill = "Cell Type"
        ) +
        publication_theme +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          legend.position = "right",
          legend.key.size = unit(0.5, "cm")
        )
      
      ggsave(file.path(plot_dir, "07_Batch_Effects.png"), p_batch_dist,
             width = 14, height = 10, dpi = 300, bg = "white")
      
      message("✓ 批次效应图已保存")
    }
  }
  
  # 图8: 稳健性分析（修复版）
  if (!is.null(assessment_results$robustness$sampling_stability)) {
    
    stab_data <- assessment_results$robustness$sampling_stability
    
    if (length(stab_data) > 0) {
      stab_df <- data.frame()
      
      for (rate in names(stab_data)) {
        scores <- stab_data[[rate]]$scores
        if (!is.null(scores) && length(scores) > 0) {
          stab_df <- rbind(stab_df, data.frame(
            SamplingRate = rep(paste0(as.numeric(rate) * 100, "%"), length(scores)),
            ARI = scores
          ))
        }
      }
      
      if (nrow(stab_df) > 0) {
        
        stab_df$SamplingRate <- factor(stab_df$SamplingRate, 
                                       levels = unique(stab_df$SamplingRate))
        
        p_stability <- ggplot(stab_df, aes(x = SamplingRate, y = ARI)) +
          geom_violin(aes(fill = SamplingRate), alpha = 0.6, trim = FALSE) +
          geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", alpha = 0.8) +
          geom_jitter(width = 0.15, alpha = 0.4, size = 2.5, color = "#313695") +
          scale_fill_brewer(palette = "Set2") +
          geom_hline(yintercept = 0.7, linetype = "dashed", 
                    color = "#D73027", size = 1.2) +
          annotate("text", x = Inf, y = 0.7, label = "Good Threshold", 
                  vjust = -0.5, hjust = 1.1, size = 4.5, color = "#D73027", fontface = "bold") +
          labs(
            title = "Subsampling Stability Analysis",
            subtitle = "ARI between original and subsampled annotations",
            x = "Sampling Rate",
            y = "Adjusted Rand Index (ARI)"
          ) +
          publication_theme +
          theme(
            legend.position = "none"
          )
        
        ggsave(file.path(plot_dir, "08_Sampling_Stability.png"), p_stability,
               width = 12, height = 8, dpi = 300, bg = "white")
        
        message("✓ 采样稳定性图已保存")
      }
    }
  }
  
  # 图9: 参数敏感性分析（新增）
  if (!is.null(assessment_results$robustness$parameter_sensitivity)) {
    
    param_data <- assessment_results$robustness$parameter_sensitivity
    
    if (length(param_data) > 0) {
      param_df <- data.frame(
        Resolution = sapply(param_data, function(x) x$resolution),
        ARI = sapply(param_data, function(x) x$ari_with_annotation),
        NClusters = sapply(param_data, function(x) x$n_clusters)
      )
      
      p_param1 <- ggplot(param_df, aes(x = Resolution, y = ARI)) +
        geom_line(color = "#313695", size = 1.5) +
        geom_point(size = 4, color = "#D73027") +
        geom_hline(yintercept = 0.7, linetype = "dashed", color = "#1A9850", size = 1) +
        labs(
          title = "Parameter Sensitivity Analysis",
          subtitle = "Effect of resolution on annotation quality",
          x = "Clustering Resolution",
          y = "Adjusted Rand Index (ARI)"
        ) +
        publication_theme
      
      p_param2 <- ggplot(param_df, aes(x = Resolution, y = NClusters)) +
        geom_line(color = "#4575B4", size = 1.5) +
        geom_point(size = 4, color = "#F46D43") +
        labs(
          title = "Number of Clusters vs Resolution",
          x = "Clustering Resolution",
          y = "Number of Clusters"
        ) +
        publication_theme
      
      p_combined_param <- p_param1 / p_param2
      
      ggsave(file.path(plot_dir, "09_Parameter_Sensitivity.png"), p_combined_param,
             width = 12, height = 10, dpi = 300, bg = "white")
      
      message("✓ 参数敏感性图已保存")
    }
  }
  
  # 图10: 细胞类型比例分布（新增）
  if (!is.null(assessment_results$biological_validity$proportion_analysis)) {
    
    prop_data <- assessment_results$biological_validity$proportion_analysis
    
    if (!is.null(prop_data$type_proportions)) {
      prop_df <- data.frame(
        CellType = names(prop_data$type_proportions),
        Proportion = as.numeric(prop_data$type_proportions)
      )
      prop_df <- prop_df[order(-prop_df$Proportion), ]
      prop_df <- head(prop_df, 30)
      
      p_prop <- ggplot(prop_df, aes(x = reorder(CellType, Proportion), y = Proportion)) +
        geom_bar(stat = "identity", fill = "#4575B4", width = 0.7) +
        geom_text(aes(label = paste0(round(Proportion * 100, 1), "%")), 
                 hjust = -0.1, size = 3.5, fontface = "bold") +
        coord_flip() +
        labs(
          title = paste(cancer_type, "Cell Type Proportion Distribution"),
          subtitle = paste0("Shannon Diversity: ", round(prop_data$normalized_shannon, 3),
                          " | Gini Coefficient: ", round(prop_data$gini_coefficient, 3)),
          x = "Cell Type",
          y = "Proportion"
        ) +
        publication_theme
      
      ggsave(file.path(plot_dir, "10_Cell_Type_Proportions.png"), p_prop,
             width = 12, height = max(8, nrow(prop_df) * 0.3), dpi = 300, bg = "white")
      
      message("✓ 细胞类型比例图已保存")
    }
  }
  
  # 图11: 综合质量评估雷达图（详细版，新增）
  if (nrow(quality_metrics) > 0) {
    
    # 创建更详细的质量指标
    detailed_metrics <- quality_metrics
    
    # 添加颜色分级
    detailed_metrics$Grade <- cut(detailed_metrics$Score, 
                                  breaks = c(0, 0.5, 0.7, 0.85, 1),
                                  labels = c("Poor", "Moderate", "Good", "Excellent"))
    
    p_detailed_radar <- ggplot(detailed_metrics, aes(x = Metric, y = Score, group = 1)) +
      geom_polygon(fill = "#74ADD1", alpha = 0.3, color = "#313695", size = 1.5) +
      geom_point(aes(color = Grade), size = 6, shape = 16) +
      scale_color_manual(values = c("Poor" = "#D73027", 
                                    "Moderate" = "#FC8D59",
                                    "Good" = "#91CF60",
                                    "Excellent" = "#1A9850")) +
      geom_hline(yintercept = c(0.5, 0.7, 0.85), 
                linetype = "dashed", color = "gray60", alpha = 0.6, size = 0.8) +
      coord_polar() +
      ylim(0, 1) +
      labs(
        title = paste(cancer_type, "Detailed Quality Assessment"),
        subtitle = paste("Overall Score:", round(mean(detailed_metrics$Score), 3)),
        color = "Quality Grade"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray40"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom"
      )
    
    ggsave(file.path(plot_dir, "11_Detailed_Quality_Radar.png"), p_detailed_radar,
           width = 12, height = 12, dpi = 300, bg = "white")
    
    message("✓ 详细质量雷达图已保存")
  }
  
  rm(plot_obj)
  gc()
  
  message("✓ 所有可视化完成")
  
  return(invisible(NULL))
}

# ============================================================================
# 5. 报告生成模块（修复版）
# ============================================================================

generate_comprehensive_report <- function(cancer_type, 
                                         assessment_results,
                                         seurat_obj,
                                         output_dir) {
  message("\n=== 生成综合报告 ===")
  
  report_file <- file.path(output_dir, "reports", 
                          paste0(cancer_type, "_综合评估报告.txt"))
  
  dir.create(dirname(report_file), showWarnings = FALSE, recursive = TRUE)
  
  sink(report_file)
  
  cat("================================================================================\n")
  cat("           单细胞RNA测序注释质量综合评估报告\n")
  cat("           COMPREHENSIVE QUALITY ASSESSMENT REPORT\n")
  cat("================================================================================\n\n")
  
  cat("癌种类型: ", cancer_type, "\n")
  cat("分析日期: ", as.character(Sys.time()), "\n")
  cat("总细胞数: ", format(ncol(seurat_obj), big.mark = ","), "\n")
  
  if ("fine_cell_type" %in% colnames(seurat_obj@meta.data)) {
    cat("已注释细胞类型数: ", 
        length(unique(seurat_obj$fine_cell_type)), "\n")
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("第一部分: 注释一致性评估\n")
  cat("================================================================================\n\n")
  
  if (!is.null(assessment_results$consistency)) {
    
    if (!is.null(assessment_results$consistency$overall_consistency)) {
      consistency_score <- assessment_results$consistency$overall_consistency
      
      if (!is.na(consistency_score)) {
        cat("总体一致性得分 (ARI): ", round(consistency_score, 4), "\n")
        
        consistency_level <- if (consistency_score > 0.8) {
          "优秀 (EXCELLENT)"
        } else if (consistency_score > 0.6) {
          "良好 (GOOD)"
        } else if (consistency_score > 0.4) {
          "中等 (MODERATE)"
        } else {
          "较差 (POOR)"
        }
        
        cat("一致性等级: ", consistency_level, "\n\n")
      } else {
        cat("总体一致性得分: 无数据\n\n")
      }
    }
    
    if (!is.null(assessment_results$consistency$consistency_matrices)) {
      cat("跨方法一致性指标:\n")
      cat("----------------------------------\n")
      
      matrices <- assessment_results$consistency$consistency_matrices
      
      for (metric_name in names(matrices)) {
        mat <- matrices[[metric_name]]
        if (!is.null(mat) && !all(is.na(mat))) {
          avg_score <- mean(mat[upper.tri(mat)], na.rm = TRUE)
          if (!is.na(avg_score)) {
            cat(sprintf("  %-25s: %.4f\n", metric_name, avg_score))
          }
        }
      }
      cat("\n")
    }
    
    if (!is.null(assessment_results$consistency$hierarchy_consistency)) {
      hc <- assessment_results$consistency$hierarchy_consistency
      cat("层级一致性:\n")
      cat("-------------------------\n")
      cat("  平均纯度: ", round(hc$mean_purity, 4), "\n")
      cat("  中位纯度: ", round(hc$median_purity, 4), "\n\n")
    }
    
    if (!is.null(assessment_results$consistency$mapping_quality)) {
      mq <- assessment_results$consistency$mapping_quality
      cat("注释映射质量:\n")
      cat("---------------------------\n")
      cat("  细胞类型总数: ", mq$total_types, "\n")
      cat("  单细胞类型数: ", mq$singleton_types, "\n")
      cat("  微小类型数 (<10细胞): ", mq$tiny_types, "\n")
      cat("  未知/未分配比例: ", round(mq$unknown_proportion * 100, 2), "%\n\n")
    }
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("第二部分: 生物学合理性评估\n")
  cat("================================================================================\n\n")
  
  if (!is.null(assessment_results$biological_validity)) {
    
    if (!is.null(assessment_results$biological_validity$marker_validation_summary)) {
      mvs <- assessment_results$biological_validity$marker_validation_summary
      cat("Marker基因验证:\n")
      cat("-----------------------\n")
      cat("  验证细胞类型数: ", mvs$n_types_validated, "\n")
      cat("  平均 Log2FC: ", round(mvs$avg_log2fc, 4), "\n")
      cat("  平均百分比差异: ", round(mvs$avg_pct_difference, 4), "\n\n")
    }
    
    if (!is.null(assessment_results$biological_validity$cluster_annotation_alignment)) {
      caa <- assessment_results$biological_validity$cluster_annotation_alignment
      cat("聚类-注释对齐度:\n")
      cat("-----------------------------\n")
      cat("  平均聚类纯度: ", round(caa$mean_cluster_purity, 4), "\n")
      cat("  平均聚类熵: ", round(caa$mean_cluster_entropy, 4), "\n")
      cat("  平均每聚类细胞类型数: ", 
          round(mean(caa$cluster_purity$n_types), 2), "\n\n")
    }
    
    if (!is.null(assessment_results$biological_validity$silhouette_analysis)) {
      sil <- assessment_results$biological_validity$silhouette_analysis
      cat("轮廓系数分析:\n")
      cat("--------------------\n")
      cat("  平均轮廓系数: ", round(sil$mean_silhouette, 4), "\n")
      
      sil_interpretation <- if (sil$mean_silhouette > 0.7) {
        "强结构 (Strong structure)"
      } else if (sil$mean_silhouette > 0.5) {
        "合理结构 (Reasonable structure)"
      } else if (sil$mean_silhouette > 0.25) {
        "弱结构 (Weak structure)"
      } else {
        "无明显结构 (No substantial structure)"
      }
      cat("  解释: ", sil_interpretation, "\n\n")
    }
    
    if (!is.null(assessment_results$biological_validity$spatial_coherence)) {
      sc <- assessment_results$biological_validity$spatial_coherence
      cat("空间连贯性:\n")
      cat("------------------\n")
      if (!is.na(sc$separation_index)) {
        cat("  分离指数: ", round(sc$separation_index, 4), "\n")
      }
      if (!is.na(sc$morans_i)) {
        cat("  Moran's I: ", round(sc$morans_i, 4), "\n")
      }
      cat("\n")
    }
    
    if (!is.null(assessment_results$biological_validity$proportion_analysis)) {
      pa <- assessment_results$biological_validity$proportion_analysis
      cat("细胞类型多样性:\n")
      cat("--------------------\n")
      cat("  类型数量: ", pa$n_types, "\n")
      cat("  Shannon多样性: ", round(pa$shannon_diversity, 4), "\n")
      cat("  标准化Shannon: ", round(pa$normalized_shannon, 4), "\n")
      cat("  Simpson多样性: ", round(pa$simpson_diversity, 4), "\n")
      cat("  基尼系数: ", round(pa$gini_coefficient, 4), "\n\n")
    }
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("第三部分: 统计显著性评估\n")
  cat("================================================================================\n\n")
  
  if (!is.null(assessment_results$statistical_significance)) {
    
    if (!is.null(assessment_results$statistical_significance$deg_analysis)) {
      deg <- assessment_results$statistical_significance$deg_analysis
      cat("差异表达分析:\n")
      cat("----------------------------------\n")
      cat("  分析细胞类型数: ", length(deg), "\n")
      
      if (length(deg) > 0) {
        total_sig <- sum(sapply(deg, function(x) x$n_significant))
        avg_sig <- mean(sapply(deg, function(x) x$n_significant))
        cat("  显著基因总数: ", total_sig, "\n")
        cat("  平均每细胞类型: ", round(avg_sig, 2), "\n\n")
      }
    }
    
    if (!is.null(assessment_results$statistical_significance$gsva_analysis)) {
      cat("基因集变异分析:\n")
      cat("----------------------------\n")
      gsva <- assessment_results$statistical_significance$gsva_analysis
      cat("  分析基因集数: ", nrow(gsva$scores_matrix), "\n")
      cat("  分析细胞类型数: ", ncol(gsva$scores_matrix), "\n\n")
    }
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("第四部分: 预测置信度评估\n")
  cat("================================================================================\n\n")
  
  if (!is.null(assessment_results$confidence)) {
    
    score_fields <- grep("_scores$", names(assessment_results$confidence), value = TRUE)
    
    if (length(score_fields) > 0) {
      cat("预测得分分析:\n")
      cat("--------------------------\n")
      
      for (score_field in score_fields) {
        score_data <- assessment_results$confidence[[score_field]]
        cat("\n", score_field, ":\n")
        cat("  平均最大得分: ", round(score_data$mean_max_score, 4), "\n")
        cat("  中位最大得分: ", round(score_data$median_max_score, 4), "\n")
        cat("  高置信度 (>0.8): ", 
            round(score_data$high_confidence_pct * 100, 2), "%\n")
        cat("  低置信度 (<0.5): ", 
            round(score_data$low_confidence_pct * 100, 2), "%\n")
      }
      cat("\n")
    }
    
    if (!is.null(assessment_results$confidence$method_consensus)) {
      mc <- assessment_results$confidence$method_consensus
      cat("方法一致性:\n")
      cat("-----------------\n")
      cat("  总体一致率: ", round(mc$overall_consensus_rate * 100, 2), "%\n")
      cat("  有效一致率: ", round(mc$valid_consensus_rate * 100, 2), "%\n\n")
    }
    
    if (!is.null(assessment_results$confidence$neighborhood_confidence)) {
      nc <- assessment_results$confidence$neighborhood_confidence
      cat("邻域纯度:\n")
      cat("--------------------\n")
      cat("  平均纯度: ", round(nc$mean_purity, 4), "\n")
      cat("  中位纯度: ", round(nc$median_purity, 4), "\n")
      cat("  低纯度 (<0.5): ", round(nc$low_purity_pct * 100, 2), "%\n\n")
    }
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("第五部分: 批次效应评估\n")
  cat("================================================================================\n\n")
  
  if (!is.null(assessment_results$batch_effects)) {
    
    if (!is.null(assessment_results$batch_effects$batch_independence)) {
      bi <- assessment_results$batch_effects$batch_independence
      cat("批次独立性检验:\n")
      cat("------------------------\n")
      cat("  卡方统计量: ", round(bi$chi_square, 2), "\n")
      cat("  P值: ", format(bi$p_value, scientific = TRUE), "\n")
      cat("  Cramér's V: ", round(bi$cramers_v, 4), "\n\n")
      
      batch_effect_level <- if (bi$cramers_v < 0.1) {
        "可忽略 (NEGLIGIBLE)"
      } else if (bi$cramers_v < 0.3) {
        "较小 (SMALL)"
      } else if (bi$cramers_v < 0.5) {
        "中等 (MEDIUM)"
      } else {
        "较大 (LARGE)"
      }
      cat("  批次效应等级: ", batch_effect_level, "\n\n")
    }
    
    if (!is.null(assessment_results$batch_effects$batch_mixing)) {
      bm <- assessment_results$batch_effects$batch_mixing
      if (!is.na(bm$mixing_score)) {
        cat("批次混合分析:\n")
        cat("----------------------\n")
        cat("  混合得分: ", round(bm$mixing_score, 4), "\n")
        cat("  低混合 (<0.5): ", round(bm$low_mixing_pct * 100, 2), "%\n\n")
      }
    }
  } else {
    cat("无批次信息可用\n\n")
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("第六部分: 稳健性评估\n")
  cat("================================================================================\n\n")
  
  if (!is.null(assessment_results$robustness)) {
    
    if (!is.null(assessment_results$robustness$sampling_stability)) {
      ss <- assessment_results$robustness$sampling_stability
      cat("采样稳定性:\n")
      cat("-------------------\n")
      
      for (rate in names(ss)) {
        stab_data <- ss[[rate]]
        if (!is.null(stab_data$mean_ari) && !is.na(stab_data$mean_ari)) {
          cat(sprintf("  采样率 %.0f%%: ARI = %.4f (SD = %.4f)\n",
                     as.numeric(rate) * 100,
                     stab_data$mean_ari,
                     stab_data$sd_ari))
        }
      }
      cat("\n")
    }
    
    if (!is.null(assessment_results$robustness$parameter_sensitivity)) {
      ps <- assessment_results$robustness$parameter_sensitivity
      cat("参数敏感性:\n")
      cat("----------------------\n")
      
      for (param in names(ps)) {
        param_data <- ps[[param]]
        if (!is.null(param_data$ari_with_annotation) && 
            !is.na(param_data$ari_with_annotation)) {
          cat(sprintf("  分辨率 %.1f: ARI = %.4f (聚类数 = %d)\n",
                     param_data$resolution,
                     param_data$ari_with_annotation,
                     param_data$n_clusters))
        }
      }
      cat("\n")
    }
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("综合评估总结\n")
  cat("================================================================================\n\n")
  
  # 收集所有评分
  scores <- c()
  weights <- c()
  score_names <- c()
  
  # 1. 一致性分数
  if (!is.null(assessment_results$consistency$overall_consistency)) {
    cons_score <- assessment_results$consistency$overall_consistency
    if (!is.na(cons_score)) {
      scores <- c(scores, cons_score)
      weights <- c(weights, 0.25)
      score_names <- c(score_names, "注释一致性")
      cat("✓ 注释一致性得分: ", round(cons_score, 3), "\n")
    }
  }
  
  # 2. 生物学有效性分数
  if (!is.null(assessment_results$biological_validity$cluster_annotation_alignment)) {
    bio_score <- assessment_results$biological_validity$cluster_annotation_alignment$mean_cluster_purity
    if (!is.na(bio_score)) {
      scores <- c(scores, bio_score)
      weights <- c(weights, 0.25)
      score_names <- c(score_names, "生物学合理性")
      cat("✓ 生物学合理性得分: ", round(bio_score, 3), "\n")
    }
  }
  
  # 3. Marker验证分数
  if (!is.null(assessment_results$biological_validity$marker_validation_summary)) {
    marker_fc <- assessment_results$biological_validity$marker_validation_summary$avg_log2fc
    if (!is.na(marker_fc)) {
      marker_score <- min(1, max(0, marker_fc / 2))
      scores <- c(scores, marker_score)
      weights <- c(weights, 0.20)
      score_names <- c(score_names, "Marker验证")
      cat("✓ Marker验证得分: ", round(marker_score, 3), "\n")
    }
  }
  
  # 4. 稳健性分数
  if (!is.null(assessment_results$robustness$sampling_stability)) {
    stab_scores <- sapply(assessment_results$robustness$sampling_stability, 
                         function(x) ifelse(is.null(x$mean_ari) || is.na(x$mean_ari), 
                                          NA, x$mean_ari))
    stab_scores <- stab_scores[!is.na(stab_scores)]
    
    if (length(stab_scores) > 0) {
      robustness_score <- mean(stab_scores)
      scores <- c(scores, robustness_score)
      weights <- c(weights, 0.15)
      score_names <- c(score_names, "稳健性")
      cat("✓ 稳健性得分: ", round(robustness_score, 3), "\n")
    }
  }
  
  # 5. 置信度分数
  conf_scores <- c()
  for (field in names(assessment_results$confidence)) {
    if (grepl("_scores$", field)) {
      mean_score <- assessment_results$confidence[[field]]$mean_max_score
      if (!is.na(mean_score)) {
        conf_scores <- c(conf_scores, mean_score)
      }
    }
  }
  
  if (length(conf_scores) > 0) {
    confidence_score <- mean(conf_scores)
    scores <- c(scores, confidence_score)
    weights <- c(weights, 0.15)
    score_names <- c(score_names, "预测置信度")
    cat("✓ 预测置信度得分: ", round(confidence_score, 3), "\n")
  }
  
  cat("\n")
  
  # 计算总分
  if (length(scores) > 0) {
    # 标准化权重
    weights <- weights / sum(weights)
    overall_score <- sum(scores * weights)
    
    cat("───────────────────────────────────────\n")
    cat("总体质量得分: ", round(overall_score, 3), "\n")
    cat("───────────────────────────────────────\n\n")
    
    # 评级
    if (!is.na(overall_score)) {
      quality_grade <- if (overall_score >= 0.85) {
        "优秀 - 注释高度可靠"
      } else if (overall_score >= 0.70) {
        "良好 - 注释总体可靠"
      } else if (overall_score >= 0.55) {
        "中等 - 注释需谨慎解释"
      } else {
        "较差 - 注释需要改进"
      }
      
      cat("质量等级: ", quality_grade, "\n\n")
      
      # 详细评分表
      cat("详细评分明细:\n")
      cat("-------------------------\n")
      score_df <- data.frame(
        Component = score_names,
        Score = round(scores, 4),
        Weight = round(weights, 4),
        Contribution = round(scores * weights, 4)
      )
      
      for (i in 1:nrow(score_df)) {
        cat(sprintf("  %-20s: %.4f (权重: %.2f, 贡献: %.4f)\n",
                   score_df$Component[i],
                   score_df$Score[i],
                   score_df$Weight[i],
                   score_df$Contribution[i]))
      }
      cat("\n")
    } else {
      cat("总体质量得分: 无数据 (数据不足)\n\n")
    }
  } else {
    cat("总体质量得分: 无数据 (无有效指标)\n\n")
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("建议与改进方向\n")
  cat("================================================================================\n\n")
  
  # 基于评分给出建议
  recommendations <- c()
  
  if (!is.null(assessment_results$consistency$overall_consistency)) {
    if (!is.na(assessment_results$consistency$overall_consistency) &&
        assessment_results$consistency$overall_consistency < 0.6) {
      recommendations <- c(recommendations,
        "• 检测到较低的注释一致性。建议使用不同方法重新运行注释流程。")
    }
  }
  
  if (!is.null(assessment_results$biological_validity$cluster_annotation_alignment)) {
    if (!is.na(assessment_results$biological_validity$cluster_annotation_alignment$mean_cluster_purity) &&
        assessment_results$biological_validity$cluster_annotation_alignment$mean_cluster_purity < 0.7) {
      recommendations <- c(recommendations,
        "• 聚类纯度较低。建议调整聚类参数或进行手动优化。")
    }
  }
  
  if (!is.null(assessment_results$consistency$mapping_quality)) {
    if (assessment_results$consistency$mapping_quality$unknown_proportion > 0.2) {
      recommendations <- c(recommendations,
        paste0("• 未知/未分配细胞比例较高 (", 
               round(assessment_results$consistency$mapping_quality$unknown_proportion * 100, 1),
               "%)。建议进行额外的注释工作。"))
    }
  }
  
  if (length(recommendations) > 0) {
    cat("建议采取的措施:\n")
    for (rec in recommendations) {
      cat(rec, "\n")
    }
  } else {
    cat("• 注释质量令人满意。未检测到重大问题。\n")
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("报告结束\n")
  cat("================================================================================\n")
  cat("\n生成工具: 高级单细胞注释质量评估系统 v2.4\n")
  cat("报告日期: ", as.character(Sys.time()), "\n")
  
  sink()
  
  message("✓ 综合报告已生成: ", report_file)
  
  return(report_file)
}

# ============================================================================
# 6. 主执行函数
# ============================================================================

run_comprehensive_quality_assessment <- function(cancer_type,
                                                 annotation_dir,
                                                 output_dir,
                                                 max_cells = 100000,
                                                 sample_fraction = 1.0,
                                                 n_cores = 16) {
  
  message("\n", paste(rep("=", 80), collapse = ""))
  message(" 单细胞注释质量综合评估")
  message(" 癌种: ", cancer_type)
  message(paste(rep("=", 80), collapse = ""))
  
  initialize_environment(n_cores = n_cores)
  
  message("\n步骤 1: 加载数据...")
  
  data <- load_annotation_data_optimized(
    cancer_type, annotation_dir,
    max_cells = max_cells,
    sample_fraction = sample_fraction
  )
  
  if (is.null(data)) {
    stop("无法加载 ", cancer_type, " 的数据")
  }
  
  seurat_obj <- data$main
  
  message("处理 ", ncol(seurat_obj), " 个细胞")
  
  assessment_results <- list()
  
  message("\n步骤 2: 运行综合评估...")
  
  assessment_results$consistency <- assess_annotation_consistency(seurat_obj)
  assessment_results$biological_validity <- assess_biological_validity(seurat_obj)
  assessment_results$statistical_significance <- assess_statistical_significance(seurat_obj)
  assessment_results$confidence <- assess_prediction_confidence(seurat_obj)
  assessment_results$batch_effects <- assess_batch_effects(seurat_obj)
  assessment_results$robustness <- assess_robustness(seurat_obj)
  
  message("\n步骤 3: 生成可视化...")
  
  generate_comprehensive_visualizations(
    seurat_obj, assessment_results, output_dir, cancer_type
  )
  
  message("\n步骤 4: 生成综合报告...")
  
  report_file <- generate_comprehensive_report(
    cancer_type, assessment_results, seurat_obj, output_dir
  )
  
  message("\n步骤 5: 保存结果...")
  
  results_file <- file.path(output_dir, "tables",
                           paste0(cancer_type, "_评估结果.rds"))
  
  dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
  
  saveRDS(list(
    cancer_type = cancer_type,
    n_cells = ncol(seurat_obj),
    assessment_results = assessment_results,
    timestamp = Sys.time()
  ), results_file)
  
  message("✓ 结果已保存: ", results_file)
  
  message("\n", paste(rep("=", 80), collapse = ""))
  message(" 评估完成")
  message(paste(rep("=", 80), collapse = ""))
  
  return(assessment_results)
}

# ============================================================================
# 7. 批量处理函数
# ============================================================================

run_batch_quality_assessment <- function(annotation_dir,
                                         output_dir,
                                         max_cells = 100000,
                                         sample_fraction = 1.0,
                                         n_cores = 16,
                                         cancer_types = NULL) {
  
  message("\n", paste(rep("=", 80), collapse = ""))
  message(" 批量单细胞注释质量评估系统")
  message(paste(rep("=", 80), collapse = ""))
  
  # 安装和加载包
  install_if_missing(required_packages)
  load_packages(required_packages)
  
  # 自动检测癌种
  if (is.null(cancer_types)) {
    cancer_types <- detect_cancer_types(annotation_dir)
  }
  
  n_cancer_types <- length(cancer_types)
  message("\n将处理 ", n_cancer_types, " 个癌种")
  
  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "reports"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
  
  # 批量处理
  all_results <- list()
  
  for (i in seq_along(cancer_types)) {
    
    cancer_type <- cancer_types[i]
    
    message("\n[", i, "/", n_cancer_types, "] 正在处理: ", cancer_type)
    
    result <- tryCatch({
      
      run_comprehensive_quality_assessment(
        cancer_type = cancer_type,
        annotation_dir = annotation_dir,
        output_dir = output_dir,
        max_cells = max_cells,
        sample_fraction = sample_fraction,
        n_cores = n_cores
      )
      
    }, error = function(e) {
      message("\n✗ 处理 ", cancer_type, " 时出错: ", e$message)
      message("跳过此癌种，继续处理下一个...\n")
      return(NULL)
    })
    
    all_results[[cancer_type]] <- result
    
    gc()
  }
  
  # 生成汇总报告
  message("\n生成汇总报告...")
  
  generate_summary_report(all_results, cancer_types, output_dir)
  
  message("\n", paste(rep("=", 80), collapse = ""))
  message(" 批量评估完成")
  message(" 输出目录: ", output_dir)
  message(paste(rep("=", 80), collapse = ""))
  
  return(all_results)
}

# ============================================================================
# 8. 汇总报告生成
# ============================================================================

generate_summary_report <- function(all_results, cancer_types, output_dir) {
  
  summary_file <- file.path(output_dir, "reports", "00_综合汇总报告.txt")
  
  sink(summary_file)
  
  cat("================================================================================\n")
  cat("                      批量质量评估汇总报告\n")
  cat("                 BATCH QUALITY ASSESSMENT SUMMARY\n")
  cat("================================================================================\n\n")
  
  cat("分析日期: ", as.character(Sys.time()), "\n")
  cat("癌种数量: ", length(cancer_types), "\n")
  cat("成功处理: ", sum(!sapply(all_results, is.null)), "\n")
  cat("处理失败: ", sum(sapply(all_results, is.null)), "\n\n")
  
  cat("================================================================================\n")
  cat("各癌种质量得分汇总\n")
  cat("================================================================================\n\n")
  
  # 收集所有癌种的质量得分
  summary_table <- data.frame()
  
  for (cancer_type in cancer_types) {
    
    result <- all_results[[cancer_type]]
    
    if (is.null(result)) {
      summary_table <- rbind(summary_table, data.frame(
        CancerType = cancer_type,
        OverallScore = NA,
        ConsistencyScore = NA,
        BiologicalScore = NA,
        MarkerScore = NA,
        RobustnessScore = NA,
        ConfidenceScore = NA,
        Status = "Failed"
      ))
      next
    }
    
    # 提取各项得分
    consistency_score <- NA
    if (!is.null(result$consistency$overall_consistency)) {
      consistency_score <- result$consistency$overall_consistency
    }
    
    biological_score <- NA
    if (!is.null(result$biological_validity$cluster_annotation_alignment)) {
      biological_score <- result$biological_validity$cluster_annotation_alignment$mean_cluster_purity
    }
    
    marker_score <- NA
    if (!is.null(result$biological_validity$marker_validation_summary)) {
      marker_fc <- result$biological_validity$marker_validation_summary$avg_log2fc
      marker_score <- min(1, max(0, marker_fc / 2))
    }
    
    robustness_score <- NA
    if (!is.null(result$robustness$sampling_stability)) {
      stab_scores <- sapply(result$robustness$sampling_stability, 
                           function(x) ifelse(is.null(x$mean_ari) || is.na(x$mean_ari), 
                                            NA, x$mean_ari))
      stab_scores <- stab_scores[!is.na(stab_scores)]
      if (length(stab_scores) > 0) {
        robustness_score <- mean(stab_scores)
      }
    }
    
    confidence_score <- NA
    conf_scores <- c()
    for (field in names(result$confidence)) {
      if (grepl("_scores$", field)) {
        mean_score <- result$confidence[[field]]$mean_max_score
        if (!is.na(mean_score)) {
          conf_scores <- c(conf_scores, mean_score)
        }
      }
    }
    if (length(conf_scores) > 0) {
      confidence_score <- mean(conf_scores)
    }
    
    # 计算总分
    scores <- c(consistency_score, biological_score, marker_score, 
                robustness_score, confidence_score)
    weights <- c(0.25, 0.25, 0.20, 0.15, 0.15)
    
    valid_scores <- !is.na(scores)
    if (sum(valid_scores) > 0) {
      overall_score <- sum(scores[valid_scores] * weights[valid_scores]) / 
                      sum(weights[valid_scores])
    } else {
      overall_score <- NA
    }
    
    summary_table <- rbind(summary_table, data.frame(
      CancerType = cancer_type,
      OverallScore = round(overall_score, 4),
      ConsistencyScore = round(consistency_score, 4),
      BiologicalScore = round(biological_score, 4),
      MarkerScore = round(marker_score, 4),
      RobustnessScore = round(robustness_score, 4),
      ConfidenceScore = round(confidence_score, 4),
      Status = "Success"
    ))
  }
  
  # 按总分排序
  summary_table <- summary_table[order(-summary_table$OverallScore, na.last = TRUE), ]
  
  # 打印汇总表
  cat(sprintf("%-15s %12s %12s %12s %12s %12s %12s %10s\n",
             "Cancer Type", "Overall", "Consistency", "Biological", 
             "Marker", "Robustness", "Confidence", "Status"))
  cat(paste(rep("-", 110), collapse = ""), "\n")
  
  for (i in 1:nrow(summary_table)) {
    row <- summary_table[i, ]
    cat(sprintf("%-15s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %10s\n",
               row$CancerType,
               ifelse(is.na(row$OverallScore), 0, row$OverallScore),
               ifelse(is.na(row$ConsistencyScore), 0, row$ConsistencyScore),
               ifelse(is.na(row$BiologicalScore), 0, row$BiologicalScore),
               ifelse(is.na(row$MarkerScore), 0, row$MarkerScore),
               ifelse(is.na(row$RobustnessScore), 0, row$RobustnessScore),
               ifelse(is.na(row$ConfidenceScore), 0, row$ConfidenceScore),
               row$Status))
  }
  
  cat("\n")
  
  # 统计信息
  cat("================================================================================\n")
  cat("统计信息\n")
  cat("================================================================================\n\n")
  
  valid_scores <- summary_table$OverallScore[!is.na(summary_table$OverallScore)]
  
  if (length(valid_scores) > 0) {
    cat("总体得分统计:\n")
    cat("  平均值: ", round(mean(valid_scores), 4), "\n")
    cat("  中位数: ", round(median(valid_scores), 4), "\n")
    cat("  最小值: ", round(min(valid_scores), 4), "\n")
    cat("  最大值: ", round(max(valid_scores), 4), "\n")
    cat("  标准差: ", round(sd(valid_scores), 4), "\n\n")
    
    # 质量等级分布
    excellent <- sum(valid_scores >= 0.85)
    good <- sum(valid_scores >= 0.70 & valid_scores < 0.85)
    moderate <- sum(valid_scores >= 0.55 & valid_scores < 0.70)
    poor <- sum(valid_scores < 0.55)
    
    cat("质量等级分布:\n")
    cat("  优秀 (≥0.85): ", excellent, " (", 
        round(excellent / length(valid_scores) * 100, 1), "%)\n")
    cat("  良好 (0.70-0.84): ", good, " (", 
        round(good / length(valid_scores) * 100, 1), "%)\n")
    cat("  中等 (0.55-0.69): ", moderate, " (", 
        round(moderate / length(valid_scores) * 100, 1), "%)\n")
    cat("  较差 (<0.55): ", poor, " (", 
        round(poor / length(valid_scores) * 100, 1), "%)\n\n")
  }
  
  cat("================================================================================\n")
  cat("报告结束\n")
  cat("================================================================================\n")
  
  sink()
  
  # 保存汇总表为CSV
  csv_file <- file.path(output_dir, "tables", "summary_scores.csv")
  write.csv(summary_table, csv_file, row.names = FALSE)
  
  message("✓ 汇总报告已生成: ", summary_file)
  message("✓ 汇总表已保存: ", csv_file)
  
  # 生成汇总可视化
  generate_summary_visualizations(summary_table, output_dir)
  
  return(summary_table)
}

# ============================================================================
# 9. 汇总可视化
# ============================================================================

generate_summary_visualizations <- function(summary_table, output_dir) {
  
  message("生成汇总可视化...")
  
  plot_dir <- file.path(output_dir, "figures", "summary")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  publication_theme <- theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black", size = 1)
    )
  
  # 只保留成功的癌种
  summary_table <- summary_table[summary_table$Status == "Success", ]
  
  if (nrow(summary_table) == 0) {
    message("没有成功的评估结果，跳过汇总可视化")
    return(invisible(NULL))
  }
  
  # 图1: 总体得分排名
  valid_overall <- summary_table[!is.na(summary_table$OverallScore), ]
  
  if (nrow(valid_overall) > 0) {
    
    valid_overall$Grade <- cut(valid_overall$OverallScore,
                               breaks = c(0, 0.55, 0.70, 0.85, 1),
                               labels = c("Poor", "Moderate", "Good", "Excellent"))
    
    p_overall <- ggplot(valid_overall, 
                       aes(x = reorder(CancerType, OverallScore), 
                           y = OverallScore)) +
      geom_bar(stat = "identity", aes(fill = Grade), width = 0.7) +
      geom_text(aes(label = round(OverallScore, 3)), 
               hjust = -0.2, size = 3.5, fontface = "bold") +
      geom_hline(yintercept = c(0.55, 0.70, 0.85), 
                linetype = "dashed", alpha = 0.6) +
      scale_fill_manual(values = c("Poor" = "#D73027", 
                                   "Moderate" = "#FC8D59",
                                   "Good" = "#91CF60",
                                   "Excellent" = "#1A9850")) +
      coord_flip() +
      ylim(0, 1.15) +
      labs(
        title = "Overall Quality Score Ranking Across Cancer Types",
        subtitle = paste("Mean Score:", round(mean(valid_overall$OverallScore), 3)),
        x = "Cancer Type",
        y = "Overall Quality Score",
        fill = "Quality Grade"
      ) +
      publication_theme +
      theme(
        legend.position = "right"
      )
    
    ggsave(file.path(plot_dir, "01_Overall_Score_Ranking.png"), p_overall,
           width = 14, height = max(10, nrow(valid_overall) * 0.4), 
           dpi = 300, bg = "white")
    
    message("✓ 总体得分排名图已保存")
  }
  
  # 图2: 各维度得分热图
  score_cols <- c("ConsistencyScore", "BiologicalScore", "MarkerScore", 
                  "RobustnessScore", "ConfidenceScore")
  
  heatmap_data <- summary_table[, c("CancerType", score_cols)]
  heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]
  
  if (nrow(heatmap_data) > 0) {
    
    # 转换为长格式
    heatmap_long <- pivot_longer(heatmap_data, 
                                 cols = all_of(score_cols),
                                 names_to = "Metric",
                                 values_to = "Score")
    
    # 清理指标名称
    heatmap_long$Metric <- gsub("Score$", "", heatmap_long$Metric)
    
    p_heatmap <- ggplot(heatmap_long, aes(x = Metric, y = CancerType, fill = Score)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = round(Score, 2)), color = "black", 
               size = 3, fontface = "bold") +
      scale_fill_gradient2(low = "#D73027", mid = "#FEE090", high = "#1A9850",
                          midpoint = 0.7, limits = c(0, 1)) +
      labs(
        title = "Quality Metrics Heatmap Across Cancer Types",
        x = "Quality Metric",
        y = "Cancer Type",
        fill = "Score"
      ) +
      publication_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 10),
        legend.position = "right"
      )
    
    ggsave(file.path(plot_dir, "02_Quality_Metrics_Heatmap.png"), p_heatmap,
           width = 12, height = max(8, nrow(heatmap_data) * 0.3), 
           dpi = 300, bg = "white")
    
    message("✓ 质量指标热图已保存")
  }
  
  # 图3: 得分分布箱线图
  if (nrow(heatmap_long) > 0) {
    
    p_boxplot <- ggplot(heatmap_long, aes(x = Metric, y = Score, fill = Metric)) +
      geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
      geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
      geom_hline(yintercept = 0.7, linetype = "dashed", 
                color = "#1A9850", size = 1) +
      scale_fill_brewer(palette = "Set2") +
      labs(
        title = "Distribution of Quality Scores by Metric",
        subtitle = "Dashed line indicates good quality threshold (0.7)",
        x = "Quality Metric",
        y = "Score"
      ) +
      publication_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(file.path(plot_dir, "03_Score_Distribution_Boxplot.png"), p_boxplot,
           width = 12, height = 8, dpi = 300, bg = "white")
    
    message("✓ 得分分布箱线图已保存")
  }
  
  # 图4: 质量等级饼图
  if (nrow(valid_overall) > 0) {
    
    grade_counts <- as.data.frame(table(valid_overall$Grade))
    names(grade_counts) <- c("Grade", "Count")
    grade_counts$Percentage <- grade_counts$Count / sum(grade_counts$Count) * 100
    
    p_pie <- ggplot(grade_counts, aes(x = "", y = Count, fill = Grade)) +
      geom_bar(stat = "identity", width = 1, color = "white", size = 1.5) +
      coord_polar("y", start = 0) +
      geom_text(aes(label = paste0(Count, "\n(", round(Percentage, 1), "%)")),
               position = position_stack(vjust = 0.5),
               size = 5, fontface = "bold", color = "white") +
      scale_fill_manual(values = c("Poor" = "#D73027", 
                                   "Moderate" = "#FC8D59",
                                   "Good" = "#91CF60",
                                   "Excellent" = "#1A9850")) +
      labs(
        title = "Quality Grade Distribution",
        subtitle = paste("Total Cancer Types:", nrow(valid_overall)),
        fill = "Quality Grade"
      ) +
      theme_void(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    ggsave(file.path(plot_dir, "04_Quality_Grade_Distribution.png"), p_pie,
           width = 10, height = 10, dpi = 300, bg = "white")
    
    message("✓ 质量等级分布图已保存")
  }
  
  message("✓ 所有汇总可视化完成")
  
  return(invisible(NULL))
}

# ============================================================================
# 10. 使用示例
# ============================================================================

# 主执行代码
if (interactive()) {
  
  # 设置路径
  annotation_dir <- "/mnt/public7/pancancercol/hefeng/annotation-quality-batch/data"
  output_dir <- "/mnt/public7/pancancercol/hefeng/annotation-quality-batch"
  
  # 批量运行评估
  results <- run_batch_quality_assessment(
    annotation_dir = annotation_dir,
    output_dir = output_dir,
    max_cells = 50000,        # 每个癌种最多采样50000个细胞
    sample_fraction = 1.0,     # 使用全部细胞（如果少于max_cells）
    n_cores = 16               # 使用16个CPU核心
  )
  
  # 或者单独运行某个癌种
  # single_result <- run_comprehensive_quality_assessment(
  #   cancer_type = "COADREAD",
  #   annotation_dir = annotation_dir,
  #   output_dir = output_dir,
  #   max_cells = 50000,
  #   sample_fraction = 1.0,
  #   n_cores = 16
  # )
  
  message("\n全部分析完成!")
  message("结果保存在: ", output_dir)
  message("  - 报告: ", file.path(output_dir, "reports"))
  message("  - 图片: ", file.path(output_dir, "figures"))
  message("  - 表格: ", file.path(output_dir, "tables"))
}

