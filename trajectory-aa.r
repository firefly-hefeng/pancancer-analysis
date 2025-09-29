#!/usr/bin/env Rscript

# =============================================================================
# 修复后的轨迹分析质量评估代码
# =============================================================================

# 加载必要的库
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(cluster)
  library(stats)
})

# 设置全局变量
trajectory_output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/trajectory_analysis_2"
quality_output_dir <- file.path(trajectory_output_dir, "quality_assessment")

# 创建输出目录
if (!dir.exists(quality_output_dir)) {
  dir.create(quality_output_dir, recursive = TRUE)
}

# =============================================================================
# 辅助函数
# =============================================================================

safe_log <- function(message, file = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_message <- paste0("[", timestamp, "] ", message)
  cat(full_message, "\n")
  if (!is.null(file)) {
    cat(full_message, "\n", file = file, append = TRUE)
  }
}

# =============================================================================
# 1. Pseudotime质量评估
# =============================================================================

assess_pseudotime_quality <- function(cds, cancer_type, cell_type) {
  safe_log(paste("评估Pseudotime质量:", cancer_type, "-", cell_type))
  
  # 检查pseudotime列是否存在
  if(!"pseudotime" %in% colnames(colData(cds))) {
    return(list(
      status = "failed",
      reason = "No pseudotime found in CDS object"
    ))
  }
  
  pseudotime_values <- colData(cds)$pseudotime
  
  # 检查有效值
  valid_pt <- !is.na(pseudotime_values) & is.finite(pseudotime_values)
  
  if(sum(valid_pt) == 0) {
    return(list(
      status = "failed",
      reason = "No valid pseudotime values"
    ))
  }
  
  # 基本统计
  pt_stats <- list(
    n_cells = length(pseudotime_values),
    n_valid = sum(valid_pt),
    valid_rate = sum(valid_pt) / length(pseudotime_values),
    range = range(pseudotime_values, na.rm = TRUE),
    mean = mean(pseudotime_values, na.rm = TRUE),
    median = median(pseudotime_values, na.rm = TRUE),
    sd = sd(pseudotime_values, na.rm = TRUE)
  )
  
  # 分布分析
  valid_values <- pseudotime_values[valid_pt]
  
  # 均匀性检验
  uniformity_test <- tryCatch({
    ks.test(valid_values, "punif", min(valid_values), max(valid_values))
  }, error = function(e) {
    list(p.value = NA, statistic = NA)
  })
  
  # 单调性检验 - 与主成分的相关性
  monotonicity <- NA
  if("UMAP" %in% names(reducedDims(cds))) {
    umap_coords <- reducedDims(cds)$UMAP
    pc1_cor <- tryCatch({
      abs(cor(valid_values, umap_coords[valid_pt, 1], method = "spearman"))
    }, error = function(e) NA)
    
    pc2_cor <- tryCatch({
      abs(cor(valid_values, umap_coords[valid_pt, 2], method = "spearman"))
    }, error = function(e) NA)
    
    monotonicity <- max(c(pc1_cor, pc2_cor), na.rm = TRUE)
  }
  
  # 边界细胞检查
  boundary_threshold <- 0.05
  min_boundary <- quantile(valid_values, boundary_threshold)
  max_boundary <- quantile(valid_values, 1 - boundary_threshold)
  
  boundary_check <- list(
    min_boundary_cells = sum(valid_values <= min_boundary),
    max_boundary_cells = sum(valid_values >= max_boundary),
    boundary_proportion = (sum(valid_values <= min_boundary) + sum(valid_values >= max_boundary)) / length(valid_values)
  )
  
  return(list(
    status = "success",
    basic_stats = pt_stats,
    uniformity_test = uniformity_test,
    monotonicity = monotonicity,
    boundary_check = boundary_check
  ))
}

# =============================================================================
# 2. 修复后的生物学一致性验证
# =============================================================================

validate_biological_consistency <- function(cds, original_seurat, cancer_type, cell_type) {
  safe_log(paste("验证生物学一致性:", cancer_type, "-", cell_type))
  
  # 获取共同细胞
  common_cells <- intersect(colnames(cds), colnames(original_seurat))
  
  if(length(common_cells) == 0) {
    return(list(
      status = "failed",
      reason = "No common cells between CDS and Seurat objects"
    ))
  }
  
  # 根据细胞类型选择标记基因
  if(grepl("B_cell|B.cell|bcell", cell_type, ignore.case = TRUE)) {
    # B细胞标记基因
    marker_genes <- list(
      naive_b = c("CD19", "CD20", "MS4A1", "TCL1A", "FCER2", "IGHD", "IGHM"),
      memory_b = c("CD27", "CD38", "TNFRSF13B", "TNFRSF17", "AICDA"),
      plasma = c("CD138", "SDC1", "XBP1", "IRF4", "PRDM1", "MZB1"),
      germinal_center = c("BCL6", "AICDA", "MEF2B", "LMO2", "FOXO1"),
      activation = c("CD69", "CD25", "IL2RA", "ICOS", "CD44", "CD86"),
      proliferation = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "STMN1")
    )
  } else if(grepl("T_cell|T.cell|tcell", cell_type, ignore.case = TRUE)) {
    # T细胞标记基因
    marker_genes <- list(
      naive_memory = c("TCF7", "LEF1", "CCR7", "SELL", "IL7R"),
      activation = c("CD69", "CD25", "IL2RA", "ICOS", "CD44"),
      effector = c("GZMA", "GZMB", "PRF1", "IFNG", "TNF"),
      exhaustion = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "EOMES"),
      proliferation = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1"),
      memory = c("IL7R", "CCR7", "TCF7", "SELL", "CD62L")
    )
  } else if(grepl("NK|nk", cell_type, ignore.case = TRUE)) {
    # NK细胞标记基因
    marker_genes <- list(
      cytotoxic = c("GZMA", "GZMB", "PRF1", "GNLY", "NKG7"),
      activation = c("CD69", "IFNG", "TNF", "FASLG"),
      receptors = c("KLRB1", "KLRD1", "KLRF1", "NCR1", "NCR3"),
      proliferation = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1"),
      maturation = c("CD16", "FCGR3A", "CD57", "NCAM1"),
      inhibitory = c("KLRC1", "KLRC2", "LILRB1", "KIR2DL1")
    )
  } else if(grepl("myeloid|monocyte|macrophage|DC", cell_type, ignore.case = TRUE)) {
    # 髓系细胞标记基因
    marker_genes <- list(
      monocyte = c("CD14", "CD16", "FCGR3A", "CCR2", "CX3CR1"),
      macrophage = c("CD68", "CD163", "MRC1", "ARG1", "IL10"),
      dendritic = c("CD1C", "CLEC9A", "XCR1", "BATF3", "IRF8"),
      activation = c("CD69", "CD80", "CD86", "IL1B", "TNF"),
      proliferation = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1"),
      m1_polarization = c("NOS2", "IL1B", "IL6", "TNF", "CXCL10")
    )
  } else {
    # 通用标记基因（适用于未知细胞类型）
    marker_genes <- list(
      activation = c("CD69", "CD25", "IL2RA", "ICOS", "CD44"),
      proliferation = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1"),
      apoptosis = c("BAX", "BCL2", "TP53", "CASP3", "CASP8"),
      stress = c("HSP90AA1", "HSPA1A", "DNAJB1", "HSPH1"),
      metabolism = c("ENO1", "PKM", "LDHA", "GAPDH", "PGK1"),
      ribosomal = c("RPS6", "RPL10", "RPS27", "RPL32", "RPS29")
    )
  }
  
  # 获取表达数据
  expr_data <- tryCatch({
    GetAssayData(original_seurat[, common_cells], layer = "data")
  }, error = function(e) {
    tryCatch({
      GetAssayData(original_seurat[, common_cells], slot = "data")
    }, error = function(e2) {
      return(NULL)
    })
  })
  
  if(is.null(expr_data)) {
    return(list(
      status = "failed",
      reason = "Cannot access expression data from Seurat object"
    ))
  }
  
  # 计算基因集评分和相关性
  correlations <- list()
  gene_set_scores <- list()
  
  pseudotime_values <- colData(cds[, common_cells])$pseudotime
  
  for(set_name in names(marker_genes)) {
    genes <- marker_genes[[set_name]]
    available_genes <- intersect(genes, rownames(expr_data))
    
    if(length(available_genes) > 0) {
      # 计算基因集评分
      if(length(available_genes) == 1) {
        gene_scores <- as.numeric(expr_data[available_genes, ])
      } else {
        gene_scores <- colMeans(as.matrix(expr_data[available_genes, , drop = FALSE]))
      }
      
      gene_set_scores[[set_name]] <- gene_scores
      
      # 计算与pseudotime的相关性
      if(length(gene_scores) > 0 && !all(is.na(gene_scores)) && 
         length(pseudotime_values) > 0 && !all(is.na(pseudotime_values)) &&
         length(unique(gene_scores[!is.na(gene_scores)])) > 1) {
        
        cor_result <- tryCatch({
          cor(pseudotime_values, gene_scores, method = "spearman", use = "complete.obs")
        }, error = function(e) {
          NA
        })
        
        correlations[[set_name]] <- cor_result
      } else {
        correlations[[set_name]] <- NA
      }
    } else {
      correlations[[set_name]] <- NA
      gene_set_scores[[set_name]] <- NULL
    }
  }
  
  # 计算模式一致性
  valid_correlations <- correlations[!is.na(unlist(correlations))]
  
  if(length(valid_correlations) > 0) {
    # 检查生物学合理的模式
    if(grepl("B_cell|B.cell|bcell", cell_type, ignore.case = TRUE)) {
      # B细胞期望模式
      pattern_consistency <- c(
        naive_b_negative = if(!is.na(correlations$naive_b)) correlations$naive_b < -0.1 else NA,
        plasma_positive = if(!is.na(correlations$plasma)) correlations$plasma > 0.1 else NA,
        proliferation_detectable = if(!is.na(correlations$proliferation)) abs(correlations$proliferation) > 0.05 else NA
      )
    } else if(grepl("T_cell|T.cell|tcell", cell_type, ignore.case = TRUE)) {
      # T细胞期望模式
      pattern_consistency <- c(
        naive_memory_negative = if(!is.na(correlations$naive_memory)) correlations$naive_memory < 0 else NA,
        effector_positive = if(!is.na(correlations$effector)) correlations$effector > 0 else NA,
        activation_positive = if(!is.na(correlations$activation)) correlations$activation > 0 else NA
      )
    } else {
      # 其他细胞类型的通用模式
      pattern_consistency <- c(
        activation_detectable = if(!is.na(correlations$activation)) abs(correlations$activation) > 0.05 else NA,
        proliferation_detectable = if(!is.na(correlations$proliferation)) abs(correlations$proliferation) > 0.05 else NA
      )
    }
    
    consistency_rate <- mean(pattern_consistency, na.rm = TRUE)
    if(is.nan(consistency_rate)) consistency_rate <- 0
  } else {
    pattern_consistency <- rep(NA, length(correlations))
    names(pattern_consistency) <- names(correlations)
    consistency_rate <- 0
  }
  
  return(list(
    status = "success",
    correlations = correlations,
    gene_set_scores = gene_set_scores,
    pattern_consistency = pattern_consistency,
    consistency_rate = consistency_rate,
    available_gene_sets = length(valid_correlations),
    common_cells_used = length(common_cells)
  ))
}

# =============================================================================
# 3. 聚类质量评估
# =============================================================================

assess_clustering_quality <- function(cds, cancer_type, cell_type) {
  safe_log(paste("评估聚类质量:", cancer_type, "-", cell_type))
  
  # 检查聚类信息
  if(!"manual_clusters" %in% colnames(colData(cds))) {
    return(list(
      status = "failed",
      reason = "No manual_clusters found in CDS object"
    ))
  }
  
  clusters <- colData(cds)$manual_clusters
  
  # 检查UMAP坐标
  if(!"UMAP" %in% names(reducedDims(cds))) {
    return(list(
      status = "failed",
      reason = "No UMAP coordinates found in CDS object"
    ))
  }
  
  umap_coords <- reducedDims(cds)$UMAP
  
  # 基本聚类统计
  n_clusters <- length(unique(clusters))
  cluster_sizes <- table(clusters)
  
  # 轮廓系数分析
  if(n_clusters > 1 && n_clusters < nrow(umap_coords)) {
    dist_matrix <- dist(umap_coords)
    numeric_clusters <- as.numeric(as.factor(clusters))
    
    sil_result <- tryCatch({
      silhouette(numeric_clusters, dist_matrix)
    }, error = function(e) {
      return(NULL)
    })
    
    if(!is.null(sil_result)) {
      avg_silhouette <- mean(sil_result[, 3])
      silhouette_by_cluster <- aggregate(sil_result[, 3], 
                                       by = list(sil_result[, 1]), 
                                       FUN = mean)
    } else {
      avg_silhouette <- NA
      silhouette_by_cluster <- NA
    }
  } else {
    avg_silhouette <- NA
    silhouette_by_cluster <- NA
  }
  
  # 聚类中心距离分析
  cluster_centers <- aggregate(umap_coords, by = list(clusters), FUN = mean)
  if(nrow(cluster_centers) > 1) {
    center_distances <- dist(cluster_centers[, -1])
    min_center_distance <- min(center_distances)
    avg_center_distance <- mean(center_distances)
  } else {
    min_center_distance <- NA
    avg_center_distance <- NA
  }
  
  # 聚类紧密度分析
  cluster_compactness <- numeric(n_clusters)
  for(i in seq_along(unique(clusters))) {
    cluster_id <- unique(clusters)[i]
    cluster_cells <- clusters == cluster_id
    if(sum(cluster_cells) > 1) {
      cluster_coords <- umap_coords[cluster_cells, ]
      cluster_center <- colMeans(cluster_coords)
      cluster_compactness[i] <- mean(sqrt(rowSums((cluster_coords - matrix(cluster_center, 
                                                   nrow = nrow(cluster_coords), 
                                                   ncol = 2, byrow = TRUE))^2)))
    } else {
      cluster_compactness[i] <- 0
    }
  }
  
  return(list(
    status = "success",
    n_clusters = n_clusters,
    cluster_sizes = cluster_sizes,
    avg_silhouette = avg_silhouette,
    silhouette_by_cluster = silhouette_by_cluster,
    min_center_distance = min_center_distance,
    avg_center_distance = avg_center_distance,
    cluster_compactness = cluster_compactness,
    avg_compactness = mean(cluster_compactness)
  ))
}

# =============================================================================
# 4. 基因质量评估
# =============================================================================

assess_gene_quality <- function(trajectory_genes, cds, cancer_type, cell_type) {
  safe_log(paste("评估基因特征质量:", cancer_type, "-", cell_type))
  
  if(is.null(trajectory_genes) || nrow(trajectory_genes) == 0) {
    return(list(
      status = "failed",
      reason = "No trajectory genes provided or empty data frame"
    ))
  }
  
  # 基本基因统计
  total_genes <- nrow(trajectory_genes)
  
  # 检查必需列
  required_cols <- c("correlation", "p_value", "q_value")
  missing_cols <- setdiff(required_cols, colnames(trajectory_genes))
  
  if(length(missing_cols) > 0) {
    return(list(
      status = "failed",
      reason = paste("Missing required columns:", paste(missing_cols, collapse = ", "))
    ))
  }
  
  # 显著性分析
  significant_genes <- sum(trajectory_genes$q_value < 0.05, na.rm = TRUE)
  highly_significant <- sum(trajectory_genes$q_value < 0.01, na.rm = TRUE)
  significance_rate <- significant_genes / total_genes
  
  # 相关性分析
  correlations <- trajectory_genes$correlation
  valid_correlations <- correlations[!is.na(correlations)]
  
  correlation_stats <- list(
    mean_correlation = mean(abs(valid_correlations)),
    median_correlation = median(abs(valid_correlations)),
    max_correlation = max(abs(valid_correlations)),
    positive_genes = sum(valid_correlations > 0),
    negative_genes = sum(valid_correlations < 0),
    strong_correlation = sum(abs(valid_correlations) > 0.3)
  )
  
  # p值分布分析
  p_values <- trajectory_genes$p_value
  valid_p_values <- p_values[!is.na(p_values) & p_values > 0]
  
  if(length(valid_p_values) > 0) {
    log_p_values <- -log10(valid_p_values)
    p_value_stats <- list(
      mean_log_pvalue = mean(log_p_values),
      median_log_pvalue = median(log_p_values),
      very_significant = sum(log_p_values > 10)  # p < 1e-10
    )
  } else {
    p_value_stats <- list(
      mean_log_pvalue = NA,
      median_log_pvalue = NA,
      very_significant = 0
    )
  }
  
  # FDR控制质量
  q_values <- trajectory_genes$q_value
  valid_q_values <- q_values[!is.na(q_values)]
  
  if(length(valid_q_values) > 0 && length(valid_p_values) > 0) {
    # 计算FDR校正的质量
    fdr_ratio <- length(valid_q_values[valid_q_values < 0.05]) / 
                 length(valid_p_values[valid_p_values < 0.05])
    avg_fdr_inflation <- mean(valid_q_values / valid_p_values[seq_along(valid_q_values)], na.rm = TRUE)
  } else {
    fdr_ratio <- NA
    avg_fdr_inflation <- NA
  }
  
  return(list(
    status = "success",
    gene_statistics = list(
      total_genes = total_genes,
      significant_genes = significant_genes,
      highly_significant = highly_significant,
      significance_rate = significance_rate
    ),
    correlation_statistics = correlation_stats,
    pvalue_statistics = p_value_stats,
    fdr_control = list(
      fdr_ratio = fdr_ratio,
      avg_fdr_inflation = avg_fdr_inflation
    ),
    expression_quality = NULL  # 可以后续添加表达质量评估
  ))
}

# =============================================================================
# 5. 采样质量评估
# =============================================================================

assess_sampling_representativeness <- function(original_seurat, sampled_cds, cancer_type, cell_type) {
  safe_log(paste("评估采样代表性:", cancer_type, "-", cell_type))
  
  # 基本采样统计
  n_original <- ncol(original_seurat)
  n_sampled <- ncol(sampled_cds)
  sampling_rate <- n_sampled / n_original
  
  # 检查细胞类型分布代表性
  common_cells <- intersect(colnames(sampled_cds), colnames(original_seurat))
  
  if(length(common_cells) == 0) {
    return(list(
      status = "failed",
      reason = "No common cells between original and sampled datasets"
    ))
  }
  
  # 细胞类型分布比较
  type_distribution <- NULL
  if("fine_cell_type" %in% colnames(original_seurat@meta.data)) {
    original_types <- table(original_seurat@meta.data$fine_cell_type)
    sampled_types <- table(original_seurat@meta.data[common_cells, "fine_cell_type"])
    
    # 计算比例差异
    original_props <- original_types / sum(original_types)
    sampled_props <- sampled_types / sum(sampled_types)
    
    # 确保所有类型都被考虑
    all_types <- union(names(original_props), names(sampled_props))
    orig_props_full <- numeric(length(all_types))
    samp_props_full <- numeric(length(all_types))
    names(orig_props_full) <- all_types
    names(samp_props_full) <- all_types
    
    orig_props_full[names(original_props)] <- original_props
    samp_props_full[names(sampled_props)] <- sampled_props
    
    proportion_differences <- abs(orig_props_full - samp_props_full)
    
    # 卡方检验
    chi_square_test <- tryCatch({
      chisq.test(rbind(original_types, sampled_types))
    }, error = function(e) {
      list(p.value = NA)
    })
    
    type_distribution <- list(
      max_proportion_difference = max(proportion_differences),
      avg_proportion_difference = mean(proportion_differences),
      chi_square_pvalue = chi_square_test$p.value,
      well_represented = max(proportion_differences) < 0.1
    )
  }
  
  # 表达谱相似性评估
  expression_similarity <- NULL
  if(length(common_cells) > 50) {  # 确保有足够的细胞进行比较
    tryCatch({
      # 获取高变基因
      top_genes <- head(VariableFeatures(original_seurat), 2000)
      available_genes <- intersect(top_genes, rownames(original_seurat))
      
      if(length(available_genes) > 100) {
        # 计算原始数据的PCA
        original_expr <- GetAssayData(original_seurat, layer = "data")[available_genes, ]
        sampled_expr <- GetAssayData(original_seurat, layer = "data")[available_genes, common_cells]
        
        # PCA分析
        original_pca <- prcomp(t(as.matrix(original_expr)), center = TRUE, scale. = TRUE)
        sampled_pca <- prcomp(t(as.matrix(sampled_expr)), center = TRUE, scale. = TRUE)
        
        # 比较前10个主成分的方差解释比例
        n_pcs <- min(10, ncol(original_pca$x), ncol(sampled_pca$x))
        original_var <- original_pca$sdev[1:n_pcs]^2 / sum(original_pca$sdev^2)
        sampled_var <- sampled_pca$sdev[1:n_pcs]^2 / sum(sampled_pca$sdev^2)
        
        pca_coverage <- cor(original_var, sampled_var, method = "spearman")
        
        expression_similarity <- list(
          pca_coverage = pca_coverage,
          good_coverage = pca_coverage > 0.8
        )
      }
    }, error = function(e) {
      expression_similarity <- list(
        pca_coverage = NA,
        good_coverage = FALSE
      )
    })
  }
  
  return(list(
    status = "success",
    sampling_rate = sampling_rate,
    n_original = n_original,
    n_sampled = n_sampled,
    type_distribution = type_distribution,
    expression_similarity = expression_similarity
  ))
}

# =============================================================================
# 6. 综合质量评估
# =============================================================================

comprehensive_quality_assessment <- function(cancer_type, cell_type) {
  safe_log(paste("开始综合质量评估:", cancer_type, "-", cell_type))
  
  # 加载数据
  sample_name <- paste(cancer_type, cell_type, sep = "_")
  results_file <- file.path(trajectory_output_dir, "monocle3_results", 
                           paste0(sample_name, "_trajectory_results.rds"))
  
  if(!file.exists(results_file)) {
    stop("结果文件不存在: ", results_file)
  }
  
  results <- readRDS(results_file)
  cds <- results$cds
  seurat_obj <- results$seurat_obj
  trajectory_genes <- results$trajectory_genes
  
  # 执行各项质量评估
  assessments <- list()
  
  # 1. Pseudotime质量
  assessments$pseudotime <- assess_pseudotime_quality(cds, cancer_type, cell_type)
  
  # 2. 生物学一致性
  assessments$biological <- validate_biological_consistency(cds, seurat_obj, cancer_type, cell_type)
  
  # 3. 聚类质量
  assessments$clustering <- assess_clustering_quality(cds, cancer_type, cell_type)
  
  # 4. 基因质量
  assessments$genes <- assess_gene_quality(trajectory_genes, cds, cancer_type, cell_type)
  
  # 5. 采样代表性
  assessments$sampling <- assess_sampling_representativeness(seurat_obj, cds, cancer_type, cell_type)
  
  # 计算综合评分
  scores <- numeric(5)
  names(scores) <- c("pseudotime", "biological", "clustering", "genes", "sampling")
  
  # Pseudotime评分
  if(assessments$pseudotime$status == "success") {
    pt_score <- 0
    if(!is.na(assessments$pseudotime$monotonicity)) {
      pt_score <- pt_score + assessments$pseudotime$monotonicity * 40
    }
    if(!is.na(assessments$pseudotime$uniformity_test$p.value)) {
      # p值越小（分布越不均匀）可能越好，但这里简单处理
      pt_score <- pt_score + 30
    }
    pt_score <- pt_score + assessments$pseudotime$basic_stats$valid_rate * 30
    scores["pseudotime"] <- min(pt_score, 100)
  } else {
    scores["pseudotime"] <- 0
  }
  
  # 生物学一致性评分
  if(assessments$biological$status == "success") {
    bio_score <- assessments$biological$consistency_rate * 100
    if(is.nan(bio_score)) bio_score <- 50  # 给予中等分数如果无法计算
    scores["biological"] <- bio_score
  } else {
    scores["biological"] <- 0
  }
  
  # 聚类质量评分
  if(assessments$clustering$status == "success") {
    cluster_score <- 0
    if(!is.na(assessments$clustering$avg_silhouette)) {
      cluster_score <- cluster_score + assessments$clustering$avg_silhouette * 70
    }
    if(assessments$clustering$n_clusters >= 2 && assessments$clustering$n_clusters <= 10) {
      cluster_score <- cluster_score + 30
    }
    scores["clustering"] <- min(cluster_score, 100)
  } else {
    scores["clustering"] <- 0
  }
  
  # 基因质量评分
  if(assessments$genes$status == "success") {
    gene_score <- assessments$genes$gene_statistics$significance_rate * 50 +
                  min(assessments$genes$correlation_statistics$mean_correlation * 100, 50)
    scores["genes"] <- gene_score
  } else {
    scores["genes"] <- 0
  }
  
  # 采样代表性评分
  if(assessments$sampling$status == "success") {
    sampling_score <- 50  # 基础分
    if(!is.null(assessments$sampling$type_distribution) && 
       !is.null(assessments$sampling$type_distribution$well_represented)) {
      if(assessments$sampling$type_distribution$well_represented) {
        sampling_score <- sampling_score + 25
      }
    }
    if(!is.null(assessments$sampling$expression_similarity) && 
       !is.null(assessments$sampling$expression_similarity$good_coverage)) {
      if(assessments$sampling$expression_similarity$good_coverage) {
        sampling_score <- sampling_score + 25
      }
    }
    scores["sampling"] <- sampling_score
  } else {
    scores["sampling"] <- 0
  }
  
  # 计算加权综合评分
  weights <- c(0.25, 0.2, 0.2, 0.2, 0.15)  # pseudotime, biological, clustering, genes, sampling
  overall_score <- sum(scores * weights)
  
  # 质量等级
  if(overall_score >= 80) {
    quality_grade <- "Excellent"
  } else if(overall_score >= 65) {
    quality_grade <- "Good"
  } else if(overall_score >= 50) {
    quality_grade <- "Fair"
  } else {
    quality_grade <- "Poor"
  }
  
  # 编译最终结果
  final_result <- list(
    sample_info = list(
      cancer_type = cancer_type,
      cell_type = cell_type,
      sample_name = sample_name,
      assessment_time = Sys.time()
    ),
    individual_assessments = assessments,
    scores = scores,
    overall_score = overall_score,
    quality_grade = quality_grade,
    recommendations = generate_recommendations(assessments, scores)
  )
  
  # 保存结果
  output_file <- file.path(quality_output_dir, paste0(sample_name, "_quality_assessment.rds"))
  saveRDS(final_result, output_file)
  
  safe_log(paste("质量评估完成，综合评分:", round(overall_score, 2), "等级:", quality_grade))
  safe_log(paste("结果已保存到:", output_file))
  
  return(final_result)
}

# =============================================================================
# 7. 生成改进建议
# =============================================================================

generate_recommendations <- function(assessments, scores) {
  recommendations <- list()
  
  # Pseudotime相关建议
  if(scores["pseudotime"] < 70) {
    recommendations <- append(recommendations, list(
      category = "pseudotime",
      issue = "Pseudotime质量较低",
      suggestions = c(
        "检查轨迹推断参数设置",
        "考虑增加更多细胞或基因",
        "验证轨迹的生物学合理性",
        "尝试不同的降维方法"
      )
    ))
  }
  
  # 生物学一致性建议
  if(scores["biological"] < 60) {
    recommendations <- append(recommendations, list(
      category = "biological",
      issue = "生物学一致性较差",
      suggestions = c(
        "检查标记基因的选择是否合适",
        "验证细胞类型注释的准确性",
        "考虑使用更具体的功能基因集",
        "检查批次效应的影响"
      )
    ))
  }
  
  # 聚类质量建议
  if(scores["clustering"] < 65) {
    recommendations <- append(recommendations, list(
      category = "clustering",
      issue = "聚类质量需要改进",
      suggestions = c(
        "调整聚类分辨率参数",
        "考虑使用不同的聚类算法",
        "检查数据预处理步骤",
        "验证聚类的生物学意义"
      )
    ))
  }
  
  # 基因质量建议
  if(scores["genes"] < 60) {
    recommendations <- append(recommendations, list(
      category = "genes",
      issue = "轨迹基因质量较低",
      suggestions = c(
        "增加统计检验的严格性",
        "使用更大的细胞群体",
        "检查基因过滤标准",
        "验证差异表达分析参数"
      )
    ))
  }
  
  # 采样代表性建议
  if(scores["sampling"] < 65) {
    recommendations <- append(recommendations, list(
      category = "sampling",
      issue = "采样代表性不足",
      suggestions = c(
        "增加采样的细胞数量",
        "确保各亚群的平衡采样",
        "检查采样偏倚",
        "考虑分层采样策略"
      )
    ))
  }
  
  return(recommendations)
}

# =============================================================================
# 主执行函数
# =============================================================================

if(!interactive()) {
  # 如果作为脚本运行，执行质量评估
  args <- commandArgs(trailingOnly = TRUE)
  
  if(length(args) >= 2) {
    cancer_type <- args[1]
    cell_type <- args[2]
    
    result <- comprehensive_quality_assessment(cancer_type, cell_type)
    cat("质量评估完成!\n")
    cat("综合评分:", result$overall_score, "\n")
    cat("质量等级:", result$quality_grade, "\n")
  } else {
    cat("用法: Rscript trajectory-quality-assessment.R <cancer_type> <cell_type>\n")
    cat("例如: Rscript trajectory-quality-assessment.R Melan B_cell\n")
  }
}