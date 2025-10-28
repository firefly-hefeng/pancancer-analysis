# trajectory_quality_assessment.R
# 轨迹分析数据质量综合评估代码

library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cluster)
library(igraph)
library(viridis)

# 设置路径
trajectory_output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/trajectory_analysis_2"
quality_output_dir <- file.path(trajectory_output_dir, "quality_assessment")
dir.create(quality_output_dir, showWarnings = FALSE, recursive = TRUE)

# ==================== 质量评估核心函数 ====================

# 1. Pseudotime质量评估
assess_pseudotime_quality <- function(cds, cancer_type, cell_type) {
  message(paste("评估Pseudotime质量:", cancer_type, "-", cell_type))
  
  if(!"pseudotime" %in% colnames(colData(cds))) {
    return(list(status = "failed", reason = "No pseudotime found"))
  }
  
  pseudotime_values <- colData(cds)$pseudotime
  valid_pt <- !is.na(pseudotime_values) & is.finite(pseudotime_values)
  
  if(sum(valid_pt) == 0) {
    return(list(status = "failed", reason = "No valid pseudotime values"))
  }
  
  # 1. 基本统计
  pt_stats <- list(
    n_cells = length(pseudotime_values),
    n_valid = sum(valid_pt),
    valid_rate = sum(valid_pt) / length(pseudotime_values),
    range = range(pseudotime_values, na.rm = TRUE),
    mean = mean(pseudotime_values, na.rm = TRUE),
    median = median(pseudotime_values, na.rm = TRUE),
    sd = sd(pseudotime_values, na.rm = TRUE)
  )
  
  # 2. 分布均匀性检查
  valid_values <- pseudotime_values[valid_pt]
  uniformity_test <- tryCatch({
    ks.test(valid_values, "punif", min(valid_values), max(valid_values))
  }, error = function(e) list(p.value = NA, statistic = NA))
  
  # 3. 单调性检查（基于UMAP距离）
  monotonicity_score <- tryCatch({
    umap_coords <- reducedDims(cds)$UMAP
    if(is.null(umap_coords)) return(NA)
    
    # 计算每个细胞到轨迹起点的距离
    start_idx <- which.min(valid_values)
    distances_to_start <- sqrt(rowSums((umap_coords - umap_coords[start_idx, ])^2))
    
    # 计算距离与pseudotime的相关性
    cor(distances_to_start, pseudotime_values, method = "spearman", use = "complete.obs")
  }, error = function(e) NA)
  
  # 4. 边界细胞检查
  boundary_check <- list(
    min_boundary_cells = sum(valid_values <= quantile(valid_values, 0.05)),
    max_boundary_cells = sum(valid_values >= quantile(valid_values, 0.95)),
    boundary_proportion = (sum(valid_values <= quantile(valid_values, 0.05)) + 
                          sum(valid_values >= quantile(valid_values, 0.95))) / length(valid_values)
  )
  
  return(list(
    status = "success",
    basic_stats = pt_stats,
    uniformity_test = uniformity_test,
    monotonicity = monotonicity_score,
    boundary_check = boundary_check
  ))
}

# 2. 生物学一致性验证
validate_biological_consistency <- function(cds, seurat_obj, cancer_type, cell_type) {
  message(paste("验证生物学一致性:", cancer_type, "-", cell_type))
  
  if(!"pseudotime" %in% colnames(colData(cds))) {
    return(list(status = "failed", reason = "No pseudotime found"))
  }
  
  pseudotime_values <- colData(cds)$pseudotime
  common_cells <- intersect(colnames(cds), colnames(seurat_obj))
  
  if(length(common_cells) < 50) {
    return(list(status = "failed", reason = "Too few common cells"))
  }
  
  # T细胞分化相关基因集
  marker_genes <- list(
    naive_memory = c("TCF7", "LEF1", "CCR7", "SELL", "IL7R"),
    activation = c("CD69", "CD25", "IL2RA", "ICOS", "CD44"),
    effector = c("GZMA", "GZMB", "PRF1", "IFNG", "TNF"),
    exhaustion = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "EOMES"),
    proliferation = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1"),
    memory = c("IL7R", "CCR7", "TCF7", "SELL", "CD62L")
  )
  
  # 计算每个基因集的评分
  correlations <- list()
  gene_set_scores <- list()
  
  for(set_name in names(marker_genes)) {
    genes <- marker_genes[[set_name]]
    available_genes <- intersect(genes, rownames(seurat_obj))
    
    if(length(available_genes) >= 2) {
      # 计算基因集评分
      expr_data <- GetAssayData(seurat_obj, layer = "data")[available_genes, common_cells, drop = FALSE]
      set_score <- colMeans(expr_data)
      
      # 与pseudotime的相关性
      pt_common <- pseudotime_values[common_cells]
      correlation <- cor(set_score, pt_common, method = "spearman", use = "complete.obs")
      
      correlations[[set_name]] <- correlation
      gene_set_scores[[set_name]] <- set_score
    } else {
      correlations[[set_name]] <- NA
    }
  }
  
  # 预期模式检查
  expected_patterns <- list(
    naive_memory = "negative",  # 预期负相关
    activation = "positive",    # 预期正相关
    effector = "positive",      # 预期正相关
    exhaustion = "positive",    # 预期正相关
    proliferation = "variable", # 可变
    memory = "negative"         # 预期负相关
  )
  
  pattern_consistency <- sapply(names(expected_patterns), function(set) {
    if(is.na(correlations[[set]])) return(NA)
    
    expected <- expected_patterns[[set]]
    observed <- correlations[[set]]
    
    if(expected == "positive") {
      return(observed > 0.1)
    } else if(expected == "negative") {
      return(observed < -0.1)
    } else {
      return(TRUE)  # variable pattern
    }
  })
  
  return(list(
    status = "success",
    correlations = correlations,
    gene_set_scores = gene_set_scores,
    pattern_consistency = pattern_consistency,
    consistency_rate = mean(pattern_consistency, na.rm = TRUE)
  ))
}

# 3. 聚类质量评估
assess_clustering_quality <- function(cds, cancer_type, cell_type) {
  message(paste("评估聚类质量:", cancer_type, "-", cell_type))
  
  if(!"manual_clusters" %in% colnames(colData(cds))) {
    return(list(status = "failed", reason = "No clusters found"))
  }
  
  clusters <- colData(cds)$manual_clusters
  umap_coords <- reducedDims(cds)$UMAP
  
  if(is.null(umap_coords) || nrow(umap_coords) == 0) {
    return(list(status = "failed", reason = "No UMAP coordinates"))
  }
  
  # 1. 聚类基本信息
  cluster_info <- table(clusters)
  n_clusters <- length(unique(clusters))
  
  # 2. 轮廓系数
  if(n_clusters > 1 && n_clusters < nrow(umap_coords)) {
    silhouette_result <- tryCatch({
      dist_matrix <- dist(umap_coords)
      silhouette(as.numeric(as.factor(clusters)), dist_matrix)
    }, error = function(e) NULL)
    
    if(!is.null(silhouette_result)) {
      avg_silhouette <- mean(silhouette_result[, 3])
      silhouette_by_cluster <- aggregate(silhouette_result[, 3], 
                                       by = list(silhouette_result[, 1]), 
                                       FUN = mean)
    } else {
      avg_silhouette <- NA
      silhouette_by_cluster <- NA
    }
  } else {
    avg_silhouette <- NA
    silhouette_by_cluster <- NA
  }
  
  # 3. 聚类间距离
  cluster_centers <- aggregate(umap_coords, by = list(clusters), FUN = mean)
  if(nrow(cluster_centers) > 1) {
    center_distances <- dist(cluster_centers[, -1])
    min_center_distance <- min(center_distances)
    avg_center_distance <- mean(center_distances)
  } else {
    min_center_distance <- NA
    avg_center_distance <- NA
  }
  
  # 4. 聚类内紧密度
  cluster_compactness <- sapply(unique(clusters), function(cl) {
    cl_cells <- which(clusters == cl)
    if(length(cl_cells) > 1) {
      cl_coords <- umap_coords[cl_cells, ]
      cl_center <- colMeans(cl_coords)
      mean(sqrt(rowSums((cl_coords - cl_center)^2)))
    } else {
      0
    }
  })
  
  return(list(
    status = "success",
    n_clusters = n_clusters,
    cluster_sizes = cluster_info,
    avg_silhouette = avg_silhouette,
    silhouette_by_cluster = silhouette_by_cluster,
    min_center_distance = min_center_distance,
    avg_center_distance = avg_center_distance,
    cluster_compactness = cluster_compactness,
    avg_compactness = mean(cluster_compactness, na.rm = TRUE)
  ))
}

# 4. 基因特征质量评估
assess_gene_quality <- function(trajectory_genes, cds, cancer_type, cell_type) {
  message(paste("评估基因特征质量:", cancer_type, "-", cell_type))
  
  if(is.null(trajectory_genes) || nrow(trajectory_genes) == 0) {
    return(list(status = "failed", reason = "No trajectory genes found"))
  }
  
  # 1. 基因数量和显著性
  gene_stats <- list(
    total_genes = nrow(trajectory_genes),
    significant_genes = sum(trajectory_genes$q_value < 0.05, na.rm = TRUE),
    highly_significant = sum(trajectory_genes$q_value < 0.01, na.rm = TRUE),
    significance_rate = sum(trajectory_genes$q_value < 0.05, na.rm = TRUE) / nrow(trajectory_genes)
  )
  
  # 2. 相关性分布
  correlation_stats <- list(
    mean_correlation = mean(abs(trajectory_genes$correlation), na.rm = TRUE),
    median_correlation = median(abs(trajectory_genes$correlation), na.rm = TRUE),
    max_correlation = max(abs(trajectory_genes$correlation), na.rm = TRUE),
    positive_genes = sum(trajectory_genes$correlation > 0, na.rm = TRUE),
    negative_genes = sum(trajectory_genes$correlation < 0, na.rm = TRUE),
    strong_correlation = sum(abs(trajectory_genes$correlation) > 0.3, na.rm = TRUE)
  )
  
  # 3. p值分布检查
  pvalue_stats <- list(
    mean_log_pvalue = mean(-log10(trajectory_genes$p_value), na.rm = TRUE),
    median_log_pvalue = median(-log10(trajectory_genes$p_value), na.rm = TRUE),
    very_significant = sum(trajectory_genes$p_value < 1e-10, na.rm = TRUE)
  )
  
  # 4. 多重检验校正效果
  if(all(c("p_value", "q_value") %in% colnames(trajectory_genes))) {
    fdr_control <- list(
      fdr_ratio = sum(trajectory_genes$q_value < 0.05, na.rm = TRUE) / 
                  sum(trajectory_genes$p_value < 0.05, na.rm = TRUE),
      avg_fdr_inflation = mean(trajectory_genes$q_value / trajectory_genes$p_value, na.rm = TRUE)
    )
  } else {
    fdr_control <- list(fdr_ratio = NA, avg_fdr_inflation = NA)
  }
  
  # 5. 基因表达质量检查（如果有CDS对象）
  expr_quality <- NULL
  if(!is.null(cds)) {
    expr_matrix <- tryCatch({
      get_cds_logcounts(cds)
    }, error = function(e) NULL)
    
    if(!is.null(expr_matrix)) {
      available_genes <- intersect(trajectory_genes$gene_short_name, rownames(expr_matrix))
      
      if(length(available_genes) > 0) {
        gene_expr <- expr_matrix[available_genes, ]
        
        expr_quality <- list(
          available_genes = length(available_genes),
          availability_rate = length(available_genes) / nrow(trajectory_genes),
          mean_expression = mean(rowMeans(gene_expr)),
          mean_detection_rate = mean(rowMeans(gene_expr > 0)),
          expression_variability = mean(apply(gene_expr, 1, cv))
        )
      }
    }
  }
  
  return(list(
    status = "success",
    gene_statistics = gene_stats,
    correlation_statistics = correlation_stats,
    pvalue_statistics = pvalue_stats,
    fdr_control = fdr_control,
    expression_quality = expr_quality
  ))
}

# 5. 采样代表性评估
assess_sampling_representativeness <- function(original_seurat, sampled_cds, cancer_type, cell_type) {
  message(paste("评估采样代表性:", cancer_type, "-", cell_type))
  
  common_cells <- intersect(colnames(original_seurat), colnames(sampled_cds))
  
  if(length(common_cells) == 0) {
    return(list(status = "failed", reason = "No common cells"))
  }
  
  # 1. 采样率
  sampling_rate <- ncol(sampled_cds) / ncol(original_seurat)
  
  # 2. 细胞类型分布比较
  if("target_cell_type" %in% colnames(original_seurat@meta.data)) {
    original_types <- table(original_seurat$target_cell_type)
    sampled_types <- table(original_seurat$target_cell_type[common_cells])
    
    # 计算比例差异
    original_props <- original_types / sum(original_types)
    sampled_props <- sampled_types / sum(sampled_types)
    
    # 只比较共同的细胞类型
    common_types <- intersect(names(original_props), names(sampled_props))
    
    if(length(common_types) > 0) {
      prop_differences <- abs(original_props[common_types] - sampled_props[common_types])
      max_prop_diff <- max(prop_differences)
      avg_prop_diff <- mean(prop_differences)
      
      # 卡方检验
      chi_square_test <- tryCatch({
        chisq.test(rbind(original_types[common_types], sampled_types[common_types]))
      }, error = function(e) list(p.value = NA))
      
      type_distribution <- list(
        max_proportion_difference = max_prop_diff,
        avg_proportion_difference = avg_prop_diff,
        chi_square_pvalue = chi_square_test$p.value,
        well_represented = max_prop_diff < 0.1
      )
    } else {
      type_distribution <- list(status = "no_common_types")
    }
  } else {
    type_distribution <- list(status = "no_type_info")
  }
  
  # 3. 表达谱相似性（PCA空间）
  if("pca" %in% names(original_seurat@reductions)) {
    original_pca <- Embeddings(original_seurat, "pca")[, 1:min(10, ncol(Embeddings(original_seurat, "pca")))]
    sampled_pca <- original_pca[common_cells, ]
    
    # 计算PCA空间覆盖度
    pca_coverage <- tryCatch({
      # 计算凸包体积比例
      original_range <- apply(original_pca, 2, range)
      sampled_range <- apply(sampled_pca, 2, range)
      
      coverage_per_pc <- sapply(1:ncol(original_pca), function(i) {
        orig_span <- diff(original_range[, i])
        samp_span <- diff(sampled_range[, i])
        min(samp_span / orig_span, 1)
      })
      
      mean(coverage_per_pc)
    }, error = function(e) NA)
    
    expression_similarity <- list(
      pca_coverage = pca_coverage,
      good_coverage = !is.na(pca_coverage) && pca_coverage > 0.8
    )
  } else {
    expression_similarity <- list(status = "no_pca")
  }
  
  return(list(
    status = "success",
    sampling_rate = sampling_rate,
    n_original = ncol(original_seurat),
    n_sampled = ncol(sampled_cds),
    type_distribution = type_distribution,
    expression_similarity = expression_similarity
  ))
}

# ==================== 综合评估函数 ====================

# 计算综合质量评分
calculate_overall_quality_score <- function(quality_results) {
  scores <- c()
  weights <- c()
  
  # 1. Pseudotime质量 (权重: 25%)
  if(quality_results$pseudotime$status == "success") {
    pt_score <- 0
    pt_weight <- 25
    
    # 有效率评分
    valid_rate <- quality_results$pseudotime$basic_stats$valid_rate
    pt_score <- pt_score + valid_rate * 40
    
    # 单调性评分
    if(!is.na(quality_results$pseudotime$monotonicity)) {
      monotonicity <- abs(quality_results$pseudotime$monotonicity)
      pt_score <- pt_score + min(monotonicity * 60, 60)
    }
    
    scores <- c(scores, min(pt_score, 100))
    weights <- c(weights, pt_weight)
  }
  
  # 2. 生物学一致性 (权重: 30%)
  if(quality_results$biological$status == "success") {
    bio_score <- quality_results$biological$consistency_rate * 100
    scores <- c(scores, bio_score)
    weights <- c(weights, 30)
  }
  
  # 3. 聚类质量 (权重: 20%)
  if(quality_results$clustering$status == "success") {
    cluster_score <- 0
    
    # 轮廓系数评分
    if(!is.na(quality_results$clustering$avg_silhouette)) {
      silhouette_score <- (quality_results$clustering$avg_silhouette + 1) * 50  # 转换到0-100
      cluster_score <- cluster_score + silhouette_score * 0.6
    }
    
    # 聚类数合理性评分
    n_clusters <- quality_results$clustering$n_clusters
    if(n_clusters >= 2 && n_clusters <= 8) {
      cluster_score <- cluster_score + 40 * 0.4
    }
    
    scores <- c(scores, min(cluster_score, 100))
    weights <- c(weights, 20)
  }
  
  # 4. 基因质量 (权重: 15%)
  if(quality_results$genes$status == "success") {
    gene_score <- 0
    
    # 显著基因比例
    sig_rate <- quality_results$genes$gene_statistics$significance_rate
    gene_score <- gene_score + sig_rate * 60
    
    # 相关性强度
    mean_corr <- quality_results$genes$correlation_statistics$mean_correlation
    gene_score <- gene_score + min(mean_corr * 200, 40)
    
    scores <- c(scores, min(gene_score, 100))
    weights <- c(weights, 15)
  }
  
  # 5. 采样质量 (权重: 10%)
  if(quality_results$sampling$status == "success") {
    sampling_score <- 50  # 基础分
    
    # 类型分布保持
    if(quality_results$sampling$type_distribution$status != "no_type_info") {
      if(quality_results$sampling$type_distribution$well_represented) {
        sampling_score <- sampling_score + 30
      }
    }
    
    # 表达谱覆盖
    if(quality_results$sampling$expression_similarity$status != "no_pca") {
      if(quality_results$sampling$expression_similarity$good_coverage) {
        sampling_score <- sampling_score + 20
      }
    }
    
    scores <- c(scores, min(sampling_score, 100))
    weights <- c(weights, 10)
  }
  
  # 计算加权平均
  if(length(scores) > 0) {
    overall_score <- sum(scores * weights) / sum(weights)
    return(overall_score / 10)  # 转换为10分制
  } else {
    return(0)
  }
}

# 主要质量评估函数
comprehensive_quality_assessment <- function(cancer_type, cell_type, 
                                           original_seurat = NULL, 
                                           force_reload = FALSE) {
  
  message(paste("\n=== 开始质量评估:", cancer_type, "-", cell_type, "==="))
  
  # 加载结果文件
  results_file <- file.path(trajectory_output_dir, "monocle3_results",
                           paste0(cancer_type, "_", cell_type, "_trajectory_results.rds"))
  
  if(!file.exists(results_file)) {
    stop(paste("结果文件不存在:", results_file))
  }
  
  results <- readRDS(results_file)
  
  if(is.null(results) || is.null(results$cds)) {
    stop("无效的结果对象")
  }
  
  # 加载原始数据（如果未提供）
  if(is.null(original_seurat)) {
    original_seurat <- results$seurat_obj
  }
  
  # 加载轨迹相关基因
  gene_file <- file.path(trajectory_output_dir, "gene_analysis",
                        paste0(cancer_type, "_", cell_type, "_trajectory_genes.csv"))
  
  trajectory_genes <- NULL
  if(file.exists(gene_file)) {
    trajectory_genes <- read.csv(gene_file, stringsAsFactors = FALSE)
  }
  
  # 执行各项质量评估
  quality_results <- list()
  
  # 1. Pseudotime质量
  quality_results$pseudotime <- assess_pseudotime_quality(results$cds, cancer_type, cell_type)
  
  # 2. 生物学一致性
  quality_results$biological <- validate_biological_consistency(results$cds, original_seurat, 
                                                              cancer_type, cell_type)
  
  # 3. 聚类质量
  quality_results$clustering <- assess_clustering_quality(results$cds, cancer_type, cell_type)
  
  # 4. 基因质量
  quality_results$genes <- assess_gene_quality(trajectory_genes, results$cds, 
                                              cancer_type, cell_type)
  
  # 5. 采样质量
  quality_results$sampling <- assess_sampling_representativeness(original_seurat, results$cds, 
                                                               cancer_type, cell_type)
  
  # 6. 计算综合评分
  quality_results$overall_score <- calculate_overall_quality_score(quality_results)
  
  # 保存质量评估结果
  quality_file <- file.path(quality_output_dir, 
                           paste0(cancer_type, "_", cell_type, "_quality_assessment.rds"))
  saveRDS(quality_results, quality_file)
  
  message(paste("质量评估完成，综合评分:", round(quality_results$overall_score, 2), "/10"))
  
  return(quality_results)
}

# 生成质量评估报告
generate_quality_report <- function(quality_results, cancer_type, cell_type) {
  
  report_file <- file.path(quality_output_dir, 
                          paste0(cancer_type, "_", cell_type, "_quality_report.txt"))
  
  sink(report_file)
  
  cat("=== 轨迹分析数据质量评估报告 ===\n")
  cat("癌种:", cancer_type, "\n")
  cat("细胞类型:", cell_type, "\n")
  cat("评估时间:", as.character(Sys.time()), "\n\n")
  
  # 1. Pseudotime质量
  cat("1. PSEUDOTIME质量评估\n")
  cat("========================\n")
  if(quality_results$pseudotime$status == "success") {
    pt <- quality_results$pseudotime
    cat(sprintf("有效细胞数: %d/%d (%.1f%%)\n", 
                pt$basic_stats$n_valid, pt$basic_stats$n_cells, 
                pt$basic_stats$valid_rate * 100))
    cat(sprintf("Pseudotime范围: %.2f - %.2f\n", 
                pt$basic_stats$range[1], pt$basic_stats$range[2]))
    cat(sprintf("单调性评分: %.3f\n", pt$monotonicity))
    cat(sprintf("分布均匀性p值: %.3e\n", pt$uniformity_test$p.value))
    cat(sprintf("边界细胞比例: %.3f\n", pt$boundary_check$boundary_proportion))
  } else {
    cat("状态: 失败 -", quality_results$pseudotime$reason, "\n")
  }
  
  # 2. 生物学一致性
  cat("\n2. 生物学一致性验证\n")
  cat("====================\n")
  if(quality_results$biological$status == "success") {
    bio <- quality_results$biological
    cat(sprintf("一致性总体评分: %.3f\n", bio$consistency_rate))
    cat("各基因集相关性:\n")
    for(set_name in names(bio$correlations)) {
      if(!is.na(bio$correlations[[set_name]])) {
        cat(sprintf("  %s: %.3f\n", set_name, bio$correlations[[set_name]]))
      }
    }
  } else {
    cat("状态: 失败 -", quality_results$biological$reason, "\n")
  }
  
  # 3. 聚类质量
  cat("\n3. 聚类质量评估\n")
  cat("================\n")
  if(quality_results$clustering$status == "success") {
    cl <- quality_results$clustering
    cat(sprintf("聚类数量: %d\n", cl$n_clusters))
    cat(sprintf("平均轮廓系数: %.3f\n", cl$avg_silhouette))
    cat(sprintf("聚类间平均距离: %.3f\n", cl$avg_center_distance))
    cat(sprintf("平均聚类紧密度: %.3f\n", cl$avg_compactness))
    cat("各聚类大小:", paste(cl$cluster_sizes, collapse = ", "), "\n")
  } else {
    cat("状态: 失败 -", quality_results$clustering$reason, "\n")
  }
  
  # 4. 基因质量
  cat("\n4. 基因特征质量\n")
  cat("================\n")
  if(quality_results$genes$status == "success") {
    gene <- quality_results$genes
    cat(sprintf("总基因数: %d\n", gene$gene_statistics$total_genes))
    cat(sprintf("显著基因数: %d (%.1f%%)\n", 
                gene$gene_statistics$significant_genes,
                gene$gene_statistics$significance_rate * 100))
    cat(sprintf("平均相关性强度: %.3f\n", gene$correlation_statistics$mean_correlation))
    cat(sprintf("最大相关性: %.3f\n", gene$correlation_statistics$max_correlation))
    cat(sprintf("FDR控制比例: %.3f\n", gene$fdr_control$fdr_ratio))
  } else {
    cat("状态: 失败 -", quality_results$genes$reason, "\n")
  }
  
  # 5. 采样质量
  cat("\n5. 采样代表性\n")
  cat("==============\n")
  if(quality_results$sampling$status == "success") {
    samp <- quality_results$sampling
    cat(sprintf("采样率: %.4f (%d/%d)\n", 
                samp$sampling_rate, samp$n_sampled, samp$n_original))
    
    if(samp$type_distribution$status != "no_type_info") {
      cat(sprintf("最大类型比例差异: %.3f\n", samp$type_distribution$max_proportion_difference))
      cat(sprintf("类型分布良好: %s\n", samp$type_distribution$well_represented))
    }
    
    if(samp$expression_similarity$status != "no_pca") {
      cat(sprintf("PCA空间覆盖度: %.3f\n", samp$expression_similarity$pca_coverage))
    }
  } else {
    cat("状态: 失败 -", quality_results$sampling$reason, "\n")
  }
  
  # 6. 综合评分
  cat("\n6. 综合质量评分\n")
  cat("================\n")
  cat(sprintf("总体评分: %.2f/10\n", quality_results$overall_score))
  
  if(quality_results$overall_score >= 8) {
    cat("质量等级: 优秀 ⭐⭐⭐⭐⭐\n")
    cat("建议: 数据质量优秀，适合进行所有类型的下游分析\n")
  } else if(quality_results$overall_score >= 6) {
    cat("质量等级: 良好 ⭐⭐⭐⭐\n")
    cat("建议: 数据质量良好，适合大部分下游分析\n")
  } else if(quality_results$overall_score >= 4) {
    cat("质量等级: 一般 ⭐⭐⭐\n")
    cat("建议: 数据质量一般，建议谨慎解释结果\n")
  } else if(quality_results$overall_score >= 2) {
    cat("质量等级: 较差 ⭐⭐\n")
    cat("建议: 数据质量较差，建议优化分析参数\n")
  } else {
    cat("质量等级: 差 ⭐\n")
    cat("建议: 数据质量较差，建议重新进行分析\n")
  }
  
  cat("\n=== 报告结束 ===\n")
  sink()
  
  message(paste("质量报告已保存:", report_file))
}

# 批量质量评估函数
batch_quality_assessment <- function(cancer_types = NULL, cell_types = NULL) {
  
  if(is.null(cancer_types)) {
    # 自动检测可用的癌种
    results_dir <- file.path(trajectory_output_dir, "monocle3_results")
    result_files <- list.files(results_dir, pattern = "_trajectory_results\\.rds$")
    
    if(length(result_files) == 0) {
      stop("未找到轨迹分析结果文件")
    }
    
    # 解析文件名
    file_info <- strsplit(gsub("_trajectory_results\\.rds$", "", result_files), "_")
    cancer_cell_combinations <- lapply(file_info, function(x) {
      if(length(x) >= 2) {
        list(cancer = x[1], cell = paste(x[-1], collapse = "_"))
      } else {
        NULL
      }
    })
    
    cancer_cell_combinations <- cancer_cell_combinations[!sapply(cancer_cell_combinations, is.null)]
  }
  
  # 执行批量评估
  all_results <- list()
  
  for(i in seq_along(cancer_cell_combinations)) {
    combo <- cancer_cell_combinations[[i]]
    cancer_type <- combo$cancer
    cell_type <- combo$cell
    
    message(paste("\n处理:", cancer_type, "-", cell_type, 
                  paste0("(", i, "/", length(cancer_cell_combinations), ")")))
    
    tryCatch({
      quality_result <- comprehensive_quality_assessment(cancer_type, cell_type)
      generate_quality_report(quality_result, cancer_type, cell_type)
      
      result_key <- paste(cancer_type, cell_type, sep = "_")
      all_results[[result_key]] <- quality_result
      
    }, error = function(e) {
      message(paste("ERROR:", cancer_type, "-", cell_type, ":", e$message))
    })
  }
  
  # 生成汇总报告
  generate_summary_report(all_results)
  
  return(all_results)
}

# 生成汇总报告
generate_summary_report <- function(all_results) {
  
  summary_file <- file.path(quality_output_dir, "quality_summary_report.txt")
  
  sink(summary_file)
  
  cat("=== 轨迹分析质量评估汇总报告 ===\n")
  cat("生成时间:", as.character(Sys.time()), "\n")
  cat("评估样本数:", length(all_results), "\n\n")
  
  # 提取所有评分
  scores <- sapply(all_results, function(x) x$overall_score)
  valid_scores <- scores[!is.na(scores)]
  
  if(length(valid_scores) > 0) {
    cat("综合质量评分统计:\n")
    cat("=================\n")
    cat(sprintf("平均评分: %.2f ± %.2f\n", mean(valid_scores), sd(valid_scores)))
    cat(sprintf("中位数评分: %.2f\n", median(valid_scores)))
    cat(sprintf("评分范围: %.2f - %.2f\n", min(valid_scores), max(valid_scores)))
    
    # 质量等级分布
    excellent <- sum(valid_scores >= 8)
    good <- sum(valid_scores >= 6 & valid_scores < 8)
    fair <- sum(valid_scores >= 4 & valid_scores < 6)
    poor <- sum(valid_scores < 4)
    
    cat("\n质量等级分布:\n")
    cat("=============\n")
    cat(sprintf("优秀 (≥8分): %d (%.1f%%)\n", excellent, excellent/length(valid_scores)*100))
    cat(sprintf("良好 (6-8分): %d (%.1f%%)\n", good, good/length(valid_scores)*100))
    cat(sprintf("一般 (4-6分): %d (%.1f%%)\n", fair, fair/length(valid_scores)*100))
    cat(sprintf("较差 (<4分): %d (%.1f%%)\n", poor, poor/length(valid_scores)*100))
    
    # 详细结果列表
    cat("\n详细评分列表:\n")
    cat("=============\n")
    sorted_results <- sort(valid_scores, decreasing = TRUE)
    for(i in 1:length(sorted_results)) {
      name <- names(sorted_results)[i]
      score <- sorted_results[i]
      cat(sprintf("%2d. %-30s: %.2f\n", i, name, score))
    }
  }
  
  cat("\n=== 汇总报告结束 ===\n")
  sink()
  
  message(paste("汇总报告已保存:", summary_file))
}

# ==================== 可视化函数 ====================

# 生成质量评估可视化
generate_quality_plots <- function(quality_results, cancer_type, cell_type) {
  
  library(ggplot2)
  library(patchwork)
  
  plots_dir <- file.path(quality_output_dir, "plots")
  dir.create(plots_dir, showWarnings = FALSE)
  
  plot_list <- list()
  
  # 1. 综合评分雷达图（如果数据足够）
  if(length(quality_results) > 1) {
    # 这里可以添加雷达图代码
  }
  
  # 2. Pseudotime分布图
  if(quality_results$pseudotime$status == "success") {
    # 这里可以添加pseudotime分布图
  }
  
  # 保存图形
  plot_file <- file.path(plots_dir, paste0(cancer_type, "_", cell_type, "_quality_plots.pdf"))
  
  if(length(plot_list) > 0) {
    pdf(plot_file, width = 12, height = 8)
    for(p in plot_list) {
      print(p)
    }
    dev.off()
    
    message(paste("质量评估图形已保存:", plot_file))
  }
}

# ==================== 主执行函数 ====================
if (FALSE){
# 使用示例
if(!interactive()) {
  message("开始批量质量评估...")
  
  # 执行批量评估
  quality_results <- batch_quality_assessment()
  
  message("质量评估完成!")
  message(paste("结果保存在:", quality_output_dir))
  
} else {
  message("质量评估代码加载完成!")
  message("使用方法:")
  message("1. 单个评估: comprehensive_quality_assessment('BC', 'T_cells')")
  message("2. 批量评估: batch_quality_assessment()")
  message("3. 生成报告: generate_quality_report(quality_results, 'BC', 'T_cells')")
}

# 辅助函数：计算变异系数
cv <- function(x) {
  if(length(x) == 0 || all(is.na(x))) return(NA)
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

# 辅助函数：获取CDS logcounts数据
get_cds_logcounts <- function(cds) {
  tryCatch({
    logcounts(cds)
  }, error = function(e) {
    tryCatch({
      assay(cds, "logcounts")
    }, error = function(e2) {
      tryCatch({
        assay(cds, "data")
      }, error = function(e3) {
        log2(counts(cds) + 1)
      })
    })
  })
}

message("轨迹分析数据质量评估代码加载完成!") 
}