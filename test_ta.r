#!/usr/bin/env Rscript

cat("正在加载 trajectory-aa.r...\n")
source("/mnt/public7/pancancercol/hefeng/cluster-annotation/trajectory-aa.r")
cat("加载完成\n")

cat("assess_pseudotime_quality 是否存在:", exists("assess_pseudotime_quality"), "\n")
# 或者如果在同一目录下
# source("trajectory-aa.r")

# 测试函数
test_individual_functions <- function() {
  
  # 设置路径（确保与原文件一致）
  trajectory_output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/trajectory_analysis_2"
  
  # 加载数据
  results_file <- file.path(trajectory_output_dir, "monocle3_results", "Melan_B_cell_trajectory_results.rds")
  results <- readRDS(results_file)
  cds <- results$cds
  seurat_obj <- results$seurat_obj
  
  # 测试1: Pseudotime质量评估
  cat("=== 测试 Pseudotime 质量评估 ===\n")
  pt_result <- tryCatch({
    assess_pseudotime_quality(cds, "Melan", "B_cell")
  }, error = function(e) {
    cat(paste("ERROR in pseudotime:", e$message, "\n"))
    print(traceback())
    return(list(status = "error", message = e$message))
  })
  
  cat("Pseudotime结果:\n")
  str(pt_result, max.level = 2)
  
  # 测试2: 生物学一致性验证
  cat("\n=== 测试生物学一致性验证 ===\n")
  bio_result <- tryCatch({
    validate_biological_consistency(cds, seurat_obj, "Melan", "B_cell")
  }, error = function(e) {
    cat(paste("ERROR in biological consistency:", e$message, "\n"))
    print(traceback())
    return(list(status = "error", message = e$message))
  })
  
  cat("生物学一致性结果:\n")
  str(bio_result, max.level = 2)
  
  # 测试3: 聚类质量评估
  cat("\n=== 测试聚类质量评估 ===\n")
  cluster_result <- tryCatch({
    assess_clustering_quality(cds, "Melan", "B_cell")
  }, error = function(e) {
    cat(paste("ERROR in clustering:", e$message, "\n"))
    print(traceback())
    return(list(status = "error", message = e$message))
  })
  
  cat("聚类质量结果:\n")
  str(cluster_result, max.level = 2)
  
  # 测试4: 基因质量评估
  cat("\n=== 测试基因质量评估 ===\n")
  trajectory_genes <- results$trajectory_genes
  gene_result <- tryCatch({
    assess_gene_quality(trajectory_genes, cds, "Melan", "B_cell")
  }, error = function(e) {
    cat(paste("ERROR in gene quality:", e$message, "\n"))
    print(traceback())
    return(list(status = "error", message = e$message))
  })
  
  cat("基因质量结果:\n")
  str(gene_result, max.level = 2)
  
  # 测试5: 采样质量评估
  cat("\n=== 测试采样质量评估 ===\n")
  sampling_result <- tryCatch({
    assess_sampling_representativeness(seurat_obj, cds, "Melan", "B_cell")
  }, error = function(e) {
    cat(paste("ERROR in sampling:", e$message, "\n"))
    print(traceback())
    return(list(status = "error", message = e$message))
  })
  
  cat("采样质量结果:\n")
  str(sampling_result, max.level = 2)
  
  return(list(
    pseudotime = pt_result,
    biological = bio_result,
    clustering = cluster_result,
    genes = gene_result,
    sampling = sampling_result
  ))
}

# 调试聚类函数
debug_clustering <- function() {
  trajectory_output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/trajectory_analysis_2"
  results_file <- file.path(trajectory_output_dir, "monocle3_results", "Melan_B_cell_trajectory_results.rds")
  results <- readRDS(results_file)
  cds <- results$cds
  
  clusters <- colData(cds)$manual_clusters
  umap_coords <- reducedDims(cds)$UMAP
  
  cat("=== 聚类转换调试 ===\n")
  cat(paste("clusters类型:", class(clusters), "\n"))
  cat(paste("clusters长度:", length(clusters), "\n"))
  cat(paste("unique clusters:", paste(unique(clusters), collapse = ", "), "\n"))
  
  # 检查转换
  cat("尝试转换为数值...\n")
  numeric_clusters <- tryCatch({
    as.numeric(as.factor(clusters))
  }, error = function(e) {
    cat(paste("转换错误:", e$message, "\n"))
    return(NULL)
  })
  
  if(!is.null(numeric_clusters)) {
    cat(paste("数值clusters长度:", length(numeric_clusters), "\n"))
    cat(paste("数值clusters范围:", paste(range(numeric_clusters), collapse = " - "), "\n"))
  }
  
  # 检查聚类数量
  n_clusters <- length(unique(clusters))
  n_cells <- nrow(umap_coords)
  cat(paste("聚类数量:", n_clusters, "\n"))
  cat(paste("细胞数量:", n_cells, "\n"))
  
  # 检查轮廓系数计算的前提条件
  if(n_clusters > 1 && n_clusters < n_cells) {
    cat("满足轮廓系数计算条件\n")
    
    cat("尝试计算距离矩阵...\n")
    dist_matrix <- tryCatch({
      dist(umap_coords)
    }, error = function(e) {
      cat(paste("距离计算错误:", e$message, "\n"))
      return(NULL)
    })
    
    if(!is.null(dist_matrix)) {
      cat("距离矩阵计算成功\n")
      cat(paste("距离矩阵长度:", length(dist_matrix), "\n"))
      
      if(!is.null(numeric_clusters)) {
        cat("尝试计算轮廓系数...\n")
        sil_result <- tryCatch({
          silhouette(numeric_clusters, dist_matrix)
        }, error = function(e) {
          cat(paste("轮廓系数错误:", e$message, "\n"))
          return(NULL)
        })
        
        if(!is.null(sil_result)) {
          cat("轮廓系数计算成功!\n")
        }
      }
    }
  } else {
    cat("不满足轮廓系数计算条件\n")
  }
}

# 主执行部分
cat("开始调试 Melan B_cell 样本...\n")

# 首先运行聚类调试
debug_clustering()

# 然后运行完整测试
test_results <- test_individual_functions()

cat("\n调试完成!\n")