# trajectory_analysis_v18.R
# 修复tradeSeq函数加载问题的最终版本

library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)
library(DelayedArray)
library(BiocParallel)
library(stringr)
library(viridis)
library(RColorBrewer)
library(S4Vectors)
library(Matrix)
library(SingleCellExperiment)

# 尝试加载tradeSeq和CytoTRACE，如果失败则跳过
tradeseq_available <- FALSE
cytotrace_available <- FALSE

# 安全加载tradeSeq
tryCatch({
  library(tradeSeq)
  tradeseq_available <- TRUE
  message("tradeSeq加载成功")
}, error = function(e) {
  message("tradeSeq加载失败，将跳过tradeSeq分析")
})

# 安全加载CytoTRACE
tryCatch({
  library(CytoTRACE)
  cytotrace_available <- TRUE
  message("CytoTRACE加载成功")
}, error = function(e) {
  message("CytoTRACE加载失败，将使用简化多样性分析")
})

# 设置并行计算
future::plan(multicore, workers = 4)
DelayedArray::setAutoBPPARAM(MulticoreParam(workers = 4))
register(MulticoreParam(workers = 4))

# 参数设置
output_base_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation"
trajectory_output_dir <- file.path(output_base_dir, "trajectory_analysis")
dir.create(trajectory_output_dir, showWarnings = FALSE, recursive = TRUE)

# 创建各子目录
sub_dirs <- c("cytotrace_results", "monocle3_results", "tradeseq_results", 
              "plots", "gene_analysis", "pan_cancer_analysis", "qc_reports", "error_logs", "debug_info")
for (dir in sub_dirs) {
  dir.create(file.path(trajectory_output_dir, dir), showWarnings = FALSE, recursive = TRUE)
}

# 全局错误收集器
error_collector <- list()

# 调试信息保存函数
save_debug_info <- function(obj, step, cancer_type, cell_type, extra_info = NULL) {
  debug_file <- file.path(trajectory_output_dir, "debug_info",
                         paste0(cancer_type, "_", cell_type, "_", step, "_debug.txt"))
  
  sink(debug_file)
  cat("=== 调试信息 ===\n")
  cat("步骤:", step, "\n")
  cat("癌种:", cancer_type, "\n")
  cat("细胞类型:", cell_type, "\n")
  cat("时间:", as.character(Sys.time()), "\n\n")
  
  if (!is.null(extra_info)) {
    cat("额外信息:\n")
    print(extra_info)
  }
  
  cat("==================\n")
  sink()
}

# 安全读取RDS文件的函数
safe_read_rds <- function(file_path) {
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  tryCatch({
    file_size <- file.info(file_path)$size
    if (is.na(file_size) || file_size < 100) {
      return(NULL)
    }
    
    result <- readRDS(file_path)
    if (is.null(result)) {
      return(NULL)
    }
    
    return(result)
    
  }, error = function(e) {
    message(paste("读取文件失败:", file_path, "-", e$message))
    backup_path <- paste0(file_path, ".corrupted.", Sys.Date())
    if (file.exists(file_path)) {
      file.rename(file_path, backup_path)
    }
    return(NULL)
  })
}

# 错误记录函数
log_error <- function(step, cancer_type, cell_type, error_msg, details = NULL) {
  error_id <- paste(cancer_type, cell_type, step, sep = "_")
  
  error_info <- list(
    timestamp = Sys.time(),
    cancer_type = cancer_type,
    cell_type = cell_type,
    step = step,
    error_message = as.character(error_msg),
    details = details
  )
  
  error_collector[[error_id]] <<- error_info
  
  error_file <- file.path(trajectory_output_dir, "error_logs", 
                         paste0(cancer_type, "_", cell_type, "_", step, "_error.txt"))
  
  sink(error_file)
  cat("=== 错误详情 ===\n")
  cat("时间:", as.character(error_info$timestamp), "\n")
  cat("癌种:", cancer_type, "\n")
  cat("细胞类型:", cell_type, "\n")
  cat("步骤:", step, "\n")
  cat("错误信息:", error_info$error_message, "\n")
  if (!is.null(details)) {
    cat("详细信息:\n")
    print(details)
  }
  cat("===================\n")
  sink()
  
  message(paste("ERROR [", step, "]:", cancer_type, "-", cell_type, ":", error_msg))
}

# 安全执行函数
safe_execute <- function(expr, step, cancer_type, cell_type, default_return = NULL) {
  tryCatch({
    result <- expr
    message(paste("SUCCESS [", step, "]:", cancer_type, "-", cell_type))
    return(result)
  }, error = function(e) {
    log_error(step, cancer_type, cell_type, e$message, list(call = deparse(substitute(expr))))
    return(default_return)
  }, warning = function(w) {
    message(paste("WARNING [", step, "]:", cancer_type, "-", cell_type, ":", w$message))
    suppressWarnings({
      result <- expr
    })
    return(result)
  })
}

# 安全保存RDS文件的函数
safe_save_rds <- function(object, file_path, compress = "xz") {
  tryCatch({
    temp_file <- paste0(file_path, ".tmp")
    saveRDS(object, temp_file, compress = compress)
    test_read <- readRDS(temp_file)
    if (is.null(test_read)) {
      stop("保存验证失败")
    }
    file.rename(temp_file, file_path)
    message(paste("文件保存成功:", file_path))
    return(TRUE)
  }, error = function(e) {
    message(paste("保存文件失败:", file_path, "-", e$message))
    temp_file <- paste0(file_path, ".tmp")
    if (file.exists(temp_file)) {
      file.remove(temp_file)
    }
    return(FALSE)
  })
}

# CDS对象数据访问函数
get_cds_counts <- function(cds) {
  tryCatch({
    counts(cds)
  }, error = function(e) {
    tryCatch({
      assay(cds, "counts")
    }, error = function(e2) {
      as.matrix(assays(cds)[[1]])
    })
  })
}

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

# 读取已修复数据的函数
load_fixed_data <- function(cancer_type, cell_type) {
  message(paste("加载已修复数据:", cancer_type, "-", cell_type))
  
  quickfix_file <- file.path(output_base_dir, cancer_type, 
                            paste0(cancer_type, "_fine_annotated_quickfix.rds"))
  main_file <- file.path(output_base_dir, cancer_type, 
                        paste0(cancer_type, "_fine_annotated.rds"))
  detailed_file <- file.path(output_base_dir, cancer_type, 
                           paste0(cancer_type, "_", cell_type, "_detailed.rds"))
  
  if (file.exists(quickfix_file)) {
    use_file <- quickfix_file
    message("  使用快速修复版本")
  } else if (file.exists(main_file)) {
    use_file <- main_file
    message("  使用原始版本")
  } else {
    message(paste("  主文件不存在:", cancer_type))
    return(NULL)
  }
  
  if (!file.exists(detailed_file)) {
    message(paste("  详细文件不存在:", detailed_file))
    return(NULL)
  }
  
  return(safe_execute({
    main_obj <- safe_read_rds(use_file)
    if (is.null(main_obj)) {
      stop("主文件读取失败")
    }
    
    if ("fine_cell_type" %in% colnames(main_obj@meta.data)) {
      cell_type_col <- "fine_cell_type"
    } else if ("fine_cell_type_fixed" %in% colnames(main_obj@meta.data)) {
      cell_type_col <- "fine_cell_type_fixed"
    } else {
      stop("未找到细胞类型注释列")
    }
    
    target_pattern <- paste0("^", cell_type, "_")
    target_cells <- colnames(main_obj)[grepl(target_pattern, main_obj@meta.data[[cell_type_col]])]
    
    message(paste("  目标细胞数:", length(target_cells)))
    
    if (length(target_cells) < 50) {
      stop(paste("目标细胞数量过少:", length(target_cells)))
    }
    
    detailed_obj <- safe_read_rds(detailed_file)
    if (is.null(detailed_obj)) {
      stop("详细文件读取失败")
    }
    
    common_cells <- intersect(target_cells, colnames(detailed_obj))
    
    if (length(common_cells) < 30) {
      stop(paste("详细文件中目标细胞过少:", length(common_cells)))
    }
    
    clean_subset <- subset(detailed_obj, cells = common_cells)
    clean_subset$target_cell_type <- main_obj@meta.data[common_cells, cell_type_col]
    
    message(paste("  最终分析细胞数:", ncol(clean_subset)))
    
    min_genes <- 200
    min_cells_per_gene <- 3
    
    gene_counts <- Matrix::rowSums(GetAssayData(clean_subset, layer = "counts") > 0)
    keep_genes <- gene_counts >= min_cells_per_gene
    
    cell_gene_counts <- Matrix::colSums(GetAssayData(clean_subset, layer = "counts") > 0)
    keep_cells <- cell_gene_counts >= min_genes
    
    message(paste("  质量过滤: 保留", sum(keep_genes), "基因,", sum(keep_cells), "细胞"))
    
    if (sum(keep_cells) < 30) {
      stop("质量过滤后细胞数过少")
    }
    
    clean_subset <- subset(clean_subset, features = rownames(clean_subset)[keep_genes], 
                          cells = colnames(clean_subset)[keep_cells])
    
    return(clean_subset)
    
  }, "data_loading", cancer_type, cell_type, NULL))
}

# 优化的基因多样性分析（减少内存使用）
calculate_gene_diversity_score <- function(seurat_obj, cancer_type, cell_type) {
  message(paste("计算基因表达多样性评分:", cancer_type, "-", cell_type))
  
  return(safe_execute({
    # 首先进行细胞采样以减少内存使用
    n_cells <- ncol(seurat_obj)
    max_cells_for_diversity <- 3000  # 进一步减少
    
    if (n_cells > max_cells_for_diversity) {
      message(paste("  对", n_cells, "个细胞采样至", max_cells_for_diversity, "进行多样性计算"))
      set.seed(42)
      sample_cells <- sample(colnames(seurat_obj), max_cells_for_diversity)
      subset_obj <- subset(seurat_obj, cells = sample_cells)
    } else {
      subset_obj <- seurat_obj
    }
    
    expr_matrix <- GetAssayData(subset_obj, assay = "RNA", layer = "data")
    
    n_cells_calc <- ncol(expr_matrix)
    n_genes <- nrow(expr_matrix)
    
    message(paste("  计算矩阵大小:", n_genes, "基因,", n_cells_calc, "细胞"))
    
    # 进一步限制基因数量
    if (n_genes > 500) {
      if (length(VariableFeatures(subset_obj)) > 0) {
        var_genes <- VariableFeatures(subset_obj)[1:min(500, length(VariableFeatures(subset_obj)))]
        expr_matrix <- expr_matrix[var_genes, ]
      } else {
        gene_means <- Matrix::rowMeans(expr_matrix)
        top_genes <- names(sort(gene_means, decreasing = TRUE))[1:500]
        expr_matrix <- expr_matrix[top_genes, ]
      }
      message(paste("  限制基因数量:", nrow(expr_matrix)))
    }
    
    # 分批计算多样性以节省内存
    batch_size <- 200
    n_batches <- ceiling(ncol(expr_matrix) / batch_size)
    diversity_scores_all <- numeric(ncol(expr_matrix))
    names(diversity_scores_all) <- colnames(expr_matrix)
    
    message(paste("  分", n_batches, "批计算多样性"))
    
    for (i in 1:n_batches) {
      start_idx <- (i - 1) * batch_size + 1
      end_idx <- min(i * batch_size, ncol(expr_matrix))
      
      batch_matrix <- expr_matrix[, start_idx:end_idx]
      
      batch_diversity <- apply(batch_matrix, 2, function(x) {
        x <- x[x > 0]
        if (length(x) < 3) return(NA)
        
        p <- x / sum(x)
        -sum(p * log2(p + 1e-10))
      })
      
      diversity_scores_all[start_idx:end_idx] <- batch_diversity
    }
    
    # 将多样性评分应用到原始对象的所有细胞
    diversity_scores_full <- rep(NA, ncol(seurat_obj))
    names(diversity_scores_full) <- colnames(seurat_obj)
    
    # 对于计算过的细胞，使用实际值
    computed_cells <- intersect(names(diversity_scores_all), names(diversity_scores_full))
    diversity_scores_full[computed_cells] <- diversity_scores_all[computed_cells]
    
    # 对于未计算的细胞，使用插值或平均值
    if (sum(is.na(diversity_scores_full)) > 0) {
      mean_diversity <- mean(diversity_scores_all, na.rm = TRUE)
      diversity_scores_full[is.na(diversity_scores_full)] <- mean_diversity
    }
    
    # 添加到Seurat对象
    seurat_obj$diversity_score <- diversity_scores_full
    seurat_obj$diversity_rank <- rank(-diversity_scores_full, na.last = "keep")
    
    # 归一化到0-1范围
    valid_scores <- diversity_scores_full[!is.na(diversity_scores_full)]
    if (length(valid_scores) > 0) {
      min_score <- min(valid_scores)
      max_score <- max(valid_scores)
      seurat_obj$normalized_diversity <- (diversity_scores_full - min_score) / (max_score - min_score)
    } else {
      seurat_obj$normalized_diversity <- NA
    }
    
    message(paste("  计算完成，有效细胞数:", sum(!is.na(diversity_scores_full))))
    
    return(list(
      seurat_obj = seurat_obj,
      diversity_results = list(
        diversity_scores = diversity_scores_full,
        summary = summary(valid_scores)
      )
    ))
    
  }, "diversity_calculation", cancer_type, cell_type, NULL))
}

# 强制聚类合并函数（按就近原则）
force_merge_to_target_clusters <- function(umap_coords, initial_clusters, target_n_clusters) {
  message(paste("    强制合并", length(unique(initial_clusters)), "个聚类到", target_n_clusters, "个"))
  
  if (length(unique(initial_clusters)) <= target_n_clusters) {
    message("    当前聚类数已满足要求")
    return(initial_clusters)
  }
  
  # 计算每个聚类的中心和大小
  cluster_info <- data.frame()
  unique_clusters <- unique(initial_clusters)
  
  for (cluster_id in unique_clusters) {
    cluster_cells <- which(initial_clusters == cluster_id)
    center <- colMeans(umap_coords[cluster_cells, , drop = FALSE])
    
    cluster_info <- rbind(cluster_info, data.frame(
      cluster = cluster_id,
      x = center[1],
      y = center[2],
      size = length(cluster_cells),
      stringsAsFactors = FALSE
    ))
  }
  
  # 计算聚类间距离矩阵
  suppressWarnings({
    dist_matrix <- as.matrix(dist(cluster_info[, c("x", "y")]))
  })
  rownames(dist_matrix) <- colnames(dist_matrix) <- cluster_info$cluster
  
  # 使用层次聚类方法进行合并
  suppressWarnings({
    hc <- hclust(as.dist(dist_matrix), method = "ward.D2")
  })
  
  # 切割树形结构获得目标聚类数
  new_cluster_assignments <- cutree(hc, k = target_n_clusters)
  
  # 创建映射关系
  cluster_mapping <- data.frame(
    old_cluster = cluster_info$cluster,
    new_cluster = new_cluster_assignments,
    stringsAsFactors = FALSE
  )
  
  # 应用映射到所有细胞
  merged_clusters <- initial_clusters
  for (i in 1:nrow(cluster_mapping)) {
    old_id <- cluster_mapping$old_cluster[i]
    new_id <- cluster_mapping$new_cluster[i]
    merged_clusters[initial_clusters == old_id] <- new_id
  }
  
  # 重新编号确保连续性
  unique_new <- unique(merged_clusters)
  final_clusters <- merged_clusters
  for (i in seq_along(unique_new)) {
    final_clusters[merged_clusters == unique_new[i]] <- i
  }
  
  message(paste("    成功合并到", length(unique(final_clusters)), "个聚类"))
  
  # 打印合并统计
  merge_stats <- cluster_mapping %>%
    group_by(new_cluster) %>%
    summarise(
      n_merged = n(),
      total_cells = sum(cluster_info$size[cluster_info$cluster %in% old_cluster]),
      .groups = 'drop'
    )
  
  for (i in 1:nrow(merge_stats)) {
    message(paste("      新聚类", i, ":", merge_stats$n_merged[i], "个原聚类合并,", 
                  merge_stats$total_cells[i], "个细胞"))
  }
  
  return(as.character(final_clusters))
}

# 简化的Monocle3分析函数
run_monocle3_analysis <- function(seurat_obj, cancer_type, cell_type) {
  message(paste("运行Monocle3分析:", cancer_type, "-", cell_type))
  
  return(safe_execute({
    # 数据采样
    n_cells <- ncol(seurat_obj)
    max_cells_for_monocle <- 800  # 进一步减少
    
    if (n_cells > max_cells_for_monocle) {
      message(paste("  对", n_cells, "个细胞进行采样，保留", max_cells_for_monocle, "个"))
      set.seed(42)
      
      # 分层采样确保每个亚型都有代表
      if ("target_cell_type" %in% colnames(seurat_obj@meta.data)) {
        cell_types <- unique(seurat_obj$target_cell_type)
        sample_cells <- c()
        
        for (ct in cell_types) {
          ct_cells <- colnames(seurat_obj)[seurat_obj$target_cell_type == ct]
          n_sample <- min(length(ct_cells), ceiling(max_cells_for_monocle / length(cell_types)))
          if (n_sample > 0) {
            sample_cells <- c(sample_cells, sample(ct_cells, n_sample))
          }
        }
        
        # 如果还有剩余配额，随机补充
        if (length(sample_cells) < max_cells_for_monocle) {
          remaining_cells <- setdiff(colnames(seurat_obj), sample_cells)
          additional_n <- min(length(remaining_cells), max_cells_for_monocle - length(sample_cells))
          if (additional_n > 0) {
            sample_cells <- c(sample_cells, sample(remaining_cells, additional_n))
          }
        }
      } else {
        sample_cells <- sample(colnames(seurat_obj), max_cells_for_monocle)
      }
      
      seurat_obj <- subset(seurat_obj, cells = sample_cells)
      message(paste("  实际采样细胞数:", ncol(seurat_obj)))
    }
    
    # 确保有必要的降维结果
    message("  准备降维结果...")
    if (!"pca" %in% names(seurat_obj@reductions)) {
      seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
    }
    
    if (!"umap" %in% names(seurat_obj@reductions)) {
      seurat_obj <- RunUMAP(seurat_obj, dims = 1:8, verbose = FALSE)
    }
    
    # 创建初始聚类（使用k-means）
    message("  创建初始聚类...")
    umap_coords <- Embeddings(seurat_obj, "umap")
    
    # 检查UMAP坐标是否包含NA或无限值
    if (any(is.na(umap_coords)) || any(!is.finite(umap_coords))) {
      message("    检测到UMAP坐标中有异常值，进行清理...")
      
      # 移除含有NA或无限值的细胞
      valid_cells <- complete.cases(umap_coords) & 
                    apply(umap_coords, 1, function(x) all(is.finite(x)))
      
      if (sum(valid_cells) < 50) {
        stop("有效UMAP坐标细胞数过少")
      }
      
      seurat_obj <- subset(seurat_obj, cells = rownames(umap_coords)[valid_cells])
      umap_coords <- umap_coords[valid_cells, ]
      message(paste("    清理后细胞数:", nrow(umap_coords)))
    }
    
    # 使用k-means创建非常少的初始聚类
    set.seed(42)
    initial_k <- min(6, floor(ncol(seurat_obj) / 120))  # 每120个细胞一个聚类
    initial_k <- max(initial_k, 3)  # 至少3个聚类
    
    kmeans_result <- kmeans(umap_coords, centers = initial_k, nstart = 10, iter.max = 50)
    initial_clusters <- kmeans_result$cluster
    
    message(paste("  k-means创建了", length(unique(initial_clusters)), "个初始聚类"))
    
    # 强制合并到目标聚类数
    target_clusters <- 3  # 减少到3个聚类
    final_clusters <- force_merge_to_target_clusters(umap_coords, initial_clusters, target_clusters)
    
    # 确保聚类标签是字符串且无异常值
    final_clusters <- as.character(final_clusters)
    final_clusters[is.na(final_clusters) | final_clusters == ""] <- "1"
    
    # 添加聚类信息到Seurat对象
    seurat_obj$final_clusters <- final_clusters
    
    # 创建CDS对象
    message("  创建CDS对象...")
    expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
    
    # 进一步简化基因集
    if (nrow(expression_matrix) > 1500) {
      if (length(VariableFeatures(seurat_obj)) > 0) {
        var_genes <- VariableFeatures(seurat_obj)[1:min(1500, length(VariableFeatures(seurat_obj)))]
        expression_matrix <- expression_matrix[var_genes, ]
      } else {
        gene_means <- Matrix::rowMeans(expression_matrix)
        top_genes <- names(sort(gene_means, decreasing = TRUE))[1:1500]
        expression_matrix <- expression_matrix[top_genes, ]
      }
      message(paste("  CDS基因数限制为:", nrow(expression_matrix)))
    }
    
    cell_metadata <- seurat_obj@meta.data
    gene_metadata <- data.frame(
      gene_short_name = rownames(expression_matrix),
      row.names = rownames(expression_matrix)
    )
    
    cds <- new_cell_data_set(expression_matrix,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_metadata)
    
    message("  预处理CDS对象...")
    cds <- preprocess_cds(cds, num_dim = 6, verbose = FALSE, method = "PCA")
    
    if (is.null(reducedDims(cds)$PCA)) {
      stop("PCA降维失败")
    }
    
    message("  设置UMAP降维...")
    # 直接使用Seurat的UMAP结果
    reducedDims(cds)$UMAP <- umap_coords
    
    message("  简化聚类流程...")
    # 完全不使用Monocle3聚类，直接创建最简单的数据结构
    
    # 方法：直接在colData中添加聚类信息，避免使用clusters()函数
    colData(cds)$manual_clusters <- factor(final_clusters)
    colData(cds)$cluster_labels <- paste0("Cluster_", final_clusters)
    
    # 验证聚类设置
    n_clusters <- length(unique(final_clusters))
    message(paste("  设置聚类数:", n_clusters))
    
    # 跳过正式的learn_graph，使用最简单的轨迹构建
    message("  构建简化轨迹...")
    
    # 创建基于聚类中心的简单轨迹
    cluster_centers <- data.frame()
    for (i in unique(final_clusters)) {
      cluster_cells <- which(final_clusters == i)
      center <- colMeans(umap_coords[cluster_cells, , drop = FALSE])
      cluster_centers <- rbind(cluster_centers, data.frame(
        cluster = i,
        UMAP_1 = center[1],
        UMAP_2 = center[2],
        size = length(cluster_cells)
      ))
    }
    
    # 计算聚类间的连接（最小生成树）
    if (nrow(cluster_centers) > 1) {
      dist_matrix <- as.matrix(dist(cluster_centers[, c("UMAP_1", "UMAP_2")]))
      
      # 创建最小生成树
      mst_edges <- data.frame()
      connected <- c(1)  # 从第一个聚类开始
      
      while (length(connected) < nrow(cluster_centers)) {
        min_dist <- Inf
        best_from <- NULL
        best_to <- NULL
        
        for (from in connected) {
          for (to in 1:nrow(cluster_centers)) {
            if (!(to %in% connected) && dist_matrix[from, to] < min_dist) {
              min_dist <- dist_matrix[from, to]
              best_from <- from
              best_to <- to
            }
          }
        }
        
        if (!is.null(best_to)) {
          mst_edges <- rbind(mst_edges, data.frame(
            from = best_from,
            to = best_to,
            weight = min_dist
          ))
          connected <- c(connected, best_to)
        }
      }
      
      message(paste("  构建了", nrow(mst_edges), "条轨迹边"))
    }
    
    # 计算简化的pseudotime
    message("  计算简化的pseudotime...")
    
    # 选择根聚类（最大的聚类）
    root_cluster <- cluster_centers$cluster[which.max(cluster_centers$size)]
    message(paste("    选择根聚类:", root_cluster))
    
    # 基于到根聚类中心的距离计算pseudotime
    root_center <- cluster_centers[cluster_centers$cluster == root_cluster, c("UMAP_1", "UMAP_2")]
    
    cell_pseudotime <- numeric(nrow(umap_coords))
    for (i in 1:nrow(umap_coords)) {
      dist_to_root <- sqrt(sum((umap_coords[i, ] - root_center)^2))
      cell_pseudotime[i] <- dist_to_root
    }
    
    # 归一化pseudotime到0-100范围
    cell_pseudotime <- (cell_pseudotime - min(cell_pseudotime)) / 
                      (max(cell_pseudotime) - min(cell_pseudotime)) * 100
    
    # 添加到CDS对象
    colData(cds)$pseudotime <- cell_pseudotime
    
    # 验证pseudotime计算
    valid_pseudotime <- sum(!is.na(cell_pseudotime) & is.finite(cell_pseudotime))
    
    if (valid_pseudotime == 0) {
      stop("所有pseudotime值无效")
    }
    
    message(paste("  简化Monocle3分析完成，有效pseudotime:", valid_pseudotime, "/", ncol(cds)))
    return(cds)
    
  }, "monocle3", cancer_type, cell_type, NULL))
}

# 简化的tradeSeq分析函数（如果包可用）
run_tradeseq_analysis <- function(cds, cancer_type, cell_type) {
  message(paste("运行tradeSeq分析:", cancer_type, "-", cell_type))
  
  if (!tradeseq_available) {
    message("  跳过tradeSeq分析（包不可用）")
    return(NULL)
  }
  
  if (is.null(cds)) {
    log_error("tradeseq", cancer_type, cell_type, "CDS对象为空")
    return(NULL)
  }
  
  return(safe_execute({
    if (!"pseudotime" %in% colnames(colData(cds))) {
      stop("pseudotime未计算")
    }
    
    pseudotime_values <- colData(cds)$pseudotime
    valid_cells <- !is.na(pseudotime_values) & is.finite(pseudotime_values)
    
    if (sum(valid_cells) < 30) {
      stop(paste("有效pseudotime细胞数过少:", sum(valid_cells)))
    }
    
    cds_filtered <- cds[, valid_cells]
    pseudotime_filtered <- pseudotime_values[valid_cells]
    cell_weights <- rep(1, length(pseudotime_filtered))
    
    # 使用正确的CDS数据访问方法
    counts <- get_cds_counts(cds_filtered)
    min_cells <- ceiling(ncol(counts) * 0.05)
    gene_filter <- Matrix::rowSums(counts > 0) >= min_cells
    counts_filtered <- counts[gene_filter, ]
    
    message(paste("  tradeSeq输入:", nrow(counts_filtered), "基因,", ncol(counts_filtered), "细胞"))
    
    if (nrow(counts_filtered) < 20) {
      stop("可分析基因数过少")
    }
    
    if (nrow(counts_filtered) > 50) {
      set.seed(42)
      sample_genes <- sample(rownames(counts_filtered), 50)
      counts_filtered <- counts_filtered[sample_genes, ]
      message(paste("  随机采样", nrow(counts_filtered), "基因进行分析"))
    }
    
    message("  拟合GAM模型...")
    sce <- fitGAM(counts = as.matrix(counts_filtered), 
                  pseudotime = pseudotime_filtered,
                  cellWeights = cell_weights,
                  nknots = 3,
                  verbose = FALSE)
    
    message("  进行Wald检验...")
    # 检查waldTest函数是否可用
    if (!exists("waldTest", mode = "function")) {
      stop("waldTest函数不可用")
    }
    
    wald_test <- waldTest(sce, name = "pseudotime")
    
    if (is.null(wald_test) || nrow(wald_test) == 0) {
      stop("Wald检验失败")
    }
    
    sig_genes <- rownames(wald_test)[wald_test$waldStat_pvalue < 0.05]
    
    message(paste("  发现", length(sig_genes), "个轨迹相关基因"))
    
    return(list(
      sce = sce, 
      wald_test = wald_test, 
      sig_genes = sig_genes,
      n_cells = ncol(counts_filtered),
      n_genes = nrow(counts_filtered)
    ))
    
  }, "tradeseq", cancer_type, cell_type, NULL))
}

# 基因分析函数（简化版，修复数据访问）
analyze_trajectory_genes_simple <- function(cds, cancer_type, cell_type) {
  message("分析轨迹相关基因（简化版）...")
  
  if (is.null(cds)) {
    log_error("gene_analysis", cancer_type, cell_type, "CDS对象为空")
    return(data.frame())
  }
  
  return(safe_execute({
    if (!"pseudotime" %in% colnames(colData(cds))) {
      stop("pseudotime未计算")
    }
    
    pseudotime_values <- colData(cds)$pseudotime
    valid_cells <- !is.na(pseudotime_values) & is.finite(pseudotime_values)
    
    if (sum(valid_cells) < 15) {
      stop("有效pseudotime细胞数过少")
    }
    
    # 使用基于相关性的简单分析替代graph_test
    cds_filtered <- cds[, valid_cells]
    
    # 使用正确的CDS数据访问方法
    expr_matrix <- get_cds_logcounts(cds_filtered)
    pseudotime_filtered <- pseudotime_values[valid_cells]
    
    # 限制基因数量
    if (nrow(expr_matrix) > 50) {
      gene_means <- Matrix::rowMeans(expr_matrix)
      top_genes <- names(sort(gene_means, decreasing = TRUE))[1:50]
      expr_matrix <- expr_matrix[top_genes, ]
    }
    
    message(paste("  分析", nrow(expr_matrix), "个基因与pseudotime的相关性"))
    
    # 计算与pseudotime的相关性
    correlation_results <- data.frame()
    
    for (gene in rownames(expr_matrix)) {
      gene_expr <- as.numeric(expr_matrix[gene, ])
      
      if (sum(gene_expr > 0) >= 3) {  # 至少3个细胞表达
        cor_test <- tryCatch({
          cor.test(gene_expr, pseudotime_filtered, method = "spearman")
        }, error = function(e) {
          NULL
        })
        
        if (!is.null(cor_test) && !is.na(cor_test$p.value)) {
          correlation_results <- rbind(correlation_results, data.frame(
            gene_short_name = gene,
            correlation = cor_test$estimate,
            p_value = cor_test$p.value,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    if (nrow(correlation_results) > 0) {
      # 计算q值
      correlation_results$q_value <- p.adjust(correlation_results$p_value, method = "BH")
      correlation_results <- correlation_results[order(correlation_results$q_value), ]
      
      sig_genes <- subset(correlation_results, q_value < 0.05)
      
      if (nrow(sig_genes) > 0) {
        gene_file <- file.path(trajectory_output_dir, "gene_analysis",
                              paste0(cancer_type, "_", cell_type, "_trajectory_genes.csv"))
        write.csv(sig_genes, gene_file, row.names = FALSE)
        message(paste("  发现", nrow(sig_genes), "个显著轨迹基因"))
      } else {
        message("  未发现显著轨迹基因")
      }
      
      return(correlation_results)
    } else {
      message("  无法分析基因相关性")
      return(data.frame())
    }
    
  }, "gene_analysis", cancer_type, cell_type, data.frame()))
}

# trajectory_analysis_v19.R
# 改进可视化：添加轨迹方向和细胞类群标注的版本

# [前面的代码保持不变，只修改可视化相关函数]

# 获取细胞类型到聚类的映射（增强版）
get_cluster_cell_type_mapping_enhanced <- function(seurat_obj, cds) {
  if (is.null(cds)) {
    return(list(mapping_df = data.frame(), cluster_summary = data.frame()))
  }
  
  if (!"manual_clusters" %in% colnames(colData(cds))) {
    return(list(mapping_df = data.frame(), cluster_summary = data.frame()))
  }
  
  clusters_vec <- colData(cds)$manual_clusters
  cell_types <- seurat_obj$target_cell_type[colnames(cds)]
  
  mapping_df <- data.frame(
    cell_id = colnames(cds),
    cluster = clusters_vec,
    cell_type = cell_types,
    stringsAsFactors = FALSE
  )
  
  # 计算每个聚类的主要细胞类型
  cluster_mapping <- mapping_df %>%
    group_by(cluster) %>%
    summarise(
      dominant_cell_type = names(sort(table(cell_type), decreasing = TRUE))[1],
      cell_count = n(),
      type_proportion = max(table(cell_type)) / n(),
      all_types = paste(names(table(cell_type)), collapse = ", "),
      .groups = 'drop'
    )
  
  # 创建更详细的聚类标签
  cluster_mapping$detailed_label <- apply(cluster_mapping, 1, function(row) {
    dominant_type <- row[["dominant_cell_type"]]
    proportion <- as.numeric(row[["type_proportion"]])
    cell_count <- as.numeric(row[["cell_count"]])
    
    # 简化细胞类型名称
    simplified_type <- gsub(".*_", "", dominant_type)
    simplified_type <- gsub("cells?", "", simplified_type, ignore.case = TRUE)
    
    # 根据比例确定标签
    if (proportion >= 0.8) {
      label <- simplified_type
    } else if (proportion >= 0.5) {
      label <- paste0(simplified_type, "-like")
    } else {
      label <- "Mixed"
    }
    
    return(paste0(label, "\n(n=", cell_count, ")"))
  })
  
  # 创建简短标签用于图例
  cluster_mapping$short_label <- apply(cluster_mapping, 1, function(row) {
    dominant_type <- row[["dominant_cell_type"]]
    simplified_type <- gsub(".*_", "", dominant_type)
    simplified_type <- gsub("cells?", "", simplified_type, ignore.case = TRUE)
    return(simplified_type)
  })
  
  return(list(
    mapping_df = mapping_df,
    cluster_summary = cluster_mapping
  ))
}

# 计算轨迹路径和方向
calculate_trajectory_path <- function(cds, cluster_centers) {
  if (is.null(cds) || nrow(cluster_centers) < 2) {
    return(NULL)
  }
  
  # 获取pseudotime值
  pseudotime_values <- colData(cds)$pseudotime
  clusters_vec <- colData(cds)$manual_clusters
  
  # 计算每个聚类的平均pseudotime
  cluster_pseudotime <- data.frame()
  for (i in unique(clusters_vec)) {
    cluster_cells <- which(clusters_vec == i)
    avg_pseudotime <- mean(pseudotime_values[cluster_cells], na.rm = TRUE)
    
    cluster_pseudotime <- rbind(cluster_pseudotime, data.frame(
      cluster = i,
      avg_pseudotime = avg_pseudotime,
      stringsAsFactors = FALSE
    ))
  }
  
  # 合并聚类中心和pseudotime信息
  cluster_info <- merge(cluster_centers, cluster_pseudotime, by = "cluster")
  cluster_info <- cluster_info[order(cluster_info$avg_pseudotime), ]
  
  # 创建轨迹路径
  if (nrow(cluster_info) >= 2) {
    trajectory_path <- data.frame(
      x = cluster_info$UMAP_1,
      y = cluster_info$UMAP_2,
      pseudotime_order = 1:nrow(cluster_info),
      cluster = cluster_info$cluster
    )
    
    return(trajectory_path)
  }
  
  return(NULL)
}

# 改进的可视化函数
generate_trajectory_plots_enhanced <- function(seurat_obj, cds, diversity_results, cancer_type, cell_type) {
  message("生成增强轨迹可视化图...")
  
  plot_base <- file.path(trajectory_output_dir, "plots", 
                        paste0(cancer_type, "_", cell_type, "_trajectory_plots"))
  
  return(safe_execute({
    base_theme <- theme_void() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    # 1. 多样性评分可视化
    if (!is.null(diversity_results) && "normalized_diversity" %in% colnames(seurat_obj@meta.data)) {
      p1 <- FeaturePlot(seurat_obj, features = "normalized_diversity", reduction = "umap") +
        scale_color_viridis_c(name = "Gene\nDiversity", option = "plasma") +
        ggtitle(paste(cancer_type, "-", gsub("_", " ", cell_type)),
                subtitle = "Gene Expression Diversity Score") +
        base_theme
      
      ggsave(paste0(plot_base, "_diversity.png"), p1, width = 10, height = 8, dpi = 300, bg = "white")
      ggsave(paste0(plot_base, "_diversity.pdf"), p1, width = 10, height = 8, bg = "white")
    }
    
    # 2. 增强的轨迹图
    if (!is.null(cds)) {
      umap_coords <- reducedDims(cds)$UMAP
      
      # 获取聚类信息
      cluster_info <- get_cluster_cell_type_mapping_enhanced(seurat_obj, cds)
      
      # 计算聚类中心
      cluster_centers <- data.frame()
      if ("manual_clusters" %in% colnames(colData(cds))) {
        clusters_vec <- colData(cds)$manual_clusters
        for (i in unique(clusters_vec)) {
          cluster_cells <- which(clusters_vec == i)
          center <- colMeans(umap_coords[cluster_cells, , drop = FALSE])
          cluster_centers <- rbind(cluster_centers, data.frame(
            cluster = i,
            UMAP_1 = center[1],
            UMAP_2 = center[2],
            size = length(cluster_cells)
          ))
        }
      }
      
      # 2a. 增强的Pseudotime轨迹图（带轨迹和方向）
      if ("pseudotime" %in% colnames(colData(cds))) {
        pt_values <- colData(cds)$pseudotime
        valid_pt <- !is.na(pt_values) & is.finite(pt_values)
        
        if (sum(valid_pt) > 0) {
          # 创建基础数据
          plot_data <- data.frame(
            UMAP_1 = umap_coords[, 1],
            UMAP_2 = umap_coords[, 2],
            Pseudotime = pt_values,
            Cell_Type = seurat_obj$target_cell_type[colnames(cds)]
          )
          
          # 计算轨迹路径
          trajectory_path <- calculate_trajectory_path(cds, cluster_centers)
          
          # 创建pseudotime图
          p2 <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Pseudotime)) +
            geom_point(size = 0.8, alpha = 0.7) +
            scale_color_viridis_c(name = "Pseudotime", option = "viridis") +
            ggtitle(paste(cancer_type, "-", gsub("_", " ", cell_type)),
                    subtitle = "Cell State Transition Trajectory with Direction") +
            base_theme
          
          # 添加轨迹路径和方向
          if (!is.null(trajectory_path) && nrow(trajectory_path) >= 2) {
            # 添加轨迹线
            p2 <- p2 + 
              geom_path(data = trajectory_path, 
                       aes(x = x, y = y), 
                       color = "black", 
                       size = 2, 
                       alpha = 0.8,
                       inherit.aes = FALSE) +
              geom_path(data = trajectory_path, 
                       aes(x = x, y = y), 
                       color = "white", 
                       size = 1, 
                       alpha = 0.9,
                       inherit.aes = FALSE)
            
            # 添加方向箭头
            for (i in 1:(nrow(trajectory_path) - 1)) {
              start_point <- trajectory_path[i, ]
              end_point <- trajectory_path[i + 1, ]
              
              # 计算箭头位置（路径中点）
              arrow_x <- (start_point$x + end_point$x) / 2
              arrow_y <- (start_point$y + end_point$y) / 2
              
              # 计算方向
              dx <- end_point$x - start_point$x
              dy <- end_point$y - start_point$y
              
              p2 <- p2 + 
                geom_segment(aes(x = arrow_x - dx*0.1, y = arrow_y - dy*0.1,
                                xend = arrow_x + dx*0.1, yend = arrow_y + dy*0.1),
                            arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                            color = "red", size = 1.2, inherit.aes = FALSE)
            }
            
            # 添加起点和终点标注
            start_point <- trajectory_path[1, ]
            end_point <- trajectory_path[nrow(trajectory_path), ]
            
            p2 <- p2 + 
              annotate("text", x = start_point$x, y = start_point$y, 
                      label = "START", color = "darkgreen", size = 4, 
                      fontface = "bold", vjust = -1) +
              annotate("text", x = end_point$x, y = end_point$y, 
                      label = "END", color = "darkred", size = 4, 
                      fontface = "bold", vjust = -1)
          }
          
          ggsave(paste0(plot_base, "_pseudotime.png"), p2, width = 14, height = 10, dpi = 300, bg = "white")
          ggsave(paste0(plot_base, "_pseudotime.pdf"), p2, width = 14, height = 10, bg = "white")
        }
      }
      
      # 2b. 增强的聚类图（带具体细胞类群名称）
      if ("manual_clusters" %in% colnames(colData(cds)) && nrow(cluster_info$cluster_summary) > 0) {
        plot_data <- data.frame(
          UMAP_1 = umap_coords[, 1],
          UMAP_2 = umap_coords[, 2],
          Cluster = colData(cds)$manual_clusters,
          Cell_Type = seurat_obj$target_cell_type[colnames(cds)]
        )
        
        # 创建聚类颜色映射
        cluster_colors <- RColorBrewer::brewer.pal(min(nrow(cluster_info$cluster_summary), 11), "Spectral")
        names(cluster_colors) <- cluster_info$cluster_summary$cluster
        
        # 创建聚类图
        p3 <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = factor(Cluster))) +
          geom_point(size = 0.8, alpha = 0.7) +
          scale_color_manual(name = "Cluster", 
                           values = cluster_colors,
                           labels = cluster_info$cluster_summary$short_label) +
          ggtitle(paste(cancer_type, "-", gsub("_", " ", cell_type)),
                  subtitle = "Cluster Annotation with Cell Type Identity") +
          base_theme
        
        # 添加聚类中心标注
        if (nrow(cluster_centers) > 0) {
          # 合并聚类标签信息
          cluster_centers_labeled <- merge(cluster_centers, 
                                          cluster_info$cluster_summary[, c("cluster", "detailed_label")], 
                                          by = "cluster", all.x = TRUE)
          
          # 添加聚类中心点
          p3 <- p3 + 
            geom_point(data = cluster_centers_labeled, 
                      aes(x = UMAP_1, y = UMAP_2), 
                      color = "black", size = 4, shape = 21, 
                      fill = "white", stroke = 2, inherit.aes = FALSE)
          
          # 添加聚类标签
          p3 <- p3 + 
            geom_text(data = cluster_centers_labeled, 
                     aes(x = UMAP_1, y = UMAP_2, label = detailed_label), 
                     color = "black", size = 3.5, fontface = "bold", 
                     vjust = -1.5, hjust = 0.5, inherit.aes = FALSE,
                     bg.colour = "white", bg.r = 0.1)
        }
        
        ggsave(paste0(plot_base, "_clusters.png"), p3, width = 16, height = 12, dpi = 300, bg = "white")
        ggsave(paste0(plot_base, "_clusters.pdf"), p3, width = 16, height = 12, bg = "white")
        
        # 保存详细的聚类映射信息
        mapping_file <- file.path(trajectory_output_dir, "plots",
                                 paste0(cancer_type, "_", cell_type, "_cluster_mapping.csv"))
        write.csv(cluster_info$cluster_summary, mapping_file, row.names = FALSE)
      }
    }
    
    # 3. 细胞类型分布
    if ("target_cell_type" %in% colnames(seurat_obj@meta.data)) {
      # 计算细胞类型分布统计
      cell_type_counts <- table(seurat_obj$target_cell_type)
      cell_type_labels <- paste0(names(cell_type_counts), "\n(n=", cell_type_counts, ")")
      names(cell_type_labels) <- names(cell_type_counts)
      
      p4 <- DimPlot(seurat_obj, group.by = "target_cell_type", reduction = "umap") +
        scale_color_brewer(name = "Cell Type", palette = "Set2",
                          labels = cell_type_labels[levels(factor(seurat_obj$target_cell_type))]) +
        ggtitle(paste(cancer_type, "-", gsub("_", " ", cell_type)),
                subtitle = "Original Cell Subtype Distribution") +
        base_theme +
        theme(legend.text = element_text(size = 9))
      
      ggsave(paste0(plot_base, "_celltypes.png"), p4, width = 14, height = 10, dpi = 300, bg = "white")
      ggsave(paste0(plot_base, "_celltypes.pdf"), p4, width = 14, height = 10, bg = "white")
    }
    
    # 4. 轨迹分析总结图
    if (!is.null(cds) && "pseudotime" %in% colnames(colData(cds)) && "manual_clusters" %in% colnames(colData(cds))) {
      # 创建四面板图
      pt_values <- colData(cds)$pseudotime
      clusters_vec <- colData(cds)$manual_clusters
      
      # 面板数据
      summary_data <- data.frame(
        UMAP_1 = umap_coords[, 1],
        UMAP_2 = umap_coords[, 2],
        Pseudotime = pt_values,
        Cluster = clusters_vec,
        Diversity = seurat_obj$normalized_diversity[colnames(cds)],
        Cell_Type = seurat_obj$target_cell_type[colnames(cds)]
      )
      
      # 创建组合图
      p_pt <- ggplot(summary_data, aes(x = UMAP_1, y = UMAP_2, color = Pseudotime)) +
        geom_point(size = 0.5, alpha = 0.8) +
        scale_color_viridis_c(option = "viridis") +
        ggtitle("Pseudotime") +
        theme_void() + theme(legend.position = "bottom")
      
      p_cl <- ggplot(summary_data, aes(x = UMAP_1, y = UMAP_2, color = factor(Cluster))) +
        geom_point(size = 0.5, alpha = 0.8) +
        scale_color_brewer(palette = "Set3") +
        ggtitle("Clusters") +
        theme_void() + theme(legend.position = "bottom")
      
      p_div <- ggplot(summary_data, aes(x = UMAP_1, y = UMAP_2, color = Diversity)) +
        geom_point(size = 0.5, alpha = 0.8) +
        scale_color_viridis_c(option = "plasma") +
        ggtitle("Diversity") +
        theme_void() + theme(legend.position = "bottom")
      
      p_ct <- ggplot(summary_data, aes(x = UMAP_1, y = UMAP_2, color = Cell_Type)) +
        geom_point(size = 0.5, alpha = 0.8) +
        scale_color_brewer(palette = "Set2") +
        ggtitle("Cell Types") +
        theme_void() + theme(legend.position = "bottom", legend.text = element_text(size = 6))
      
      # 组合图
      combined_plot <- (p_pt + p_cl) / (p_div + p_ct) +
        plot_annotation(
          title = paste(cancer_type, "-", gsub("_", " ", cell_type), "Trajectory Analysis Summary"),
          theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
        )
      
      ggsave(paste0(plot_base, "_summary.png"), combined_plot, width = 16, height = 12, dpi = 300, bg = "white")
      ggsave(paste0(plot_base, "_summary.pdf"), combined_plot, width = 16, height = 12, bg = "white")
    }
    
    message(paste("  增强轨迹图已保存:", paste0(plot_base, "*")))
    
  }, "visualization", cancer_type, cell_type, NULL))
}

# 更新主分析函数中的可视化步骤
analyze_single_cell_type_enhanced <- function(cancer_type, cell_type, force_rerun = FALSE) {
  message(paste("\n", paste(rep("=", 60), collapse = "")))
  message(paste("开始分析:", cancer_type, "-", cell_type))
  message(paste(rep("=", 60), collapse = ""))
  
  results_file <- file.path(trajectory_output_dir, "monocle3_results",
                           paste0(cancer_type, "_", cell_type, "_trajectory_results.rds"))
  
  if (!force_rerun && file.exists(results_file)) {
    existing_result <- safe_read_rds(results_file)
    if (!is.null(existing_result)) {
      message(paste("成功加载已有结果:", cancer_type, "-", cell_type))
      return(existing_result)
    }
  }
  
  # 1. 数据加载
  message("步骤1: 数据加载")
  seurat_obj <- load_fixed_data(cancer_type, cell_type)
  
  if (is.null(seurat_obj)) {
    log_error("complete_analysis", cancer_type, cell_type, "数据加载失败")
    return(NULL)
  }
  
  # 2. 基因多样性分析
  message("步骤2: 基因多样性分析")
  diversity_results <- calculate_gene_diversity_score(seurat_obj, cancer_type, cell_type)
  
  if (!is.null(diversity_results)) {
    seurat_obj <- diversity_results$seurat_obj
    diversity_data <- diversity_results$diversity_results
    
    diversity_file <- file.path(trajectory_output_dir, "cytotrace_results",
                               paste0(cancer_type, "_", cell_type, "_diversity.rds"))
    safe_execute({
      safe_save_rds(diversity_data, diversity_file)
    }, "save_diversity", cancer_type, cell_type, NULL)
  } else {
    diversity_data <- NULL
  }
  
  # 3. 简化Monocle3分析
  message("步骤3: 简化Monocle3分析")
  cds <- run_monocle3_analysis(seurat_obj, cancer_type, cell_type)
  
  # 4. tradeSeq分析
  message("步骤4: tradeSeq分析")
  tradeseq_results <- run_tradeseq_analysis(cds, cancer_type, cell_type)
  
  # 5. 简化基因分析
  message("步骤5: 简化基因分析")
  trajectory_genes <- analyze_trajectory_genes_simple(cds, cancer_type, cell_type)
  
  # 6. 增强可视化
  message("步骤6: 生成增强可视化")
  generate_trajectory_plots_enhanced(seurat_obj, cds, diversity_data, cancer_type, cell_type)
  
  # 7. 质控报告
  message("步骤7: 生成质控报告")
  generate_qc_report(cancer_type, cell_type, seurat_obj, diversity_data, cds, tradeseq_results, trajectory_genes)
  
  # 整合结果
  results <- list(
    cancer_type = cancer_type,
    cell_type = cell_type,
    seurat_obj = seurat_obj,
    cds = cds,
    diversity_results = diversity_data,
    tradeseq_results = tradeseq_results,
    trajectory_genes = trajectory_genes,
    n_cells = ncol(seurat_obj),
    analysis_time = Sys.time(),
    success = list(
      diversity_analysis = !is.null(diversity_data),
      monocle3_simplified = !is.null(cds),
      tradeseq = !is.null(tradeseq_results),
      gene_analysis = !is.null(trajectory_genes) && nrow(trajectory_genes) > 0
    )
  )
  
  # 保存结果
  save_success <- safe_execute({
    safe_save_rds(results, results_file)
  }, "save_results", cancer_type, cell_type, FALSE)
  
  message(paste("分析完成:", cancer_type, "-", cell_type))
  return(results)
}

# 更新主分析函数
main_trajectory_analysis_enhanced <- function(force_rerun = FALSE) {
  message("开始增强细胞状态转换轨迹分析...")
  message(paste("输出目录:", trajectory_output_dir))
  
  analysis_summary_file <- file.path(output_base_dir, "analysis_summary.csv")
  if (!file.exists(analysis_summary_file)) {
    stop("分析汇总文件不存在:", analysis_summary_file)
  }
  
  analysis_summary <- read.csv(analysis_summary_file, stringsAsFactors = FALSE)
  cancer_types <- analysis_summary$Cancer_Type
  
  message(paste("发现癌种:", length(cancer_types)))
  
  all_results <- list()
  success_count <- 0
  total_count <- 0
  
  for (cancer_type in cancer_types) {
    message(paste("\n处理癌种:", cancer_type))
    
    cancer_dir <- file.path(output_base_dir, cancer_type)
    if (!dir.exists(cancer_dir)) {
      message(paste("癌种目录不存在:", cancer_dir))
      next
    }
    
    detailed_files <- list.files(cancer_dir, pattern = "_detailed\\.rds$", full.names = TRUE)
    if (length(detailed_files) == 0) {
      message(paste("未找到详细文件:", cancer_type))
      next
    }
    
    available_cell_types <- gsub(paste0(".*", cancer_type, "_(.+)_detailed\\.rds$"), "\\1", basename(detailed_files))
    cell_types_to_analyze <- intersect(available_cell_types, trajectory_cell_types)
    
    message(paste("将分析细胞类型:", paste(cell_types_to_analyze, collapse = ", ")))
    
    for (cell_type in cell_types_to_analyze) {
      total_count <- total_count + 1
      
      result <- analyze_single_cell_type_enhanced(cancer_type, cell_type, force_rerun)
      
      if (!is.null(result)) {
        result_key <- paste(cancer_type, cell_type, sep = "_")
        all_results[[result_key]] <- result
        success_count <- success_count + 1
      }
      
      gc()
    }
  }
  
  message(paste("\n分析完成! 成功:", success_count, "/ 总计:", total_count))
  
  if (length(error_collector) > 0) {
    error_summary_file <- file.path(trajectory_output_dir, "error_logs", "error_summary.rds")
    safe_save_rds(error_collector, error_summary_file)
    
    message(paste("发现", length(error_collector), "个错误"))
    
    cat("\n=== 错误摘要 ===\n")
    error_by_step <- table(sapply(error_collector, function(x) x$step))
    print(error_by_step)
  }
  
  if (length(all_results) > 0) {
    final_file <- file.path(trajectory_output_dir, "trajectory_analysis_results.rds")
    save_success <- safe_save_rds(all_results, final_file)
    
    if (save_success) {
      message(paste("成功结果保存在:", final_file))
    }
  }
  
  return(list(
    results = all_results,
    errors = error_collector,
    success_rate = success_count / total_count
  ))
}

# 执行增强分析
if (!interactive()) {
  message("开始执行增强轨迹分析...")
  analysis_output <- main_trajectory_analysis_enhanced(force_rerun = FALSE)
  
  message("\n=== 最终统计 ===")
  message(paste("成功分析:", length(analysis_output$results)))
  message(paste("错误数量:", length(analysis_output$errors)))
  message(paste("成功率:", round(analysis_output$success_rate * 100, 1), "%"))
  
} else {
  message("脚本在交互模式下运行")
  message("请运行: analysis_output <- main_trajectory_analysis_enhanced()")
}

message("\n轨迹分析脚本v19加载完成!")
message("增强可视化功能:")
message("1. Pseudotime图: 添加轨迹路径、方向箭头、起点终点标注")
message("2. Cluster图: 添加细胞类群名称、聚类中心标注、详细标签")
message("3. 新增四面板总结图")
message("4. 改进图例和标签显示")
message("5. 保存更高质量的图片")