# conserved_specific_analysis.R
# 跨癌种T细胞轨迹保守和特异模式分析

library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)

# 设置路径
trajectory_output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/trajectory_analysis_2"
conserved_output_dir <- file.path(trajectory_output_dir, "conserved_specific_analysis")
dir.create(conserved_output_dir, showWarnings = FALSE, recursive = TRUE)

# 创建子目录
sub_dirs <- c("conserved_modules", "specific_modules", "comparative_plots", 
              "network_analysis", "functional_analysis")
for(dir in sub_dirs) {
  dir.create(file.path(conserved_output_dir, dir), showWarnings = FALSE)
}

# ==================== 数据加载和预处理 ====================

# 加载所有可用的轨迹分析结果
load_all_trajectory_results <- function() {
  message("加载所有轨迹分析结果...")
  
  results_dir <- file.path(trajectory_output_dir, "monocle3_results")
  result_files <- list.files(results_dir, pattern = "_trajectory_results\\.rds$", full.names = TRUE)
  
  all_results <- list()
  
  for(file in result_files) {
    result_name <- gsub("_trajectory_results\\.rds$", "", basename(file))
    
    tryCatch({
      result <- readRDS(file)
      if(!is.null(result) && !is.null(result$cds)) {
        all_results[[result_name]] <- result
        message(paste("加载成功:", result_name))
      }
    }, error = function(e) {
      message(paste("加载失败:", result_name, "-", e$message))
    })
  }
  
  message(paste("共加载", length(all_results), "个结果"))
  return(all_results)
}

# 加载所有轨迹相关基因
load_all_trajectory_genes <- function() {
  message("加载所有轨迹相关基因...")
  
  gene_dir <- file.path(trajectory_output_dir, "gene_analysis")
  gene_files <- list.files(gene_dir, pattern = "_trajectory_genes\\.csv$", full.names = TRUE)
  
  all_genes <- list()
  
  for(file in gene_files) {
    gene_name <- gsub("_trajectory_genes\\.csv$", "", basename(file))
    
    tryCatch({
      genes <- read.csv(file, stringsAsFactors = FALSE)
      if(nrow(genes) > 0) {
        all_genes[[gene_name]] <- genes
        message(paste("加载基因:", gene_name, "-", nrow(genes), "个基因"))
      }
    }, error = function(e) {
      message(paste("基因加载失败:", gene_name, "-", e$message))
    })
  }
  
  return(all_genes)
}

# 标准化基因名称（ENSEMBL -> SYMBOL）
convert_ensembl_to_symbol <- function(ensembl_ids) {
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
  tryCatch({
    symbols <- mapIds(org.Hs.eg.db, 
                     keys = ensembl_ids,
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
    
    # 保留成功转换的基因
    valid_symbols <- symbols[!is.na(symbols)]
    return(valid_symbols)
    
  }, error = function(e) {
    message("基因名转换失败，保持原名")
    names(ensembl_ids) <- ensembl_ids
    return(ensembl_ids)
  })
}

# ==================== 保守模块识别 ====================

# 识别跨癌种保守的轨迹基因
identify_conserved_trajectory_genes <- function(all_trajectory_genes, 
                                              min_cancer_types = 2,
                                              p_threshold = 0.05,
                                              correlation_threshold = 0.15) {
  
  message("识别保守轨迹基因...")
  
  if(length(all_trajectory_genes) < min_cancer_types) {
    stop("可用癌种数量不足")
  }
  
  # 1. 提取所有显著基因
  all_significant_genes <- list()
  
  for(cancer_cell in names(all_trajectory_genes)) {
    genes <- all_trajectory_genes[[cancer_cell]]
    
    # 筛选显著基因
    significant <- genes[genes$q_value < p_threshold & 
                        abs(genes$correlation) > correlation_threshold, ]
    
    if(nrow(significant) > 0) {
      # 转换基因名
      gene_symbols <- convert_ensembl_to_symbol(significant$gene_short_name)
      
      if(length(gene_symbols) > 0) {
        significant$gene_symbol <- gene_symbols[significant$gene_short_name]
        significant <- significant[!is.na(significant$gene_symbol), ]
        
        all_significant_genes[[cancer_cell]] <- significant
      }
    }
  }
  
  if(length(all_significant_genes) == 0) {
    stop("未找到显著基因")
  }
  
  # 2. 计算基因出现频率
  all_gene_symbols <- unlist(lapply(all_significant_genes, function(x) x$gene_symbol))
  gene_frequency <- table(all_gene_symbols)
  
  # 3. 识别保守基因（在至少min_cancer_types个癌种中出现）
  conserved_genes <- names(gene_frequency)[gene_frequency >= min_cancer_types]
  
  message(paste("发现", length(conserved_genes), "个保守基因"))
  
  # 4. 详细分析保守基因
  conserved_analysis <- data.frame()
  
  for(gene in conserved_genes) {
    gene_info <- data.frame()
    
    for(cancer_cell in names(all_significant_genes)) {
      gene_data <- all_significant_genes[[cancer_cell]]
      gene_row <- gene_data[gene_data$gene_symbol == gene, ]
      
      if(nrow(gene_row) > 0) {
        parts <- strsplit(cancer_cell, "_")[[1]]
        cancer_type <- parts[1]
        cell_type <- paste(parts[-1], collapse = "_")
        
        gene_info <- rbind(gene_info, data.frame(
          gene_symbol = gene,
          cancer_type = cancer_type,
          cell_type = cell_type,
          correlation = gene_row$correlation[1],
          p_value = gene_row$p_value[1],
          q_value = gene_row$q_value[1],
          ensembl_id = gene_row$gene_short_name[1]
        ))
      }
    }
    
    if(nrow(gene_info) > 0) {
      # 计算保守性统计
      conservation_stats <- data.frame(
        gene_symbol = gene,
        n_cancer_types = nrow(gene_info),
        mean_correlation = mean(gene_info$correlation),
        sd_correlation = sd(gene_info$correlation),
        consistent_direction = length(unique(sign(gene_info$correlation))) == 1,
        mean_significance = mean(-log10(gene_info$p_value)),
        conservation_score = nrow(gene_info) * abs(mean(gene_info$correlation))
      )
      
      conserved_analysis <- rbind(conserved_analysis, conservation_stats)
    }
  }
  
  # 按保守性评分排序
  conserved_analysis <- conserved_analysis[order(conserved_analysis$conservation_score, 
                                               decreasing = TRUE), ]
  
  return(list(
    conserved_genes = conserved_genes,
    conserved_analysis = conserved_analysis,
    detailed_data = all_significant_genes,
    gene_frequency = gene_frequency
  ))
}

# 识别癌种特异性基因
identify_cancer_specific_genes <- function(all_trajectory_genes, 
                                         conserved_genes,
                                         p_threshold = 0.05,
                                         correlation_threshold = 0.2) {
  
  message("识别癌种特异性基因...")
  
  specific_genes <- list()
  
  for(cancer_cell in names(all_trajectory_genes)) {
    genes <- all_trajectory_genes[[cancer_cell]]
    
    # 筛选显著且非保守的基因
    significant <- genes[genes$q_value < p_threshold & 
                        abs(genes$correlation) > correlation_threshold, ]
    
    if(nrow(significant) > 0) {
      # 转换基因名
      gene_symbols <- convert_ensembl_to_symbol(significant$gene_short_name)
      valid_genes <- gene_symbols[!is.na(gene_symbols)]
      
      # 筛选非保守基因
      specific_to_cancer <- setdiff(valid_genes, conserved_genes)
      
      if(length(specific_to_cancer) > 0) {
        specific_data <- significant[significant$gene_short_name %in% 
                                   names(gene_symbols)[gene_symbols %in% specific_to_cancer], ]
        specific_data$gene_symbol <- gene_symbols[specific_data$gene_short_name]
        
        specific_genes[[cancer_cell]] <- specific_data
        
        message(paste(cancer_cell, "特异基因:", length(specific_to_cancer), "个"))
      }
    }
  }
  
  return(specific_genes)
}

# ==================== 轨迹模式比较分析 ====================

# 比较不同癌种的轨迹特征
compare_trajectory_patterns <- function(all_results) {
  message("比较轨迹模式...")
  
  trajectory_comparison <- data.frame()
  
  for(result_name in names(all_results)) {
    result <- all_results[[result_name]]
    cds <- result$cds
    
    if(!is.null(cds) && "pseudotime" %in% colnames(colData(cds))) {
      parts <- strsplit(result_name, "_")[[1]]
      cancer_type <- parts[1]
      cell_type <- paste(parts[-1], collapse = "_")
      
      # 提取轨迹特征
      pseudotime_values <- colData(cds)$pseudotime
      valid_pt <- !is.na(pseudotime_values) & is.finite(pseudotime_values)
      
      if(sum(valid_pt) > 0) {
        pt_valid <- pseudotime_values[valid_pt]
        
        # 计算轨迹特征
        features <- data.frame(
          cancer_type = cancer_type,
          cell_type = cell_type,
          n_cells = ncol(cds),
          n_valid_pseudotime = sum(valid_pt),
          pseudotime_range = diff(range(pt_valid)),
          pseudotime_mean = mean(pt_valid),
          pseudotime_sd = sd(pt_valid),
          pseudotime_skewness = calculate_skewness(pt_valid),
          n_clusters = length(unique(colData(cds)$manual_clusters)),
          stringsAsFactors = FALSE
        )
        
        # 计算轨迹复杂度
        if("diversity_score" %in% colnames(colData(cds))) {
          diversity_scores <- colData(cds)$diversity_score[valid_pt]
          features$diversity_mean <- mean(diversity_scores, na.rm = TRUE)
          features$diversity_range <- diff(range(diversity_scores, na.rm = TRUE))
        } else {
          features$diversity_mean <- NA
          features$diversity_range <- NA
        }
        
        trajectory_comparison <- rbind(trajectory_comparison, features)
      }
    }
  }
  
  return(trajectory_comparison)
}

# 计算偏度
calculate_skewness <- function(x) {
  if(length(x) < 3) return(NA)
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- sd(x)
  skew <- (n / ((n-1) * (n-2))) * sum(((x - mean_x) / sd_x)^3)
  return(skew)
}

# ==================== 功能富集分析 ====================

# 对保守基因进行功能富集分析
analyze_conserved_gene_functions <- function(conserved_genes) {
  message("分析保守基因功能...")
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  if(length(conserved_genes) < 5) {
    message("保守基因数量不足，跳过功能分析")
    return(NULL)
  }
  
  # 转换为ENTREZ ID
  entrez_ids <- tryCatch({
    bitr(conserved_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("基因ID转换失败")
    return(NULL)
  })
  
  if(is.null(entrez_ids) || nrow(entrez_ids) == 0) {
    return(NULL)
  }
  
  enrichment_results <- list()
  
  # GO富集分析
  tryCatch({
    go_bp <- enrichGO(gene = entrez_ids$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
    
    if(!is.null(go_bp) && nrow(go_bp@result) > 0) {
      enrichment_results$GO_BP <- go_bp
    }
  }, error = function(e) {
    message("GO BP富集失败:", e$message)
  })
  
  # KEGG富集分析
  tryCatch({
    kegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                      organism = 'hsa',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.1)
    
    if(!is.null(kegg) && nrow(kegg@result) > 0) {
      enrichment_results$KEGG <- kegg
    }
  }, error = function(e) {
    message("KEGG富集失败:", e$message)
  })
  
  # Reactome富集分析
  tryCatch({
    library(ReactomePA)
    reactome <- enrichPathway(gene = entrez_ids$ENTREZID,
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.1,
                            readable = TRUE)
    
    if(!is.null(reactome) && nrow(reactome@result) > 0) {
      enrichment_results$Reactome <- reactome
    }
  }, error = function(e) {
    message("Reactome富集失败:", e$message)
  })
  
  return(enrichment_results)
}

# ==================== 可视化函数 ====================

# 生成保守性热图
generate_conservation_heatmap <- function(conserved_analysis, all_trajectory_genes) {
  message("生成保守性热图...")
  
  # 选择top保守基因
  top_genes <- head(conserved_analysis$gene_symbol, 50)
  
  # 构建相关性矩阵
  correlation_matrix <- matrix(NA, nrow = length(top_genes), 
                              ncol = length(all_trajectory_genes))
  rownames(correlation_matrix) <- top_genes
  colnames(correlation_matrix) <- names(all_trajectory_genes)
  
  for(gene in top_genes) {
    for(cancer_cell in names(all_trajectory_genes)) {
      gene_data <- all_trajectory_genes[[cancer_cell]]
      
      # 转换基因名
      gene_symbols <- convert_ensembl_to_symbol(gene_data$gene_short_name)
      
      if(gene %in% gene_symbols) {
        ensembl_id <- names(gene_symbols)[gene_symbols == gene][1]
        gene_row <- gene_data[gene_data$gene_short_name == ensembl_id, ]
        
        if(nrow(gene_row) > 0) {
          correlation_matrix[gene, cancer_cell] <- gene_row$correlation[1]
        }
      }
    }
  }
  
  # 移除全NA的行和列
  valid_rows <- rowSums(!is.na(correlation_matrix)) > 0
  valid_cols <- colSums(!is.na(correlation_matrix)) > 0
  
  if(sum(valid_rows) > 0 && sum(valid_cols) > 0) {
    correlation_matrix <- correlation_matrix[valid_rows, valid_cols, drop = FALSE]
    
    # 生成热图
    library(ComplexHeatmap)
    library(circlize)
    
    # 颜色设置
    col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
    
    # 行注释（保守性评分）
    gene_scores <- conserved_analysis$conservation_score[match(rownames(correlation_matrix), 
                                                             conserved_analysis$gene_symbol)]
    row_ha <- rowAnnotation(
      Conservation_Score = anno_barplot(gene_scores, 
                                      gp = gpar(fill = "darkgreen", col = "darkgreen")),
      annotation_name_side = "top"
    )
    
    # 列注释（癌种）
    cancer_types <- sapply(strsplit(colnames(correlation_matrix), "_"), function(x) x[1])
    col_ha <- HeatmapAnnotation(
      Cancer_Type = cancer_types,
      col = list(Cancer_Type = rainbow(length(unique(cancer_types))))
    )
    
    # 创建热图
    ht <- Heatmap(
      correlation_matrix,
      name = "Correlation",
      col = col_fun,
      top_annotation = col_ha,
      right_annotation = row_ha,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8),
      heatmap_legend_param = list(title = "Spearman\nCorrelation")
    )
    
    return(ht)
  } else {
    message("没有足够的数据生成热图")
    return(NULL)
  }
}

# 生成保守基因韦恩图
generate_conservation_venn <- function(all_trajectory_genes, p_threshold = 0.05) {
  message("生成保守基因韦恩图...")
  
  # 提取每个癌种的显著基因
  cancer_genes <- list()
  
  for(cancer_cell in names(all_trajectory_genes)) {
    genes <- all_trajectory_genes[[cancer_cell]]
    significant_genes <- genes$gene_short_name[genes$q_value < p_threshold]
    
    # 转换基因名
    gene_symbols <- convert_ensembl_to_symbol(significant_genes)
    valid_symbols <- gene_symbols[!is.na(gene_symbols)]
    
    if(length(valid_symbols) > 0) {
      cancer_type <- strsplit(cancer_cell, "_")[[1]][1]
      
      if(cancer_type %in% names(cancer_genes)) {
        cancer_genes[[cancer_type]] <- union(cancer_genes[[cancer_type]], valid_symbols)
      } else {
        cancer_genes[[cancer_type]] <- valid_symbols
      }
    }
  }
  
  # 生成韦恩图（最多5个集合）
  if(length(cancer_genes) >= 2 && length(cancer_genes) <= 5) {
    library(VennDiagram)
    
    venn_file <- file.path(conserved_output_dir, "comparative_plots", "conserved_genes_venn.png")
    
    venn.diagram(
      x = cancer_genes,
      category.names = names(cancer_genes),
      filename = venn_file,
      output = TRUE,
      imagetype = "png",
      height = 2000,
      width = 2000,
      resolution = 300,
      col = "transparent",
      fill = rainbow(length(cancer_genes), alpha = 0.3),
      cex = 1.5,
      fontfamily = "serif",
      cat.cex = 1.5,
      cat.fontfamily = "serif"
    )
    
    message(paste("韦恩图已保存:", venn_file))
    return(venn_file)
  } else {
    message("癌种数量不适合生成韦恩图")
    return(NULL)
  }
}

# ==================== 主分析函数 ====================

# 综合保守性和特异性分析
comprehensive_conservation_analysis <- function(min_cancer_types = 2,
                                              p_threshold = 0.05,
                                              correlation_threshold = 0.15) {
  
  message("开始综合保守性和特异性分析...")
  
  # 1. 加载数据
  all_results <- load_all_trajectory_results()
  all_trajectory_genes <- load_all_trajectory_genes()
  
  if(length(all_trajectory_genes) < min_cancer_types) {
    stop(paste("可用数据不足，需要至少", min_cancer_types, "个癌种"))
  }
  
  # 2. 识别保守基因
  conserved_results <- identify_conserved_trajectory_genes(
    all_trajectory_genes, 
    min_cancer_types = min_cancer_types,
    p_threshold = p_threshold,
    correlation_threshold = correlation_threshold
  )
  
  # 3. 识别特异基因
  specific_results <- identify_cancer_specific_genes(
    all_trajectory_genes,
    conserved_results$conserved_genes,
    p_threshold = p_threshold,
    correlation_threshold = correlation_threshold
  )
  
  # 4. 轨迹模式比较
  trajectory_patterns <- compare_trajectory_patterns(all_results)
  
  # 5. 功能富集分析
  conserved_functions <- analyze_conserved_gene_functions(conserved_results$conserved_genes)
  
  # 6. 生成可视化
  conservation_heatmap <- generate_conservation_heatmap(
    conserved_results$conserved_analysis, 
    all_trajectory_genes
  )
  
  conservation_venn <- generate_conservation_venn(all_trajectory_genes, p_threshold)
  
  # 7. 保存结果
  final_results <- list(
    conserved_analysis = conserved_results,
    specific_analysis = specific_results,
    trajectory_patterns = trajectory_patterns,
    functional_enrichment = conserved_functions,
    visualization = list(
      heatmap = conservation_heatmap,
      venn_diagram = conservation_venn
    ),
    parameters = list(
      min_cancer_types = min_cancer_types,
      p_threshold = p_threshold,
      correlation_threshold = correlation_threshold
    )
  )
  
  # 保存完整结果
  results_file <- file.path(conserved_output_dir, "comprehensive_conservation_analysis.rds")
  saveRDS(final_results, results_file)
  
  # 保存各部分结果
  write.csv(conserved_results$conserved_analysis, 
            file.path(conserved_output_dir, "conserved_modules", "conserved_genes_analysis.csv"),
            row.names = FALSE)
  
  write.csv(trajectory_patterns,
            file.path(conserved_output_dir, "comparative_plots", "trajectory_patterns_comparison.csv"),
            row.names = FALSE)
  
  # 保存热图
  if(!is.null(conservation_heatmap)) {
    pdf(file.path(conserved_output_dir, "comparative_plots", "conservation_heatmap.pdf"),
        width = 12, height = 10)
    draw(conservation_heatmap)
    dev.off()
  }
  
  # 生成报告
  generate_conservation_report(final_results)
  
  message("保守性和特异性分析完成!")
  message(paste("结果保存在:", conserved_output_dir))
  
  return(final_results)
}

# 生成保守性分析报告
generate_conservation_report <- function(results) {
  
  report_file <- file.path(conserved_output_dir, "conservation_analysis_report.txt")
  
  sink(report_file)
  
  cat("=== 跨癌种T细胞轨迹保守性和特异性分析报告 ===\n")
  cat("分析时间:", as.character(Sys.time()), "\n\n")
  
  # 基本统计
  cat("1. 基本统计信息\n")
  cat("================\n")
  cat("分析的癌种-细胞类型组合数:", length(results$specific_analysis) + 
      length(results$conserved_analysis$detailed_data), "\n")
  cat("发现的保守基因数:", nrow(results$conserved_analysis$conserved_analysis), "\n")
  
  # 保守基因top 10
  cat("\n2. Top 10保守基因\n")
  cat("=================\n")
  top_conserved <- head(results$conserved_analysis$conserved_analysis, 10)
  for(i in 1:nrow(top_conserved)) {
    cat(sprintf("%2d. %s (保守性评分: %.3f, 出现在%d个癌种)\n",
                i, top_conserved$gene_symbol[i], 
                top_conserved$conservation_score[i],
                top_conserved$n_cancer_types[i]))
  }
  
  # 轨迹模式比较
  if(!is.null(results$trajectory_patterns) && nrow(results$trajectory_patterns) > 0) {
    cat("\n3. 轨迹模式比较\n")
    cat("================\n")
    patterns <- results$trajectory_patterns
    
    # 按癌种分组统计
    cancer_summary <- patterns %>%
      group_by(cancer_type) %>%
      summarise(
        n_cell_types = n(),
        mean_cells = mean(n_cells, na.rm = TRUE),
        mean_pseudotime_range = mean(pseudotime_range, na.rm = TRUE),
        mean_clusters = mean(n_clusters, na.rm = TRUE),
        .groups = 'drop'
      )
    
    for(i in 1:nrow(cancer_summary)) {
      cancer <- cancer_summary$cancer_type[i]
      cat(sprintf("%s: %d个细胞类型, 平均%.0f个细胞, 平均轨迹长度%.2f\n",
                  cancer, cancer_summary$n_cell_types[i],
                  cancer_summary$mean_cells[i], 
                  cancer_summary$mean_pseudotime_range[i]))
    }
  }
  
  # 功能富集结果
  if(!is.null(results$functional_enrichment)) {
    cat("\n4. 保守基因功能富集\n")
    cat("====================\n")
    
    if("GO_BP" %in% names(results$functional_enrichment)) {
      go_results <- results$functional_enrichment$GO_BP@result
      if(nrow(go_results) > 0) {
        cat("主要GO生物过程:\n")
        top_go <- head(go_results[go_results$p.adjust < 0.05, ], 5)
        for(i in 1:nrow(top_go)) {
          cat(sprintf("  - %s (p=%.2e)\n", top_go$Description[i], top_go$p.adjust[i]))
        }
      }
    }
    
    if("KEGG" %in% names(results$functional_enrichment)) {
      kegg_results <- results$functional_enrichment$KEGG@result
      if(nrow(kegg_results) > 0) {
        cat("\n主要KEGG通路:\n")
        top_kegg <- head(kegg_results[kegg_results$p.adjust < 0.05, ], 5)
        for(i in 1:nrow(top_kegg)) {
          cat(sprintf("  - %s (p=%.2e)\n", top_kegg$Description[i], top_kegg$p.adjust[i]))
        }
      }
    }
  }
  
  # 特异性基因统计
  cat("\n5. 癌种特异性基因\n")
  cat("==================\n")
  if(length(results$specific_analysis) > 0) {
    for(cancer_cell in names(results$specific_analysis)) {
      n_specific <- nrow(results$specific_analysis[[cancer_cell]])
      cat(sprintf("%s: %d个特异基因\n", cancer_cell, n_specific))
    }
  } else {
    cat("未发现显著的癌种特异性基因\n")
  }
  
  cat("\n=== 报告结束 ===\n")
  sink()
  
  message(paste("分析报告已保存:", report_file))
}

# ==================== 执行分析 ====================

if(!interactive()) {
  message("开始执行保守性和特异性分析...")
  
  # 执行分析
  conservation_results <- comprehensive_conservation_analysis(
    min_cancer_types = 2,
    p_threshold = 0.05,
    correlation_threshold = 0.15
  )
  
  message("保守性和特异性分析完成!")
  
} else {
  message("保守性和特异性分析代码加载完成!")
  message("使用方法:")
  message("conservation_results <- comprehensive_conservation_analysis()")
}

message("跨癌种保守性和特异性分析代码加载完成!")