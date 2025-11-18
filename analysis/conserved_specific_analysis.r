# conserved_specific_analysis_enhanced.R
# 跨癌种T细胞轨迹保守和特异模式分析（增强可视化版 v3.0）

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
library(reshape2)
library(ggrepel)
library(scales)
library(gridExtra)
library(viridis)

# 设置路径
trajectory_output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/trajectory_analysis_2"
conserved_output_dir <- file.path(trajectory_output_dir, "conserved_specific_analysis")
dir.create(conserved_output_dir, showWarnings = FALSE, recursive = TRUE)

# 本地数据库路径
enrich_db_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/enrich-database"

# 创建子目录
sub_dirs <- c("conserved_modules", "specific_modules", "comparative_plots", 
              "network_analysis", "functional_analysis", "detailed_tables",
              "enrichment_visualization", "gene_expression_patterns",
              "interactive_reports")
for(dir in sub_dirs) {
  dir.create(file.path(conserved_output_dir, dir), showWarnings = FALSE)
}

# ==================== 数据加载 ====================

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
        message(paste("✓", result_name))
      }
    }, error = function(e) {
      message(paste("✗", result_name, "-", e$message))
    })
  }
  
  message(paste("共加载", length(all_results), "个结果"))
  return(all_results)
}

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
      }
    }, error = function(e) {
      message(paste("基因加载失败:", gene_name))
    })
  }
  
  message(paste("加载了", length(all_genes), "个基因集"))
  return(all_genes)
}

convert_ensembl_to_symbol <- function(ensembl_ids) {
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
  tryCatch({
    symbols <- suppressMessages(
      mapIds(org.Hs.eg.db, 
             keys = ensembl_ids,
             column = "SYMBOL",
             keytype = "ENSEMBL",
             multiVals = "first")
    )
    
    valid_symbols <- symbols[!is.na(symbols)]
    return(valid_symbols)
    
  }, error = function(e) {
    names(ensembl_ids) <- ensembl_ids
    return(ensembl_ids)
  })
}

# ==================== 保守模块识别（增强版） ====================

identify_conserved_trajectory_genes <- function(all_trajectory_genes, 
                                              min_cancer_types = 2,
                                              p_threshold = 0.05,
                                              correlation_threshold = 0.15) {
  
  message("识别保守轨迹基因（增强版）...")
  
  all_significant_genes <- list()
  
  # 收集所有显著基因
  for(cancer_cell in names(all_trajectory_genes)) {
    genes <- all_trajectory_genes[[cancer_cell]]
    
    significant <- genes[genes$q_value < p_threshold & 
                        abs(genes$correlation) > correlation_threshold, ]
    
    if(nrow(significant) > 0) {
      gene_symbols <- convert_ensembl_to_symbol(significant$gene_short_name)
      
      if(length(gene_symbols) > 0) {
        significant$gene_symbol <- gene_symbols[significant$gene_short_name]
        significant <- significant[!is.na(significant$gene_symbol), ]
        
        all_significant_genes[[cancer_cell]] <- significant
      }
    }
  }
  
  all_gene_symbols <- unlist(lapply(all_significant_genes, function(x) x$gene_symbol))
  gene_frequency <- table(all_gene_symbols)
  
  conserved_genes <- names(gene_frequency)[gene_frequency >= min_cancer_types]
  
  message(paste("发现", length(conserved_genes), "个保守基因"))
  
  # 🆕 详细的保守性分析
  conserved_analysis <- data.frame()
  conserved_details_list <- list()
  
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
          cancer_cell = cancer_cell,
          correlation = gene_row$correlation[1],
          p_value = gene_row$p_value[1],
          q_value = gene_row$q_value[1],
          ensembl_id = gene_row$gene_short_name[1],
          morans_I = ifelse("morans_I" %in% names(gene_row), gene_row$morans_I[1], NA),
          morans_test_statistic = ifelse("morans_test_statistic" %in% names(gene_row), 
                                        gene_row$morans_test_statistic[1], NA)
        ))
      }
    }
    
    if(nrow(gene_info) > 0) {
      # 保存每个基因的详细信息
      conserved_details_list[[gene]] <- gene_info
      
      # 计算多维度保守性指标
      conservation_stats <- data.frame(
        gene_symbol = gene,
        n_cancer_types = length(unique(gene_info$cancer_type)),
        n_cell_contexts = nrow(gene_info),
        
        # 相关性统计
        mean_correlation = mean(gene_info$correlation),
        median_correlation = median(gene_info$correlation),
        sd_correlation = sd(gene_info$correlation),
        min_correlation = min(gene_info$correlation),
        max_correlation = max(gene_info$correlation),
        
        # 方向一致性
        consistent_direction = length(unique(sign(gene_info$correlation))) == 1,
        positive_ratio = sum(gene_info$correlation > 0) / nrow(gene_info),
        
        # 显著性统计
        mean_neg_log10_pval = mean(-log10(gene_info$p_value)),
        mean_neg_log10_qval = mean(-log10(gene_info$q_value)),
        
        # 空间自相关（如果有）
        mean_morans_I = mean(gene_info$morans_I, na.rm = TRUE),
        
        # 综合评分
        conservation_score = nrow(gene_info) * abs(mean(gene_info$correlation)) * 
                           mean(-log10(gene_info$q_value)),
        
        # 稳定性评分
        stability_score = 1 / (1 + sd(gene_info$correlation))
      )
      
      conserved_analysis <- rbind(conserved_analysis, conservation_stats)
    }
  }
  
  conserved_analysis <- conserved_analysis[order(conserved_analysis$conservation_score, 
                                               decreasing = TRUE), ]
  
  return(list(
    conserved_genes = conserved_genes,
    conserved_analysis = conserved_analysis,
    detailed_data = all_significant_genes,
    conserved_details_list = conserved_details_list,
    gene_frequency = gene_frequency
  ))
}

# ==================== 特异性基因识别（增强版） ====================

identify_cancer_specific_genes <- function(all_trajectory_genes, 
                                         conserved_genes,
                                         p_threshold = 0.05,
                                         correlation_threshold = 0.2) {
  
  message("识别癌种特异性基因（增强版）...")
  
  specific_genes <- list()
  specific_summary <- data.frame()
  specific_functional_categories <- list()
  
  for(cancer_cell in names(all_trajectory_genes)) {
    genes <- all_trajectory_genes[[cancer_cell]]
    
    significant <- genes[genes$q_value < p_threshold & 
                        abs(genes$correlation) > correlation_threshold, ]
    
    if(nrow(significant) > 0) {
      gene_symbols <- convert_ensembl_to_symbol(significant$gene_short_name)
      valid_genes <- gene_symbols[!is.na(gene_symbols)]
      
      specific_to_cancer <- setdiff(valid_genes, conserved_genes)
      
      if(length(specific_to_cancer) > 0) {
        specific_data <- significant[significant$gene_short_name %in% 
                                   names(gene_symbols)[gene_symbols %in% specific_to_cancer], ]
        specific_data$gene_symbol <- gene_symbols[specific_data$gene_short_name]
        
        # 🆕 按相关性分类
        specific_data$regulation <- ifelse(specific_data$correlation > 0, "Upregulated", "Downregulated")
        specific_data$strength <- cut(abs(specific_data$correlation), 
                                     breaks = c(0, 0.3, 0.5, 1),
                                     labels = c("Weak", "Moderate", "Strong"))
        
        specific_genes[[cancer_cell]] <- specific_data
        
        parts <- strsplit(cancer_cell, "_")[[1]]
        cancer_type <- parts[1]
        cell_type <- paste(parts[-1], collapse = "_")
        
        # 🆕 更详细的汇总
        summary_row <- data.frame(
          cancer_cell = cancer_cell,
          cancer_type = cancer_type,
          cell_type = cell_type,
          
          n_specific_genes = length(specific_to_cancer),
          n_upregulated = sum(specific_data$correlation > 0),
          n_downregulated = sum(specific_data$correlation < 0),
          
          mean_correlation = mean(specific_data$correlation),
          median_correlation = median(specific_data$correlation),
          max_abs_correlation = max(abs(specific_data$correlation)),
          
          top_upregulated_gene = ifelse(sum(specific_data$correlation > 0) > 0,
                                       specific_data$gene_symbol[which.max(specific_data$correlation)], NA),
          top_up_correlation = ifelse(sum(specific_data$correlation > 0) > 0,
                                     max(specific_data$correlation), NA),
          
          top_downregulated_gene = ifelse(sum(specific_data$correlation < 0) > 0,
                                         specific_data$gene_symbol[which.min(specific_data$correlation)], NA),
          top_down_correlation = ifelse(sum(specific_data$correlation < 0) > 0,
                                       min(specific_data$correlation), NA),
          
          mean_significance = mean(-log10(specific_data$q_value))
        )
        
        specific_summary <- rbind(specific_summary, summary_row)
        
        message(paste(cancer_cell, ":", length(specific_to_cancer), "个特异基因",
                     "(上调:", sum(specific_data$correlation > 0), 
                     "下调:", sum(specific_data$correlation < 0), ")"))
      }
    }
  }
  
  return(list(
    specific_genes = specific_genes,
    specific_summary = specific_summary
  ))
}

# ==================== 功能富集分析（本地版） ====================

load_pathway_databases <- function() {
  message("加载本地pathway数据库...")
  
  db_list <- list()
  
  msigdb_file <- file.path(enrich_db_dir, "msigdb_kegg_genesets.RData")
  if(file.exists(msigdb_file)) {
    load(msigdb_file)
    if(exists("msigdb_data")) {
      db_list$msigdb_kegg <- msigdb_data
      message("✓ MSigDB KEGG")
    }
  }
  
  reactome_file <- file.path(enrich_db_dir, "reactome_genesets.RData")
  if(file.exists(reactome_file)) {
    load(reactome_file)
    if(exists("reactome_data")) {
      db_list$reactome <- reactome_data
      message("✓ Reactome")
    }
  }
  
  return(db_list)
}

manual_go_enrichment <- function(gene_list, ont = "BP", pval_cutoff = 0.05, qval_cutoff = 0.2) {
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  
  go_all <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
      columns = c("ENTREZID", "GO", "ONTOLOGY"),
      keytype = "ENTREZID"
    )
  )
  
  go_all <- go_all[go_all$ONTOLOGY == ont & !is.na(go_all$GO), ]
  
  go_term_genes <- split(go_all$ENTREZID, go_all$GO)
  go_term_genes <- lapply(go_term_genes, unique)
  
  go_sizes <- sapply(go_term_genes, length)
  go_term_genes <- go_term_genes[go_sizes >= 10 & go_sizes <= 500]
  
  if (length(go_term_genes) == 0) return(NULL)
  
  universe <- unique(go_all$ENTREZID)
  N <- length(universe)
  n <- length(gene_list)
  
  results_list <- lapply(names(go_term_genes), function(go_id) {
    term_genes <- go_term_genes[[go_id]]
    M <- length(term_genes)
    overlap <- intersect(gene_list, term_genes)
    k <- length(overlap)
    
    if (k == 0) return(NULL)
    
    pval <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    
    # 🆕 添加富集倍数
    expected <- (n * M) / N
    fold_enrichment <- k / expected
    
    data.frame(
      GO_ID = go_id,
      Count = k,
      Total = M,
      GeneRatio = paste0(k, "/", n),
      BgRatio = paste0(M, "/", N),
      pvalue = pval,
      FoldEnrichment = fold_enrichment,
      geneID = paste(overlap, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results_list)
  if (is.null(results_df) || nrow(results_df) == 0) return(NULL)
  
  results_df$p.adjust <- p.adjust(results_df$pvalue, method = "BH")
  
  go_terms <- suppressMessages(
    AnnotationDbi::select(
      GO.db::GO.db,
      keys = results_df$GO_ID,
      columns = c("GOID", "TERM"),
      keytype = "GOID"
    )
  )
  
  results_df <- merge(results_df, go_terms, by.x = "GO_ID", by.y = "GOID", all.x = TRUE)
  results_df <- results_df[results_df$pvalue < pval_cutoff & results_df$p.adjust < qval_cutoff, ]
  results_df <- results_df[order(results_df$p.adjust), ]
  
  return(results_df)
}

manual_pathway_enrichment <- function(gene_list, term2gene, term2name, pval_cutoff = 0.05, qval_cutoff = 0.2) {
  
  if (!is.data.frame(term2gene) || ncol(term2gene) < 2) return(NULL)
  
  colnames(term2gene)[1:2] <- c("pathway", "gene")
  
  pathway_genes <- split(term2gene$gene, term2gene$pathway)
  pathway_genes <- lapply(pathway_genes, unique)
  
  pathway_sizes <- sapply(pathway_genes, length)
  pathway_genes <- pathway_genes[pathway_sizes >= 10 & pathway_sizes <= 500]
  
  if (length(pathway_genes) == 0) return(NULL)
  
  universe <- unique(term2gene$gene)
  N <- length(universe)
  n <- length(gene_list)
  
  results_list <- lapply(names(pathway_genes), function(pathway_id) {
    pathway_genes_set <- pathway_genes[[pathway_id]]
    M <- length(pathway_genes_set)
    overlap <- intersect(gene_list, pathway_genes_set)
    k <- length(overlap)
    
    if (k == 0) return(NULL)
    
    pval <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    
    # 🆕 添加富集倍数
    expected <- (n * M) / N
    fold_enrichment <- k / expected
    
    data.frame(
      Pathway_ID = pathway_id,
      Count = k,
      Total = M,
      GeneRatio = paste0(k, "/", n),
      BgRatio = paste0(M, "/", N),
      pvalue = pval,
      FoldEnrichment = fold_enrichment,
      geneID = paste(overlap, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results_list)
  if (is.null(results_df) || nrow(results_df) == 0) return(NULL)
  
  results_df$p.adjust <- p.adjust(results_df$pvalue, method = "BH")
  
  if (!is.null(term2name) && is.data.frame(term2name) && ncol(term2name) >= 2) {
    colnames(term2name)[1:2] <- c("pathway", "description")
    results_df <- merge(results_df, term2name, by.x = "Pathway_ID", by.y = "pathway", all.x = TRUE)
  } else {
    results_df$description <- results_df$Pathway_ID
  }
  
  results_df <- results_df[results_df$pvalue < pval_cutoff & results_df$p.adjust < qval_cutoff, ]
  results_df <- results_df[order(results_df$p.adjust), ]
  
  return(results_df)
}

analyze_conserved_gene_functions <- function(conserved_genes) {
  message("分析保守基因功能...")
  
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
  if(length(conserved_genes) < 5) {
    message("保守基因数量不足")
    return(NULL)
  }
  
  entrez_ids <- tryCatch({
    suppressMessages(
      mapIds(
        org.Hs.eg.db,
        keys = conserved_genes,
        keytype = "SYMBOL",
        column = "ENTREZID",
        multiVals = "first"
      )
    )
  }, error = function(e) NULL)
  
  if(is.null(entrez_ids)) return(NULL)
  
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  entrez_ids <- unique(as.character(entrez_ids))
  
  if(length(entrez_ids) < 5) return(NULL)
  
  message(paste("成功转换", length(entrez_ids), "个基因ID"))
  
  enrichment_results <- list()
  
  # GO分析
  for(ont in c("BP", "MF", "CC")) {
    tryCatch({
      go_result <- manual_go_enrichment(entrez_ids, ont = ont)
      if(!is.null(go_result) && nrow(go_result) > 0) {
        enrichment_results[[paste0("GO_", ont)]] <- go_result
        message(paste("GO", ont, ":", nrow(go_result), "个通路"))
      }
    }, error = function(e) message(paste("GO", ont, "失败")))
  }
  
  # 本地数据库富集
  local_dbs <- load_pathway_databases()
  
  if("msigdb_kegg" %in% names(local_dbs)) {
    tryCatch({
      kegg_result <- manual_pathway_enrichment(
        gene_list = entrez_ids,
        term2gene = local_dbs$msigdb_kegg$term2gene,
        term2name = local_dbs$msigdb_kegg$term2name
      )
      
      if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
        enrichment_results$KEGG <- kegg_result
        message(paste("KEGG:", nrow(kegg_result), "个通路"))
      }
    }, error = function(e) message("KEGG失败"))
  }
  
  if("reactome" %in% names(local_dbs)) {
    tryCatch({
      reactome_result <- manual_pathway_enrichment(
        gene_list = entrez_ids,
        term2gene = local_dbs$reactome$term2gene,
        term2name = local_dbs$reactome$term2name
      )
      
      if(!is.null(reactome_result) && nrow(reactome_result) > 0) {
        enrichment_results$Reactome <- reactome_result
        message(paste("Reactome:", nrow(reactome_result), "个通路"))
      }
    }, error = function(e) message("Reactome失败"))
  }
  
  return(enrichment_results)
}

# 🆕 特异性基因功能分析
analyze_specific_gene_functions <- function(specific_genes_list) {
  message("分析特异性基因功能...")
  
  library(org.Hs.eg.db)
  
  specific_enrichment <- list()
  
  for(cancer_cell in names(specific_genes_list)) {
    gene_data <- specific_genes_list[[cancer_cell]]
    genes <- unique(gene_data$gene_symbol)
    
    if(length(genes) < 5) next
    
    entrez_ids <- tryCatch({
      suppressMessages(
        mapIds(org.Hs.eg.db,
               keys = genes,
               keytype = "SYMBOL",
               column = "ENTREZID",
               multiVals = "first")
      )
    }, error = function(e) NULL)
    
    if(is.null(entrez_ids)) next
    
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    entrez_ids <- unique(as.character(entrez_ids))
    
    if(length(entrez_ids) < 5) next
    
    # GO BP分析
    go_bp <- tryCatch({
      manual_go_enrichment(entrez_ids, ont = "BP", qval_cutoff = 0.25)
    }, error = function(e) NULL)
    
    if(!is.null(go_bp) && nrow(go_bp) > 0) {
      specific_enrichment[[cancer_cell]] <- list(GO_BP = go_bp)
      message(paste(cancer_cell, ":", nrow(go_bp), "个GO通路"))
    }
  }
  
  return(specific_enrichment)
}

# ==================== 🎨 高级可视化函数 ====================

# 1. 富集分析气泡图
plot_enrichment_bubble <- function(enrichment_results, database = "GO_BP", top_n = 20) {
  
  if(is.null(enrichment_results) || !database %in% names(enrichment_results)) {
    message(paste("未找到", database, "数据"))
    return(NULL)
  }
  
  enrich_data <- enrichment_results[[database]]
  
  if(nrow(enrich_data) == 0) return(NULL)
  
  plot_data <- head(enrich_data, top_n)
  
  # 准备数据
  if("TERM" %in% names(plot_data)) {
    plot_data$Description <- plot_data$TERM
  } else if("description" %in% names(plot_data)) {
    plot_data$Description <- plot_data$description
  } else {
    plot_data$Description <- ifelse("GO_ID" %in% names(plot_data), 
                                   plot_data$GO_ID, plot_data$Pathway_ID)
  }
  
  # 截断过长的描述
  plot_data$Description <- ifelse(nchar(plot_data$Description) > 60,
                                 paste0(substr(plot_data$Description, 1, 57), "..."),
                                 plot_data$Description)
  
  plot_data$neg_log10_pval <- -log10(plot_data$p.adjust)
  
  # 提取GeneRatio数值
  if("GeneRatio" %in% names(plot_data)) {
    gene_ratio_parts <- strsplit(as.character(plot_data$GeneRatio), "/")
    plot_data$GeneRatio_numeric <- sapply(gene_ratio_parts, function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
  } else {
    plot_data$GeneRatio_numeric <- plot_data$Count / sum(plot_data$Count)
  }
  
  # 气泡图
  p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = reorder(Description, GeneRatio_numeric))) +
    geom_point(aes(size = Count, color = neg_log10_pval)) +
    scale_color_gradient(low = "blue", high = "red", 
                        name = "-log10(FDR)") +
    scale_size_continuous(range = c(3, 10), name = "Gene Count") +
    labs(title = paste(database, "Enrichment Analysis"),
         x = "Gene Ratio",
         y = "") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  return(p)
}

# 2. 富集分析条形图（带富集倍数）
plot_enrichment_barplot <- function(enrichment_results, database = "GO_BP", top_n = 20) {
  
  if(is.null(enrichment_results) || !database %in% names(enrichment_results)) {
    return(NULL)
  }
  
  enrich_data <- enrichment_results[[database]]
  
  if(nrow(enrich_data) == 0) return(NULL)
  
  plot_data <- head(enrich_data, top_n)
  
  if("TERM" %in% names(plot_data)) {
    plot_data$Description <- plot_data$TERM
  } else if("description" %in% names(plot_data)) {
    plot_data$Description <- plot_data$description
  } else {
    plot_data$Description <- ifelse("GO_ID" %in% names(plot_data), 
                                   plot_data$GO_ID, plot_data$Pathway_ID)
  }
  
  plot_data$Description <- ifelse(nchar(plot_data$Description) > 60,
                                 paste0(substr(plot_data$Description, 1, 57), "..."),
                                 plot_data$Description)
  
  plot_data$neg_log10_pval <- -log10(plot_data$p.adjust)
  
  p <- ggplot(plot_data, aes(x = reorder(Description, neg_log10_pval), 
                            y = neg_log10_pval)) +
    geom_bar(stat = "identity", aes(fill = FoldEnrichment)) +
    scale_fill_viridis_c(option = "plasma", direction = -1, 
                        name = "Fold Enrichment") +
    coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(title = paste(database, "Enrichment"),
         x = "",
         y = "-log10(Adjusted P-value)") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  return(p)
}

# 3. GO层级树状图
plot_go_hierarchy <- function(enrichment_results, top_n = 15) {
  
  if(is.null(enrichment_results) || !"GO_BP" %in% names(enrichment_results)) {
    return(NULL)
  }
  
  go_data <- head(enrichment_results$GO_BP, top_n)
  
  if(nrow(go_data) == 0) return(NULL)
  
  library(igraph)
  library(ggraph)
  
  # 简化的层级关系（基于p值排序）
  go_data$Level <- cut(1:nrow(go_data), breaks = 3, labels = c("High", "Medium", "Low"))
  go_data$Size <- -log10(go_data$p.adjust)
  
  # 创建边
  edges <- data.frame()
  for(i in 2:nrow(go_data)) {
    edges <- rbind(edges, data.frame(
      from = "Root",
      to = go_data$TERM[i],
      weight = go_data$Size[i]
    ))
  }
  
  if(nrow(edges) == 0) return(NULL)
  
  # 创建图
  g <- graph_from_data_frame(edges, directed = TRUE)
  
  # 顶点属性
  V(g)$size <- c(10, go_data$Size[-1])
  V(g)$color <- c("gray", as.numeric(factor(go_data$Level[-1])))
  
  p <- ggraph(g, layout = 'tree') +
    geom_edge_link(aes(width = weight), alpha = 0.5, color = "gray50") +
    geom_node_point(aes(size = size, color = factor(color))) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_size_continuous(range = c(3, 10)) +
    scale_color_viridis_d() +
    theme_void() +
    labs(title = "GO Term Hierarchy") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  return(p)
}

# 4. 富集通路网络图
plot_pathway_network <- function(enrichment_results, database = "KEGG", top_n = 20) {
  
  if(is.null(enrichment_results) || !database %in% names(enrichment_results)) {
    return(NULL)
  }
  
  pathway_data <- head(enrichment_results[[database]], top_n)
  
  if(nrow(pathway_data) == 0) return(NULL)
  
  library(igraph)
  library(ggraph)
  
  # 创建基因-通路关系
  edges <- data.frame()
  
  for(i in 1:nrow(pathway_data)) {
    pathway_name <- ifelse(!is.na(pathway_data$description[i]), 
                          pathway_data$description[i], 
                          pathway_data$Pathway_ID[i])
    
    # 截断名称
    pathway_name <- ifelse(nchar(pathway_name) > 40,
                          paste0(substr(pathway_name, 1, 37), "..."),
                          pathway_name)
    
    genes <- strsplit(pathway_data$geneID[i], "/")[[1]]
    
    if(length(genes) > 0) {
      for(gene in head(genes, 5)) {  # 每个通路最多显示5个基因
        edges <- rbind(edges, data.frame(
          from = pathway_name,
          to = gene,
          weight = -log10(pathway_data$p.adjust[i])
        ))
      }
    }
  }
  
  if(nrow(edges) == 0) return(NULL)
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # 顶点类型
  V(g)$type <- ifelse(V(g)$name %in% unique(edges$from), "Pathway", "Gene")
  V(g)$size <- ifelse(V(g)$type == "Pathway", 8, 4)
  V(g)$color <- ifelse(V(g)$type == "Pathway", "coral", "lightblue")
  
  p <- ggraph(g, layout = 'fr') +
    geom_edge_link(aes(width = weight), alpha = 0.3, color = "gray70") +
    geom_node_point(aes(size = size, color = color)) +
    geom_node_text(aes(label = name, filter = type == "Pathway"), 
                  repel = TRUE, size = 3, fontface = "bold") +
    scale_size_identity() +
    scale_color_identity() +
    theme_void() +
    labs(title = paste(database, "Pathway-Gene Network")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  return(p)
}

# 5. 热图：多数据库富集比较
plot_enrichment_heatmap <- function(enrichment_results, top_n = 15) {
  
  # 收集所有数据库的top通路
  all_pathways <- list()
  
  for(db_name in names(enrichment_results)) {
    db_data <- enrichment_results[[db_name]]
    
    if(nrow(db_data) == 0) next
    
    top_pathways <- head(db_data, top_n)
    
    pathway_names <- if("TERM" %in% names(top_pathways)) {
      top_pathways$TERM
    } else if("description" %in% names(top_pathways)) {
      top_pathways$description
    } else {
      top_pathways$Pathway_ID
    }
    
    pathway_names <- ifelse(nchar(pathway_names) > 50,
                           paste0(substr(pathway_names, 1, 47), "..."),
                           pathway_names)
    
    all_pathways[[db_name]] <- data.frame(
      Pathway = pathway_names,
      Database = db_name,
      NegLog10Pval = -log10(top_pathways$p.adjust),
      Count = top_pathways$Count
    )
  }
  
  if(length(all_pathways) == 0) return(NULL)
  
  combined_data <- do.call(rbind, all_pathways)
  
  p <- ggplot(combined_data, aes(x = Database, y = Pathway, fill = NegLog10Pval)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "white", mid = "yellow", high = "red",
                        midpoint = 2, name = "-log10(FDR)") +
    geom_text(aes(label = Count), size = 3) +
    theme_minimal() +
    labs(title = "Enrichment Comparison Across Databases",
         x = "Database",
         y = "Pathway") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  return(p)
}

# 6. 基因表达模式聚类热图
plot_gene_expression_clusters <- function(conserved_details_list, top_genes = 50) {
  
  message("生成基因表达模式聚类热图...")
  
  # 选择top基因
  if(length(conserved_details_list) == 0) return(NULL)
  
  gene_names <- names(conserved_details_list)
  if(length(gene_names) > top_genes) {
    gene_names <- gene_names[1:top_genes]
  }
  
  # 构建表达矩阵
  all_contexts <- unique(unlist(lapply(conserved_details_list, function(x) x$cancer_cell)))
  
  expr_matrix <- matrix(NA, nrow = length(gene_names), ncol = length(all_contexts))
  rownames(expr_matrix) <- gene_names
  colnames(expr_matrix) <- all_contexts
  
  for(gene in gene_names) {
    gene_info <- conserved_details_list[[gene]]
    for(i in 1:nrow(gene_info)) {
      expr_matrix[gene, gene_info$cancer_cell[i]] <- gene_info$correlation[i]
    }
  }
  
  # 移除全NA行列
  valid_rows <- rowSums(!is.na(expr_matrix)) > 0
  valid_cols <- colSums(!is.na(expr_matrix)) > 0
  
  expr_matrix <- expr_matrix[valid_rows, valid_cols, drop = FALSE]
  expr_matrix[is.na(expr_matrix)] <- 0
  
  if(nrow(expr_matrix) < 2 || ncol(expr_matrix) < 2) return(NULL)
  
  # 癌种注释
  cancer_types <- sapply(strsplit(colnames(expr_matrix), "_"), function(x) x[1])
  annotation_col <- data.frame(
    Cancer = cancer_types,
    row.names = colnames(expr_matrix)
  )
  
  cancer_colors <- setNames(rainbow(length(unique(cancer_types))), unique(cancer_types))
  
  p <- pheatmap(
    expr_matrix,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-0.5, 0.5, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    annotation_col = annotation_col,
    annotation_colors = list(Cancer = cancer_colors),
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    main = "Conserved Gene Expression Patterns"
  )
  
  return(p)
}

# 7. 保守基因雷达图（多维度评估）
plot_conservation_radar <- function(conserved_analysis, top_n = 10) {
  
  top_genes <- head(conserved_analysis, top_n)
  
  # 标准化指标到0-1
  normalize <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  
  radar_data <- data.frame(
    Gene = top_genes$gene_symbol,
    Prevalence = normalize(top_genes$n_cell_contexts),
    Correlation = normalize(abs(top_genes$mean_correlation)),
    Significance = normalize(top_genes$mean_neg_log10_qval),
    Stability = normalize(top_genes$stability_score),
    Conservation = normalize(top_genes$conservation_score)
  )
  
  # 转换为长格式
  radar_long <- reshape2::melt(radar_data, id.vars = "Gene")
  
  library(ggradar)
  
  # 如果没有ggradar，用基础ggplot绘制
  p <- ggplot(radar_long, aes(x = variable, y = value, group = Gene, color = Gene)) +
    geom_polygon(fill = NA, size = 1) +
    geom_point(size = 2) +
    coord_polar() +
    ylim(0, 1) +
    labs(title = "Multi-dimensional Conservation Assessment",
         x = "", y = "") +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  return(p)
}

# 8. 特异性基因比较（癌种间）
plot_specific_genes_comparison <- function(specific_summary) {
  
  if(is.null(specific_summary) || nrow(specific_summary) == 0) return(NULL)
  
  # 重塑数据
  plot_data <- specific_summary[, c("cancer_cell", "n_upregulated", "n_downregulated")]
  plot_long <- reshape2::melt(plot_data, id.vars = "cancer_cell")
  
  plot_long$variable <- factor(plot_long$variable, 
                               levels = c("n_upregulated", "n_downregulated"),
                               labels = c("Upregulated", "Downregulated"))
  
  p <- ggplot(plot_long, aes(x = reorder(cancer_cell, value), y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    coord_flip() +
    scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
    labs(title = "Cancer-Specific Gene Regulation",
         x = "Cancer Type - Cell Type",
         y = "Number of Genes",
         fill = "Regulation") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  return(p)
}

# 9. 保守性分布小提琴图
plot_conservation_distribution <- function(conserved_analysis) {
  
  # 准备数据
  dist_data <- data.frame(
    Metric = rep(c("Correlation", "Significance", "Stability"), each = nrow(conserved_analysis)),
    Value = c(
      abs(conserved_analysis$mean_correlation),
      conserved_analysis$mean_neg_log10_qval,
      conserved_analysis$stability_score
    ),
    Gene = rep(conserved_analysis$gene_symbol, 3)
  )
  
  # 标准化
  dist_data <- dist_data %>%
    group_by(Metric) %>%
    mutate(Value_scaled = (Value - min(Value)) / (max(Value) - min(Value))) %>%
    ungroup()
  
  p <- ggplot(dist_data, aes(x = Metric, y = Value_scaled, fill = Metric)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Distribution of Conservation Metrics",
         x = "Metric",
         y = "Scaled Value") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  return(p)
}

# 10. 基因频率柱状图
plot_gene_frequency <- function(gene_frequency, top_n = 30) {
  
  freq_df <- data.frame(
    Gene = names(gene_frequency),
    Frequency = as.numeric(gene_frequency)
  )
  
  freq_df <- freq_df[order(freq_df$Frequency, decreasing = TRUE), ]
  top_freq <- head(freq_df, top_n)
  
  p <- ggplot(top_freq, aes(x = reorder(Gene, Frequency), y = Frequency)) +
    geom_bar(stat = "identity", aes(fill = Frequency)) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    coord_flip() +
    labs(title = "Gene Frequency Across Cancer Types",
         x = "Gene",
         y = "Number of Cancer-Cell Contexts") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  return(p)
}

# ==================== 🆕 综合可视化面板 ====================

generate_comprehensive_visualization <- function(results) {
  
  message("生成综合可视化...")
  
  viz_dir <- file.path(conserved_output_dir, "enrichment_visualization")
  
  # 1. 富集分析可视化
  if(!is.null(results$functional_enrichment)) {
    
    for(db_name in names(results$functional_enrichment)) {
      
      # 气泡图
      p_bubble <- plot_enrichment_bubble(results$functional_enrichment, db_name, top_n = 20)
      if(!is.null(p_bubble)) {
        ggsave(file.path(viz_dir, paste0(db_name, "_bubble.pdf")),
               p_bubble, width = 10, height = 8)
      }
      
      # 条形图
      p_bar <- plot_enrichment_barplot(results$functional_enrichment, db_name, top_n = 20)
      if(!is.null(p_bar)) {
        ggsave(file.path(viz_dir, paste0(db_name, "_barplot.pdf")),
               p_bar, width = 10, height = 8)
      }
    }
    
    # 多数据库热图
    p_heatmap <- plot_enrichment_heatmap(results$functional_enrichment, top_n = 15)
    if(!is.null(p_heatmap)) {
      ggsave(file.path(viz_dir, "enrichment_comparison_heatmap.pdf"),
             p_heatmap, width = 12, height = 10)
    }
    
    # GO层级图
    p_hierarchy <- plot_go_hierarchy(results$functional_enrichment, top_n = 15)
    if(!is.null(p_hierarchy)) {
      ggsave(file.path(viz_dir, "GO_hierarchy.pdf"),
             p_hierarchy, width = 12, height = 10)
    }
    
    # 通路网络图
    for(db in c("KEGG", "Reactome")) {
      p_network <- plot_pathway_network(results$functional_enrichment, db, top_n = 15)
      if(!is.null(p_network)) {
        ggsave(file.path(viz_dir, paste0(db, "_network.pdf")),
               p_network, width = 12, height = 10)
      }
    }
  }
  
  # 2. 基因表达模式
  expr_dir <- file.path(conserved_output_dir, "gene_expression_patterns")
  
  p_cluster <- plot_gene_expression_clusters(
    results$conserved_analysis$conserved_details_list, 
    top_genes = 50
  )
  if(!is.null(p_cluster)) {
    pdf(file.path(expr_dir, "gene_expression_clusters.pdf"), width = 14, height = 12)
    print(p_cluster)
    dev.off()
  }
  
  # 3. 保守性评估
  p_radar <- plot_conservation_radar(results$conserved_analysis$conserved_analysis, top_n = 10)
  if(!is.null(p_radar)) {
    ggsave(file.path(expr_dir, "conservation_radar.pdf"),
           p_radar, width = 10, height = 8)
  }
  
  p_dist <- plot_conservation_distribution(results$conserved_analysis$conserved_analysis)
  if(!is.null(p_dist)) {
    ggsave(file.path(expr_dir, "conservation_distribution.pdf"),
           p_dist, width = 10, height = 6)
  }
  
  # 4. 基因频率
  p_freq <- plot_gene_frequency(results$conserved_analysis$gene_frequency, top_n = 30)
  if(!is.null(p_freq)) {
    ggsave(file.path(expr_dir, "gene_frequency.pdf"),
           p_freq, width = 10, height = 8)
  }
  
  # 5. 特异性基因比较
  p_specific <- plot_specific_genes_comparison(results$specific_analysis$specific_summary)
  if(!is.null(p_specific)) {
    ggsave(file.path(viz_dir, "specific_genes_comparison.pdf"),
           p_specific, width = 10, height = 8)
  }
  
  message("✓ 综合可视化完成")
}

# ==================== 保存详细结果 ====================

save_detailed_results <- function(results) {
  message("保存详细结果...")
  
  detail_dir <- file.path(conserved_output_dir, "detailed_tables")
  
  # 1. 保守基因详细信息（合并所有）
  conserved_detail_file <- file.path(detail_dir, "conserved_genes_detailed.csv")
  
  all_conserved_details <- do.call(rbind, lapply(names(results$conserved_analysis$conserved_details_list), function(gene) {
    gene_data <- results$conserved_analysis$conserved_details_list[[gene]]
    gene_data$gene_symbol <- gene
    return(gene_data)
  }))
  
  write.csv(all_conserved_details, conserved_detail_file, row.names = FALSE)
  message(paste("✓ 保守基因详细:", nrow(all_conserved_details), "条记录"))
  
  # 2. 每个保守基因单独文件
  conserved_genes_dir <- file.path(conserved_output_dir, "conserved_modules", "individual_genes")
  dir.create(conserved_genes_dir, showWarnings = FALSE)
  
  for(gene in names(results$conserved_analysis$conserved_details_list)) {
    gene_file <- file.path(conserved_genes_dir, paste0(gene, "_details.csv"))
    write.csv(results$conserved_analysis$conserved_details_list[[gene]], 
             gene_file, row.names = FALSE)
  }
  message(paste("✓", length(results$conserved_analysis$conserved_details_list), "个基因详细文件"))
  
  # 3. 特异性基因（每个癌种）
  specific_dir <- file.path(conserved_output_dir, "specific_modules")
  
  for(cancer_cell in names(results$specific_analysis$specific_genes)) {
    specific_file <- file.path(specific_dir, paste0(cancer_cell, "_specific_genes.csv"))
    write.csv(results$specific_analysis$specific_genes[[cancer_cell]], 
             specific_file, row.names = FALSE)
  }
  message(paste("✓", length(results$specific_analysis$specific_genes), "个特异性文件"))
  
  # 4. 特异性汇总
  if(!is.null(results$specific_analysis$specific_summary)) {
    summary_file <- file.path(specific_dir, "specific_genes_summary.csv")
    write.csv(results$specific_analysis$specific_summary, summary_file, row.names = FALSE)
  }
  
  # 5. 基因频率
  freq_file <- file.path(detail_dir, "gene_frequency.csv")
  freq_df <- data.frame(
    gene = names(results$conserved_analysis$gene_frequency),
    frequency = as.numeric(results$conserved_analysis$gene_frequency)
  )
  freq_df <- freq_df[order(freq_df$frequency, decreasing = TRUE), ]
  write.csv(freq_df, freq_file, row.names = FALSE)
  
  # 6. 保守性统计
  conserved_stats_file <- file.path(detail_dir, "conservation_statistics.csv")
  write.csv(results$conserved_analysis$conserved_analysis, 
           conserved_stats_file, row.names = FALSE)
  
  message("所有详细结果已保存")
}

# ==================== 生成HTML报告 ====================

generate_html_report <- function(results) {
  
  html_file <- file.path(conserved_output_dir, "interactive_reports", "analysis_report.html")
  
  html_content <- paste0('
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Trajectory Conservation Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }
        .container { max-width: 1200px; margin: auto; background: white; padding: 30px; border-radius: 10px; }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        .summary-box { background: #ecf0f1; padding: 15px; border-radius: 5px; margin: 20px 0; }
        table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        th, td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background-color: #3498db; color: white; }
        tr:hover { background-color: #f5f5f5; }
        .metric { display: inline-block; margin: 10px 20px; }
        .metric-value { font-size: 24px; font-weight: bold; color: #e74c3c; }
        .metric-label { color: #7f8c8d; }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Trajectory Conservation Analysis Report</h1>
        
        <div class="summary-box">
            <h2>📊 Summary Statistics</h2>
            <div class="metric">
                <div class="metric-value">', nrow(results$conserved_analysis$conserved_analysis), '</div>
                <div class="metric-label">Conserved Genes</div>
            </div>
            <div class="metric">
                <div class="metric-value">', length(results$specific_analysis$specific_genes), '</div>
                <div class="metric-label">Cancer-Specific Contexts</div>
            </div>
            <div class="metric">
                <div class="metric-value">', 
                ifelse(!is.null(results$functional_enrichment), 
                      sum(sapply(results$functional_enrichment, nrow)), 0), '</div>
                <div class="metric-label">Enriched Pathways</div>
            </div>
        </div>
        
        <h2>🔝 Top 10 Conserved Genes</h2>
        <table>
            <tr>
                <th>Rank</th>
                <th>Gene</th>
                <th>Contexts</th>
                <th>Mean Correlation</th>
                <th>Conservation Score</th>
            </tr>')
  
  top_genes <- head(results$conserved_analysis$conserved_analysis, 10)
  for(i in 1:nrow(top_genes)) {
    html_content <- paste0(html_content, '
            <tr>
                <td>', i, '</td>
                <td><b>', top_genes$gene_symbol[i], '</b></td>
                <td>', top_genes$n_cell_contexts[i], '</td>
                <td>', sprintf("%.3f", top_genes$mean_correlation[i]), '</td>
                <td>', sprintf("%.2f", top_genes$conservation_score[i]), '</td>
            </tr>')
  }
  
  html_content <- paste0(html_content, '
        </table>
        
        <h2>📈 Enrichment Analysis</h2>')
  
  if(!is.null(results$functional_enrichment)) {
    for(db_name in names(results$functional_enrichment)) {
      db_data <- results$functional_enrichment[[db_name]]
      if(nrow(db_data) > 0) {
        html_content <- paste0(html_content, '
        <h3>', db_name, '</h3>
        <table>
            <tr>
                <th>Pathway</th>
                <th>Count</th>
                <th>FDR</th>
                <th>Fold Enrichment</th>
            </tr>')
        
        top_pathways <- head(db_data, 5)
        for(i in 1:nrow(top_pathways)) {
          pathway_name <- if("TERM" %in% names(top_pathways)) {
            top_pathways$TERM[i]
          } else if("description" %in% names(top_pathways)) {
            top_pathways$description[i]
          } else {
            top_pathways$Pathway_ID[i]
          }
          
          html_content <- paste0(html_content, '
            <tr>
                <td>', pathway_name, '</td>
                <td>', top_pathways$Count[i], '</td>
                <td>', sprintf("%.2e", top_pathways$p.adjust[i]), '</td>
                <td>', sprintf("%.2f", top_pathways$FoldEnrichment[i]), '</td>
            </tr>')
        }
        
        html_content <- paste0(html_content, '
        </table>')
      }
    }
  }
  
  html_content <- paste0(html_content, '
        
        <h2>📁 Output Files</h2>
        <ul>
            <li>Conserved genes: <code>conserved_modules/conserved_genes_analysis.csv</code></li>
            <li>Specific genes: <code>specific_modules/*.csv</code></li>
            <li>Detailed tables: <code>detailed_tables/</code></li>
            <li>Visualizations: <code>enrichment_visualization/</code></li>
            <li>Gene patterns: <code>gene_expression_patterns/</code></li>
        </ul>
        
        <hr>
        <p style="text-align: center; color: #7f8c8d;">
            Generated: ', as.character(Sys.time()), '
        </p>
    </div>
</body>
</html>')
  
  writeLines(html_content, html_file)
  message(paste("✓ HTML报告:", html_file))
}

# ==================== 主分析函数（完整版） ====================

comprehensive_conservation_analysis <- function(min_cancer_types = 2,
                                              p_threshold = 0.05,
                                              correlation_threshold = 0.15) {
  
  message("\n========== 综合保守性和特异性分析（增强版） ==========\n")
  
  # 1. 加载数据
  all_results <- load_all_trajectory_results()
  all_trajectory_genes <- load_all_trajectory_genes()
  
  if(length(all_trajectory_genes) < min_cancer_types) {
    stop(paste("数据不足，需要至少", min_cancer_types, "个癌种"))
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
  
  # 4. 功能富集分析
  conserved_functions <- analyze_conserved_gene_functions(conserved_results$conserved_genes)
  
  # 🆕 5. 特异性基因功能分析
  specific_functions <- analyze_specific_gene_functions(specific_results$specific_genes)
  
  # 6. 整合结果
  final_results <- list(
    conserved_analysis = conserved_results,
    specific_analysis = specific_results,
    functional_enrichment = conserved_functions,
    specific_enrichment = specific_functions,
    parameters = list(
      min_cancer_types = min_cancer_types,
      p_threshold = p_threshold,
      correlation_threshold = correlation_threshold
    )
  )
  
  # 7. 保存RDS
  results_file <- file.path(conserved_output_dir, "comprehensive_conservation_analysis.rds")
  saveRDS(final_results, results_file)
  message(paste("✓ RDS结果:", results_file))
  
  # 8. 保存详细表格
  save_detailed_results(final_results)
  
  # 9. 保存基本CSV
  write.csv(conserved_results$conserved_analysis, 
            file.path(conserved_output_dir, "conserved_modules", "conserved_genes_analysis.csv"),
            row.names = FALSE)
  
  # 10. 保存富集结果
  if(!is.null(conserved_functions)) {
    for(db_name in names(conserved_functions)) {
      if(is.data.frame(conserved_functions[[db_name]])) {
        write.csv(conserved_functions[[db_name]],
                  file.path(conserved_output_dir, "functional_analysis", 
                           paste0(db_name, "_enrichment.csv")),
                  row.names = FALSE)
      }
    }
  }
  
  # 🆕 11. 保存特异性富集
  if(!is.null(specific_functions) && length(specific_functions) > 0) {
    specific_enrich_dir <- file.path(conserved_output_dir, "functional_analysis", "specific_enrichment")
    dir.create(specific_enrich_dir, showWarnings = FALSE)
    
    for(cancer_cell in names(specific_functions)) {
      if("GO_BP" %in% names(specific_functions[[cancer_cell]])) {
        write.csv(specific_functions[[cancer_cell]]$GO_BP,
                 file.path(specific_enrich_dir, paste0(cancer_cell, "_GO_BP.csv")),
                 row.names = FALSE)
      }
    }
  }
  
  # 🆕 12. 生成所有可视化
  generate_comprehensive_visualization(final_results)
  
  # 🆕 13. 生成HTML报告
  generate_html_report(final_results)
  
  # 14. 生成文本报告
  generate_conservation_report(final_results)
  
  message("\n========== 分析完成 ==========")
  message(paste("结果目录:", conserved_output_dir))
  message("\n主要输出:")
  message("  • conserved_modules/ - 保守基因")
  message("  • specific_modules/ - 特异性基因")
  message("  • detailed_tables/ - 详细数据表")
  message("  • enrichment_visualization/ - 富集分析图")
  message("  • gene_expression_patterns/ - 基因表达模式图")
  message("  • interactive_reports/ - HTML交互报告")
  
  return(final_results)
}

# ==================== 文本报告 ====================

generate_conservation_report <- function(results) {
  
  report_file <- file.path(conserved_output_dir, "conservation_analysis_report.txt")
  
  sink(report_file)
  
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║   跨癌种T细胞轨迹保守性和特异性分析报告（增强版）          ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  
  cat("分析时间:", as.character(Sys.time()), "\n")
  cat(rep("=", 70), "\n\n")
  
  # 1. 概览
  cat("【1】分析概览\n")
  cat(rep("-", 70), "\n")
  cat(sprintf("保守基因总数: %d\n", nrow(results$conserved_analysis$conserved_analysis)))
  cat(sprintf("特异性基因组合: %d\n", length(results$specific_analysis$specific_genes)))
  cat(sprintf("总体分析上下文: %d\n", 
             length(unique(unlist(lapply(results$conserved_analysis$detailed_data, 
                                        function(x) unique(x$cancer_cell)))))))
  cat("\n")
  
  # 2. Top保守基因
  cat("【2】Top 20 保守基因\n")
  cat(rep("-", 70), "\n")
  top_conserved <- head(results$conserved_analysis$conserved_analysis, 20)
  
  cat(sprintf("%-4s %-12s %-10s %-12s %-12s %-10s\n",
             "排名", "基因", "出现次数", "平均相关性", "保守评分", "稳定性"))
  cat(rep("-", 70), "\n")
  
  for(i in 1:nrow(top_conserved)) {
    cat(sprintf("%-4d %-12s %-10d %-12.3f %-12.2f %-10.3f\n",
                i, 
                top_conserved$gene_symbol[i],
                top_conserved$n_cell_contexts[i],
                top_conserved$mean_correlation[i],
                top_conserved$conservation_score[i],
                top_conserved$stability_score[i]))
  }
  cat("\n")
  
  # 3. 特异性基因Top 20
  cat("【3】特异性基因 Top 20（按基因数量）\n")
  cat(rep("-", 70), "\n")
  
  if(!is.null(results$specific_analysis$specific_summary)) {
    top_specific <- head(results$specific_analysis$specific_summary[
      order(results$specific_analysis$specific_summary$n_specific_genes, decreasing = TRUE), ], 20)
    
    cat(sprintf("%-4s %-25s %-10s %-10s %-10s\n",
               "排名", "癌种-细胞", "总数", "上调", "下调"))
    cat(rep("-", 70), "\n")
    
    for(i in 1:nrow(top_specific)) {
      cat(sprintf("%-4d %-25s %-10d %-10d %-10d\n",
                  i,
                  top_specific$cancer_cell[i],
                  top_specific$n_specific_genes[i],
                  top_specific$n_upregulated[i],
                  top_specific$n_downregulated[i]))
    }
  }
  cat("\n")
  
  # 4. 功能富集摘要
  cat("【4】保守基因功能富集摘要\n")
  cat(rep("-", 70), "\n")
  
  if(!is.null(results$functional_enrichment)) {
    
    for(db_name in names(results$functional_enrichment)) {
      db_data <- results$functional_enrichment[[db_name]]
      
      if(nrow(db_data) == 0) next
      
      cat(sprintf("\n%s (%d 条通路):\n", db_name, nrow(db_data)))
      cat(rep("·", 70), "\n")
      
      top_pathways <- head(db_data, 10)
      
      for(i in 1:nrow(top_pathways)) {
        pathway_name <- if("TERM" %in% names(top_pathways)) {
          top_pathways$TERM[i]
        } else if("description" %in% names(top_pathways)) {
          top_pathways$description[i]
        } else {
          top_pathways$Pathway_ID[i]
        }
        
        # 截断名称
        if(nchar(pathway_name) > 55) {
          pathway_name <- paste0(substr(pathway_name, 1, 52), "...")
        }
        
        cat(sprintf("  %2d. %s\n      基因数: %d, FDR: %.2e, 富集倍数: %.2f\n",
                    i, pathway_name,
                    top_pathways$Count[i],
                    top_pathways$p.adjust[i],
                    top_pathways$FoldEnrichment[i]))
      }
    }
  } else {
    cat("未进行富集分析或无显著结果\n")
  }
  cat("\n")
  
  # 5. 特异性富集摘要
  cat("【5】特异性基因功能富集摘要\n")
  cat(rep("-", 70), "\n")
  
  if(!is.null(results$specific_enrichment) && length(results$specific_enrichment) > 0) {
    cat(sprintf("共%d个癌种-细胞组合有显著富集\n\n", length(results$specific_enrichment)))
    
    for(cancer_cell in head(names(results$specific_enrichment), 5)) {
      cat(sprintf("\n%s:\n", cancer_cell))
      
      if("GO_BP" %in% names(results$specific_enrichment[[cancer_cell]])) {
        go_data <- head(results$specific_enrichment[[cancer_cell]]$GO_BP, 3)
        
        for(i in 1:nrow(go_data)) {
          term <- ifelse(!is.na(go_data$TERM[i]), go_data$TERM[i], go_data$GO_ID[i])
          cat(sprintf("  • %s (FDR: %.2e)\n", term, go_data$p.adjust[i]))
        }
      }
    }
  } else {
    cat("无特异性基因功能富集结果\n")
  }
  cat("\n")
  
  # 6. 输出文件
  cat("【6】生成的输出文件\n")
  cat(rep("-", 70), "\n")
  cat("主要目录结构:\n\n")
  cat("conserved_modules/\n")
  cat("  ├── conserved_genes_analysis.csv (保守基因统计)\n")
  cat("  └── individual_genes/ (每个基因详细信息)\n\n")
  
  cat("specific_modules/\n")
  cat("  ├── *_specific_genes.csv (每个癌种特异基因)\n")
  cat("  └── specific_genes_summary.csv (特异性汇总)\n\n")
  
  cat("detailed_tables/\n")
  cat("  ├── conserved_genes_detailed.csv (保守基因详细信息)\n")
  cat("  ├── gene_frequency.csv (基因频率)\n")
  cat("  └── conservation_statistics.csv (保守性统计)\n\n")
  
  cat("enrichment_visualization/\n")
  cat("  ├── *_bubble.pdf (富集气泡图)\n")
  cat("  ├── *_barplot.pdf (富集条形图)\n")
  cat("  ├── enrichment_comparison_heatmap.pdf (多数据库比较)\n")
  cat("  ├── GO_hierarchy.pdf (GO层级图)\n")
  cat("  ├── *_network.pdf (通路网络图)\n")
  cat("  └── specific_genes_comparison.pdf (特异性比较)\n\n")
  
  cat("gene_expression_patterns/\n")
  cat("  ├── gene_expression_clusters.pdf (基因表达聚类)\n")
  cat("  ├── conservation_radar.pdf (保守性雷达图)\n")
  cat("  ├── conservation_distribution.pdf (保守性分布)\n")
  cat("  └── gene_frequency.pdf (基因频率)\n\n")
  
  cat("functional_analysis/\n")
  cat("  ├── *_enrichment.csv (富集结果表)\n")
  cat("  └── specific_enrichment/ (特异性富集)\n\n")
  
  cat("interactive_reports/\n")
  cat("  └── analysis_report.html (HTML交互报告)\n\n")
  
  # 7. 分析参数
  cat("【7】分析参数\n")
  cat(rep("-", 70), "\n")
  cat(sprintf("最小癌种数: %d\n", results$parameters$min_cancer_types))
  cat(sprintf("P值阈值: %.3f\n", results$parameters$p_threshold))
  cat(sprintf("相关性阈值: %.3f\n", results$parameters$correlation_threshold))
  cat("\n")
  
  cat(rep("=", 70), "\n")
  cat("报告结束\n")
  cat(rep("=", 70), "\n")
  
  sink()
  
  message(paste("✓ 文本报告:", report_file))
}

# ==================== 执行分析 ====================

if(!interactive()) {
  message("开始执行增强版分析...")
  
  conservation_results <- comprehensive_conservation_analysis(
    min_cancer_types = 2,
    p_threshold = 0.05,
    correlation_threshold = 0.15
  )
  
  message("\n✅ 所有分析完成!")
  message("请查看输出目录的各个子文件夹")
  
} else {
  message("增强版代码已加载!")
  message("运行: conservation_results <- comprehensive_conservation_analysis()")
}

message("跨癌种保守性和特异性分析（增强可视化版）加载完成!")