# functional_enrichment.R
# 功能富集分析

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(ReactomePA)
library(msigdbr)
library(dplyr)
library(ggplot2)

# 参数设置
output_base_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation"
enrichment_output_dir <- file.path(output_base_dir, "functional_enrichment")
dir.create(enrichment_output_dir, showWarnings = FALSE, recursive = TRUE)

# 获取基因集数据库
get_genesets <- function() {
  # 获取MSigDB基因集
  h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") # Hallmark
  c2_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
  c5_gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  
  return(list(
    hallmark = h_gene_sets,
    kegg = c2_gene_sets,
    gobp = c5_gene_sets
  ))
}

# 基因ID转换函数
convert_gene_ids <- function(genes, from = "SYMBOL", to = "ENTREZID") {
  converted <- bitr(genes, fromType = from, toType = to, OrgDb = org.Hs.eg.db)
  return(converted)
}

# 运行多种富集分析
run_enrichment_analysis <- function(genes, analysis_name) {
  print(paste("运行富集分析:", analysis_name))
  
  if (length(genes) < 10) {
    warning(paste("基因数量过少:", length(genes)))
    return(NULL)
  }
  
  # 基因ID转换
  gene_df <- convert_gene_ids(genes)
  
  if (nrow(gene_df) < 5) {
    warning(paste("转换后基因数量过少:", nrow(gene_df)))
    return(NULL)
  }
  
  entrez_genes <- gene_df$ENTREZID
  
  results <- list()
  
  # 1. GO富集分析
  tryCatch({
    go_bp <- enrichGO(gene = entrez_genes,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
    
    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
      results$GO_BP <- go_bp
    }
  }, error = function(e) {
    warning(paste("GO分析失败:", e$message))
  })
  
  # 2. KEGG富集分析
  tryCatch({
    kegg <- enrichKEGG(gene = entrez_genes,
                      organism = "hsa",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH")
    
    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      results$KEGG <- kegg
    }
  }, error = function(e) {
    warning(paste("KEGG分析失败:", e$message))
  })
  
  # 3. Reactome富集分析
  tryCatch({
    reactome <- enrichPathway(gene = entrez_genes,
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH",
                             readable = TRUE)
    
    if (!is.null(reactome) && nrow(reactome@result) > 0) {
      results$Reactome <- reactome
    }
  }, error = function(e) {
    warning(paste("Reactome分析失败:", e$message))
  })
  
  # 4. Hallmark基因集分析
  tryCatch({
    genesets <- get_genesets()
    
    # 准备Hallmark基因集
    hallmark_list <- split(genesets$hallmark$gene_symbol, genesets$hallmark$gs_name)
    
    hallmark <- enricher(gene = genes,
                        TERM2GENE = genesets$hallmark[, c("gs_name", "gene_symbol")],
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
    
    if (!is.null(hallmark) && nrow(hallmark@result) > 0) {
      results$Hallmark <- hallmark
    }
  }, error = function(e) {
    warning(paste("Hallmark分析失败:", e$message))
  })
  
  return(results)
}

# 可视化富集结果
plot_enrichment_results <- function(enrichment_results, analysis_name) {
  print(paste("生成富集分析图表:", analysis_name))
  
  plot_dir <- file.path(enrichment_output_dir, "plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  pdf(file.path(plot_dir, paste0(analysis_name, "_enrichment.pdf")),
      width = 14, height = 10)
  
  for (db_name in names(enrichment_results)) {
    result <- enrichment_results[[db_name]]
    
    if (is.null(result) || nrow(result@result) == 0) next
    
    # 柱状图
    p1 <- barplot(result, showCategory = 15) +
      ggtitle(paste(analysis_name, "-", db_name, "Enrichment"))
    
    # 点图
    p2 <- dotplot(result, showCategory = 15) +
      ggtitle(paste(analysis_name, "-", db_name, "Dotplot"))
    
    print(p1)
    print(p2)
    
    # 如果有足够的条目，绘制网络图
    if (nrow(result@result) >= 5) {
      tryCatch({
        p3 <- emapplot(pairwise_termsim(result), showCategory = 15) +
          ggtitle(paste(analysis_name, "-", db_name, "Network"))
        print(p3)
      }, error = function(e) {
        warning("网络图生成失败")
      })
    }
  }
  
  dev.off()
}

# 保存富集结果
save_enrichment_results <- function(enrichment_results, analysis_name) {
  for (db_name in names(enrichment_results)) {
    result <- enrichment_results[[db_name]]
    
    if (!is.null(result) && nrow(result@result) > 0) {
      write.csv(result@result, 
                file.path(enrichment_output_dir, 
                         paste0(analysis_name, "_", db_name, "_enrichment.csv")),
                row.names = FALSE)
    }
  }
}

# 主富集分析函数
main_enrichment_analysis <- function() {
  print("开始功能富集分析...")
  
  # 1. 分析CellChat保守相互作用相关基因
  cellchat_file <- file.path(output_base_dir, "cellchat_analysis", "conserved_interactions.csv")
  if (file.exists(cellchat_file)) {
    print("分析CellChat保守相互作用...")
    
    conserved_interactions <- read.csv(cellchat_file, stringsAsFactors = FALSE)
    
    # 提取配体和受体基因
    ligand_receptor_genes <- unique(c(conserved_interactions$Var1, conserved_interactions$Var2))
    ligand_receptor_genes <- ligand_receptor_genes[!is.na(ligand_receptor_genes)]
    
    if (length(ligand_receptor_genes) > 10) {
      cellchat_enrichment <- run_enrichment_analysis(ligand_receptor_genes, "CellChat_conserved")
      
      if (!is.null(cellchat_enrichment)) {
        plot_enrichment_results(cellchat_enrichment, "CellChat_conserved")
        save_enrichment_results(cellchat_enrichment, "CellChat_conserved")
      }
    }
  }
  
  # 2. 分析轨迹保守基因
  trajectory_file <- file.path(output_base_dir, "trajectory_analysis", "conserved_trajectory_genes.csv")
  if (file.exists(trajectory_file)) {
    print("分析轨迹保守基因...")
    
    conserved_trajectory <- read.csv(trajectory_file, stringsAsFactors = FALSE)
    
    # 选择在至少3个癌种中保守的基因
    highly_conserved_genes <- conserved_trajectory$gene[conserved_trajectory$frequency >= 3]
    
    if (length(highly_conserved_genes) > 10) {
      trajectory_enrichment <- run_enrichment_analysis(highly_conserved_genes, "Trajectory_conserved")
      
      if (!is.null(trajectory_enrichment)) {
        plot_enrichment_results(trajectory_enrichment, "Trajectory_conserved")
        save_enrichment_results(trajectory_enrichment, "Trajectory_conserved")
      }
    }
  }
  
  # 3. 分析每个癌种的差异表达基因（如果有的话）
  analyze_cancer_specific_genes()
  
  # 4. 整合分析报告
  generate_integrated_report()
}

# 分析癌种特异性基因
analyze_cancer_specific_genes <- function() {
  print("分析癌种特异性差异表达基因...")
  
  valid_samples <- read.csv(file.path(output_base_dir, "valid_samples.csv"), 
                           stringsAsFactors = FALSE)
  cancer_types <- unique(valid_samples$cancer_type)
  
  for (cancer_type in cancer_types) {
    print(paste("分析癌种:", cancer_type))
    
    # 检查是否已经分析过
    enrichment_file <- file.path(enrichment_output_dir, 
                                paste0(cancer_type, "_specific_enrichment.rds"))
    
    if (file.exists(enrichment_file)) {
      next
    }
    
    # 读取整合的Seurat对象
    integrated_file <- file.path(output_base_dir, cancer_type, 
                                paste0(cancer_type, "_integrated_annotated.rds"))
    
    if (!file.exists(integrated_file)) {
      next
    }
    
    tryCatch({
      seurat_obj <- readRDS(integrated_file)
      
      # 寻找标志基因
      Idents(seurat_obj) <- "singleR_labels"
      markers <- FindAllMarkers(seurat_obj, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
      
      if (nrow(markers) > 0) {
        # 对每种细胞类型的标志基因进行富集分析
        cell_types <- unique(markers$cluster)
        
        for (cell_type in cell_types) {
          cell_markers <- markers$gene[markers$cluster == cell_type]
          
          if (length(cell_markers) >= 10) {
            enrichment_results <- run_enrichment_analysis(
              cell_markers, 
              paste0(cancer_type, "_", gsub("[^A-Za-z0-9]", "_", cell_type))
            )
            
            if (!is.null(enrichment_results)) {
              save_enrichment_results(
                enrichment_results, 
                paste0(cancer_type, "_", gsub("[^A-Za-z0-9]", "_", cell_type))
              )
            }
          }
        }
        
        # 保存标志基因
        write.csv(markers, 
                  file.path(enrichment_output_dir, paste0(cancer_type, "_markers.csv")),
                  row.names = FALSE)
      }
      
      # 标记为已处理
      saveRDS(TRUE, enrichment_file)
      
    }, error = function(e) {
      warning(paste("癌种", cancer_type, "富集分析失败:", e$message))
    })
  }
}

# 生成整合分析报告
generate_integrated_report <- function() {
  print("生成整合分析报告...")
  
  # 收集所有富集分析结果
  enrichment_files <- list.files(enrichment_output_dir, 
                                pattern = "*_enrichment.csv$", 
                                full.names = TRUE)
  
  if (length(enrichment_files) == 0) {
    warning("没有找到富集分析结果文件")
    return(NULL)
  }
  
  # 统计分析
  summary_stats <- data.frame(
    Analysis = character(),
    Database = character(),
    Significant_Terms = integer(),
    Top_Term = character(),
    Top_Pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (file in enrichment_files) {
    file_name <- basename(file)
    parts <- strsplit(gsub("_enrichment.csv$", "", file_name), "_")[[1]]
    
    analysis_name <- paste(parts[-length(parts)], collapse = "_")
    database <- parts[length(parts)]
    
    enrichment_data <- read.csv(file, stringsAsFactors = FALSE)
    
    if (nrow(enrichment_data) > 0) {
      summary_stats <- rbind(summary_stats, data.frame(
        Analysis = analysis_name,
        Database = database,
        Significant_Terms = nrow(enrichment_data),
        Top_Term = enrichment_data$Description[1],
        Top_Pvalue = enrichment_data$pvalue[1],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 保存汇总统计
  write.csv(summary_stats, 
            file.path(enrichment_output_dir, "enrichment_summary.csv"),
            row.names = FALSE)
  
  # 生成汇总图表
  if (nrow(summary_stats) > 0) {
    pdf(file.path(enrichment_output_dir, "enrichment_summary.pdf"),
        width = 12, height = 8)
    
    # 显著条目数量比较
    p1 <- ggplot(summary_stats, aes(x = Analysis, y = Significant_Terms, fill = Database)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "富集分析结果汇总", y = "显著条目数量")
    
    print(p1)
    
    dev.off()
  }
  
  print("整合分析报告生成完成")
  return(summary_stats)
}

# 执行富集分析
main_enrichment_analysis()

print("功能富集分析完成！")
print(paste("结果保存在:", enrichment_output_dir))