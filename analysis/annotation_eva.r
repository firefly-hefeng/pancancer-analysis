library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(stringr)
library(future)
library(future.apply)

# ============================================================================
# 并行计算设置
# ============================================================================
# 根据你的服务器核心数调整workers数量
# 建议设置为服务器核心数的70-80%
plan(multisession, workers = 10)  # 使用8个核心并行计算
options(future.globals.maxSize = 8000 * 1024^2)  # 增加全局变量大小限制到8GB

# ============================================================================
# 参数设置
# ============================================================================
input_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/fine_annotation"
output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/annotation_evaluation_results"

# 创建输出目录
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "umap_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "enrichment_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "marker_heatmaps"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "summary_stats"), showWarnings = FALSE, recursive = TRUE)

# 优化参数
OPTIMIZATION_PARAMS <- list(
  max_cells_per_ident = 500,      # FindMarkers每组最多使用500个细胞
  max_markers_per_type = 100,     # 每个细胞类型保留top 100 markers
  enrichment_gene_number = 50,    # 富集分析使用top 50基因
  skip_heatmap = TRUE,            # 跳过耗时的热图绘制
  draw_major_umaps_only = TRUE,   # 只绘制主要细胞类型的UMAP
  min_cells_for_analysis = 10     # 细胞数少于10的类型跳过
)

# ============================================================================
# 1. Marker基因列表
# ============================================================================
get_marker_genes <- function() {
  message("加载marker基因列表...")
  
  # 免疫细胞marker
  immune_markers <- list(
    "Immune" = c("PTPRC"),
    "T" = c("CD3D", "CD3E", "CD3G", "CD2", "THEMIS"),
    "CD8T" = c("CD8A", "CD8B"),
    "CD8 Tn" = c("CCR7", "SELL", "IL7R"),
    "CD8 Tex" = c("GZMK", "TCF7", "OXPHOS"),
    "CD8 Trm" = c("ZNF683", "CXCR6"),
    "CD8 Tm" = c("IL7R", "ZNF683"),
    "CD8 Tem" = c("GZMK"),
    "CD8 Temra" = c("CX3CR1"),
    "CD4T" = c("CD4"),
    "CD4 Tn" = c("ADSL", "IL7R"),
    "CD4 Tm" = c("CAPG", "TIMP1", "CCL5", "CREM", "AREG"),
    "CD4 Tem" = c("GZMK"),
    "CD4 Temra" = c("CX3CR1"),
    "TH1" = c("IFNG", "IFNGR1", "IFNGR2"),
    "TH2" = c("GATA3", "IL4", "IL10"),
    "TH9" = c("IL9", "TGFBR2", "IL17RB"),
    "TH17" = c("IL17A", "CCR6", "BATF", "IRF4", "IL23R", "RORA", "STAT3", "IL26"),
    "TH22" = c("AHR", "IL22"),
    "Tfh" = c("BCL6", "CXCR5", "IL21"),
    "γδT" = c("TRDC"),
    "Treg" = c("FOXP3", "TNFRSF9"),
    "B" = c("CD79A", "MS4A1", "CD19"),
    "Plasma" = c("MZB1", "IGHG1", "IGHG4", "JCHAIN"),
    "NK" = c("NKG7", "NCR1"),
    "NKT" = c("NCAM1"),
    "Myeloid" = c("ITGAM"),
    "Monocyte" = c("FCN1", "CD14"),
    "Macrophage" = c("CD68", "MARCO"),
    "M1" = c("CD86", "CD80", "INOS"),
    "M2" = c("CD163"),
    "cDC" = c("CLEC9A", "CD1A", "CD1C", "XCR1", "ITGAX"),
    "cDC1" = c("XCR1", "CLEC9A"),
    "cDC2" = c("CD1C", "CLEC10A", "SIRPA", "HLA-DQA1"),
    "LAM3_cDC" = c("LAMP3", "CCR7", "FSCN1"),
    "pDC" = c("LILRA4", "IL3RA", "CLEC4C"),
    "Neutrophil" = c("CSF3R", "S100A8", "S100A9", "FCGR3B"),
    "Mast" = c("KIT", "CPA3", "TPSB2", "TPSAB1", "FCER1A"),
    "Epithelial" = c("EPCAM", "TSPAN1", "ELF3", "MUC1", "CEACAM5"),
    "Endothelial" = c("PECAM1", "CLDN5", "ECSCR", "ENG", "VWF", "FLT1", "PLVAP"),
    "Fibroblast" = c("DCN", "LUM", "COL3A1", "COL1A1", "COL6A2", "MMP2", "C1S", "PDGFRB"),
    "myCAF" = c("RGS5", "ACTA2", "TAGLN"),
    "iCAF" = c("IL6", "CXCL12"),
    "apCAF" = c("CD74", "HLA")
  )
  
  # 癌症特异性marker
  cancer_markers <- list(
    "OVCA" = c("WFDC2"),
    "UCEC" = c("WFDC2"),
    "BRCA" = c("CLU", "MGP"),
    "PDAC" = c("LAMC2", "TM4SF1"),
    "COAD" = c("CEACAM5"),
    "LIHC" = c("APOH"),
    "PRAD" = c("TMPRSS2", "CLDN4"),
    "KIRC" = c("CA9"),
    "GIST" = c("PDGFRA", "KCNK3", "ANO1")
  )
  
  message("✓ 免疫细胞marker: ", length(immune_markers), " 个细胞类型")
  message("✓ 癌症marker: ", length(cancer_markers), " 个癌种")
  
  return(list(immune = immune_markers, cancer = cancer_markers))
}

marker_genes <- get_marker_genes()

# ============================================================================
# 2. 优化版 Marker基因UMAP可视化
# ============================================================================
plot_marker_umap_optimized <- function(seurat_obj, markers_list, cancer_type, plot_type = "immune") {
  message("生成", plot_type, "marker基因UMAP图（优化版）...")
  
  # 如果启用了只绘制主要类型
  if (OPTIMIZATION_PARAMS$draw_major_umaps_only && plot_type == "immune") {
    major_cell_types <- c("T", "CD8T", "CD4T", "B", "NK", "Myeloid", 
                          "Macrophage", "Epithelial", "Fibroblast", "Endothelial")
    markers_list <- markers_list[names(markers_list) %in% major_cell_types]
    message("  仅绘制主要细胞类型 (", length(markers_list), " 个)")
  }
  
  plot_list <- list()
  
  for (cell_type in names(markers_list)) {
    genes <- markers_list[[cell_type]]
    
    # 转换为ENSEMBL ID（如果需要）
    if (!all(genes %in% rownames(seurat_obj))) {
      genes_ensembl <- mapIds(
        org.Hs.eg.db,
        keys = genes,
        keytype = "SYMBOL",
        column = "ENSEMBL",
        multiVals = "first"
      )
      genes <- genes_ensembl[!is.na(genes_ensembl)]
    }
    
    available_genes <- intersect(genes, rownames(seurat_obj))
    
    if (length(available_genes) == 0) {
      next
    }
    
    message("  绘制 ", cell_type, " (", length(available_genes), " 个基因)")
    
    # 计算module score
    seurat_obj <- AddModuleScore(
      seurat_obj,
      features = list(available_genes),
      name = paste0(cell_type, "_score"),
      assay = "RNA"
    )
    
    score_col <- paste0(cell_type, "_score1")
    
    # 绘制UMAP
    p <- FeaturePlot(
      seurat_obj,
      features = score_col,
      reduction = "umap",
      pt.size = 0.3  # 减小点大小以加快绘图
    ) + 
      scale_color_viridis(option = "magma") +
      ggtitle(paste0(cell_type, " Marker Score")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right"
      )
    
    plot_list[[cell_type]] <- p
    
    # 保存单个图
    filename <- file.path(
      output_dir, 
      "umap_plots",
      paste0(cancer_type, "_", plot_type, "_", cell_type, "_umap.pdf")
    )
    
    ggsave(filename, p, width = 8, height = 6, device = "pdf")
  }
  
  # 创建组合图
  if (length(plot_list) > 0) {
    n_plots <- length(plot_list)
    ncol <- min(3, ceiling(sqrt(n_plots)))  # 最多3列
    nrow <- ceiling(n_plots / ncol)
    
    combined_plot <- wrap_plots(plot_list, ncol = ncol)
    
    combined_filename <- file.path(
      output_dir,
      "umap_plots",
      paste0(cancer_type, "_", plot_type, "_markers_combined.pdf")
    )
    
    ggsave(
      combined_filename,
      combined_plot,
      width = ncol * 6,
      height = nrow * 5,
      device = "pdf",
      limitsize = FALSE
    )
    
    message("✓ 组合图已保存: ", combined_filename)
  }
  
  return(plot_list)
}

# ============================================================================
# 3. 优化版 FindMarkers分析（并行）
# ============================================================================
perform_findmarkers_optimized <- function(seurat_obj, cancer_type) {
  message("对细胞类型进行FindMarkers分析（并行优化版）...")
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  # 获取所有细胞类型并过滤
  cell_types <- unique(seurat_obj$fine_cell_type)
  cell_types <- cell_types[!is.na(cell_types) & cell_types != "Other"]
  
  # 统计每个细胞类型的细胞数
  cell_counts <- table(seurat_obj$fine_cell_type)
  valid_types <- names(cell_counts[cell_counts >= OPTIMIZATION_PARAMS$min_cells_for_analysis])
  cell_types <- intersect(cell_types, valid_types)
  
  message("发现 ", length(cell_types), " 个细胞类型（过滤后，>=", 
          OPTIMIZATION_PARAMS$min_cells_for_analysis, " 个细胞）")
  
  # 并行执行FindMarkers
  message("开始并行计算（", 10 , " 个核心）...")
  
  all_markers <- future_lapply(cell_types, function(cell_type) {
    
    n_cells <- sum(seurat_obj$fine_cell_type == cell_type)
    message("  [", Sys.time(), "] 处理: ", cell_type, " (", n_cells, " 个细胞)")
    
    tryCatch({
      markers <- FindMarkers(
        seurat_obj,
        ident.1 = cell_type,
        group.by = "fine_cell_type",
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.25,
        test.use = "wilcox",
        max.cells.per.ident = OPTIMIZATION_PARAMS$max_cells_per_ident,  # 限制细胞数
        verbose = FALSE
      )
      
      if (nrow(markers) > 0) {
        markers$gene <- rownames(markers)
        markers$cell_type <- cell_type
        
        # 按log2FC排序并只保留top markers
        markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), ]
        markers <- head(markers, OPTIMIZATION_PARAMS$max_markers_per_type)
        
        message("    ✓ ", cell_type, ": 发现 ", nrow(markers), " 个marker基因")
        return(markers)
      }
      
      message("    - ", cell_type, ": 未发现显著marker")
      return(NULL)
      
    }, error = function(e) {
      message("    ✗ ", cell_type, " 失败: ", e$message)
      return(NULL)
    })
    
  }, future.seed = TRUE)
  
  # 过滤NULL结果并合并
  all_markers <- all_markers[!sapply(all_markers, is.null)]
  
  if (length(all_markers) > 0) {
    all_markers_df <- do.call(rbind, all_markers)
    
    output_file <- file.path(
      output_dir,
      "summary_stats",
      paste0(cancer_type, "_findmarkers_results.csv")
    )
    
    write.csv(all_markers_df, output_file, row.names = FALSE)
    message("✓ FindMarkers结果已保存: ", output_file)
    message("  总计: ", nrow(all_markers_df), " 个marker基因，", 
            length(unique(all_markers_df$cell_type)), " 个细胞类型")
    
    return(all_markers_df)
  } else {
    message("✗ 未生成FindMarkers结果")
    return(NULL)
  }
}

# ============================================================================
# 4. 优化版 GO和KEGG富集分析（并行）
# ============================================================================
perform_enrichment_analysis_optimized <- function(markers_df, cancer_type) {
  message("进行GO和KEGG富集分析（并行优化版）...")
  
  if (is.null(markers_df) || nrow(markers_df) == 0) {
    message("  没有marker数据，跳过富集分析")
    return(NULL)
  }
  
  cell_types <- unique(markers_df$cell_type)
  message("对 ", length(cell_types), " 个细胞类型进行富集分析")
  
  # 并行执行富集分析
  enrichment_results <- future_lapply(cell_types, function(cell_type) {
    
    message("  [", Sys.time(), "] 分析: ", cell_type)
    
    # 提取该细胞类型的marker基因
    markers_subset <- markers_df[markers_df$cell_type == cell_type, ]
    markers_subset <- markers_subset[order(markers_subset$avg_log2FC, decreasing = TRUE), ]
    
    # 只取top N基因
    top_n <- min(OPTIMIZATION_PARAMS$enrichment_gene_number, nrow(markers_subset))
    top_markers <- head(markers_subset, top_n)
    
    gene_ids <- top_markers$gene
    
    # 转换ENSEMBL ID到ENTREZ ID
    if (all(grepl("^ENS", gene_ids))) {
      entrez_ids <- mapIds(
        org.Hs.eg.db,
        keys = gene_ids,
        keytype = "ENSEMBL",
        column = "ENTREZID",
        multiVals = "first"
      )
    } else {
      entrez_ids <- mapIds(
        org.Hs.eg.db,
        keys = gene_ids,
        keytype = "SYMBOL",
        column = "ENTREZID",
        multiVals = "first"
      )
    }
    
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    
    if (length(entrez_ids) < 5) {
      message("    - 可用基因太少 (", length(entrez_ids), ")，跳过")
      return(NULL)
    }
    
    message("    使用 ", length(entrez_ids), " 个基因")
    
    result <- list()
    
    # GO富集分析
    tryCatch({
      go_bp <- enrichGO(
        gene = entrez_ids,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE,
        minGSSize = 10,
        maxGSSize = 500
      )
      
      if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
        result$GO_BP <- go_bp
        message("    ✓ GO BP: ", nrow(go_bp@result), " 个通路")
      }
    }, error = function(e) {
      message("    - GO BP分析失败: ", e$message)
    })
    
    # KEGG富集分析
    tryCatch({
      kegg <- enrichKEGG(
        gene = entrez_ids,
        organism = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        minGSSize = 5,
        maxGSSize = 500
      )
      
      if (!is.null(kegg) && nrow(kegg@result) > 0) {
        result$KEGG <- kegg
        message("    ✓ KEGG: ", nrow(kegg@result), " 个通路")
      }
    }, error = function(e) {
      message("    - KEGG分析失败: ", e$message)
    })
    
    if (length(result) > 0) {
      return(list(cell_type = cell_type, results = result))
    }
    
    return(NULL)
    
  }, future.seed = TRUE)
  
  # 整理结果
  enrichment_results <- enrichment_results[!sapply(enrichment_results, is.null)]
  
  if (length(enrichment_results) > 0) {
    # 转换为命名列表
    named_results <- list()
    for (item in enrichment_results) {
      named_results[[item$cell_type]] <- item$results
    }
    
    message("✓ 完成富集分析，", length(named_results), " 个细胞类型有结果")
    
    save_enrichment_results(named_results, cancer_type)
    plot_enrichment_results(named_results, cancer_type)
    
    return(named_results)
  } else {
    message("✗ 未生成富集分析结果")
    return(NULL)
  }
}

# ============================================================================
# 5. 保存富集分析结果
# ============================================================================
save_enrichment_results <- function(enrichment_results, cancer_type) {
  message("保存富集分析结果...")
  
  for (cell_type in names(enrichment_results)) {
    # 保存GO结果
    if (!is.null(enrichment_results[[cell_type]]$GO_BP)) {
      go_file <- file.path(
        output_dir,
        "summary_stats",
        paste0(cancer_type, "_", cell_type, "_GO_BP.csv")
      )
      
      write.csv(
        enrichment_results[[cell_type]]$GO_BP@result,
        go_file,
        row.names = FALSE
      )
    }
    
    # 保存KEGG结果
    if (!is.null(enrichment_results[[cell_type]]$KEGG)) {
      kegg_file <- file.path(
        output_dir,
        "summary_stats",
        paste0(cancer_type, "_", cell_type, "_KEGG.csv")
      )
      
      write.csv(
        enrichment_results[[cell_type]]$KEGG@result,
        kegg_file,
        row.names = FALSE
      )
    }
  }
  
  message("✓ 富集分析结果已保存")
}

# ============================================================================
# 6. 绘制富集分析图
# ============================================================================
plot_enrichment_results <- function(enrichment_results, cancer_type) {
  message("绘制富集分析图...")
  
  for (cell_type in names(enrichment_results)) {
    # GO BP图
    if (!is.null(enrichment_results[[cell_type]]$GO_BP)) {
      tryCatch({
        # Dotplot
        p_go <- dotplot(enrichment_results[[cell_type]]$GO_BP, showCategory = 15) +
          ggtitle(paste0(cell_type, " - GO BP")) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
        go_plot_file <- file.path(
          output_dir,
          "enrichment_plots",
          paste0(cancer_type, "_", cell_type, "_GO_BP_dotplot.pdf")
        )
        
        ggsave(go_plot_file, p_go, width = 10, height = 8, device = "pdf")
        
      }, error = function(e) {
        message("  GO图绘制失败: ", e$message)
      })
    }
    
    # KEGG图
    if (!is.null(enrichment_results[[cell_type]]$KEGG)) {
      tryCatch({
        # Dotplot
        p_kegg <- dotplot(enrichment_results[[cell_type]]$KEGG, showCategory = 15) +
          ggtitle(paste0(cell_type, " - KEGG")) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
        kegg_plot_file <- file.path(
          output_dir,
          "enrichment_plots",
          paste0(cancer_type, "_", cell_type, "_KEGG_dotplot.pdf")
        )
        
        ggsave(kegg_plot_file, p_kegg, width = 10, height = 8, device = "pdf")
        
      }, error = function(e) {
        message("  KEGG图绘制失败: ", e$message)
      })
    }
  }
  
  message("✓ 富集分析图已保存")
}

# ============================================================================
# 7. 生成评估报告
# ============================================================================
generate_evaluation_report <- function(cancer_type, markers_df, enrichment_results) {
  message("生成评估报告...")
  
  report_file <- file.path(
    output_dir,
    "summary_stats",
    paste0(cancer_type, "_evaluation_report.txt")
  )
  
  sink(report_file)
  
  cat("=", rep("=", 70), "\n", sep = "")
  cat("细胞注释质量评估报告（优化版）\n")
  cat("癌种: ", cancer_type, "\n")
  cat("生成时间: ", as.character(Sys.time()), "\n")
  cat("=", rep("=", 70), "\n\n", sep = "")
  
  # 优化参数
  cat("优化参数:\n")
  cat("-", rep("-", 70), "\n", sep = "")
  cat("  每组最大细胞数: ", OPTIMIZATION_PARAMS$max_cells_per_ident, "\n")
  cat("  每类型最大marker数: ", OPTIMIZATION_PARAMS$max_markers_per_type, "\n")
  cat("  富集分析基因数: ", OPTIMIZATION_PARAMS$enrichment_gene_number, "\n")
  cat("  最小细胞数阈值: ", OPTIMIZATION_PARAMS$min_cells_for_analysis, "\n\n")
  
  # FindMarkers统计
  cat("1. FINDMARKERS分析结果\n")
  cat("-", rep("-", 70), "\n", sep = "")
  
  if (!is.null(markers_df)) {
    cell_types <- unique(markers_df$cell_type)
    cat("细胞类型数量: ", length(cell_types), "\n\n")
    
    for (ct in cell_types) {
      ct_markers <- markers_df[markers_df$cell_type == ct, ]
      cat("  ", ct, ":\n")
      cat("    Marker基因数: ", nrow(ct_markers), "\n")
      cat("    平均log2FC: ", round(mean(ct_markers$avg_log2FC), 3), "\n")
      cat("    最大log2FC: ", round(max(ct_markers$avg_log2FC), 3), "\n")
      cat("    Top 5 markers: ", paste(head(ct_markers$gene, 5), collapse = ", "), "\n\n")
    }
  } else {
    cat("  无FindMarkers结果\n\n")
  }
  
  # 富集分析统计
  cat("\n2. 富集分析结果\n")
  cat("-", rep("-", 70), "\n", sep = "")
  
  if (!is.null(enrichment_results) && length(enrichment_results) > 0) {
    for (ct in names(enrichment_results)) {
      cat("  ", ct, ":\n")
      
      if (!is.null(enrichment_results[[ct]]$GO_BP)) {
        n_go <- nrow(enrichment_results[[ct]]$GO_BP@result)
        cat("    GO BP通路数: ", n_go, "\n")
        if (n_go > 0) {
          top_go <- head(enrichment_results[[ct]]$GO_BP@result$Description, 3)
          cat("    Top 3 GO: ", paste(top_go, collapse = "; "), "\n")
        }
      }
      
      if (!is.null(enrichment_results[[ct]]$KEGG)) {
        n_kegg <- nrow(enrichment_results[[ct]]$KEGG@result)
        cat("    KEGG通路数: ", n_kegg, "\n")
        if (n_kegg > 0) {
          top_kegg <- head(enrichment_results[[ct]]$KEGG@result$Description, 3)
          cat("    Top 3 KEGG: ", paste(top_kegg, collapse = "; "), "\n")
        }
      }
      
      cat("\n")
    }
  } else {
    cat("  无富集分析结果\n\n")
  }
  
  cat("\n", "=", rep("=", 70), "\n", sep = "")
  cat("评估完成\n")
  cat("=", rep("=", 70), "\n", sep = "")
  
  sink()
  
  message("✓ 评估报告已保存: ", report_file)
}

# ============================================================================
# 8. 主评估函数（优化版）
# ============================================================================
evaluate_cancer_annotation <- function(cancer_type) {
  message("\n", rep("=", 80))
  message("开始评估癌种: ", cancer_type)
  message("时间: ", Sys.time())
  message(rep("=", 80))
  
  start_time <- Sys.time()
  
  # 读取数据
  seurat_file <- file.path(
    input_dir,
    cancer_type,
    paste0(cancer_type, "_fine_annotated.rds")
  )
  
  if (!file.exists(seurat_file)) {
    warning("文件不存在: ", seurat_file)
    return(NULL)
  }
  
  message("读取Seurat对象...")
  seurat_obj <- readRDS(seurat_file)
  message("✓ 细胞数: ", ncol(seurat_obj))
  message("✓ 细胞类型数: ", length(unique(seurat_obj$fine_cell_type)))
  
  # 检查是否有UMAP
  if (!"umap" %in% names(seurat_obj@reductions)) {
    message("重新计算UMAP...")
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, verbose = FALSE)
  }
  
  # Step 1: Marker基因UMAP可视化
  message("\n--- Step 1: Marker基因UMAP可视化 ---")
  step1_start <- Sys.time()
  plot_marker_umap_optimized(seurat_obj, marker_genes$immune, cancer_type, "immune")
  
  if (cancer_type %in% names(marker_genes$cancer)) {
    cancer_specific_markers <- list()
    cancer_specific_markers[[cancer_type]] <- marker_genes$cancer[[cancer_type]]
    plot_marker_umap_optimized(seurat_obj, cancer_specific_markers, cancer_type, "cancer_specific")
  }
  step1_time <- as.numeric(difftime(Sys.time(), step1_start, units = "mins"))
  message("Step 1 完成，耗时: ", round(step1_time, 2), " 分钟")
  
  # Step 2: 跳过热图（可选）
  if (!OPTIMIZATION_PARAMS$skip_heatmap) {
    message("\n--- Step 2: Marker基因热图 ---")
    step2_start <- Sys.time()
    # plot_marker_heatmap函数保持不变
    step2_time <- as.numeric(difftime(Sys.time(), step2_start, units = "mins"))
    message("Step 2 完成，耗时: ", round(step2_time, 2), " 分钟")
  } else {
    message("\n--- Step 2: 跳过热图绘制（优化选项） ---")
  }
  
  # Step 3: FindMarkers分析
  message("\n--- Step 3: FindMarkers分析（并行） ---")
  step3_start <- Sys.time()
  markers_df <- perform_findmarkers_optimized(seurat_obj, cancer_type)
  step3_time <- as.numeric(difftime(Sys.time(), step3_start, units = "mins"))
  message("Step 3 完成，耗时: ", round(step3_time, 2), " 分钟")
  
  # Step 4: 富集分析
  message("\n--- Step 4: GO和KEGG富集分析（并行） ---")
  step4_start <- Sys.time()
  enrichment_results <- perform_enrichment_analysis_optimized(markers_df, cancer_type)
  step4_time <- as.numeric(difftime(Sys.time(), step4_start, units = "mins"))
  message("Step 4 完成，耗时: ", round(step4_time, 2), " 分钟")
  
  # Step 5: 生成评估报告
  message("\n--- Step 5: 生成评估报告 ---")
  generate_evaluation_report(cancer_type, markers_df, enrichment_results)
  
  # 总结
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message("\n", rep("=", 80))
  message("✓ 癌种 ", cancer_type, " 评估完成")
  message("总耗时: ", round(total_time, 2), " 分钟")
  message(rep("=", 80))
  
  # 清理内存
  gc()
  
  return(list(
    markers = markers_df,
    enrichment = enrichment_results,
    time_usage = list(
      step1 = step1_time,
      step3 = step3_time,
      step4 = step4_time,
      total = total_time
    )
  ))
}

# ============================================================================
# 9. 批量评估所有癌种（优化版）
# ============================================================================
evaluate_all_cancers <- function() {
  message("开始批量评估所有癌种（优化版）...")
  message("并行核心数: ", 10)
  
  overall_start <- Sys.time()
  
  # 获取所有癌种目录
  cancer_dirs <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
  cancer_dirs <- cancer_dirs[cancer_dirs != "reference" & 
                             cancer_dirs != "evaluation_results" &
                             !grepl("^\\.", cancer_dirs)]
  
  message("发现 ", length(cancer_dirs), " 个癌种目录")
  
  results_summary <- list()
  
  for (i in seq_along(cancer_dirs)) {
    cancer_type <- cancer_dirs[i]
    message("\n", rep("#", 80))
    message("进度: [", i, "/", length(cancer_dirs), "] ", cancer_type)
    message(rep("#", 80))
    
    tryCatch({
      result <- evaluate_cancer_annotation(cancer_type)
      
      if (!is.null(result)) {
        results_summary[[cancer_type]] <- list(
          n_cell_types = if(!is.null(result$markers)) length(unique(result$markers$cell_type)) else 0,
          n_markers = if(!is.null(result$markers)) nrow(result$markers) else 0,
          n_enriched_types = if(!is.null(result$enrichment)) length(result$enrichment) else 0,
          time_minutes = result$time_usage$total
        )
      }
      
    }, error = function(e) {
      warning("评估 ", cancer_type, " 时出错: ", e$message)
      results_summary[[cancer_type]] <- list(
        n_cell_types = NA,
        n_markers = NA,
        n_enriched_types = NA,
        time_minutes = NA,
        error = e$message
      )
    })
    
    # 每个癌种完成后清理内存
    gc()
  }
  
  # 保存总体汇总
  if (length(results_summary) > 0) {
    summary_df <- do.call(rbind, lapply(names(results_summary), function(ct) {
      data.frame(
        Cancer_Type = ct,
        N_CellTypes = results_summary[[ct]]$n_cell_types,
        N_Markers = results_summary[[ct]]$n_markers,
        N_Enriched = results_summary[[ct]]$n_enriched_types,
        Time_Minutes = round(results_summary[[ct]]$time_minutes, 2),
        stringsAsFactors = FALSE
      )
    }))
    
    summary_file <- file.path(output_dir, "evaluation_summary.csv")
    write.csv(summary_df, summary_file, row.names = FALSE)
    
    # 计算总时间
    overall_time <- as.numeric(difftime(Sys.time(), overall_start, units = "hours"))
    
    message("\n", rep("=", 80))
    message("✓ 所有癌种评估完成")
    message("总耗时: ", round(overall_time, 2), " 小时")
    message("总体评估汇总已保存: ", summary_file)
    message(rep("=", 80))
    
    print(summary_df)
  }
  
  return(results_summary)
}

# ============================================================================
# 10. 执行评估
# ============================================================================
message("\n")
message("╔", rep("═", 78), "╗")
message("║", sprintf("%-78s", "  细胞注释质量评估系统 - 优化版"), "║")
message("╠", rep("═", 78), "╣")
message("║", sprintf("%-78s", paste0("  输出目录: ", output_dir)), "║")
message("╚", rep("═", 78), "╝")
message("\n")

# ============================================================================
# 选择运行模式
# ============================================================================

# 模式1: 测试单个癌种（推荐先测试）
# result <- evaluate_cancer_annotation("BRCA")

# 模式2: 批量评估所有癌种
evaluation_summary <- evaluate_all_cancers()

message("\n所有评估任务完成!")
message("结果保存在: ", output_dir)

# 关闭并行计算
plan(sequential)