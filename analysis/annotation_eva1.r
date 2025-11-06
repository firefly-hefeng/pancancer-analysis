library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
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
# 设置并行worker数量（根据你的CPU核心数调整）
# 对于大型Seurat对象，建议减少worker数量以避免内存问题
n_workers <- min(parallel::detectCores() - 2, 12)  # 保留2个核心，最多使用8个
options(future.globals.maxSize = 8 * 1024^3)  # 增加到4GB
plan(multisession, workers = n_workers)
message("并行计算已启用，使用 ", n_workers, " 个workers")

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

# ============================================================================
# 1. Marker基因列表（直接定义）
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
  
  message("免疫细胞marker: ", length(immune_markers), " 个细胞类型")
  message("癌症marker: ", length(cancer_markers), " 个癌种")
  
  return(list(immune = immune_markers, cancer = cancer_markers))
}

# 加载marker基因
marker_genes <- get_marker_genes()

# ============================================================================
# 1.5 绘制细胞类型分布UMAP图
# ============================================================================
plot_celltype_umap <- function(seurat_obj, cancer_type) {
  message("生成细胞类型分布UMAP图...")
  
  # 获取细胞类型数量
  n_celltypes <- length(unique(seurat_obj$fine_cell_type))
  
  # 生成颜色方案
  if (n_celltypes <= 12) {
    colors <- brewer.pal(min(n_celltypes, 12), "Set3")
  } else if (n_celltypes <= 20) {
    colors <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Set2"))
  } else {
    # 使用多个调色板组合
    colors <- c(
      brewer.pal(12, "Set3"),
      brewer.pal(8, "Set2"),
      brewer.pal(9, "Set1"),
      brewer.pal(8, "Dark2")
    )
  }
  
  # 如果颜色还不够，使用viridis补充
  if (n_celltypes > length(colors)) {
    additional_colors <- viridis(n_celltypes - length(colors))
    colors <- c(colors, additional_colors)
  }
  
  # 绘制UMAP
  p <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "fine_cell_type",
    label = TRUE,
    label.size = 3,
    repel = TRUE,
    pt.size = 0.5,
    cols = colors[1:n_celltypes]
  ) +
    ggtitle(paste0(cancer_type, " - Cell Type Distribution")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.text = element_text(size = 8)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
  
  # 保存单独的图
  filename <- file.path(
    output_dir,
    "umap_plots",
    paste0(cancer_type, "_celltype_distribution_umap.pdf")
  )
  
  ggsave(
    filename,
    p,
    width = 10,
    height = 8,
    device = "pdf"
  )
  
  message("✓ 细胞类型分布UMAP图已保存: ", filename)
  
  return(p)
}

# ============================================================================
# 2. Marker基因UMAP可视化函数
# ============================================================================
plot_marker_umap <- function(seurat_obj, markers_list, cancer_type, plot_type = "immune") {
  message("生成", plot_type, "marker基因UMAP图...")
  
  # 首先生成细胞类型分布图
  celltype_plot <- plot_celltype_umap(seurat_obj, cancer_type)
  
  # 串行处理marker基因（避免内存问题）
  plot_list <- list()
  
  for (cell_type in names(markers_list)) {
    genes <- markers_list[[cell_type]]
    
    # 转换为ENSEMBL ID（如果需要）
    if (!all(genes %in% rownames(seurat_obj))) {
      genes_ensembl <- tryCatch({
        mapIds(
          org.Hs.eg.db,
          keys = genes,
          keytype = "SYMBOL",
          column = "ENSEMBL",
          multiVals = "first"
        )
      }, error = function(e) {
        return(NULL)
      })
      
      if (!is.null(genes_ensembl)) {
        genes <- genes_ensembl[!is.na(genes_ensembl)]
      }
    }
    
    # 找到在数据中存在的基因
    available_genes <- intersect(genes, rownames(seurat_obj))
    
    if (length(available_genes) == 0) {
      message("  跳过 ", cell_type, ": 无可用marker基因")
      next
    }
    
    message("  绘制 ", cell_type, " (", length(available_genes), "/", length(markers_list[[cell_type]]), " 个基因)")
    
    tryCatch({
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
        pt.size = 0.5
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
      
      ggsave(
        filename,
        p,
        width = 8,
        height = 6,
        device = "pdf"
      )
      
    }, error = function(e) {
      message("  ✗ ", cell_type, " 绘图失败: ", e$message)
    })
    
    # 清理内存
    gc()
  }
  
  # 创建组合图（包含细胞类型分布图）
  if (length(plot_list) > 0) {
    # 提取所有plot
    all_plots <- c(list(celltype_plot), plot_list)
    
    n_plots <- length(all_plots)
    ncol <- ceiling(sqrt(n_plots))
    nrow <- ceiling(n_plots / ncol)
    
    combined_plot <- wrap_plots(all_plots, ncol = ncol)
    
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
    
    message("✓ 组合图已保存（包含细胞类型分布图）: ", combined_filename)
  }
  
  return(plot_list)
}

# ============================================================================
# 3. Marker基因热图可视化
# ============================================================================
plot_marker_heatmap <- function(seurat_obj, markers_list, cancer_type, plot_type = "immune") {
  message("生成", plot_type, "marker基因热图...")
  
  # 准备所有marker基因
  all_genes <- unique(unlist(markers_list))
  
  # 转换为ENSEMBL ID
  if (!all(all_genes %in% rownames(seurat_obj))) {
    genes_ensembl <- mapIds(
      org.Hs.eg.db,
      keys = all_genes,
      keytype = "SYMBOL",
      column = "ENSEMBL",
      multiVals = "first"
    )
    
    # 创建映射关系
    gene_mapping <- data.frame(
      SYMBOL = names(genes_ensembl),
      ENSEMBL = as.character(genes_ensembl),
      stringsAsFactors = FALSE
    )
    gene_mapping <- gene_mapping[!is.na(gene_mapping$ENSEMBL), ]
  }
  
  # 对每个细胞类型绘制热图
  for (cell_type in names(markers_list)) {
    genes_symbol <- markers_list[[cell_type]]
    
    if (!all(genes_symbol %in% rownames(seurat_obj))) {
      genes_ensembl <- mapIds(
        org.Hs.eg.db,
        keys = genes_symbol,
        keytype = "SYMBOL",
        column = "ENSEMBL",
        multiVals = "first"
      )
      genes_to_plot <- genes_ensembl[!is.na(genes_ensembl)]
    } else {
      genes_to_plot <- genes_symbol
    }
    
    available_genes <- intersect(genes_to_plot, rownames(seurat_obj))
    
    if (length(available_genes) < 2) {
      message("  跳过 ", cell_type, ": marker基因不足")
      next
    }
    
    # 绘制热图
    tryCatch({
      p <- DoHeatmap(
        seurat_obj,
        features = available_genes,
        group.by = "fine_cell_type",
        size = 3,
        angle = 90
      ) +
        ggtitle(paste0(cell_type, " Markers")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      filename <- file.path(
        output_dir,
        "marker_heatmaps",
        paste0(cancer_type, "_", plot_type, "_", cell_type, "_heatmap.pdf")
      )
      
      ggsave(filename, p, width = 12, height = 8, device = "pdf")
      message("  ✓ ", cell_type, " 热图已保存")
      
    }, error = function(e) {
      message("  ✗ ", cell_type, " 热图绘制失败: ", e$message)
    })
    
    gc()
  }
}

# ============================================================================
# 4. FindMarkers分析
# ============================================================================
perform_findmarkers <- function(seurat_obj, cancer_type) {
  message("对细胞类型进行FindMarkers分析...")
  
  # 设置默认assay
  DefaultAssay(seurat_obj) <- "RNA"
  
  # 获取所有细胞类型
  cell_types <- unique(seurat_obj$fine_cell_type)
  cell_types <- cell_types[!is.na(cell_types) & cell_types != "Other"]
  
  message("发现 ", length(cell_types), " 个细胞类型")
  
  all_markers <- list()
  
  # 串行处理（避免并行时的内存问题）
  for (cell_type in cell_types) {
    message("  处理: ", cell_type)
    
    # 检查该细胞类型的细胞数
    n_cells <- sum(seurat_obj$fine_cell_type == cell_type)
    
    if (n_cells < 3) {
      message("    跳过 (细胞数不足): ", n_cells)
      next
    }
    
    tryCatch({
      markers <- FindMarkers(
        seurat_obj,
        ident.1 = cell_type,
        group.by = "fine_cell_type",
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.25,
        test.use = "wilcox"
      )
      
      if (nrow(markers) > 0) {
        markers$gene <- rownames(markers)
        markers$cell_type <- cell_type
        all_markers[[cell_type]] <- markers
        message("    ✓ 发现 ", nrow(markers), " 个marker基因")
      } else {
        message("    未发现显著marker")
      }
      
    }, error = function(e) {
      message("    ✗ FindMarkers失败: ", e$message)
    })
    
    gc()
  }
  
  if (length(all_markers) > 0) {
    # 合并结果
    all_markers_df <- do.call(rbind, all_markers)
    
    # 保存结果
    output_file <- file.path(
      output_dir,
      "summary_stats",
      paste0(cancer_type, "_findmarkers_results.csv")
    )
    
    write.csv(all_markers_df, output_file, row.names = FALSE)
    message("✓ FindMarkers结果已保存: ", output_file)
    
    return(all_markers_df)
  } else {
    message("✗ 未生成FindMarkers结果")
    return(NULL)
  }
}

# ============================================================================
# 5. GO和KEGG富集分析
# ============================================================================
perform_enrichment_analysis <- function(markers_df, cancer_type) {
  message("进行GO和KEGG富集分析...")
  
  if (is.null(markers_df) || nrow(markers_df) == 0) {
    message("  没有marker数据，跳过富集分析")
    return(NULL)
  }
  
  cell_types <- unique(markers_df$cell_type)
  enrichment_results <- list()
  
  # 串行处理富集分析
  for (cell_type in cell_types) {
    message("  分析: ", cell_type)
    
    # 提取该细胞类型的top marker基因
    markers_subset <- markers_df[markers_df$cell_type == cell_type, ]
    markers_subset <- markers_subset[order(markers_subset$avg_log2FC, decreasing = TRUE), ]
    
    # 取top 200或全部（如果少于200）
    top_n <- min(200, nrow(markers_subset))
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
    
    if (length(entrez_ids) < 10) {
      message("    可用基因太少 (", length(entrez_ids), ")，跳过")
      next
    }
    
    message("    使用 ", length(entrez_ids), " 个基因进行富集分析")
    
    # GO富集分析
    tryCatch({
      go_bp <- enrichGO(
        gene = entrez_ids,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
      
      if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
        enrichment_results[[cell_type]]$GO_BP <- go_bp
        message("    ✓ GO BP: ", nrow(go_bp@result), " 个富集通路")
      }
    }, error = function(e) {
      message("    GO BP分析失败: ", e$message)
    })
    
    # KEGG富集分析
    tryCatch({
      kegg <- enrichKEGG(
        gene = entrez_ids,
        organism = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
      
      if (!is.null(kegg) && nrow(kegg@result) > 0) {
        enrichment_results[[cell_type]]$KEGG <- kegg
        message("    ✓ KEGG: ", nrow(kegg@result), " 个富集通路")
      }
    }, error = function(e) {
      message("    KEGG分析失败: ", e$message)
    })
    
    gc()
  }
  
  # 保存富集结果
  if (length(enrichment_results) > 0) {
    save_enrichment_results(enrichment_results, cancer_type)
    plot_enrichment_results(enrichment_results, cancer_type)
  }
  
  return(enrichment_results)
}

# ============================================================================
# 6. 保存富集分析结果
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
# 7. 绘制富集分析图
# ============================================================================
plot_enrichment_results <- function(enrichment_results, cancer_type) {
  message("绘制富集分析图...")
  
  for (cell_type in names(enrichment_results)) {
    # GO BP dotplot
    if (!is.null(enrichment_results[[cell_type]]$GO_BP)) {
      tryCatch({
        p_go <- dotplot(enrichment_results[[cell_type]]$GO_BP, showCategory = 20) +
          ggtitle(paste0(cell_type, " - GO Biological Process")) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
        go_plot_file <- file.path(
          output_dir,
          "enrichment_plots",
          paste0(cancer_type, "_", cell_type, "_GO_BP_dotplot.pdf")
        )
        
        ggsave(go_plot_file, p_go, width = 10, height = 8, device = "pdf")
        
        # GO barplot
        p_go_bar <- barplot(enrichment_results[[cell_type]]$GO_BP, showCategory = 20) +
          ggtitle(paste0(cell_type, " - GO Biological Process")) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
        go_bar_file <- file.path(
          output_dir,
          "enrichment_plots",
          paste0(cancer_type, "_", cell_type, "_GO_BP_barplot.pdf")
        )
        
        ggsave(go_bar_file, p_go_bar, width = 10, height = 8, device = "pdf")
        
      }, error = function(e) {
        message("  GO图绘制失败: ", e$message)
      })
    }
    
    # KEGG dotplot
    if (!is.null(enrichment_results[[cell_type]]$KEGG)) {
      tryCatch({
        p_kegg <- dotplot(enrichment_results[[cell_type]]$KEGG, showCategory = 20) +
          ggtitle(paste0(cell_type, " - KEGG Pathway")) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
        kegg_plot_file <- file.path(
          output_dir,
          "enrichment_plots",
          paste0(cancer_type, "_", cell_type, "_KEGG_dotplot.pdf")
        )
        
        ggsave(kegg_plot_file, p_kegg, width = 10, height = 8, device = "pdf")
        
        # KEGG barplot
        p_kegg_bar <- barplot(enrichment_results[[cell_type]]$KEGG, showCategory = 20) +
          ggtitle(paste0(cell_type, " - KEGG Pathway")) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
        kegg_bar_file <- file.path(
          output_dir,
          "enrichment_plots",
          paste0(cancer_type, "_", cell_type, "_KEGG_barplot.pdf")
        )
        
        ggsave(kegg_bar_file, p_kegg_bar, width = 10, height = 8, device = "pdf")
        
      }, error = function(e) {
        message("  KEGG图绘制失败: ", e$message)
      })
    }
    
    gc()
  }
  
  message("✓ 富集分析图已保存")
}

# ============================================================================
# 8. 生成评估报告
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
  cat("细胞注释质量评估报告\n")
  cat("癌种: ", cancer_type, "\n")
  cat("生成时间: ", as.character(Sys.time()), "\n")
  cat("=", rep("=", 70), "\n\n", sep = "")
  
  # 1. FindMarkers统计
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
      cat("    Top 5 markers: ", paste(head(ct_markers$gene, 5), collapse = ", "), "\n\n")
    }
  } else {
    cat("  无FindMarkers结果\n\n")
  }
  
  # 2. 富集分析统计
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
# 9. 主评估函数
# ============================================================================
evaluate_cancer_annotation <- function(cancer_type) {
  message("\n", rep("=", 80))
  message("开始评估癌种: ", cancer_type)
  message(rep("=", 80))
  
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
  message("细胞数: ", ncol(seurat_obj))
  message("细胞类型数: ", length(unique(seurat_obj$fine_cell_type)))
  
  # 检查是否有UMAP
  if (!"umap" %in% names(seurat_obj@reductions)) {
    message("重新计算UMAP...")
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, verbose = FALSE)
  }
  
  # 1. Marker基因UMAP可视化（包含细胞类型分布图）
  message("\n--- Step 1: Marker基因UMAP可视化 ---")
  plot_marker_umap(seurat_obj, marker_genes$immune, cancer_type, "immune")
  
  # 如果该癌种有特异性marker，也绘制
  if (cancer_type %in% names(marker_genes$cancer)) {
    cancer_specific_markers <- list()
    cancer_specific_markers[[cancer_type]] <- marker_genes$cancer[[cancer_type]]
    plot_marker_umap(seurat_obj, cancer_specific_markers, cancer_type, "cancer_specific")
  }
  
  # 2. Marker基因热图
  message("\n--- Step 2: Marker基因热图 ---")
  plot_marker_heatmap(seurat_obj, marker_genes$immune, cancer_type, "immune")
  
  # 3. FindMarkers分析
  message("\n--- Step 3: FindMarkers分析 ---")
  markers_df <- perform_findmarkers(seurat_obj, cancer_type)
  
  # 4. 富集分析
  message("\n--- Step 4: GO和KEGG富集分析 ---")
  enrichment_results <- perform_enrichment_analysis(markers_df, cancer_type)
  
  # 5. 生成评估报告
  message("\n--- Step 5: 生成评估报告 ---")
  generate_evaluation_report(cancer_type, markers_df, enrichment_results)
  
  message("\n✓ 癌种 ", cancer_type, " 评估完成")
  
  # 清理内存
  gc()
  
  return(list(
    markers = markers_df,
    enrichment = enrichment_results
  ))
}

# ============================================================================
# 10. 批量评估所有癌种（改用串行，避免并行问题）
# ============================================================================
evaluate_all_cancers <- function() {
  message("开始批量评估所有癌种...")
  
  # 获取所有癌种目录
  cancer_dirs <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
  cancer_dirs <- cancer_dirs[cancer_dirs != "reference" & 
                             cancer_dirs != "evaluation_results" &
                             !grepl("^\\.", cancer_dirs)]
  
  message("发现 ", length(cancer_dirs), " 个癌种目录")
  
  results_summary <- list()
  
  # 串行处理所有癌种
  for (cancer_type in cancer_dirs) {
    tryCatch({
      result <- evaluate_cancer_annotation(cancer_type)
      
      if (!is.null(result)) {
        results_summary[[cancer_type]] <- list(
          n_cell_types = length(unique(result$markers$cell_type)),
          n_markers = nrow(result$markers),
          n_enriched_types = length(result$enrichment)
        )
      }
      
    }, error = function(e) {
      warning("评估 ", cancer_type, " 时出错: ", e$message)
    })
    
    # 每个癌种处理完后清理内存
    gc()
  }
  
  # 保存总体汇总
  if (length(results_summary) > 0) {
    summary_df <- do.call(rbind, lapply(names(results_summary), function(ct) {
      data.frame(
        Cancer_Type = ct,
        N_CellTypes = results_summary[[ct]]$n_cell_types,
        N_Markers = results_summary[[ct]]$n_markers,
        N_Enriched = results_summary[[ct]]$n_enriched_types
      )
    }))
    
    summary_file <- file.path(output_dir, "evaluation_summary.csv")
    write.csv(summary_df, summary_file, row.names = FALSE)
    
    message("\n✓ 总体评估汇总已保存: ", summary_file)
    print(summary_df)
  }
  
  return(results_summary)
}

# ============================================================================
# 执行评估
# ============================================================================
message("细胞注释质量评估系统")
message("=", rep("=", 80))

# 选择评估模式
# 模式1: 评估单个癌种
# evaluate_cancer_annotation("BRCA")

# 模式2: 评估所有癌种（串行，稳定但较慢）
evaluation_summary <- evaluate_all_cancers()

# 关闭并行计算
plan(sequential)

message("\n所有评估任务完成!")
message("结果保存在: ", output_dir)