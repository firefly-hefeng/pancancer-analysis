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
library(ggsci)

# ============================================================================
# 并行计算设置
# ============================================================================
n_workers <- min(parallel::detectCores() - 2, 25)
options(future.globals.maxSize = 250 * 1024^3)
plan(multisession, workers = n_workers)
message("并行计算已启用，使用 ", n_workers, " 个workers")

# ============================================================================
# 参数设置
# ============================================================================
input_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/fine_annotation"
output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/annotation_evaluation_results_11_9"

# 创建输出目录
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "umap_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "enrichment_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "marker_heatmaps"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "summary_stats"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "processed_rds"), showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 🔧 改进的保存函数 - 同时保存PDF和PNG
# ============================================================================
save_plot_both_formats <- function(plot, filename_base, width, height, dpi = 300) {
  # 保存PDF
  pdf_file <- paste0(filename_base, ".pdf")
  ggsave(pdf_file, plot, width = width, height = height, device = "pdf")
  
  # 保存PNG
  png_file <- paste0(filename_base, ".png")
  ggsave(png_file, plot, width = width, height = height, device = "png", dpi = dpi)
  
  message("  ✓ 已保存: ", basename(pdf_file), " 和 ", basename(png_file))
}

# ============================================================================
# 🎨 专业配色方案生成函数
# ============================================================================
generate_professional_colors <- function(cell_types) {
  cell_category_colors <- list(
    "T_cell" = c("CD3", "CD4", "CD8", "Treg", "Tfh", "TH", "Tn", "Tm", "Tem", "Temra", "Tex", "Trm", "γδT", "NKT"),
    "B_cell" = c("B cell", "Plasma", "Plasmablast"),
    "NK_cell" = c("NK"),
    "Myeloid" = c("Monocyte", "Macrophage", "M1", "M2", "DC", "cDC", "pDC", "Dendritic", "LAM3"),
    "Neutrophil" = c("Neutrophil", "Granulocyte"),
    "Mast" = c("Mast"),
    "Epithelial" = c("Epithelial", "Cancer", "Tumor", "Malignant"),
    "Endothelial" = c("Endothelial", "Vascular"),
    "Fibroblast" = c("Fibroblast", "CAF", "myCAF", "iCAF", "apCAF"),
    "Other" = c("Other", "Unknown")
  )
  
  color_palettes <- list(
    "T_cell" = colorRampPalette(c("#E64B35", "#F39B7F", "#FFDC91"))(20),
    "B_cell" = colorRampPalette(c("#4DBBD5", "#91D1C2", "#B4E7CE"))(10),
    "NK_cell" = colorRampPalette(c("#00A087", "#3CB371"))(5),
    "Myeloid" = colorRampPalette(c("#3C5488", "#6A8AC7", "#9BB7D4"))(15),
    "Neutrophil" = colorRampPalette(c("#F39B7F", "#FDB462"))(5),
    "Mast" = colorRampPalette(c("#8491B4", "#B09CC1"))(5),
    "Epithelial" = colorRampPalette(c("#E64B35", "#DC0000", "#B24745"))(10),
    "Endothelial" = colorRampPalette(c("#7E6148", "#B39C7D"))(10),
    "Fibroblast" = colorRampPalette(c("#00A087", "#469E7A", "#7CB992"))(10),
    "Other" = c("#CCCCCC", "#999999", "#666666")
  )
  
  color_assignment <- rep(NA, length(cell_types))
  names(color_assignment) <- cell_types
  
  for (i in seq_along(cell_types)) {
    ct <- cell_types[i]
    assigned <- FALSE
    
    for (category in names(cell_category_colors)) {
      keywords <- cell_category_colors[[category]]
      
      if (any(sapply(keywords, function(kw) grepl(kw, ct, ignore.case = TRUE)))) {
        category_cells <- cell_types[sapply(cell_types, function(x) {
          any(sapply(keywords, function(kw) grepl(kw, x, ignore.case = TRUE)))
        })]
        
        idx <- which(category_cells == ct)
        palette <- color_palettes[[category]]
        
        color_idx <- ((idx - 1) %% length(palette)) + 1
        color_assignment[ct] <- palette[color_idx]
        assigned <- TRUE
        break
      }
    }
    
    if (!assigned) {
      color_assignment[ct] <- "#CCCCCC"
    }
  }
  
  return(color_assignment)
}

# ============================================================================
# 🎨 优化的细胞类型UMAP图 - 双格式保存
# ============================================================================
plot_celltype_umap_enhanced <- function(seurat_obj, cancer_type) {
  message("生成专业版细胞类型分布UMAP图...")
  
  cell_types <- unique(seurat_obj$fine_cell_type)
  cell_types <- cell_types[!is.na(cell_types)]
  n_celltypes <- length(cell_types)
  
  message("细胞类型数量: ", n_celltypes)
  
  colors <- generate_professional_colors(cell_types)
  
  umap_coords <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  umap_coords$cell_type <- seurat_obj$fine_cell_type
  
  label_data <- umap_coords %>%
    group_by(cell_type) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = 'drop'
    ) %>%
    filter(!is.na(cell_type))
  
  # 主UMAP图
  p_main <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = colors, name = "Cell Type") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    ) +
    labs(
      title = paste0(cancer_type, " - Cell Type Distribution"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    geom_text(
      data = label_data,
      aes(x = UMAP_1, y = UMAP_2, label = cell_type),
      size = 2.5,
      fontface = "bold",
      color = "black",
      bg.color = "white",
      bg.r = 0.1,
      check_overlap = TRUE
    )
  
  filename_main <- file.path(
    output_dir,
    "umap_plots",
    paste0(cancer_type, "_celltype_umap_main")
  )
  
  save_plot_both_formats(p_main, filename_main, width = 14, height = 12)
  
  # 带图例版本
  max_per_col <- 15
  n_cols <- ceiling(n_celltypes / max_per_col)
  
  p_with_legend <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = colors, name = "Cell Type") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.title = element_text(face = "bold", size = 11),
      legend.key.size = unit(0.4, "cm"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    ) +
    labs(
      title = paste0(cancer_type, " - Cell Type Distribution"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 4, alpha = 1),
        ncol = n_cols,
        byrow = FALSE
      )
    )
  
  filename_legend <- file.path(
    output_dir,
    "umap_plots",
    paste0(cancer_type, "_celltype_umap_with_legend")
  )
  
  legend_width <- 12 + (n_cols - 1) * 3
  
  save_plot_both_formats(p_with_legend, filename_legend, width = legend_width, height = 12)
  
  # 分面图（如果细胞类型较多）
  if (n_celltypes > 30) {
    message("细胞类型较多，生成分面图...")
    
    top_types <- umap_coords %>%
      count(cell_type, sort = TRUE) %>%
      head(30) %>%
      pull(cell_type)
    
    facet_data <- umap_coords %>%
      filter(cell_type %in% top_types)
    
    p_facet <- ggplot(facet_data, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(
        data = umap_coords,
        color = "gray90",
        size = 0.1,
        alpha = 0.3
      ) +
      geom_point(
        aes(color = cell_type),
        size = 0.5,
        alpha = 0.8
      ) +
      scale_color_manual(values = colors) +
      facet_wrap(~ cell_type, ncol = 6) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_rect(fill = "gray95", color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        axis.text = element_text(size = 7)
      ) +
      labs(
        title = paste0(cancer_type, " - Top 30 Cell Types (Faceted)"),
        x = "UMAP 1",
        y = "UMAP 2"
      )
    
    filename_facet <- file.path(
      output_dir,
      "umap_plots",
      paste0(cancer_type, "_celltype_umap_faceted")
    )
    
    save_plot_both_formats(p_facet, filename_facet, width = 20, height = 16)
  }
  
  return(p_main)
}

# ============================================================================
# 🎨 优化的Marker基因UMAP图 - 双格式保存
# ============================================================================
plot_marker_umap_enhanced <- function(seurat_obj, markers_list, cancer_type, plot_type = "immune") {
  message("生成", plot_type, "marker基因UMAP图（专业版）...")
  
  plot_list <- list()
  
  for (cell_type in names(markers_list)) {
    genes <- markers_list[[cell_type]]
    
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
    
    available_genes <- intersect(genes, rownames(seurat_obj))
    
    if (length(available_genes) == 0) {
      message("  跳过 ", cell_type, ": 无可用marker基因")
      next
    }
    
    message("  绘制 ", cell_type, " (", length(available_genes), "/", length(markers_list[[cell_type]]), " 个基因)")
    
    tryCatch({
      seurat_obj <- AddModuleScore(
        seurat_obj,
        features = list(available_genes),
        name = paste0(cell_type, "_score"),
        assay = "RNA"
      )
      
      score_col <- paste0(cell_type, "_score1")
      
      p <- FeaturePlot(
        seurat_obj,
        features = score_col,
        reduction = "umap",
        pt.size = 0.5,
        order = TRUE
      ) + 
        scale_color_gradientn(
          colors = c("#F0F0F0", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C", "#08306B"),
          name = "Score",
          guide = guide_colorbar(
            barwidth = 1,
            barheight = 8,
            title.position = "top",
            title.hjust = 0.5
          )
        ) +
        ggtitle(paste0(cell_type, " Marker Expression")) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "#2C3E50"),
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold", size = 12)
        ) +
        labs(x = "UMAP 1", y = "UMAP 2")
      
      plot_list[[cell_type]] <- p
      
      filename <- file.path(
        output_dir, 
        "umap_plots",
        paste0(cancer_type, "_", plot_type, "_", cell_type, "_umap")
      )
      
      save_plot_both_formats(p, filename, width = 10, height = 8)
      
    }, error = function(e) {
      message("  ✗ ", cell_type, " 绘图失败: ", e$message)
    })
    
    gc()
  }
  
  # 组合图
  if (length(plot_list) > 1) {
    n_plots <- length(plot_list)
    ncol <- min(3, ceiling(sqrt(n_plots)))
    nrow <- ceiling(n_plots / ncol)
    
    combined_plot <- wrap_plots(plot_list, ncol = ncol) +
      plot_annotation(
        title = paste0(cancer_type, " - ", toupper(plot_type), " Marker Expression"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 20, color = "#2C3E50")
        )
      )
    
    combined_filename <- file.path(
      output_dir,
      "umap_plots",
      paste0(cancer_type, "_", plot_type, "_markers_combined")
    )
    
    save_plot_both_formats(
      combined_plot,
      combined_filename,
      width = ncol * 8,
      height = nrow * 7
    )
  }
  
  return(plot_list)
}

# ============================================================================
# 加载预下载的数据库
# ============================================================================
load_pathway_databases <- function() {
  message("加载预下载的pathway数据库...")
  
  db_list <- list()
  
  msigdb_file <- file.path(output_dir, "msigdb_kegg_genesets.RData")
  if (file.exists(msigdb_file)) {
    load(msigdb_file)
    db_list$msigdb_kegg <- msigdb_data
    message("✓ 加载MSigDB KEGG数据库")
    message("  - Pathways: ", nrow(msigdb_data$term2name))
  } else {
    message("✗ MSigDB KEGG数据库未找到")
  }
  
  reactome_file <- file.path(output_dir, "reactome_genesets.RData")
  if (file.exists(reactome_file)) {
    load(reactome_file)
    db_list$reactome <- reactome_data
    message("✓ 加载Reactome数据库")
    message("  - Pathways: ", nrow(reactome_data$term2name))
  } else {
    message("✗ Reactome数据库未找到")
  }
  
  kegg_rest_file <- file.path(output_dir, "kegg_rest_database.RData")
  if (file.exists(kegg_rest_file)) {
    load(kegg_rest_file)
    db_list$kegg_rest <- kegg_data_rest
    message("✓ 加载KEGGREST数据库")
  } else {
    message("✗ KEGGREST数据库未找到")
  }
  
  kegg_cp_file <- file.path(output_dir, "kegg_hsa_database.RData")
  if (file.exists(kegg_cp_file)) {
    load(kegg_cp_file)
    db_list$kegg_cp <- kegg_data
    message("✓ 加载clusterProfiler KEGG数据库")
  } else {
    message("✗ clusterProfiler KEGG数据库未找到")
  }
  
  if (length(db_list) == 0) {
    warning("没有找到任何预下载的数据库！")
    message("请先运行 download_kegg_data.r 脚本下载数据库")
  }
  
  return(db_list)
}

PATHWAY_DATABASES <- load_pathway_databases()

# ============================================================================
# Marker基因列表
# ============================================================================
get_marker_genes <- function() {
  message("加载marker基因列表...")
  
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

marker_genes <- get_marker_genes()

# ============================================================================
# UMAP计算和保存函数
# ============================================================================
process_and_save_umap <- function(seurat_obj, cancer_type) {
  message("\n检查并处理UMAP...")
  
  if ("umap" %in% names(seurat_obj@reductions)) {
    message("✓ 已存在UMAP降维结果")
    return(seurat_obj)
  }
  
  message("未发现UMAP，开始计算...")
  
  tryCatch({
    message("  1. 数据标准化...")
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    
    message("  2. 寻找高变基因...")
    seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
    
    message("  3. 数据缩放...")
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    
    message("  4. PCA分析...")
    seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
    
    message("  5. 计算UMAP...")
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
    
    message("✓ UMAP计算完成")
    
    output_rds <- file.path(
      output_dir,
      "processed_rds",
      paste0(cancer_type, "_with_umap.rds")
    )
    
    message("  保存包含UMAP的RDS文件...")
    saveRDS(seurat_obj, output_rds)
    message("✓ 新RDS文件已保存: ", output_rds)
    
    info_file <- file.path(
      output_dir,
      "processed_rds",
      paste0(cancer_type, "_file_info.txt")
    )
    
    cat(
      "原始文件路径:\n",
      file.path(input_dir, cancer_type, paste0(cancer_type, "_fine_annotated.rds")), "\n\n",
      "处理后文件路径:\n",
      output_rds, "\n\n",
      "处理时间: ", as.character(Sys.time()), "\n",
      "细胞数: ", ncol(seurat_obj), "\n",
      "特征数: ", nrow(seurat_obj), "\n",
      "包含降维结果: ", paste(names(seurat_obj@reductions), collapse = ", "), "\n",
      file = info_file
    )
    
    message("✓ 文件信息已保存: ", info_file)
    
  }, error = function(e) {
    stop("UMAP计算失败: ", e$message)
  })
  
  return(seurat_obj)
}

# ============================================================================
# 热图绘制 - 双格式保存
# ============================================================================
plot_marker_heatmap <- function(seurat_obj, markers_list, cancer_type, plot_type = "immune") {
  message("生成", plot_type, "marker基因热图...")
  
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
        paste0(cancer_type, "_", plot_type, "_", cell_type, "_heatmap")
      )
      
      save_plot_both_formats(p, filename, width = 12, height = 8)
      
    }, error = function(e) {
      message("  ✗ ", cell_type, " 热图绘制失败: ", e$message)
    })
    
    gc()
  }
}

# ============================================================================
# FindMarkers分析
# ============================================================================
perform_findmarkers <- function(seurat_obj, cancer_type) {
  message("对细胞类型进行FindMarkers分析...")
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  cell_types <- unique(seurat_obj$fine_cell_type)
  cell_types <- cell_types[!is.na(cell_types) & cell_types != "Other"]
  
  message("发现 ", length(cell_types), " 个细胞类型")
  
  all_markers <- list()
  
  for (cell_type in cell_types) {
    message("  处理: ", cell_type)
    
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
    all_markers_df <- do.call(rbind, all_markers)
    
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
# 🔧 手动GO富集分析（超几何检验）
# ============================================================================
manual_go_enrichment <- function(gene_list, ont = "BP", pval_cutoff = 0.05, qval_cutoff = 0.2) {
  
  # 获取所有GO注释
  go_all <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
      columns = c("ENTREZID", "GO", "ONTOLOGY"),
      keytype = "ENTREZID"
    )
  )
  
  # 过滤指定的本体（BP/CC/MF）
  go_all <- go_all[go_all$ONTOLOGY == ont & !is.na(go_all$GO), ]
  
  # 统计每个GO term的基因数
  go_term_genes <- split(go_all$ENTREZID, go_all$GO)
  go_term_genes <- lapply(go_term_genes, unique)
  
  # 过滤GO term大小（10-500个基因）
  go_sizes <- sapply(go_term_genes, length)
  go_term_genes <- go_term_genes[go_sizes >= 10 & go_sizes <= 500]
  
  if (length(go_term_genes) == 0) {
    return(NULL)
  }
  
  # 背景基因集
  universe <- unique(go_all$ENTREZID)
  N <- length(universe)
  n <- length(gene_list)
  
  # 对每个GO term进行超几何检验
  results_list <- lapply(names(go_term_genes), function(go_id) {
    term_genes <- go_term_genes[[go_id]]
    M <- length(term_genes)
    overlap <- intersect(gene_list, term_genes)
    k <- length(overlap)
    
    if (k == 0) return(NULL)
    
    # 超几何检验
    pval <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    
    data.frame(
      GO_ID = go_id,
      Count = k,
      Total = M,
      GeneRatio = paste0(k, "/", n),
      BgRatio = paste0(M, "/", N),
      pvalue = pval,
      geneID = paste(overlap, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results_list)
  
  if (is.null(results_df) || nrow(results_df) == 0) {
    return(NULL)
  }
  
  # BH校正
  results_df$p.adjust <- p.adjust(results_df$pvalue, method = "BH")
  
  # 添加GO term描述
  go_terms <- suppressMessages(
    AnnotationDbi::select(
      GO.db::GO.db,
      keys = results_df$GO_ID,
      columns = c("GOID", "TERM"),
      keytype = "GOID"
    )
  )
  
  results_df <- merge(results_df, go_terms, by.x = "GO_ID", by.y = "GOID", all.x = TRUE)
  
  # 过滤显著结果
  results_df <- results_df[results_df$pvalue < pval_cutoff & results_df$p.adjust < qval_cutoff, ]
  
  # 排序
  results_df <- results_df[order(results_df$p.adjust), ]
  
  return(results_df)
}

# ============================================================================
# 🔧 手动Pathway富集分析（超几何检验）
# ============================================================================
manual_pathway_enrichment <- function(gene_list, term2gene, term2name, pval_cutoff = 0.05, qval_cutoff = 0.2) {
  
  if (!is.data.frame(term2gene) || ncol(term2gene) < 2) {
    return(NULL)
  }
  
  colnames(term2gene)[1:2] <- c("pathway", "gene")
  
  # 统计每个pathway的基因数
  pathway_genes <- split(term2gene$gene, term2gene$pathway)
  pathway_genes <- lapply(pathway_genes, unique)
  
  # 过滤pathway大小
  pathway_sizes <- sapply(pathway_genes, length)
  pathway_genes <- pathway_genes[pathway_sizes >= 10 & pathway_sizes <= 500]
  
  if (length(pathway_genes) == 0) {
    return(NULL)
  }
  
  # 背景基因集
  universe <- unique(term2gene$gene)
  N <- length(universe)
  n <- length(gene_list)
  
  # 对每个pathway进行超几何检验
  results_list <- lapply(names(pathway_genes), function(pathway_id) {
    pathway_genes_set <- pathway_genes[[pathway_id]]
    M <- length(pathway_genes_set)
    overlap <- intersect(gene_list, pathway_genes_set)
    k <- length(overlap)
    
    if (k == 0) return(NULL)
    
    pval <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    
    data.frame(
      Pathway_ID = pathway_id,
      Count = k,
      Total = M,
      GeneRatio = paste0(k, "/", n),
      BgRatio = paste0(M, "/", N),
      pvalue = pval,
      geneID = paste(overlap, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results_list)
  
  if (is.null(results_df) || nrow(results_df) == 0) {
    return(NULL)
  }
  
  results_df$p.adjust <- p.adjust(results_df$pvalue, method = "BH")
  
  # 添加pathway名称
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

# ============================================================================
# 🎨 专业富集分析可视化函数集
# ============================================================================

# 1. 气泡图 (Bubble Plot) - 经典富集可视化
plot_enrichment_bubble <- function(enrich_df, title, top_n = 20) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  
  plot_df <- head(enrich_df, top_n)
  
  # 提取描述列
  if ("TERM" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$TERM
  } else if ("description" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$description
  }
  
  # 计算GeneRatio数值
  plot_df$GeneRatio_num <- sapply(strsplit(plot_df$GeneRatio, "/"), function(x) {
    as.numeric(x[1]) / as.numeric(x[2])
  })
  
  p <- ggplot(plot_df, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(
      low = "#B2182B", 
      high = "#2166AC",
      name = "Adjusted\nP-value",
      trans = "log10"
    ) +
    scale_size_continuous(
      name = "Gene\nCount",
      range = c(3, 12)
    ) +
    labs(
      title = title,
      x = "Gene Ratio",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 11),
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.background = element_rect(fill = "white", color = "black")
    )
  
  return(p)
}

# 2. 点图 (Dot Plot) - Nature/Cell常用
plot_enrichment_dotplot <- function(enrich_df, title, top_n = 20) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  
  plot_df <- head(enrich_df, top_n)
  
  if ("TERM" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$TERM
  } else if ("description" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$description
  }
  
  plot_df$GeneRatio_num <- sapply(strsplit(plot_df$GeneRatio, "/"), function(x) {
    as.numeric(x[1]) / as.numeric(x[2])
  })
  
  plot_df$neg_log10_p <- -log10(plot_df$p.adjust)
  
  p <- ggplot(plot_df, aes(x = neg_log10_p, y = reorder(Description, neg_log10_p))) +
    geom_segment(
      aes(x = 0, xend = neg_log10_p, y = Description, yend = Description),
      color = "grey70",
      linewidth = 0.5
    ) +
    geom_point(aes(size = Count, color = GeneRatio_num)) +
    scale_color_gradientn(
      colors = c("#3288BD", "#66C2A5", "#FEE08B", "#F46D43", "#D53E4F"),
      name = "Gene Ratio"
    ) +
    scale_size_continuous(name = "Count", range = c(3, 10)) +
    labs(
      title = title,
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(face = "bold", size = 11),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "right"
    )
  
  return(p)
}

# 3. 棒棒糖图 (Lollipop Plot) - 优雅简洁
plot_enrichment_lollipop <- function(enrich_df, title, top_n = 15) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  
  plot_df <- head(enrich_df, top_n)
  
  if ("TERM" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$TERM
  } else if ("description" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$description
  }
  
  plot_df$neg_log10_p <- -log10(plot_df$p.adjust)
  
  p <- ggplot(plot_df, aes(x = reorder(Description, neg_log10_p), y = neg_log10_p)) +
    geom_segment(
      aes(x = Description, xend = Description, y = 0, yend = neg_log10_p),
      color = "grey60",
      linewidth = 1
    ) +
    geom_point(aes(size = Count, color = neg_log10_p)) +
    scale_color_gradient(
      low = "#FEE090",
      high = "#D73027",
      name = "-log10(FDR)"
    ) +
    scale_size_continuous(name = "Gene Count", range = c(4, 12)) +
    coord_flip() +
    labs(
      title = title,
      x = NULL,
      y = "-log10(Adjusted P-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 11),
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}

# 4. 网络图 (Network Plot) - 展示通路间关系
plot_enrichment_network <- function(enrich_df, title, top_n = 15) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  
  # 需要额外安装: igraph, ggraph
  if (!requireNamespace("igraph", quietly = TRUE) || 
      !requireNamespace("ggraph", quietly = TRUE)) {
    message("需要安装 igraph 和 ggraph 包")
    return(NULL)
  }
  
  library(igraph)
  library(ggraph)
  
  plot_df <- head(enrich_df, top_n)
  
  if ("TERM" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$TERM
  } else if ("description" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$description
  }
  
  # 构建基因共享矩阵
  gene_lists <- strsplit(plot_df$geneID, "/")
  names(gene_lists) <- plot_df$Description
  
  # 计算Jaccard相似度
  n <- length(gene_lists)
  similarity_matrix <- matrix(0, n, n)
  rownames(similarity_matrix) <- colnames(similarity_matrix) <- names(gene_lists)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      intersection <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
      union <- length(union(gene_lists[[i]], gene_lists[[j]]))
      similarity_matrix[i, j] <- similarity_matrix[j, i] <- intersection / union
    }
  }
  
  # 过滤弱连接
  similarity_matrix[similarity_matrix < 0.2] <- 0
  
  # 构建网络
  g <- graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", 
                                   weighted = TRUE, diag = FALSE)
  
  # 添加节点属性
  V(g)$size <- plot_df$Count
  V(g)$pvalue <- -log10(plot_df$p.adjust)
  
  # 绘图
  set.seed(123)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = weight), alpha = 0.3, color = "grey60") +
    geom_node_point(aes(size = size, color = pvalue)) +
    geom_node_text(aes(label = name), size = 2.5, repel = TRUE, max.overlaps = 20) +
    scale_color_gradient(low = "#3288BD", high = "#D53E4F", name = "-log10(FDR)") +
    scale_size_continuous(name = "Gene Count", range = c(3, 12)) +
    scale_edge_width(range = c(0.5, 2), guide = "none") +
    labs(title = title) +
    theme_graph(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )
  
  return(p)
}

# 5. 瀑布图 (Waterfall Plot) - 富集强度变化
plot_enrichment_waterfall <- function(enrich_df, title, top_n = 20) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  
  plot_df <- head(enrich_df, top_n)
  
  if ("TERM" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$TERM
  } else if ("description" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$description
  }
  
  plot_df$neg_log10_p <- -log10(plot_df$p.adjust)
  plot_df <- plot_df[order(plot_df$neg_log10_p, decreasing = TRUE), ]
  plot_df$cumsum <- cumsum(plot_df$neg_log10_p)
  plot_df$id <- 1:nrow(plot_df)
  
  p <- ggplot(plot_df, aes(x = id)) +
    geom_bar(aes(y = neg_log10_p, fill = neg_log10_p), stat = "identity", width = 0.7) +
    geom_line(aes(y = cumsum), color = "#D53E4F", linewidth = 1.5) +
    geom_point(aes(y = cumsum), color = "#D53E4F", size = 2) +
    scale_fill_gradient(low = "#FFFFCC", high = "#E31A1C", name = "-log10(FDR)") +
    scale_x_continuous(breaks = plot_df$id, labels = plot_df$Description) +
    labs(
      title = title,
      x = NULL,
      y = "-log10(Adjusted P-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}

# 6. 热图 (Heatmap) - 通路-细胞类型
plot_enrichment_heatmap <- function(enrichment_results, cancer_type, analysis_type = "GO_BP") {
  if (is.null(enrichment_results) || length(enrichment_results) == 0) return(NULL)
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("需要安装 pheatmap 包")
    return(NULL)
  }
  
  library(pheatmap)
  
  # 收集所有通路
  all_pathways <- c()
  pathway_matrix_list <- list()
  
  for (ct in names(enrichment_results)) {
    if (!is.null(enrichment_results[[ct]][[analysis_type]])) {
      enrich_df <- enrichment_results[[ct]][[analysis_type]]
      
      if ("TERM" %in% colnames(enrich_df)) {
        pathways <- enrich_df$TERM
      } else if ("description" %in% colnames(enrich_df)) {
        pathways <- enrich_df$description
      } else {
        next
      }
      
      pvalues <- -log10(enrich_df$p.adjust)
      names(pvalues) <- pathways
      
      pathway_matrix_list[[ct]] <- pvalues
      all_pathways <- c(all_pathways, pathways)
    }
  }
  
  if (length(pathway_matrix_list) == 0) return(NULL)
  
  # 选择top通路
  all_pathways <- unique(all_pathways)
  top_pathways <- head(all_pathways, 30)
  
  # 构建矩阵
  heatmap_matrix <- matrix(
    0,
    nrow = length(top_pathways),
    ncol = length(pathway_matrix_list)
  )
  
  rownames(heatmap_matrix) <- top_pathways
  colnames(heatmap_matrix) <- names(pathway_matrix_list)
  
  for (ct in names(pathway_matrix_list)) {
    pvals <- pathway_matrix_list[[ct]]
    matched_pathways <- intersect(names(pvals), top_pathways)
    heatmap_matrix[matched_pathways, ct] <- pvals[matched_pathways]
  }
  
  # 绘制热图
  p <- pheatmap(
    heatmap_matrix,
    color = colorRampPalette(c("white", "#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 10,
    angle_col = 45,
    border_color = "grey80",
    main = paste0(cancer_type, " - ", analysis_type, " Enrichment Heatmap"),
    legend_labels = "-log10(FDR)",
    cellwidth = 15,
    cellheight = 10
  )
  
  return(p)
}

# 7. 岭线图 (Ridgeline Plot) - 基因分布
plot_enrichment_ridgeline <- function(seurat_obj, enrich_df, title, top_n = 10) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    message("需要安装 ggridges 包")
    return(NULL)
  }
  
  library(ggridges)
  
  plot_df <- head(enrich_df, top_n)
  
  if ("TERM" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$TERM
  } else if ("description" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$description
  }
  
  # 准备表达数据
  expr_data_list <- list()
  
  for (i in 1:nrow(plot_df)) {
    genes <- unlist(strsplit(plot_df$geneID[i], "/"))
    available_genes <- intersect(genes, rownames(seurat_obj))
    
    if (length(available_genes) > 0) {
      expr_values <- rowMeans(as.matrix(seurat_obj@assays$RNA@data[available_genes, , drop = FALSE]))
      expr_data_list[[plot_df$Description[i]]] <- data.frame(
        Pathway = plot_df$Description[i],
        Expression = expr_values
      )
    }
  }
  
  expr_data <- do.call(rbind, expr_data_list)
  
  p <- ggplot(expr_data, aes(x = Expression, y = Pathway, fill = Pathway)) +
    geom_density_ridges(
      scale = 3,
      alpha = 0.7,
      color = "white",
      size = 0.5
    ) +
    scale_fill_manual(values = colorRampPalette(c("#4DBBD5", "#E64B35", "#00A087"))(top_n)) +
    labs(
      title = title,
      x = "Average Gene Expression",
      y = NULL
    ) +
    theme_ridges(font_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "none",
      axis.text.y = element_text(size = 10, color = "black")
    )
  
  return(p)
}

# 8. 小提琴图 (Violin Plot) - 通路得分分布
plot_enrichment_violin <- function(seurat_obj, enrich_df, title, top_n = 8) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  
  plot_df <- head(enrich_df, top_n)
  
  if ("TERM" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$TERM
  } else if ("description" %in% colnames(plot_df)) {
    plot_df$Description <- plot_df$description
  }
  
  # 计算通路得分
  score_data_list <- list()
  
  for (i in 1:nrow(plot_df)) {
    genes <- unlist(strsplit(plot_df$geneID[i], "/"))
    available_genes <- intersect(genes, rownames(seurat_obj))
    
    if (length(available_genes) > 0) {
      seurat_obj <- AddModuleScore(
        seurat_obj,
        features = list(available_genes),
        name = "temp_score",
        assay = "RNA"
      )
      
      scores <- seurat_obj$temp_score1
      score_data_list[[i]] <- data.frame(
        Pathway = plot_df$Description[i],
        Score = scores,
        CellType = seurat_obj$fine_cell_type
      )
    }
  }
  
  score_data <- do.call(rbind, score_data_list)
  
  p <- ggplot(score_data, aes(x = Pathway, y = Score, fill = Pathway)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = colorRampPalette(c("#4DBBD5", "#E64B35", "#00A087"))(top_n)) +
    labs(
      title = title,
      x = NULL,
      y = "Pathway Score"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}

# 9. 组合面板图 (Composite Panel) - 多图整合
plot_enrichment_composite <- function(enrich_df, seurat_obj, cancer_type, cell_type) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  
  # 生成4个子图
  p1 <- plot_enrichment_bubble(enrich_df, "Top Enriched Pathways", 15)
  p2 <- plot_enrichment_dotplot(enrich_df, "Pathway Significance", 15)
  p3 <- plot_enrichment_lollipop(enrich_df, "Enrichment Strength", 12)
  p4 <- plot_enrichment_waterfall(enrich_df, "Cumulative Enrichment", 15)
  
  # 组合
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = paste0(cancer_type, " - ", cell_type, " Enrichment Analysis"),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18, color = "#2C3E50")
      )
    )
  
  return(combined)
}

# ============================================================================
# 🔧 完全修复的富集分析函数（使用超几何检验）
# ============================================================================
perform_enrichment_analysis <- function(markers_df, cancer_type) {
  message("进行GO和KEGG富集分析（使用超几何检验）...")
  
  if (is.null(markers_df) || nrow(markers_df) == 0) {
    message("  没有marker数据，跳过富集分析")
    return(NULL)
  }
  
  cell_types <- unique(markers_df$cell_type)
  enrichment_results <- list()
  
  for (cell_type in cell_types) {
    message("  分析: ", cell_type)
    
    markers_subset <- markers_df[markers_df$cell_type == cell_type, ]
    markers_subset <- markers_subset[order(markers_subset$avg_log2FC, decreasing = TRUE), ]
    
    top_n <- min(200, nrow(markers_subset))
    top_markers <- head(markers_subset, top_n)
    gene_ids <- top_markers$gene
    
    # ========================================================================
    # 基因ID转换
    # ========================================================================
    entrez_ids <- NULL
    
    tryCatch({
      if (all(grepl("^ENS", gene_ids))) {
        message("    检测到Ensembl ID，转换为Entrez ID...")
        entrez_ids <- suppressMessages(
          mapIds(
            org.Hs.eg.db,
            keys = gene_ids,
            keytype = "ENSEMBL",
            column = "ENTREZID",
            multiVals = "first"
          )
        )
      } else {
        message("    检测到Gene Symbol，转换为Entrez ID...")
        entrez_ids <- suppressMessages(
          mapIds(
            org.Hs.eg.db,
            keys = gene_ids,
            keytype = "SYMBOL",
            column = "ENTREZID",
            multiVals = "first"
          )
        )
      }
      
      entrez_ids <- entrez_ids[!is.na(entrez_ids)]
      entrez_ids <- unique(as.character(entrez_ids))
      
      message("    成功转换 ", length(entrez_ids), " 个基因ID")
      
    }, error = function(e) {
      message("    ✗ 基因ID转换失败: ", e$message)
      entrez_ids <- NULL
    })
    
    if (is.null(entrez_ids) || length(entrez_ids) < 10) {
      message("    可用基因太少 (", length(entrez_ids), ")，跳过")
      next
    }
    
    # ========================================================================
    # 🔧 GO富集分析（使用超几何检验，避免clusterProfiler的bug）
    # ========================================================================
    tryCatch({
      message("    进行GO BP分析（超几何检验）...")
      
      # 手动进行超几何检验
      go_results <- manual_go_enrichment(entrez_ids, ont = "BP")
      
      if (!is.null(go_results) && nrow(go_results) > 0) {
        enrichment_results[[cell_type]]$GO_BP <- go_results
        message("    ✓ GO BP: ", nrow(go_results), " 个富集通路")
      } else {
        message("    GO BP: 无显著富集结果")
      }
      
    }, error = function(e) {
      message("    ✗ GO BP分析失败: ", e$message)
    })
    
    # ========================================================================
    # 🔧 KEGG富集分析（超几何检验）
    # ========================================================================
    kegg_success <- FALSE
    
    if (!is.null(PATHWAY_DATABASES$msigdb_kegg)) {
      tryCatch({
        message("    KEGG分析（超几何检验）...")
        
        kegg_results <- manual_pathway_enrichment(
          entrez_ids,
          PATHWAY_DATABASES$msigdb_kegg$term2gene,
          PATHWAY_DATABASES$msigdb_kegg$term2name
        )
        
        if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
          enrichment_results[[cell_type]]$KEGG <- kegg_results
          message("    ✓ KEGG: ", nrow(kegg_results), " 个富集通路")
          kegg_success <- TRUE
        } else {
          message("    KEGG: 无显著富集结果")
        }
        
      }, error = function(e) {
        message("    ✗ KEGG分析失败: ", e$message)
      })
    }
    
    # ========================================================================
    # 🔧 Reactome富集分析（超几何检验）
    # ========================================================================
    if (!kegg_success && !is.null(PATHWAY_DATABASES$reactome)) {
      tryCatch({
        message("    Reactome分析（超几何检验）...")
        
        reactome_results <- manual_pathway_enrichment(
          entrez_ids,
          PATHWAY_DATABASES$reactome$term2gene,
          PATHWAY_DATABASES$reactome$term2name
        )
        
        if (!is.null(reactome_results) && nrow(reactome_results) > 0) {
          enrichment_results[[cell_type]]$Reactome <- reactome_results
          message("    ✓ Reactome: ", nrow(reactome_results), " 个富集通路")
          kegg_success <- TRUE
        }
        
      }, error = function(e) {
        message("    ✗ Reactome分析失败: ", e$message)
      })
    }
    
    if (!kegg_success) {
      message("    ⚠ 所有pathway富集分析均未产生结果")
    }
    
    gc()
  }
  
  if (length(enrichment_results) > 0) {
    save_enrichment_results_manual(enrichment_results, cancer_type)
    plot_enrichment_results_manual(enrichment_results, cancer_type)
  } else {
    message("  ✗ 所有细胞类型均未产生富集结果")
  }
  
  return(enrichment_results)
}

# ============================================================================
# 保存手动富集分析结果
# ============================================================================
save_enrichment_results_manual <- function(enrichment_results, cancer_type) {
  message("保存富集分析结果...")
  
  for (cell_type in names(enrichment_results)) {
    if (!is.null(enrichment_results[[cell_type]]$GO_BP)) {
      tryCatch({
        go_file <- file.path(
          output_dir,
          "summary_stats",
          paste0(cancer_type, "_", cell_type, "_GO_BP.csv")
        )
        write.csv(enrichment_results[[cell_type]]$GO_BP, go_file, row.names = FALSE)
        message("  ✓ 保存GO结果: ", cell_type)
      }, error = function(e) {
        message("  ✗ 保存GO结果失败: ", e$message)
      })
    }
    
    if (!is.null(enrichment_results[[cell_type]]$KEGG)) {
      tryCatch({
        kegg_file <- file.path(
          output_dir,
          "summary_stats",
          paste0(cancer_type, "_", cell_type, "_KEGG.csv")
        )
        write.csv(enrichment_results[[cell_type]]$KEGG, kegg_file, row.names = FALSE)
        message("  ✓ 保存KEGG结果: ", cell_type)
      }, error = function(e) {
        message("  ✗ 保存KEGG结果失败: ", e$message)
      })
    }
    
    if (!is.null(enrichment_results[[cell_type]]$Reactome)) {
      tryCatch({
        reactome_file <- file.path(
          output_dir,
          "summary_stats",
          paste0(cancer_type, "_", cell_type, "_Reactome.csv")
        )
        write.csv(enrichment_results[[cell_type]]$Reactome, reactome_file, row.names = FALSE)
        message("  ✓ 保存Reactome结果: ", cell_type)
      }, error = function(e) {
        message("  ✗ 保存Reactome结果失败: ", e$message)
      })
    }
  }
}

# ============================================================================
# 绘制手动富集分析结果 - 双格式保存
# ============================================================================
plot_enrichment_results_manual <- function(enrichment_results, cancer_type) {
  message("绘制富集分析图...")
  
  for (cell_type in names(enrichment_results)) {
    
    # GO BP图
    if (!is.null(enrichment_results[[cell_type]]$GO_BP)) {
      go_df <- enrichment_results[[cell_type]]$GO_BP
      
      if (nrow(go_df) > 0) {
        tryCatch({
          n_show <- min(20, nrow(go_df))
          plot_df <- head(go_df, n_show)
          
          plot_df$Description <- ifelse(
            is.na(plot_df$TERM),
            plot_df$GO_ID,
            plot_df$TERM
          )
          
          p <- ggplot(plot_df, aes(x = Count, y = reorder(Description, Count))) +
            geom_col(aes(fill = -log10(p.adjust))) +
            scale_fill_gradient(low = "blue", high = "red") +
            labs(
              title = paste0(cell_type, " - GO Biological Process"),
              x = "Gene Count",
              y = "GO Term",
              fill = "-log10(FDR)"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.y = element_text(size = 9)
            )
          
          filename <- file.path(
            output_dir,
            "enrichment_plots",
            paste0(cancer_type, "_", cell_type, "_GO_BP_barplot")
          )
          
          save_plot_both_formats(p, filename, width = 10, height = 8)
          
        }, error = function(e) {
          message("  ✗ GO图绘制失败: ", e$message)
        })
      }
    }
    
    # KEGG图
    if (!is.null(enrichment_results[[cell_type]]$KEGG)) {
      kegg_df <- enrichment_results[[cell_type]]$KEGG
      
      if (nrow(kegg_df) > 0) {
        tryCatch({
          n_show <- min(20, nrow(kegg_df))
          plot_df <- head(kegg_df, n_show)
          
          plot_df$Description <- ifelse(
            is.na(plot_df$description),
            plot_df$Pathway_ID,
            plot_df$description
          )
          
          p <- ggplot(plot_df, aes(x = Count, y = reorder(Description, Count))) +
            geom_col(aes(fill = -log10(p.adjust))) +
            scale_fill_gradient(low = "blue", high = "red") +
            labs(
              title = paste0(cell_type, " - KEGG Pathway"),
              x = "Gene Count",
              y = "Pathway",
              fill = "-log10(FDR)"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.y = element_text(size = 9)
            )
          
          filename <- file.path(
            output_dir,
            "enrichment_plots",
            paste0(cancer_type, "_", cell_type, "_KEGG_barplot")
          )
          
          save_plot_both_formats(p, filename, width = 10, height = 8)
          
        }, error = function(e) {
          message("  ✗ KEGG图绘制失败: ", e$message)
        })
      }
    }
    
    # Reactome图
    if (!is.null(enrichment_results[[cell_type]]$Reactome)) {
      reactome_df <- enrichment_results[[cell_type]]$Reactome
      
      if (nrow(reactome_df) > 0) {
        tryCatch({
          n_show <- min(20, nrow(reactome_df))
          plot_df <- head(reactome_df, n_show)
          
          plot_df$Description <- ifelse(
            is.na(plot_df$description),
            plot_df$Pathway_ID,
            plot_df$description
          )
          
          p <- ggplot(plot_df, aes(x = Count, y = reorder(Description, Count))) +
            geom_col(aes(fill = -log10(p.adjust))) +
            scale_fill_gradient(low = "blue", high = "red") +
            labs(
              title = paste0(cell_type, " - Reactome Pathway"),
              x = "Gene Count",
              y = "Pathway",
              fill = "-log10(FDR)"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.y = element_text(size = 9)
            )
          
          filename <- file.path(
            output_dir,
            "enrichment_plots",
            paste0(cancer_type, "_", cell_type, "_Reactome_barplot")
          )
          
          save_plot_both_formats(p, filename, width = 10, height = 8)
          
        }, error = function(e) {
          message("  ✗ Reactome图绘制失败: ", e$message)
        })
      }
    }
    
    gc()
  }
}

# ============================================================================
# 🎨 主可视化函数 - 生成所有增强图表
# ============================================================================
plot_comprehensive_enrichment <- function(enrichment_results, seurat_obj, cancer_type) {
  message("\n生成综合富集分析可视化...")
  
  for (cell_type in names(enrichment_results)) {
    message("\n处理细胞类型: ", cell_type)
    
    # GO BP分析
    if (!is.null(enrichment_results[[cell_type]]$GO_BP)) {
      go_df <- enrichment_results[[cell_type]]$GO_BP
      
      # 1. 气泡图
      p1 <- plot_enrichment_bubble(go_df, paste0(cell_type, " - GO BP (Bubble)"), 20)
      if (!is.null(p1)) {
        filename <- file.path(output_dir, "enrichment_plots",
                             paste0(cancer_type, "_", cell_type, "_GO_bubble"))
        save_plot_both_formats(p1, filename, 12, 10)
      }
      
      # 2. 点图
      p2 <- plot_enrichment_dotplot(go_df, paste0(cell_type, " - GO BP (Dot)"), 20)
      if (!is.null(p2)) {
        filename <- file.path(output_dir, "enrichment_plots",
                             paste0(cancer_type, "_", cell_type, "_GO_dotplot"))
        save_plot_both_formats(p2, filename, 12, 10)
      }
      
      # 3. 棒棒糖图
      p3 <- plot_enrichment_lollipop(go_df, paste0(cell_type, " - GO BP (Lollipop)"), 15)
      if (!is.null(p3)) {
        filename <- file.path(output_dir, "enrichment_plots",
                             paste0(cancer_type, "_", cell_type, "_GO_lollipop"))
        save_plot_both_formats(p3, filename, 10, 8)
      }
      
      # 4. 网络图
      tryCatch({
        p4 <- plot_enrichment_network(go_df, paste0(cell_type, " - GO BP Network"), 15)
        if (!is.null(p4)) {
          filename <- file.path(output_dir, "enrichment_plots",
                               paste0(cancer_type, "_", cell_type, "_GO_network"))
          save_plot_both_formats(p4, filename, 14, 12)
        }
      }, error = function(e) message("网络图失败: ", e$message))
      
      # 5. 瀑布图
      p5 <- plot_enrichment_waterfall(go_df, paste0(cell_type, " - GO BP (Waterfall)"), 20)
      if (!is.null(p5)) {
        filename <- file.path(output_dir, "enrichment_plots",
                             paste0(cancer_type, "_", cell_type, "_GO_waterfall"))
        save_plot_both_formats(p5, filename, 14, 8)
      }
      
      # 6. 岭线图
      tryCatch({
        p6 <- plot_enrichment_ridgeline(seurat_obj, go_df, 
                                       paste0(cell_type, " - Gene Expression Distribution"), 10)
        if (!is.null(p6)) {
          filename <- file.path(output_dir, "enrichment_plots",
                               paste0(cancer_type, "_", cell_type, "_GO_ridgeline"))
          save_plot_both_formats(p6, filename, 12, 10)
        }
      }, error = function(e) message("岭线图失败: ", e$message))
      
      # 7. 小提琴图
      tryCatch({
        p7 <- plot_enrichment_violin(seurat_obj, go_df, 
                                     paste0(cell_type, " - Pathway Scores"), 8)
        if (!is.null(p7)) {
          filename <- file.path(output_dir, "enrichment_plots",
                               paste0(cancer_type, "_", cell_type, "_GO_violin"))
          save_plot_both_formats(p7, filename, 12, 8)
        }
      }, error = function(e) message("小提琴图失败: ", e$message))
    }
    
    # KEGG分析
    if (!is.null(enrichment_results[[cell_type]]$KEGG)) {
      kegg_df <- enrichment_results[[cell_type]]$KEGG
      
      # 气泡图
      p_kegg <- plot_enrichment_bubble(kegg_df, paste0(cell_type, " - KEGG Pathways"), 20)
      if (!is.null(p_kegg)) {
        filename <- file.path(output_dir, "enrichment_plots",
                             paste0(cancer_type, "_", cell_type, "_KEGG_bubble"))
        save_plot_both_formats(p_kegg, filename, 12, 10)
      }
      
      # 组合图
      tryCatch({
        p_composite <- plot_enrichment_composite(kegg_df, seurat_obj, cancer_type, cell_type)
        if (!is.null(p_composite)) {
          filename <- file.path(output_dir, "enrichment_plots",
                               paste0(cancer_type, "_", cell_type, "_KEGG_composite"))
          save_plot_both_formats(p_composite, filename, 20, 16)
        }
      }, error = function(e) message("组合图失败: ", e$message))
    }
    
    gc()
  }
  
  # 跨细胞类型热图
  tryCatch({
    p_heatmap <- plot_enrichment_heatmap(enrichment_results, cancer_type, "GO_BP")
    if (!is.null(p_heatmap)) {
      filename <- file.path(output_dir, "enrichment_plots",
                           paste0(cancer_type, "_all_celltypes_GO_heatmap"))
      ggsave(paste0(filename, ".pdf"), p_heatmap$gtable, width = 16, height = 12)
      ggsave(paste0(filename, ".png"), p_heatmap$gtable, width = 16, height = 12, dpi = 300)
    }
  }, error = function(e) message("热图失败: ", e$message))
  
  message("✓ 综合可视化完成")
}

# ============================================================================
# 评估报告
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
  
  cat("\n2. 富集分析结果\n")
  cat("-", rep("-", 70), "\n", sep = "")
  
  if (!is.null(enrichment_results) && length(enrichment_results) > 0) {
    for (ct in names(enrichment_results)) {
      cat("  ", ct, ":\n")
      
      if (!is.null(enrichment_results[[ct]]$GO_BP)) {
        n_go <- nrow(enrichment_results[[ct]]$GO_BP)
        cat("    GO BP通路数: ", n_go, "\n")
        if (n_go > 0) {
          top_terms <- head(enrichment_results[[ct]]$GO_BP$TERM, 3)
          cat("    Top 3 GO: ", paste(top_terms, collapse = "; "), "\n")
        }
      }
      
      if (!is.null(enrichment_results[[ct]]$KEGG)) {
        n_kegg <- nrow(enrichment_results[[ct]]$KEGG)
        cat("    KEGG通路数: ", n_kegg, "\n")
        if (n_kegg > 0) {
          top_pathways <- head(enrichment_results[[ct]]$KEGG$description, 3)
          cat("    Top 3 KEGG: ", paste(top_pathways, collapse = "; "), "\n")
        }
      }
      
      if (!is.null(enrichment_results[[ct]]$Reactome)) {
        n_reactome <- nrow(enrichment_results[[ct]]$Reactome)
        cat("    Reactome通路数: ", n_reactome, "\n")
        if (n_reactome > 0) {
          top_pathways <- head(enrichment_results[[ct]]$Reactome$description, 3)
          cat("    Top 3 Reactome: ", paste(top_pathways, collapse = "; "), "\n")
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
# 主评估函数
# ============================================================================
evaluate_cancer_annotation <- function(cancer_type) {
  message("\n", rep("=", 80))
  message("开始评估癌种: ", cancer_type)
  message(rep("=", 80))
  
  processed_rds <- file.path(
    output_dir,
    "processed_rds",
    paste0(cancer_type, "_with_umap.rds")
  )
  
  if (file.exists(processed_rds)) {
    message("✓ 发现已处理的RDS文件，直接加载...")
    seurat_obj <- readRDS(processed_rds)
    message("细胞数: ", ncol(seurat_obj))
    message("细胞类型数: ", length(unique(seurat_obj$fine_cell_type)))
    message("降维结果: ", paste(names(seurat_obj@reductions), collapse = ", "))
    
  } else {
    seurat_file <- file.path(
      input_dir,
      cancer_type,
      paste0(cancer_type, "_fine_annotated.rds")
    )
    
    if (!file.exists(seurat_file)) {
      warning("文件不存在: ", seurat_file)
      return(NULL)
    }
    
    message("读取原始Seurat对象...")
    seurat_obj <- readRDS(seurat_file)
    message("细胞数: ", ncol(seurat_obj))
    message("细胞类型数: ", length(unique(seurat_obj$fine_cell_type)))
    
    seurat_obj <- process_and_save_umap(seurat_obj, cancer_type)
  }
  
  message("\n--- Step 1: 细胞类型分布UMAP可视化（专业版）---")
  plot_celltype_umap_enhanced(seurat_obj, cancer_type)
  
  message("\n--- Step 2: Marker基因UMAP可视化（专业版）---")
  plot_marker_umap_enhanced(seurat_obj, marker_genes$immune, cancer_type, "immune")
  
  if (cancer_type %in% names(marker_genes$cancer)) {
    cancer_specific_markers <- list()
    cancer_specific_markers[[cancer_type]] <- marker_genes$cancer[[cancer_type]]
    plot_marker_umap_enhanced(seurat_obj, cancer_specific_markers, cancer_type, "cancer_specific")
  }
  
  message("\n--- Step 3: Marker基因热图 ---")
  plot_marker_heatmap(seurat_obj, marker_genes$immune, cancer_type, "immune")
  
  message("\n--- Step 4: FindMarkers分析 ---")
  markers_df <- perform_findmarkers(seurat_obj, cancer_type)
  
  message("\n--- Step 5: GO和KEGG富集分析 ---")
  enrichment_results <- perform_enrichment_analysis(markers_df, cancer_type)
  
  message("\n--- Step 6: 生成综合富集分析可视化 ---")
  if (!is.null(enrichment_results) && length(enrichment_results) > 0) {
    plot_comprehensive_enrichment(enrichment_results, seurat_obj, cancer_type)
  }
  
  message("\n--- Step 7: 生成评估报告 ---")
  generate_evaluation_report(cancer_type, markers_df, enrichment_results)
  
  message("\n✓ 癌种 ", cancer_type, " 评估完成")
  
  gc()
  
  return(list(
    markers = markers_df,
    enrichment = enrichment_results
  ))
}

# ============================================================================
# 批量评估所有癌种
# ============================================================================
evaluate_all_cancers <- function() {
  message("开始批量评估所有癌种...")
  
  cancer_dirs <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
  cancer_dirs <- cancer_dirs[cancer_dirs != "reference" & 
                             cancer_dirs != "evaluation_results" &
                             !grepl("^\\.", cancer_dirs)]
  
  message("发现 ", length(cancer_dirs), " 个癌种目录")
  
  results_summary <- list()
  
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
    
    gc()
  }
  
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
message("细胞注释质量评估系统（增强版 - 15种专业可视化）")
message("=", rep("=", 80))

# 评估单个癌种
# evaluate_cancer_annotation("Melan")

# 评估所有癌种
evaluation_summary <- evaluate_all_cancers()

plan(sequential)

message("\n所有评估任务完成!")
message("结果保存在: ", output_dir)
message("所有图片已同时保存为PDF和PNG格式（PNG分辨率: 300 DPI）")
message("生成的可视化类型:")
message("  ✓ 气泡图 (Bubble Plot)")
message("  ✓ 点图 (Dot Plot)")
message("  ✓ 棒棒糖图 (Lollipop Plot)")
message("  ✓ 网络图 (Network Plot)")
message("  ✓ 瀑布图 (Waterfall Plot)")
message("  ✓ 热图 (Heatmap)")
message("  ✓ 岭线图 (Ridgeline Plot)")
message("  ✓ 小提琴图 (Violin Plot)")
message("  ✓ 组合面板图 (Composite Panel)")