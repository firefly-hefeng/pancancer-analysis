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
# 参数设置 - 修改部分
# ============================================================================
# 降维数据源目录（包含所有 *_with_umap.rds 文件）
umap_data_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/annotation_evaluation_results_11_9/processed_rds"

# 输出目录
output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/annotation_evaluation_results_major_11_12"

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
# 🎨 Major细胞类型专用配色方案
# ============================================================================
generate_major_celltype_colors <- function(cell_types) {
  # 定义major细胞类型的标准配色
  standard_colors <- c(
    "T cell" = "#E64B35",
    "T cells" = "#E64B35",
    "CD4 T" = "#F39B7F",
    "CD8 T" = "#DC0000",
    "B cell" = "#4DBBD5",
    "B cells" = "#4DBBD5",
    "Plasma" = "#91D1C2",
    "NK cell" = "#00A087",
    "NK cells" = "#00A087",
    "Myeloid" = "#3C5488",
    "Myeloid cell" = "#3C5488",
    "Myeloid cells" = "#3C5488",
    "Monocyte" = "#6A8AC7",
    "Macrophage" = "#4A708B",
    "DC" = "#9BB7D4",
    "Dendritic" = "#9BB7D4",
    "Neutrophil" = "#F39B7F",
    "Mast" = "#8491B4",
    "Mast cell" = "#8491B4",
    "Epithelial" = "#E64B35",
    "Epithelial cell" = "#E64B35",
    "Cancer" = "#B22222",
    "Tumor" = "#B22222",
    "Malignant" = "#8B0000",
    "Endothelial" = "#7E6148",
    "Endothelial cell" = "#7E6148",
    "Fibroblast" = "#00A087",
    "Stromal" = "#469E7A",
    "CAF" = "#7CB992",
    "Other" = "#CCCCCC",
    "Unknown" = "#999999"
  )
  
  # 为所有细胞类型分配颜色
  color_assignment <- rep("#CCCCCC", length(cell_types))
  names(color_assignment) <- cell_types
  
  for (i in seq_along(cell_types)) {
    ct <- cell_types[i]
    
    # 精确匹配
    if (ct %in% names(standard_colors)) {
      color_assignment[ct] <- standard_colors[ct]
      next
    }
    
    # 模糊匹配
    matched <- FALSE
    for (std_type in names(standard_colors)) {
      if (grepl(std_type, ct, ignore.case = TRUE) || 
          grepl(ct, std_type, ignore.case = TRUE)) {
        color_assignment[ct] <- standard_colors[std_type]
        matched <- TRUE
        break
      }
    }
    
    # 如果没有匹配，使用基于索引的颜色
    if (!matched) {
      base_palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", 
                       "#F39B7F", "#8491B4", "#7E6148", "#B09C85")
      color_assignment[ct] <- base_palette[((i - 1) %% length(base_palette)) + 1]
    }
  }
  
  return(color_assignment)
}

# ============================================================================
# 🎨 优化的Major细胞类型UMAP图 - 双格式保存
# ============================================================================
plot_major_celltype_umap_enhanced <- function(seurat_obj, cancer_type) {
  message("生成Major细胞类型分布UMAP图...")
  
  cell_types <- unique(seurat_obj$major_cell_type)
  cell_types <- cell_types[!is.na(cell_types)]
  n_celltypes <- length(cell_types)
  
  message("Major细胞类型数量: ", n_celltypes)
  
  colors <- generate_major_celltype_colors(cell_types)
  
  umap_coords <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  umap_coords$cell_type <- seurat_obj$major_cell_type
  
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
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = colors, name = "Major Cell Type") +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text = element_text(color = "black", size = 14),
      axis.title = element_text(face = "bold", size = 16)
    ) +
    labs(
      title = paste0(cancer_type, " - Major Cell Type Distribution"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    geom_text(
      data = label_data,
      aes(x = UMAP_1, y = UMAP_2, label = cell_type),
      size = 4,
      fontface = "bold",
      color = "black",
      bg.color = "white",
      bg.r = 0.15,
      check_overlap = TRUE
    )
  
  filename_main <- file.path(
    output_dir,
    "umap_plots",
    paste0(cancer_type, "_major_celltype_umap_main")
  )
  
  save_plot_both_formats(p_main, filename_main, width = 14, height = 12)
  
  # 带图例版本
  max_per_col <- 20
  n_cols <- ceiling(n_celltypes / max_per_col)
  
  p_with_legend <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = colors, name = "Major Cell Type") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.title = element_text(face = "bold", size = 13),
      legend.key.size = unit(0.6, "cm"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", size = 14)
    ) +
    labs(
      title = paste0(cancer_type, " - Major Cell Type Distribution"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 5, alpha = 1),
        ncol = n_cols,
        byrow = FALSE
      )
    )
  
  filename_legend <- file.path(
    output_dir,
    "umap_plots",
    paste0(cancer_type, "_major_celltype_umap_with_legend")
  )
  
  legend_width <- 14 + (n_cols - 1) * 2.5
  
  save_plot_both_formats(p_with_legend, filename_legend, width = legend_width, height = 12)
  
  # 分面图（如果细胞类型较多）
  if (n_celltypes > 15) {
    message("Major细胞类型较多，生成分面图...")
    
    p_facet <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(
        data = umap_coords,
        color = "gray90",
        size = 0.1,
        alpha = 0.3
      ) +
      geom_point(
        aes(color = cell_type),
        size = 0.8,
        alpha = 0.8
      ) +
      scale_color_manual(values = colors) +
      facet_wrap(~ cell_type, ncol = 4) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "gray95", color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(size = 9)
      ) +
      labs(
        title = paste0(cancer_type, " - Major Cell Types (Faceted)"),
        x = "UMAP 1",
        y = "UMAP 2"
      )
    
    filename_facet <- file.path(
      output_dir,
      "umap_plots",
      paste0(cancer_type, "_major_celltype_umap_faceted")
    )
    
    save_plot_both_formats(p_facet, filename_facet, width = 18, height = 14)
  }
  
  return(p_main)
}

# ============================================================================
# 🎨 Major细胞类型Marker基因UMAP图 - 双格式保存
# ============================================================================
plot_major_marker_umap_enhanced <- function(seurat_obj, markers_list, cancer_type, plot_type = "major") {
  message("生成", plot_type, " Major marker基因UMAP图...")
  
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
        pt.size = 0.8,
        order = TRUE
      ) + 
        scale_color_gradientn(
          colors = c("#F0F0F0", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C", "#08306B"),
          name = "Score",
          guide = guide_colorbar(
            barwidth = 1.2,
            barheight = 10,
            title.position = "top",
            title.hjust = 0.5
          )
        ) +
        ggtitle(paste0(cell_type, " Marker Expression")) +
        theme_minimal(base_size = 16) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18, color = "#2C3E50"),
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(face = "bold", size = 14)
        ) +
        labs(x = "UMAP 1", y = "UMAP 2")
      
      plot_list[[cell_type]] <- p
      
      filename <- file.path(
        output_dir, 
        "umap_plots",
        paste0(cancer_type, "_major_", plot_type, "_", cell_type, "_umap")
      )
      
      save_plot_both_formats(p, filename, width = 11, height = 9)
      
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
        title = paste0(cancer_type, " - Major ", toupper(plot_type), " Marker Expression"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 22, color = "#2C3E50")
        )
      )
    
    combined_filename <- file.path(
      output_dir,
      "umap_plots",
      paste0(cancer_type, "_major_", plot_type, "_markers_combined")
    )
    
    save_plot_both_formats(
      combined_plot,
      combined_filename,
      width = ncol * 9,
      height = nrow * 8
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
  
  # 从annotation_evaluation_results_11_9目录加载
  base_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/annotation_evaluation_results_11_9"
  
  msigdb_file <- file.path(base_dir, "msigdb_kegg_genesets.RData")
  if (file.exists(msigdb_file)) {
    load(msigdb_file)
    db_list$msigdb_kegg <- msigdb_data
    message("✓ 加载MSigDB KEGG数据库")
    message("  - Pathways: ", nrow(msigdb_data$term2name))
  } else {
    message("✗ MSigDB KEGG数据库未找到")
  }
  
  reactome_file <- file.path(base_dir, "reactome_genesets.RData")
  if (file.exists(reactome_file)) {
    load(reactome_file)
    db_list$reactome <- reactome_data
    message("✓ 加载Reactome数据库")
    message("  - Pathways: ", nrow(reactome_data$term2name))
  } else {
    message("✗ Reactome数据库未找到")
  }
  
  kegg_rest_file <- file.path(base_dir, "kegg_rest_database.RData")
  if (file.exists(kegg_rest_file)) {
    load(kegg_rest_file)
    db_list$kegg_rest <- kegg_data_rest
    message("✓ 加载KEGGREST数据库")
  } else {
    message("✗ KEGGREST数据库未找到")
  }
  
  kegg_cp_file <- file.path(base_dir, "kegg_hsa_database.RData")
  if (file.exists(kegg_cp_file)) {
    load(kegg_cp_file)
    db_list$kegg_cp <- kegg_data
    message("✓ 加载clusterProfiler KEGG数据库")
  } else {
    message("✗ clusterProfiler KEGG数据库未找到")
  }
  
  if (length(db_list) == 0) {
    warning("没有找到任何预下载的数据库！")
  }
  
  return(db_list)
}

PATHWAY_DATABASES <- load_pathway_databases()

# ============================================================================
# Major细胞类型Marker基因列表
# ============================================================================
get_major_marker_genes <- function() {
  message("加载Major细胞类型marker基因列表...")
  
  major_markers <- list(
    "T cell" = c("CD3D", "CD3E", "CD3G", "CD2", "THEMIS", "CD247"),
    "CD4 T" = c("CD4", "IL7R", "FOXP3"),
    "CD8 T" = c("CD8A", "CD8B", "GZMK"),
    "B cell" = c("CD79A", "MS4A1", "CD19", "CD79B"),
    "Plasma" = c("MZB1", "IGHG1", "IGHG4", "JCHAIN", "SDC1"),
    "NK cell" = c("NKG7", "NCR1", "NCAM1", "KLRD1", "GNLY"),
    "Myeloid" = c("ITGAM", "CD68", "LYZ", "CSF1R"),
    "Monocyte" = c("FCN1", "CD14", "S100A8", "S100A9"),
    "Macrophage" = c("CD68", "MARCO", "CD163", "MSR1"),
    "DC" = c("CLEC9A", "CD1A", "CD1C", "XCR1", "ITGAX", "FCER1A"),
    "Neutrophil" = c("CSF3R", "S100A8", "S100A9", "FCGR3B"),
    "Mast" = c("KIT", "CPA3", "TPSB2", "TPSAB1", "FCER1A"),
    "Epithelial" = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
    "Endothelial" = c("PECAM1", "VWF", "CDH5", "FLT1", "PLVAP"),
    "Fibroblast" = c("DCN", "LUM", "COL1A1", "COL3A1", "PDGFRB"),
    "Stromal" = c("DCN", "LUM", "COL1A1", "COL3A1", "VIM")
  )
  
  message("Major细胞类型marker: ", length(major_markers), " 个细胞类型")
  
  return(major_markers)
}

major_marker_genes <- get_major_marker_genes()

# ============================================================================
# 🔧 简化的加载函数 - 直接从指定目录加载
# ============================================================================
load_seurat_with_umap <- function(cancer_type) {
  message("\n", rep("=", 80))
  message("加载癌种: ", cancer_type)
  message(rep("=", 80))
  
  # 构建RDS文件路径
  rds_file <- file.path(umap_data_dir, paste0(cancer_type, "_with_umap.rds"))
  
  # 检查文件是否存在
  if (!file.exists(rds_file)) {
    stop("❌ 未找到文件: ", rds_file, 
         "\n请确认文件名格式为: ", cancer_type, "_with_umap.rds")
  }
  
  message("✓ 找到RDS文件: ", rds_file)
  message("加载中...")
  
  # 加载Seurat对象
  seurat_obj <- tryCatch({
    readRDS(rds_file)
  }, error = function(e) {
    stop("❌ 加载失败: ", e$message)
  })
  
  message("✓ Seurat对象加载成功")
  
  # 验证必要的内容
  message("\n检查对象内容...")
  
  # 1. 检查UMAP降维
  if (!"umap" %in% names(seurat_obj@reductions)) {
    stop("❌ 该对象缺少UMAP降维结果！")
  }
  message("  ✓ UMAP降维: 存在")
  
  # 2. 检查major_cell_type字段
  if (!"major_cell_type" %in% colnames(seurat_obj@meta.data)) {
    stop("❌ 该对象缺少 major_cell_type 字段！")
  }
  message("  ✓ major_cell_type字段: 存在")
  
  # 3. 检查RNA assay
  if (!"RNA" %in% names(seurat_obj@assays)) {
    warning("⚠ 该对象缺少RNA assay")
  } else {
    message("  ✓ RNA assay: 存在")
  }
  
  # 显示对象信息
  message("\n对象信息:")
  message("  细胞数量: ", ncol(seurat_obj))
  message("  基因数量: ", nrow(seurat_obj))
  
  # Major细胞类型统计
  major_types <- table(seurat_obj$major_cell_type)
  message("  Major细胞类型数: ", length(major_types))
  message("\n  Major细胞类型分布:")
  for (ct_name in names(sort(major_types, decreasing = TRUE))) {
    message("    ", ct_name, ": ", major_types[ct_name], " cells")
  }
  
  # 降维信息
  message("\n  降维结果: ", paste(names(seurat_obj@reductions), collapse = ", "))
  
  # 如果存在fine_cell_type，也显示信息
  if ("fine_cell_type" %in% colnames(seurat_obj@meta.data)) {
    fine_types <- table(seurat_obj$fine_cell_type)
    message("  Fine细胞类型数: ", length(fine_types), " (本次分析不使用)")
  }
  
  message("\n✓ 数据加载和验证完成")
  message(rep("=", 80), "\n")
  
  return(seurat_obj)
}

# ============================================================================
# 热图绘制 - 双格式保存
# ============================================================================
plot_major_marker_heatmap <- function(seurat_obj, markers_list, cancer_type) {
  message("生成Major细胞类型marker基因热图...")
  
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
        group.by = "major_cell_type",
        size = 4,
        angle = 90
      ) +
        ggtitle(paste0(cell_type, " Markers (Major)")) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text.y = element_text(size = 10)
        )
      
      filename <- file.path(
        output_dir,
        "marker_heatmaps",
        paste0(cancer_type, "_major_", cell_type, "_heatmap")
      )
      
      save_plot_both_formats(p, filename, width = 14, height = 10)
      
    }, error = function(e) {
      message("  ✗ ", cell_type, " 热图绘制失败: ", e$message)
    })
    
    gc()
  }
}

# ============================================================================
# FindMarkers分析 - Major细胞类型
# ============================================================================
perform_major_findmarkers <- function(seurat_obj, cancer_type) {
  message("对Major细胞类型进行FindMarkers分析...")
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  cell_types <- unique(seurat_obj$major_cell_type)
  cell_types <- cell_types[!is.na(cell_types) & cell_types != "Other"]
  
  message("发现 ", length(cell_types), " 个Major细胞类型")
  
  all_markers <- list()
  
  for (cell_type in cell_types) {
    message("  处理: ", cell_type)
    
    n_cells <- sum(seurat_obj$major_cell_type == cell_type)
    
    if (n_cells < 3) {
      message("    跳过 (细胞数不足): ", n_cells)
      next
    }
    
    tryCatch({
      markers <- FindMarkers(
        seurat_obj,
        ident.1 = cell_type,
        group.by = "major_cell_type",
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
      paste0(cancer_type, "_major_findmarkers_results.csv")
    )
    
    write.csv(all_markers_df, output_file, row.names = FALSE)
    message("✓ Major FindMarkers结果已保存: ", output_file)
    
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
# 🔧 完全修复的富集分析函数（使用超几何检验）- Major版本
# ============================================================================
perform_major_enrichment_analysis <- function(markers_df, cancer_type) {
  message("进行Major细胞类型GO和KEGG富集分析（使用超几何检验）...")
  
  if (is.null(markers_df) || nrow(markers_df) == 0) {
    message("  没有marker数据，跳过富集分析")
    return(NULL)
  }
  
  cell_types <- unique(markers_df$cell_type)
  enrichment_results <- list()
  
  for (cell_type in cell_types) {
    message("  分析Major细胞类型: ", cell_type)
    
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
    # 🔧 GO富集分析（使用超几何检验）
    # ========================================================================
    tryCatch({
      message("    进行GO BP分析（超几何检验）...")
      
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
    save_major_enrichment_results(enrichment_results, cancer_type)
    plot_major_enrichment_results(enrichment_results, cancer_type)
  } else {
    message("  ✗ 所有Major细胞类型均未产生富集结果")
  }
  
  return(enrichment_results)
}

# ============================================================================
# 保存Major细胞类型富集分析结果
# ============================================================================
save_major_enrichment_results <- function(enrichment_results, cancer_type) {
  message("保存Major细胞类型富集分析结果...")
  
  for (cell_type in names(enrichment_results)) {
    if (!is.null(enrichment_results[[cell_type]]$GO_BP)) {
      tryCatch({
        go_file <- file.path(
          output_dir,
          "summary_stats",
          paste0(cancer_type, "_major_", cell_type, "_GO_BP.csv")
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
          paste0(cancer_type, "_major_", cell_type, "_KEGG.csv")
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
          paste0(cancer_type, "_major_", cell_type, "_Reactome.csv")
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
# 绘制Major细胞类型富集分析结果 - 双格式保存
# ============================================================================
plot_major_enrichment_results <- function(enrichment_results, cancer_type) {
  message("绘制Major细胞类型富集分析图...")
  
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
              title = paste0(cell_type, " (Major) - GO Biological Process"),
              x = "Gene Count",
              y = "GO Term",
              fill = "-log10(FDR)"
            ) +
            theme_minimal(base_size = 14) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
              axis.text.y = element_text(size = 10)
            )
          
          filename <- file.path(
            output_dir,
            "enrichment_plots",
            paste0(cancer_type, "_major_", cell_type, "_GO_BP_barplot")
          )
          
          save_plot_both_formats(p, filename, width = 12, height = 9)
          
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
              title = paste0(cell_type, " (Major) - KEGG Pathway"),
              x = "Gene Count",
              y = "Pathway",
              fill = "-log10(FDR)"
            ) +
            theme_minimal(base_size = 14) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
              axis.text.y = element_text(size = 10)
            )
          
          filename <- file.path(
            output_dir,
            "enrichment_plots",
            paste0(cancer_type, "_major_", cell_type, "_KEGG_barplot")
          )
          
          save_plot_both_formats(p, filename, width = 12, height = 9)
          
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
              title = paste0(cell_type, " (Major) - Reactome Pathway"),
              x = "Gene Count",
              y = "Pathway",
              fill = "-log10(FDR)"
            ) +
            theme_minimal(base_size = 14) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
              axis.text.y = element_text(size = 10)
            )
          
          filename <- file.path(
            output_dir,
            "enrichment_plots",
            paste0(cancer_type, "_major_", cell_type, "_Reactome_barplot")
          )
          
          save_plot_both_formats(p, filename, width = 12, height = 9)
          
        }, error = function(e) {
          message("  ✗ Reactome图绘制失败: ", e$message)
        })
      }
    }
    
    gc()
  }
}

# ============================================================================
# Major细胞类型评估报告
# ============================================================================
generate_major_evaluation_report <- function(cancer_type, markers_df, enrichment_results) {
  message("生成Major细胞类型评估报告...")
  
  report_file <- file.path(
    output_dir,
    "summary_stats",
    paste0(cancer_type, "_major_evaluation_report.txt")
  )
  
  sink(report_file)
  
  cat("=", rep("=", 70), "\n", sep = "")
  cat("Major细胞类型注释质量评估报告\n")
  cat("癌种: ", cancer_type, "\n")
  cat("生成时间: ", as.character(Sys.time()), "\n")
  cat("=", rep("=", 70), "\n\n", sep = "")
  
  cat("1. FINDMARKERS分析结果 (Major Cell Types)\n")
  cat("-", rep("-", 70), "\n", sep = "")
  
  if (!is.null(markers_df)) {
    cell_types <- unique(markers_df$cell_type)
    cat("Major细胞类型数量: ", length(cell_types), "\n\n")
    
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
  
  cat("\n2. 富集分析结果 (Major Cell Types)\n")
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
  cat("Major细胞类型评估完成\n")
  cat("=", rep("=", 70), "\n", sep = "")
  
  sink()
  
  message("✓ Major评估报告已保存: ", report_file)
}

# ============================================================================
# 主评估函数 - Major细胞类型
# ============================================================================
evaluate_major_cancer_annotation <- function(cancer_type) {
  message("\n", rep("=", 80))
  message("开始评估癌种 (Major Cell Types): ", cancer_type)
  message(rep("=", 80))
  
  # 加载Seurat对象（包含UMAP）
  seurat_obj <- load_seurat_with_umap(cancer_type)
  
  message("\n--- Step 1: Major细胞类型分布UMAP可视化 ---")
  plot_major_celltype_umap_enhanced(seurat_obj, cancer_type)
  
  message("\n--- Step 2: Major细胞类型Marker基因UMAP可视化 ---")
  plot_major_marker_umap_enhanced(seurat_obj, major_marker_genes, cancer_type, "major")
  
  message("\n--- Step 3: Major细胞类型Marker基因热图 ---")
  plot_major_marker_heatmap(seurat_obj, major_marker_genes, cancer_type)
  
  message("\n--- Step 4: Major细胞类型FindMarkers分析 ---")
  markers_df <- perform_major_findmarkers(seurat_obj, cancer_type)
  
  message("\n--- Step 5: Major细胞类型GO和KEGG富集分析 ---")
  enrichment_results <- perform_major_enrichment_analysis(markers_df, cancer_type)
  
  message("\n--- Step 6: 生成Major细胞类型评估报告 ---")
  generate_major_evaluation_report(cancer_type, markers_df, enrichment_results)
  
  message("\n✓ 癌种 ", cancer_type, " Major细胞类型评估完成")
  
  gc()
  
  return(list(
    markers = markers_df,
    enrichment = enrichment_results
  ))
}

# ============================================================================
# 🔧 扫描可用的癌种 - 从UMAP数据目录
# ============================================================================
get_available_cancer_types <- function() {
  message("扫描可用的癌种...")
  
  # 列出所有 *_with_umap.rds 文件
  rds_files <- list.files(
    umap_data_dir,
    pattern = "_with_umap\\.rds$",
    full.names = FALSE
  )
  
  if (length(rds_files) == 0) {
    stop("❌ 在 ", umap_data_dir, " 中未找到 *_with_umap.rds 文件")
  }
  
  # 提取癌种名称
  cancer_types <- gsub("_with_umap\\.rds$", "", rds_files)
  
  message("✓ 发现 ", length(cancer_types), " 个癌种:")
  for (ct in cancer_types) {
    message("  - ", ct)
  }
  
  return(cancer_types)
}

# ============================================================================
# 批量评估所有癌种 - Major细胞类型
# ============================================================================
evaluate_all_major_cancers <- function() {
  message("\n", rep("=", 80))
  message("开始批量评估所有癌种的Major细胞类型...")
  message(rep("=", 80), "\n")
  
  # 自动扫描可用的癌种
  cancer_types <- get_available_cancer_types()
  
  results_summary <- list()
  
  for (cancer_type in cancer_types) {
    tryCatch({
      result <- evaluate_major_cancer_annotation(cancer_type)
      
      if (!is.null(result)) {
        results_summary[[cancer_type]] <- list(
          n_cell_types = if(!is.null(result$markers)) length(unique(result$markers$cell_type)) else 0,
          n_markers = if(!is.null(result$markers)) nrow(result$markers) else 0,
          n_enriched_types = if(!is.null(result$enrichment)) length(result$enrichment) else 0
        )
      }
      
    }, error = function(e) {
      warning("评估 ", cancer_type, " 时出错: ", e$message)
    })
    
    gc()
  }
  
  # 生成汇总报告
  if (length(results_summary) > 0) {
    summary_df <- do.call(rbind, lapply(names(results_summary), function(ct) {
      data.frame(
        Cancer_Type = ct,
        N_Major_CellTypes = results_summary[[ct]]$n_cell_types,
        N_Markers = results_summary[[ct]]$n_markers,
        N_Enriched = results_summary[[ct]]$n_enriched_types
      )
    }))
    
    summary_file <- file.path(output_dir, "major_evaluation_summary.csv")
    write.csv(summary_df, summary_file, row.names = FALSE)
    
    message("\n", rep("=", 80))
    message("✓ Major细胞类型总体评估汇总已保存: ", summary_file)
    message(rep("=", 80))
    print(summary_df)
  }
  
  return(results_summary)
}

# ============================================================================
# 执行评估
# ============================================================================
message("\n", rep("=", 80))
message("Major细胞类型注释质量评估系统（优化版v2）")
message("数据源: ", umap_data_dir)
message("输出目录: ", output_dir)
message(rep("=", 80))

# 评估单个癌种（示例）
# evaluate_major_cancer_annotation("LIHC")

# 评估所有可用癌种
evaluation_summary <- evaluate_all_major_cancers()

plan(sequential)

message("\n", rep("=", 80))
message("✓ 所有Major细胞类型评估任务完成!")
message("结果保存在: ", output_dir)
message(rep("=", 80), "\n")