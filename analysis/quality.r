################################################################################
# 🎯 TME单细胞注释质量评估系统 v3.0 - 批量处理完整版
# 作者：Claude
# 日期：2025-11-17
# 功能：自动扫描、评估、可视化、异常处理
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(cluster)
  library(ComplexHeatmap)
  library(circlize)
  library(ggpubr)
  library(cowplot)
  library(scales)
  library(viridis)
  library(patchwork)
  library(GGally)
  library(ggrepel)
  library(fmsb)
})

#===============================================================================
# 📋 模块1：配置与参考数据库
#===============================================================================

## 1.1 细胞类型层级解析 ----
parse_cell_type_hierarchy <- function(cell_type_vector) {
  hierarchy <- str_split_fixed(cell_type_vector, "_", 2)
  
  data.frame(
    full_type = cell_type_vector,
    major_type = hierarchy[, 1],
    sub_type = ifelse(hierarchy[, 2] == "", hierarchy[, 1], hierarchy[, 2]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(sub_type = str_replace_all(sub_type, "_", " "))
}

## 1.2 生物学期望Marker库 ----
create_keyword_mapping <- function() {
  list(
    major = list(
      "B" = list(
        positive = c("CD19", "CD20", "MS4A1", "CD79A", "CD79B", "BLK", "PAX5"),
        negative = c("CD3D", "CD3E", "CD14", "CD68", "EPCAM", "PECAM1", "COL1A1")
      ),
      "T" = list(
        positive = c("CD3D", "CD3E", "CD3G", "CD2", "CD7", "CD5"),
        negative = c("CD19", "CD20", "CD14", "EPCAM", "PECAM1", "COL1A1")
      ),
      "NK" = list(
        positive = c("NCAM1", "NKG7", "GNLY", "KLRD1", "KLRF1", "KLRB1"),
        negative = c("CD3D", "CD19", "CD14", "EPCAM", "PECAM1")
      ),
      "Myeloid" = list(
        positive = c("CD14", "CD68", "LYZ", "CSF1R", "FCGR3A", "CD163"),
        negative = c("CD3D", "CD19", "EPCAM", "PECAM1", "COL1A1")
      ),
      "Endothelial" = list(
        positive = c("PECAM1", "VWF", "CDH5", "CLDN5", "CD34"),
        negative = c("CD3D", "CD19", "EPCAM", "COL1A1", "CD68")
      ),
      "Epithelial" = list(
        positive = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
        negative = c("CD3D", "CD19", "PECAM1", "COL1A1", "CD68")
      ),
      "Fibroblasts" = list(
        positive = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "PDGFRB"),
        negative = c("CD3D", "EPCAM", "PECAM1", "CD68", "CD19")
      )
    )
  )
}

## 1.3 识别不合理注释 ----
identify_problematic_annotations <- function(cell_types_df) {
  
  invalid_rules <- list(
    list(major = "B", forbidden = c(
      "Classical monocytes", "Intermediate monocytes", "Non classical monocytes",
      "Classical Monocytes", "Intermediate Monocytes", "Non classical Monocytes",
      "Myeloid dendritic cells", "Plasmacytoid dendritic cells",
      "Natural killer cells", "NK cells", "Th1 cells", "Th2 cells", "Th17 cells",
      "Progenitor cells", "CD4 T cells", "CD8 T cells"
    ), reason = "B细胞不应包含髓系/NK/T细胞亚型"),
    
    list(major = "Endothelial", forbidden = c(
      "Classical monocytes", "Intermediate monocytes", "Non classical monocytes",
      "Myeloid dendritic cells", "Plasmacytoid dendritic cells",
      "Natural killer cells", "Th1 cells", "Progenitor cells"
    ), reason = "内皮细胞不应包含免疫细胞亚型"),
    
    list(major = "Epithelial", forbidden = c(
      "Classical monocytes", "Intermediate monocytes", "Non classical monocytes",
      "Myeloid dendritic cells", "Plasmacytoid dendritic cells",
      "Natural killer cells", "Th1 cells", "Plasmablasts",
      "Low density neutrophils", "Progenitor cells"
    ), reason = "上皮细胞不应包含免疫细胞亚型"),
    
    list(major = "Fibroblasts", forbidden = c(
      "Classical monocytes", "Intermediate monocytes", "Non classical monocytes",
      "Myeloid dendritic cells", "Plasmacytoid dendritic cells",
      "Exhausted B cells", "Th1 cells", "Progenitor cells"
    ), reason = "成纤维细胞不应包含免疫细胞亚型"),
    
    list(major = "NK", forbidden = c(
      "Exhausted B cells", "Plasmablasts", "Memory B cells"
    ), reason = "NK细胞不应包含B细胞亚型"),
    
    list(major = "T", forbidden = c(
      "Exhausted B cells", "Switched memory B cells", "Memory B cells",
      "Plasmablasts", "Plasma cells"
    ), reason = "T细胞不应包含B细胞亚型"),
    
    list(major = "Myeloid", forbidden = c(
      "Exhausted B cells", "Memory B cells", "Naive B cells"
    ), reason = "髓系细胞不应包含B细胞亚型")
  )
  
  problematic_list <- list()
  
  for (rule in invalid_rules) {
    matches <- cell_types_df %>%
      filter(major_type == rule$major, sub_type %in% rule$forbidden) %>%
      mutate(issue_reason = rule$reason)
    
    if (nrow(matches) > 0) {
      problematic_list[[rule$major]] <- matches
    }
  }
  
  return(bind_rows(problematic_list))
}

#===============================================================================
# 📊 模块2：核心量化指标计算
#===============================================================================

## 2.1 Marker特异性评分（MSS）----
calculate_marker_specificity <- function(seurat_obj, cell_type_col = "fine_cell_type") {
  
  cat("  🔬 计算Marker特异性（MSS）...\n")
  
  hierarchy <- parse_cell_type_hierarchy(seurat_obj@meta.data[[cell_type_col]])
  seurat_obj$major_type <- hierarchy$major_type
  seurat_obj$sub_type <- hierarchy$sub_type
  
  keyword_map <- create_keyword_mapping()
  expr_matrix <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")
  
  results <- list()
  
  for (major in unique(hierarchy$major_type)) {
    if (!major %in% names(keyword_map$major)) next
    
    cells_in <- which(hierarchy$major_type == major)
    cells_out <- which(hierarchy$major_type != major)
    
    if (length(cells_in) < 10 || length(cells_out) < 10) next
    
    markers <- keyword_map$major[[major]]
    pos_genes <- intersect(markers$positive, rownames(expr_matrix))
    neg_genes <- intersect(markers$negative, rownames(expr_matrix))
    
    if (length(pos_genes) == 0) next
    
    pos_in <- rowMeans(expr_matrix[pos_genes, cells_in, drop = FALSE])
    pos_out <- rowMeans(expr_matrix[pos_genes, cells_out, drop = FALSE])
    neg_in <- if(length(neg_genes) > 0) rowMeans(expr_matrix[neg_genes, cells_in, drop = FALSE]) else 0
    neg_out <- if(length(neg_genes) > 0) rowMeans(expr_matrix[neg_genes, cells_out, drop = FALSE]) else 0
    
    pos_enrichment <- mean((pos_in - pos_out) / (pos_in + pos_out + 1e-6))
    neg_depletion <- mean((neg_out - neg_in) / (neg_out + neg_in + 1e-6))
    mss <- (pos_enrichment + neg_depletion) / 2
    
    pos_pval <- tryCatch(t.test(pos_in, pos_out)$p.value, error = function(e) 1)
    
    results[[major]] <- list(
      MSS = mss,
      pos_enrichment = pos_enrichment,
      neg_depletion = neg_depletion,
      pos_in_mean = mean(pos_in),
      pos_out_mean = mean(pos_out),
      pvalue = pos_pval,
      n_cells = length(cells_in),
      n_pos_genes = length(pos_genes),
      n_neg_genes = length(neg_genes)
    )
  }
  
  mss_df <- map_dfr(names(results), ~data.frame(
    cell_type = .x,
    MSS = results[[.x]]$MSS,
    pos_enrichment = results[[.x]]$pos_enrichment,
    neg_depletion = results[[.x]]$neg_depletion,
    pos_in_mean = results[[.x]]$pos_in_mean,
    pos_out_mean = results[[.x]]$pos_out_mean,
    pvalue = results[[.x]]$pvalue,
    n_cells = results[[.x]]$n_cells,
    n_markers = results[[.x]]$n_pos_genes + results[[.x]]$n_neg_genes,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("     ✓ 完成 %d 个细胞类型\n", nrow(mss_df)))
  return(mss_df)
}

## 2.2 功能一致性评分（FCS）----
calculate_functional_consistency <- function(seurat_obj, cell_type_col = "fine_cell_type") {
  
  cat("  🧬 计算功能一致性（FCS）...\n")
  
  functional_modules <- list(
    "Cytotoxicity" = c("PRF1", "GZMA", "GZMB", "GZMH", "GNLY", "NKG7"),
    "Proliferation" = c("MKI67", "TOP2A", "PCNA", "STMN1", "CCNB1", "CDK1"),
    "Stress" = c("HSPA1A", "HSPA1B", "DNAJB1", "HSP90AA1"),
    "Hypoxia" = c("HIF1A", "VEGFA", "SLC2A1", "LDHA", "PDK1"),
    "Interferon" = c("ISG15", "IFIT1", "IFIT3", "MX1", "IFI44L", "STAT1"),
    "Chemokine" = c("CCL3", "CCL4", "CCL5", "CXCL8", "CXCL10"),
    "Activation" = c("CD69", "CD44", "ICOS", "CD40LG", "IL2RA"),
    "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1"),
    "Migration" = c("CCR7", "SELL", "S1PR1", "ITGAL"),
    "Metabolism" = c("PKM", "ENO1", "GAPDH", "PFKP", "LDHA")
  )
  
  expr_matrix <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")
  hierarchy <- parse_cell_type_hierarchy(seurat_obj@meta.data[[cell_type_col]])
  
  fcs_results <- list()
  
  for (major in unique(hierarchy$major_type)) {
    cells_in <- which(hierarchy$major_type == major)
    if (length(cells_in) < 10) next
    
    module_scores <- sapply(functional_modules, function(genes) {
      genes_present <- intersect(genes, rownames(expr_matrix))
      if (length(genes_present) == 0) return(0)
      mean(expr_matrix[genes_present, cells_in, drop = FALSE])
    })
    
    top_modules <- names(sort(module_scores, decreasing = TRUE)[1:3])
    top_scores <- module_scores[top_modules]
    
    cv <- sd(top_scores) / (mean(top_scores) + 1e-6)
    fcs <- 1 / (cv + 0.1)
    
    fcs_results[[major]] <- list(
      FCS = fcs,
      CV = cv,
      dominant_modules = top_modules,
      module_scores = module_scores
    )
  }
  
  fcs_df <- map_dfr(names(fcs_results), ~data.frame(
    cell_type = .x,
    FCS = fcs_results[[.x]]$FCS,
    CV = fcs_results[[.x]]$CV,
    top_function_1 = fcs_results[[.x]]$dominant_modules[1],
    top_function_2 = fcs_results[[.x]]$dominant_modules[2],
    top_function_3 = fcs_results[[.x]]$dominant_modules[3],
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("     ✓ 完成 %d 个细胞类型\n", nrow(fcs_df)))
  return(fcs_df)
}

## 2.3 细胞纯度（Silhouette）----
calculate_cell_purity <- function(seurat_obj, cell_type_col = "fine_cell_type",
                                   reduction = "umap", dims = 1:2) {
  
  cat("  🎯 计算细胞纯度（Silhouette）...\n")
  
  if (!reduction %in% names(seurat_obj@reductions)) {
    cat("     ⚠️  未找到UMAP降维，尝试PCA...\n")
    reduction <- "pca"
    dims <- 1:30
    if (!reduction %in% names(seurat_obj@reductions)) {
      cat("     ⚠️  未找到降维，跳过\n")
      return(NULL)
    }
  }
  
  embeddings <- Embeddings(seurat_obj, reduction = reduction)[, dims]
  hierarchy <- parse_cell_type_hierarchy(seurat_obj@meta.data[[cell_type_col]])
  
  # Major类型
  major_factor <- factor(hierarchy$major_type)
  if (length(levels(major_factor)) > 1) {
    sil_major <- silhouette(as.numeric(major_factor), dist(embeddings))
    major_purity <- data.frame(
      cell_type = levels(major_factor),
      silhouette = tapply(sil_major[, 3], hierarchy$major_type, mean),
      level = "major",
      stringsAsFactors = FALSE
    )
  } else {
    major_purity <- data.frame()
  }
  
  # Subtype
  subtype_counts <- table(hierarchy$sub_type)
  valid_subtypes <- names(subtype_counts[subtype_counts > 10])
  cells_valid <- which(hierarchy$sub_type %in% valid_subtypes)
  
  if (length(cells_valid) > 100 && length(valid_subtypes) > 1) {
    sil_sub <- silhouette(
      as.numeric(factor(hierarchy$sub_type[cells_valid])),
      dist(embeddings[cells_valid, ])
    )
    sub_purity <- data.frame(
      cell_type = valid_subtypes,
      silhouette = tapply(sil_sub[, 3], hierarchy$sub_type[cells_valid], mean),
      level = "subtype",
      stringsAsFactors = FALSE
    )
  } else {
    sub_purity <- data.frame()
  }
  
  purity_df <- bind_rows(major_purity, sub_purity)
  cat(sprintf("     ✓ 完成 %d 个细胞类型\n", nrow(purity_df)))
  return(purity_df)
}

## 2.4 注释置信度 ----
calculate_annotation_confidence <- function(seurat_obj, cell_type_col = "fine_cell_type",
                                             expr_threshold = 0.5) {
  
  cat("  🔍 计算注释置信度...\n")
  
  keyword_map <- create_keyword_mapping()
  expr_matrix <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")
  hierarchy <- parse_cell_type_hierarchy(seurat_obj@meta.data[[cell_type_col]])
  
  conf_results <- list()
  
  for (major in unique(hierarchy$major_type)) {
    if (!major %in% names(keyword_map$major)) next
    
    cells_in <- which(hierarchy$major_type == major)
    markers <- keyword_map$major[[major]]
    pos_genes <- intersect(markers$positive, rownames(expr_matrix))
    
    if (length(pos_genes) == 0) next
    
    marker_coverage <- sapply(pos_genes, function(g) {
      mean(expr_matrix[g, cells_in] > expr_threshold)
    })
    
    confidence <- mean(marker_coverage)
    
    conf_results[[major]] <- list(
      confidence = confidence,
      marker_coverage = marker_coverage,
      n_markers = length(pos_genes),
      n_cells = length(cells_in)
    )
  }
  
  conf_df <- map_dfr(names(conf_results), ~data.frame(
    cell_type = .x,
    confidence = conf_results[[.x]]$confidence,
    n_markers_used = conf_results[[.x]]$n_markers,
    n_cells = conf_results[[.x]]$n_cells,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("     ✓ 完成 %d 个细胞类型\n", nrow(conf_df)))
  return(conf_df)
}

## 2.5 综合质量评分 ----
calculate_composite_score <- function(mss_df, fcs_df, purity_df, conf_df) {
  
  cat("  📊 计算综合评分...\n")
  
  composite <- mss_df %>%
    left_join(fcs_df, by = "cell_type") %>%
    left_join(
      purity_df %>% filter(level == "major") %>% select(cell_type, silhouette),
      by = "cell_type"
    ) %>%
    left_join(conf_df, by = "cell_type")
  
  # Z-score标准化
  composite <- composite %>%
    mutate(
      MSS_norm = (MSS - mean(MSS, na.rm = TRUE)) / (sd(MSS, na.rm = TRUE) + 1e-6),
      FCS_norm = (FCS - mean(FCS, na.rm = TRUE)) / (sd(FCS, na.rm = TRUE) + 1e-6),
      Purity_norm = (silhouette - mean(silhouette, na.rm = TRUE)) / (sd(silhouette, na.rm = TRUE) + 1e-6),
      Confidence_norm = (confidence - mean(confidence, na.rm = TRUE)) / (sd(confidence, na.rm = TRUE) + 1e-6)
    ) %>%
    mutate(
      MSS_norm = (MSS_norm - min(MSS_norm, na.rm = TRUE)) / 
                 (max(MSS_norm, na.rm = TRUE) - min(MSS_norm, na.rm = TRUE) + 1e-6),
      FCS_norm = (FCS_norm - min(FCS_norm, na.rm = TRUE)) / 
                 (max(FCS_norm, na.rm = TRUE) - min(FCS_norm, na.rm = TRUE) + 1e-6),
      Purity_norm = (Purity_norm - min(Purity_norm, na.rm = TRUE)) / 
                    (max(Purity_norm, na.rm = TRUE) - min(Purity_norm, na.rm = TRUE) + 1e-6),
      Confidence_norm = (Confidence_norm - min(Confidence_norm, na.rm = TRUE)) / 
                        (max(Confidence_norm, na.rm = TRUE) - min(Confidence_norm, na.rm = TRUE) + 1e-6)
    )
  
  # 加权评分
  weights <- c(MSS = 0.30, FCS = 0.25, Purity = 0.25, Confidence = 0.20)
  
  composite <- composite %>%
    mutate(
      composite_score = MSS_norm * weights["MSS"] +
                        FCS_norm * weights["FCS"] +
                        Purity_norm * weights["Purity"] +
                        Confidence_norm * weights["Confidence"],
      quality_level = case_when(
        composite_score >= 0.75 ~ "Excellent",
        composite_score >= 0.60 ~ "Good",
        composite_score >= 0.40 ~ "Fair",
        TRUE ~ "Poor"
      ),
      quality_level = factor(quality_level, 
                            levels = c("Excellent", "Good", "Fair", "Poor"))
    )
  
  cat("     ✓ 综合评分完成\n")
  return(composite)
}

#===============================================================================
# 📈 模块3：论文级可视化
#===============================================================================

## 3.1 综合质量总览图 ----
plot_quality_overview <- function(composite_df, output_dir, cancer_type) {
  
  cat("  📊 生成质量总览图...\n")
  
  # A: 质量分布饼图
  p1 <- composite_df %>%
    count(quality_level) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ggplot(aes(x = "", y = pct, fill = quality_level)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 0.5) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = sprintf("%s\n%.1f%%", quality_level, pct)),
              position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("Excellent" = "#2E7D32", "Good" = "#66BB6A",
                                   "Fair" = "#FFA726", "Poor" = "#D32F2F")) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
    labs(title = "A. Quality Distribution")
  
  # B: 综合评分排序
  p2 <- composite_df %>%
    arrange(composite_score) %>%
    mutate(cell_type = factor(cell_type, levels = cell_type)) %>%
    ggplot(aes(x = cell_type, y = composite_score, fill = quality_level)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Excellent" = "#2E7D32", "Good" = "#66BB6A",
                                   "Fair" = "#FFA726", "Poor" = "#D32F2F")) +
    geom_hline(yintercept = c(0.4, 0.6, 0.75), linetype = "dashed", 
               color = "grey40", alpha = 0.5) +
    labs(title = "B. Composite Quality Score",
         x = NULL, y = "Score", fill = "Quality") +
    theme_classic(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 7),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  # C: 指标相关性
  cor_data <- composite_df %>%
    select(MSS, FCS, silhouette, confidence) %>%
    cor(use = "complete.obs")
  
  p3 <- ggplot(reshape2::melt(cor_data), aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "C. Metric Correlation", x = NULL, y = NULL, fill = "Correlation") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  # D: 细胞数量分布
  p4 <- composite_df %>%
    arrange(desc(n_cells)) %>%
    mutate(cell_type = factor(cell_type, levels = cell_type)) %>%
    head(15) %>%
    ggplot(aes(x = cell_type, y = n_cells, fill = quality_level)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("Excellent" = "#2E7D32", "Good" = "#66BB6A",
                                   "Fair" = "#FFA726", "Poor" = "#D32F2F")) +
    scale_y_continuous(labels = comma) +
    labs(title = "D. Top 15 Cell Types by Count",
         x = NULL, y = "Cell Count", fill = "Quality") +
    theme_classic(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 7),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
  
  # 组合
  layout <- "
  AABBBB
  CCDDDD
  "
  
  combined <- p1 + p2 + p3 + p4 + 
    plot_layout(design = layout) +
    plot_annotation(
      title = sprintf("%s - Quality Assessment Overview", cancer_type),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  ggsave(file.path(output_dir, sprintf("%s_Fig1_Quality_Overview.png", cancer_type)),
         combined, width = 14, height = 10, dpi = 300, bg = "white")
  
  cat("     ✓ 已保存\n")
}

## 3.2 综合质量热图 ----
plot_quality_heatmap <- function(composite_df, output_dir, cancer_type) {
  
  cat("  📊 生成综合质量热图...\n")
  
  mat <- composite_df %>%
    select(cell_type, MSS_norm, FCS_norm, Purity_norm, Confidence_norm, composite_score) %>%
    column_to_rownames("cell_type") %>%
    as.matrix()
  
  colnames(mat) <- c("MSS", "FCS", "Purity", "Confidence", "Composite")
  
  # 行注释
  row_ha <- rowAnnotation(
    Quality = composite_df$quality_level,
    `Cell Count` = anno_barplot(
      log10(composite_df$n_cells + 1),
      gp = gpar(fill = "steelblue"),
      width = unit(1.5, "cm")
    ),
    col = list(
      Quality = c("Excellent" = "#2E7D32", "Good" = "#66BB6A",
                  "Fair" = "#FFA726", "Poor" = "#D32F2F")
    ),
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    show_legend = c(TRUE, FALSE)
  )
  
  # 颜色
  col_fun <- colorRamp2(
    seq(0, 1, length.out = 100),
    colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
  )
  
  png(file.path(output_dir, sprintf("%s_Fig2_Quality_Heatmap.png", cancer_type)),
      width = 10, height = max(8, nrow(mat) * 0.3), units = "in", res = 300)
  
  ht <- Heatmap(
    mat,
    name = "Normalized\nScore",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    right_annotation = row_ha,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(4, "cm")
    ),
    column_title = sprintf("%s - Quality Assessment Metrics", cancer_type),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_title = "Cell Types",
    row_title_gp = gpar(fontsize = 11, fontface = "bold"),
    border = TRUE
  )
  
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  cat("     ✓ 已保存\n")
}

## 3.3 MSS vs FCS关系图 ----
plot_mss_fcs_relationship <- function(mss_df, fcs_df, composite_df, output_dir, cancer_type) {
  
  cat("  📊 生成MSS-FCS关系图...\n")
  
  merged <- mss_df %>%
    left_join(fcs_df, by = "cell_type") %>%
    left_join(composite_df %>% select(cell_type, quality_level, composite_score), 
              by = "cell_type")
  
  p <- ggplot(merged, aes(x = MSS, y = FCS)) +
    geom_vline(xintercept = median(merged$MSS, na.rm = TRUE), 
               linetype = "dashed", color = "grey50", alpha = 0.7) +
    geom_hline(yintercept = median(merged$FCS, na.rm = TRUE), 
               linetype = "dashed", color = "grey50", alpha = 0.7) +
    geom_point(aes(size = n_cells, color = quality_level), alpha = 0.7) +
    geom_text_repel(
      aes(label = cell_type),
      size = 2.5,
      max.overlaps = 12,
      box.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.2
    ) +
    scale_color_manual(
      values = c("Excellent" = "#2E7D32", "Good" = "#66BB6A",
                 "Fair" = "#FFA726", "Poor" = "#D32F2F")
    ) +
    scale_size_continuous(range = c(2, 8), labels = comma) +
    labs(
      title = sprintf("%s - Marker Specificity vs Functional Consistency", cancer_type),
      subtitle = "Dashed lines indicate median values",
      x = "MSS (Marker Specificity Score)",
      y = "FCS (Functional Consistency Score)",
      size = "Cell Count",
      color = "Quality"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey40", size = 9),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90", size = 0.3)
    ) +
    guides(
      size = guide_legend(order = 1),
      color = guide_legend(order = 2, override.aes = list(size = 4))
    )
  
  ggsave(file.path(output_dir, sprintf("%s_Fig3_MSS_FCS_Relationship.png", cancer_type)),
         p, width = 11, height = 8, dpi = 300, bg = "white")
  
  cat("     ✓ 已保存\n")
}

## 3.4 问题注释可视化 ----
plot_problematic_annotations <- function(problematic_df, composite_df, output_dir, cancer_type) {
  
  cat("  📊 生成问题注释图...\n")
  
  if (nrow(problematic_df) == 0) {
    cat("     ℹ️  无问题注释，跳过\n")
    return(NULL)
  }
  
  # 合并评分信息
  problem_with_score <- problematic_df %>%
    left_join(composite_df %>% select(cell_type = cell_type, composite_score, quality_level),
              by = c("major_type" = "cell_type"))
  
  p <- problem_with_score %>%
    mutate(full_type = fct_reorder(full_type, composite_score)) %>%
    ggplot(aes(x = full_type, y = composite_score, fill = quality_level)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_hline(yintercept = 0.4, linetype = "dashed", color = "red", size = 1) +
    coord_flip() +
    scale_fill_manual(values = c("Excellent" = "#2E7D32", "Good" = "#66BB6A",
                                   "Fair" = "#FFA726", "Poor" = "#D32F2F")) +
    labs(
      title = sprintf("%s - Problematic Annotations", cancer_type),
      subtitle = sprintf("Total: %d biologically inconsistent annotations", nrow(problem_with_score)),
      x = NULL,
      y = "Composite Quality Score",
      fill = "Quality"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey40", size = 9),
      axis.text.y = element_text(size = 7),
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_dir, sprintf("%s_Fig4_Problematic_Annotations.png", cancer_type)),
         p, width = 10, height = max(6, nrow(problem_with_score) * 0.2), 
         dpi = 300, bg = "white")
  
  cat("     ✓ 已保存\n")
}

#===============================================================================
# 📝 模块4：报告生成
#===============================================================================

generate_comprehensive_report <- function(composite_df, mss_df, fcs_df, 
                                          purity_df, conf_df, problematic_df, 
                                          output_dir, cancer_type) {
  
  cat("  📝 生成综合报告...\n")
  
  report_file <- file.path(output_dir, sprintf("%s_Quality_Report.txt", cancer_type))
  
  sink(report_file)
  
  cat(strrep("=", 80), "\n")
  cat(sprintf("      %s - TME单细胞注释质量综合评估报告\n", cancer_type))
  cat("      生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat(strrep("=", 80), "\n\n")
  
  # 总览
  cat("【1】总体统计\n")
  cat(strrep("-", 80), "\n")
  cat(sprintf("总细胞类型数: %d\n", nrow(composite_df)))
  cat(sprintf("  - Excellent (≥0.75): %d (%.1f%%)\n",
              sum(composite_df$quality_level == "Excellent", na.rm = TRUE),
              mean(composite_df$quality_level == "Excellent", na.rm = TRUE) * 100))
  cat(sprintf("  - Good (≥0.60):      %d (%.1f%%)\n",
              sum(composite_df$quality_level == "Good", na.rm = TRUE),
              mean(composite_df$quality_level == "Good", na.rm = TRUE) * 100))
  cat(sprintf("  - Fair (≥0.40):      %d (%.1f%%)\n",
              sum(composite_df$quality_level == "Fair", na.rm = TRUE),
              mean(composite_df$quality_level == "Fair", na.rm = TRUE) * 100))
  cat(sprintf("  - Poor (<0.40):      %d (%.1f%%)\n",
              sum(composite_df$quality_level == "Poor", na.rm = TRUE),
              mean(composite_df$quality_level == "Poor", na.rm = TRUE) * 100))
  cat("\n")
  
  # 指标总结
  cat("【2】指标总结\n")
  cat(strrep("-", 80), "\n")
  cat(sprintf("MSS: 均值=%.3f, 中位数=%.3f, 范围=[%.3f, %.3f]\n",
              mean(composite_df$MSS, na.rm = TRUE),
              median(composite_df$MSS, na.rm = TRUE),
              min(composite_df$MSS, na.rm = TRUE),
              max(composite_df$MSS, na.rm = TRUE)))
  cat(sprintf("FCS: 均值=%.3f, 中位数=%.3f, 范围=[%.3f, %.3f]\n",
              mean(composite_df$FCS, na.rm = TRUE),
              median(composite_df$FCS, na.rm = TRUE),
              min(composite_df$FCS, na.rm = TRUE),
              max(composite_df$FCS, na.rm = TRUE)))
  cat(sprintf("Purity: 均值=%.3f, 中位数=%.3f, 范围=[%.3f, %.3f]\n",
              mean(composite_df$silhouette, na.rm = TRUE),
              median(composite_df$silhouette, na.rm = TRUE),
              min(composite_df$silhouette, na.rm = TRUE),
              max(composite_df$silhouette, na.rm = TRUE)))
  cat(sprintf("Confidence: 均值=%.3f, 中位数=%.3f, 范围=[%.3f, %.3f]\n",
              mean(composite_df$confidence, na.rm = TRUE),
              median(composite_df$confidence, na.rm = TRUE),
              min(composite_df$confidence, na.rm = TRUE),
              max(composite_df$confidence, na.rm = TRUE)))
  cat("\n")
  
  # Top 10
  cat("【3】Top 10 高质量细胞类型\n")
  cat(strrep("-", 80), "\n")
  top10 <- composite_df %>% arrange(desc(composite_score)) %>% head(10)
  for (i in 1:min(10, nrow(top10))) {
    cat(sprintf("%2d. %-35s Score: %.3f (%s)\n",
                i, top10$cell_type[i], top10$composite_score[i], top10$quality_level[i]))
  }
  cat("\n")
  
  # Bottom 10
  cat("【4】Bottom 10 低质量细胞类型\n")
  cat(strrep("-", 80), "\n")
  bottom10 <- composite_df %>% arrange(composite_score) %>% head(10)
  for (i in 1:min(10, nrow(bottom10))) {
    cat(sprintf("%2d. %-35s Score: %.3f (%s)\n",
                i, bottom10$cell_type[i], bottom10$composite_score[i], bottom10$quality_level[i]))
  }
  cat("\n")
  
  # 问题注释
  cat("【5】生物学不合理的注释\n")
  cat(strrep("-", 80), "\n")
  if (nrow(problematic_df) > 0) {
    cat(sprintf("检测到 %d 个不合理的Major-Subtype组合:\n\n", nrow(problematic_df)))
    
    by_major <- split(problematic_df, problematic_df$major_type)
    for (major in names(by_major)) {
      cat(sprintf("\n[%s] - %d 个问题:\n", major, nrow(by_major[[major]])))
      for (ft in by_major[[major]]$full_type) {
        cat(sprintf("    ❌ %s\n", ft))
      }
    }
  } else {
    cat("✅ 未检测到明显不合理的注释组合\n")
  }
  cat("\n")
  
  cat(strrep("=", 80), "\n")
  cat("报告结束\n")
  cat(strrep("=", 80), "\n")
  
  sink()
  
  cat("     ✓ 报告已保存\n")
}

#===============================================================================
# 🚀 模块5：单个样本评估主流程
#===============================================================================

run_single_sample_evaluation <- function(
  rds_file,
  cell_type_col = "fine_cell_type",
  output_base_dir = "./TME_Quality_Batch_Results",
  reduction = "umap",
  dims = 1:2
) {
  
  # 提取癌种名称
  cancer_type <- str_extract(basename(rds_file), "^[A-Z]+")
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat(sprintf("🎯 处理: %s\n", cancer_type))
  cat(strrep("=", 80), "\n\n")
  
  # 创建输出目录
  output_dir <- file.path(output_base_dir, cancer_type)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 加载数据
  cat("📦 [1/9] 加载Seurat对象...\n")
  seurat_obj <- tryCatch({
    readRDS(rds_file)
  }, error = function(e) {
    cat(sprintf("   ❌ 加载失败: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(seurat_obj)) {
    return(NULL)
  }
  
  # 检查必需列
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    cat(sprintf("   ❌ 未找到列: %s\n", cell_type_col))
    return(NULL)
  }
  
  cat(sprintf("   ✓ 细胞数: %d\n\n", ncol(seurat_obj)))
  
  # 解析层级
  cat("📋 [2/9] 解析细胞类型层级...\n")
  hierarchy <- parse_cell_type_hierarchy(seurat_obj@meta.data[[cell_type_col]])
  cat(sprintf("   ✓ Major: %d | Subtype: %d\n\n",
              length(unique(hierarchy$major_type)),
              length(unique(hierarchy$sub_type))))
  
  # 识别问题
  cat("📋 [3/9] 识别问题注释...\n")
  problematic_df <- identify_problematic_annotations(hierarchy)
  cat(sprintf("   ⚠️  检测到 %d 个问题组合\n\n", nrow(problematic_df)))
  
  # MSS
  cat("📋 [4/9] 计算MSS...\n")
  mss_df <- calculate_marker_specificity(seurat_obj, cell_type_col)
  cat("\n")
  
  # FCS
  cat("📋 [5/9] 计算FCS...\n")
  fcs_df <- calculate_functional_consistency(seurat_obj, cell_type_col)
  cat("\n")
  
  # Purity
  cat("📋 [6/9] 计算Purity...\n")
  purity_df <- calculate_cell_purity(seurat_obj, cell_type_col, reduction, dims)
  cat("\n")
  
  # Confidence
  cat("📋 [7/9] 计算Confidence...\n")
  conf_df <- calculate_annotation_confidence(seurat_obj, cell_type_col)
  cat("\n")
  
  # 综合评分
  cat("📋 [8/9] 计算综合评分...\n")
  composite_df <- calculate_composite_score(mss_df, fcs_df, purity_df, conf_df)
  cat("\n")
  
  # 保存数据
  cat("📋 [9/9] 保存结果...\n")
  write.csv(composite_df, file.path(output_dir, sprintf("%s_composite_scores.csv", cancer_type)), row.names = FALSE)
  write.csv(mss_df, file.path(output_dir, sprintf("%s_MSS_scores.csv", cancer_type)), row.names = FALSE)
  write.csv(fcs_df, file.path(output_dir, sprintf("%s_FCS_scores.csv", cancer_type)), row.names = FALSE)
  if (!is.null(purity_df)) {
    write.csv(purity_df, file.path(output_dir, sprintf("%s_purity_scores.csv", cancer_type)), row.names = FALSE)
  }
  write.csv(conf_df, file.path(output_dir, sprintf("%s_confidence_scores.csv", cancer_type)), row.names = FALSE)
  write.csv(problematic_df, file.path(output_dir, sprintf("%s_problematic_annotations.csv", cancer_type)), row.names = FALSE)
  cat("   ✓ CSV文件已保存\n\n")
  
  # 可视化
  cat("📊 生成可视化...\n")
  plot_quality_overview(composite_df, output_dir, cancer_type)
  plot_quality_heatmap(composite_df, output_dir, cancer_type)
  plot_mss_fcs_relationship(mss_df, fcs_df, composite_df, output_dir, cancer_type)
  plot_problematic_annotations(problematic_df, composite_df, output_dir, cancer_type)
  cat("\n")
  
  # 生成报告
  generate_comprehensive_report(composite_df, mss_df, fcs_df, purity_df, 
                                conf_df, problematic_df, output_dir, cancer_type)
  
  cat(sprintf("✅ %s 评估完成！\n", cancer_type))
  cat(strrep("=", 80), "\n\n")
  
  return(list(
    cancer_type = cancer_type,
    composite = composite_df,
    mss = mss_df,
    fcs = fcs_df,
    purity = purity_df,
    confidence = conf_df,
    problematic = problematic_df
  ))
}

#===============================================================================
# 🌐 模块6：批量处理所有样本
#===============================================================================

run_batch_evaluation <- function(
  input_dir = "/mnt/public7/pancancercol/hefeng/cluster-annotation/annotation_evaluation_results_11_9/processed_rds",
  pattern = "_with_umap\\.rds$",
  cell_type_col = "fine_cell_type",
  output_base_dir = "./TME_Quality_Batch_Results",
  reduction = "umap",
  dims = 1:2,
  max_parallel = 1  # 建议串行处理避免内存问题
) {
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("🌐 批量质量评估启动\n")
  cat(strrep("=", 80), "\n\n")
  
  # 扫描所有RDS文件
  rds_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  cat(sprintf("📁 输入目录: %s\n", input_dir))
  cat(sprintf("📊 检测到 %d 个RDS文件\n\n", length(rds_files)))
  
  if (length(rds_files) == 0) {
    cat("❌ 未找到匹配的RDS文件！\n")
    return(NULL)
  }
  
  # 显示文件列表
  cat("待处理文件:\n")
  for (i in seq_along(rds_files)) {
    cancer <- str_extract(basename(rds_files[i]), "^[A-Z]+")
    cat(sprintf("  [%2d] %s\n", i, cancer))
  }
  cat("\n")
  
  # 创建总输出目录
  if (!dir.exists(output_base_dir)) {
    dir.create(output_base_dir, recursive = TRUE)
  }
  
  # 批量处理
  all_results <- list()
  success_count <- 0
  failed_files <- c()
  
  for (i in seq_along(rds_files)) {
    cat(sprintf("\n进度: [%d/%d]\n", i, length(rds_files)))
    
    result <- tryCatch({
      run_single_sample_evaluation(
        rds_file = rds_files[i],
        cell_type_col = cell_type_col,
        output_base_dir = output_base_dir,
        reduction = reduction,
        dims = dims
      )
    }, error = function(e) {
      cat(sprintf("❌ 处理失败: %s\n", e$message))
      return(NULL)
    })
    
    if (!is.null(result)) {
      all_results[[result$cancer_type]] <- result
      success_count <- success_count + 1
    } else {
      failed_files <- c(failed_files, basename(rds_files[i]))
    }
    
    # 释放内存
    gc()
  }
  
  # 生成汇总报告
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("📊 生成批量汇总报告...\n")
  cat(strrep("=", 80), "\n\n")
  
  generate_batch_summary(all_results, output_base_dir, failed_files)
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat(sprintf("✅ 批量评估完成！成功: %d/%d\n", success_count, length(rds_files)))
  cat(strrep("=", 80), "\n\n")
  
  return(all_results)
}

#===============================================================================
# 📊 模块7：批量汇总报告
#===============================================================================

generate_batch_summary <- function(all_results, output_base_dir, failed_files) {
  
  # 汇总所有composite scores
  summary_df <- map_dfr(names(all_results), ~{
    all_results[[.x]]$composite %>%
      mutate(cancer_type = .x) %>%
      select(cancer_type, cell_type, composite_score, quality_level, 
             MSS, FCS, silhouette, confidence, n_cells)
  })
  
  # 保存汇总CSV
  write.csv(summary_df, 
            file.path(output_base_dir, "Batch_Summary_All_Scores.csv"), 
            row.names = FALSE)
  
  # 按癌种汇总
  cancer_summary <- summary_df %>%
    group_by(cancer_type) %>%
    summarise(
      n_cell_types = n(),
      avg_composite_score = mean(composite_score, na.rm = TRUE),
      n_excellent = sum(quality_level == "Excellent", na.rm = TRUE),
      n_good = sum(quality_level == "Good", na.rm = TRUE),
      n_fair = sum(quality_level == "Fair", na.rm = TRUE),
      n_poor = sum(quality_level == "Poor", na.rm = TRUE),
      pct_excellent = n_excellent / n_cell_types * 100,
      avg_MSS = mean(MSS, na.rm = TRUE),
      avg_FCS = mean(FCS, na.rm = TRUE),
      avg_purity = mean(silhouette, na.rm = TRUE),
      avg_confidence = mean(confidence, na.rm = TRUE),
      total_cells = sum(n_cells, na.rm = TRUE)
    ) %>%
    arrange(desc(avg_composite_score))
  
  write.csv(cancer_summary, 
            file.path(output_base_dir, "Batch_Summary_By_Cancer.csv"), 
            row.names = FALSE)
  
  # 生成汇总可视化
  
  # 1. 癌种质量对比
  p1 <- ggplot(cancer_summary, aes(x = reorder(cancer_type, avg_composite_score), 
                                    y = avg_composite_score)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    geom_hline(yintercept = c(0.4, 0.6, 0.75), linetype = "dashed", 
               color = "red", alpha = 0.5) +
    coord_flip() +
    labs(title = "Average Composite Quality Score by Cancer Type",
         x = NULL, y = "Average Score") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_base_dir, "Batch_Summary_Cancer_Quality.png"),
         p1, width = 10, height = max(6, nrow(cancer_summary) * 0.4), 
         dpi = 300, bg = "white")
  
  # 2. 质量分布堆叠柱状图
  quality_dist <- summary_df %>%
    count(cancer_type, quality_level) %>%
    group_by(cancer_type) %>%
    mutate(pct = n / sum(n) * 100)
  
  p2 <- ggplot(quality_dist, aes(x = cancer_type, y = pct, fill = quality_level)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Excellent" = "#2E7D32", "Good" = "#66BB6A",
                                   "Fair" = "#FFA726", "Poor" = "#D32F2F")) +
    coord_flip() +
    labs(title = "Quality Distribution by Cancer Type",
         x = NULL, y = "Percentage (%)", fill = "Quality") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  
  ggsave(file.path(output_base_dir, "Batch_Summary_Quality_Distribution.png"),
         p2, width = 10, height = max(6, nrow(cancer_summary) * 0.4), 
         dpi = 300, bg = "white")
  
  # 3. 四指标对比热图
  metric_mat <- cancer_summary %>%
    select(cancer_type, avg_MSS, avg_FCS, avg_purity, avg_confidence) %>%
    column_to_rownames("cancer_type") %>%
    as.matrix()
  
  colnames(metric_mat) <- c("MSS", "FCS", "Purity", "Confidence")
  
  # 归一化
  metric_mat_norm <- apply(metric_mat, 2, function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  })
  
  col_fun <- colorRamp2(
    seq(0, 1, length.out = 100),
    colorRampPalette(c("#F7F7F7", "#2166AC"))(100)
  )
  
  png(file.path(output_base_dir, "Batch_Summary_Metric_Heatmap.png"),
      width = 8, height = max(6, nrow(metric_mat) * 0.4), units = "in", res = 300)
  
  ht <- Heatmap(
    metric_mat_norm,
    name = "Normalized\nScore",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_title = "Cross-Cancer Quality Metrics Comparison",
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_title = "Cancer Types",
    row_title_gp = gpar(fontsize = 11, fontface = "bold"),
    border = TRUE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )
  
  draw(ht)
  dev.off()
  
  # 文本报告
  report_file <- file.path(output_base_dir, "Batch_Summary_Report.txt")
  
  sink(report_file)
  
  cat(strrep("=", 80), "\n")
  cat("      批量质量评估汇总报告\n")
  cat("      生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat(strrep("=", 80), "\n\n")
  
  cat("【1】总体统计\n")
  cat(strrep("-", 80), "\n")
  cat(sprintf("成功处理癌种数: %d\n", length(all_results)))
  if (length(failed_files) > 0) {
    cat(sprintf("失败文件数: %d\n", length(failed_files)))
    cat("失败文件:\n")
    for (f in failed_files) {
      cat(sprintf("  - %s\n", f))
    }
  }
  cat(sprintf("总细胞类型数: %d\n", nrow(summary_df)))
  cat(sprintf("总细胞数: %s\n", format(sum(summary_df$n_cells, na.rm = TRUE), big.mark = ",")))
  cat("\n")
  
  cat("【2】各癌种质量排名（按平均评分）\n")
  cat(strrep("-", 80), "\n")
  for (i in 1:nrow(cancer_summary)) {
    cat(sprintf("%2d. %-10s | Score: %.3f | Excellent: %2d (%.1f%%) | Good: %2d | Fair: %2d | Poor: %2d\n",
                i,
                cancer_summary$cancer_type[i],
                cancer_summary$avg_composite_score[i],
                cancer_summary$n_excellent[i],
                cancer_summary$pct_excellent[i],
                cancer_summary$n_good[i],
                cancer_summary$n_fair[i],
                cancer_summary$n_poor[i]))
  }
  cat("\n")
  
  cat("【3】跨癌种指标统计\n")
  cat(strrep("-", 80), "\n")
  cat(sprintf("MSS:        均值=%.3f, 中位数=%.3f\n",
              mean(cancer_summary$avg_MSS, na.rm = TRUE),
              median(cancer_summary$avg_MSS, na.rm = TRUE)))
  cat(sprintf("FCS:        均值=%.3f, 中位数=%.3f\n",
              mean(cancer_summary$avg_FCS, na.rm = TRUE),
              median(cancer_summary$avg_FCS, na.rm = TRUE)))
  cat(sprintf("Purity:     均值=%.3f, 中位数=%.3f\n",
              mean(cancer_summary$avg_purity, na.rm = TRUE),
              median(cancer_summary$avg_purity, na.rm = TRUE)))
  cat(sprintf("Confidence: 均值=%.3f, 中位数=%.3f\n",
              mean(cancer_summary$avg_confidence, na.rm = TRUE),
              median(cancer_summary$avg_confidence, na.rm = TRUE)))
  cat("\n")
  
  cat(strrep("=", 80), "\n")
  cat("报告结束\n")
  cat(strrep("=", 80), "\n")
  
  sink()
  
  cat("✅ 汇总报告已生成\n")
  cat(sprintf("   - %s\n", report_file))
  cat(sprintf("   - %s\n", file.path(output_base_dir, "Batch_Summary_All_Scores.csv")))
  cat(sprintf("   - %s\n", file.path(output_base_dir, "Batch_Summary_By_Cancer.csv")))
}

#===============================================================================
# 🎬 执行批量评估
#===============================================================================

## 设置参数
INPUT_DIR <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/annotation_evaluation_results_11_9/processed_rds"
OUTPUT_DIR <- "./TME_Quality_Batch_Results"
CELL_TYPE_COL <- "fine_cell_type"
REDUCTION <- "umap"
DIMS <- 1:2

## 运行批量评估
results <- run_batch_evaluation(
  input_dir = INPUT_DIR,
  pattern = "_with_umap\\.rds$",
  cell_type_col = CELL_TYPE_COL,
  output_base_dir = OUTPUT_DIR,
  reduction = REDUCTION,
  dims = DIMS,
  max_parallel = 1
)

## 查看结果
cat("\n批量评估完成！结果保存在:", OUTPUT_DIR, "\n")