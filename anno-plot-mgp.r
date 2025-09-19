# ###########################################################################
#  Marker Gene Validation Pipeline (ENSEMBL → Symbol + Cache + Report)
# ###########################################################################

# ------------------ 依赖检查 & 安装 ------------------
pkg <- c("Seurat", "dplyr", "Matrix", "biomaRt", "ggplot2", "stringr")
for (p in pkg) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(Seurat)
library(dplyr)
library(Matrix)
library(biomaRt)
library(ggplot2)
library(stringr)

# ------------------ 路径设置 ------------------
analysis_dir   <- "/mnt/public7/pancancercol/hefeng/cluster-annotation"
output_dir     <- file.path(analysis_dir, "biological_validation", "marker_validation")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cache_global <- file.path(output_dir, "global_ensg_symbol_map.rds")

# ------------------ 全局 ENSEMBL → Symbol 映射表 ------------------
build_global_map <- function() {
  if (file.exists(cache_global)) {
    message("🧬 读取已缓存的全局 ENSEMBL → Symbol 映射表")
    return(readRDS(cache_global))
  }
  message("🧬 首次运行：正在构建全局 ENSEMBL → Symbol 映射表...")
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
               mart = ensembl)
  map <- map[map$hgnc_symbol != "", ]
  saveRDS(map, cache_global)
  message("✅ 全局映射表已保存：", cache_global)
  map
}

global_map <- build_global_map()

# ------------------ Marker 列表 ------------------
load_marker_genes <- function() {
  list(
    Epithelial = c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "CLDN3", "CLDN4"),
    Luminal   = c("KRT8", "KRT18", "AR", "HOXB13", "FOLH1"),
    Basal     = c("KRT5", "KRT14", "EGFR", "ITGA6"),
    T_cell    = c("CD3E", "CD3D", "LCK", "IL7R"),
    CD4_T     = c("CD4", "S100A4"),
    CD8_T     = c("CD8A", "CD8B", "GZMK", "CCL5"),
    Treg      = c("FOXP3", "IL2RA", "CTLA4"),
    NK        = c("GNLY", "NKG7", "NCAM1"),
    B_cell    = c("MS4A1", "CD79A", "CD79B"),
    Plasma    = c("IGHG1", "MZB1", "XBP1"),
    Macrophage= c("CD68", "CD163", "MSR1"),
    M1_Macro  = c("NOS2", "IL1B", "TNF"),
    M2_Macro  = c("ARG1", "IL10", "TGFB1"),
    Dendritic = c("FCER1A", "CLEC9A", "XCR1"),
    Monocyte  = c("S100A8", "S100A9", "VCAN"),
    Fibroblast= c("COL1A1", "COL1A2", "PDGFRA"),
    CAF       = c("ACTA2", "TAGLN", "FAP"),
    Endothelial=c("PECAM1", "VWF", "CLDN5"),
    Proliferating=c("MKI67", "TOP2A", "PCNA")
  )
}

# ------------------ 主验证函数 ------------------
validate_cancer_markers <- function(cancer_type) {
  message("\n=== ", cancer_type, " ===")
  in_file <- file.path(analysis_dir, cancer_type, paste0(cancer_type, "_fine_annotated.rds"))
  if (!file.exists(in_file)) {
    message("❌ 文件不存在：", in_file)
    return(NULL)
  }

  seurat_obj <- readRDS(in_file)

  # ✅ ENSEMBL → Symbol 转换
  if (all(str_detect(rownames(seurat_obj), "^ENSG"))) {
    message("🧬 检测到 ENSEMBL ID，正在转换...")
    common <- intersect(rownames(seurat_obj), global_map$ensembl_gene_id)
    seurat_obj <- seurat_obj[common, ]
    rownames(seurat_obj) <- global_map$hgnc_symbol[match(common, global_map$ensembl_gene_id)]
    message("✅ 转换完成，剩余基因数：", nrow(seurat_obj))
  }

  markers <- load_marker_genes()
  out_dir <- file.path(output_dir, cancer_type)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  summary_df <- data.frame()

  for (ct in unique(seurat_obj$fine_cell_type)) {
    cells <- WhichCells(seurat_obj, idents = ct)
    if (length(cells) == 0) next

    for (set_name in names(markers)) {
      genes <- markers[[set_name]]
      avail <- intersect(genes, rownames(seurat_obj))
      if (length(avail) == 0) next

      expr <- seurat_obj@assays$RNA@data[avail, cells, drop = FALSE]
      mean_expr <- mean(Matrix::rowMeans(expr))
      other <- setdiff(WhichCells(seurat_obj), cells)
      other_expr <- seurat_obj@assays$RNA@data[avail, other, drop = FALSE]
      other_mean <- mean(Matrix::rowMeans(other_expr))

      specificity <- ifelse(other_mean == 0, Inf, mean_expr / other_mean)
      detection_rate <- mean(Matrix::rowMeans(expr) > 0)

      summary_df <- rbind(summary_df, data.frame(
        Cell_Type          = ct,
        Marker_Set         = set_name,
        Mean_Expression    = mean_expr,
        Specificity_Score  = specificity,
        Detection_Rate     = detection_rate,
        N_Markers_Detected = sum(Matrix::rowMeans(expr) > 0),
        N_Total_Markers    = length(genes),
        Available_Markers  = paste(avail, collapse = ",")
      ))
    }
  }

  if (nrow(summary_df) == 0) {
    message("⚠️ 无有效 marker 信息")
    return(NULL)
  }

  # 保存
  write.csv(summary_df, file.path(out_dir, "marker_validation_summary.csv"), row.names = FALSE)

  # 绘图
  p1 <- ggplot(summary_df, aes(x = reorder(paste(Cell_Type, Marker_Set, sep = "_"), Specificity_Score),
                                 y = Specificity_Score)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = paste(cancer_type, "- Marker Specificity"),
         x = "Cell_Type - Marker_Set", y = "Specificity Score") +
    theme_minimal()

  p2 <- ggplot(summary_df, aes(x = reorder(paste(Cell_Type, Marker_Set, sep = "_"), Detection_Rate),
                               y = Detection_Rate)) +
    geom_col(fill = "coral") +
    coord_flip() +
    labs(title = paste(cancer_type, "- Marker Detection Rate"),
         x = "Cell_Type - Marker_Set", y = "Detection Rate") +
    theme_minimal()

  ggsave(file.path(out_dir, "specificity_scores.png"), p1, width = 10, height = max(6, nrow(summary_df)/3))
  ggsave(file.path(out_dir, "detection_rates.png"),   p2, width = 10, height = max(6, nrow(summary_df)/3))

  message("✅ 完成：", cancer_type)
  summary_df
}

# ------------------ 批量运行 & 汇总报告 ------------------
run_all <- function() {
  dirs <- list.dirs(analysis_dir, full.names = FALSE, recursive = FALSE)
  dirs <- setdiff(dirs, c("biological_validation", "visualization_analysis"))
  all_res <- list()

  for (ct in dirs) {
    res <- validate_cancer_markers(ct)
    if (!is.null(res)) all_res[[ct]] <- res
  }

  # 合并
  combined <- bind_rows(all_res, .id = "Cancer_Type")
  write.csv(combined, file.path(output_dir, "all_cancer_marker_validation.csv"), row.names = FALSE)

  # 简单统计
  stats <- combined %>%
    group_by(Cancer_Type) %>%
    summarise(N_Cell_Types = n_distinct(Cell_Type),
              N_Validated  = sum(Marker_Set != "No_Match"),
              .groups = "drop")
  write.csv(stats, file.path(output_dir, "validation_statistics.csv"), row.names = FALSE)

  message("🎉 全部完成！结果见：", output_dir)
}

# ------------------ 运行 ------------------
run_all()