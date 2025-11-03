#!/usr/bin/env Rscript
library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)
library(optparse)

option_list <- list(
  make_option(c("--root"), type = "character",
              default = "/mnt/public7/pancancercol/hefeng/cluster-annotation/fine_annotation",
              help = "根目录"),
  make_option(c("--cancer"), type = "character",
              default = "PRAD",
              help = "癌种")
)
opt <- parse_args(OptionParser(option_list = option_list))

root_dir    <- opt$root
cancer_type <- opt$cancer
cdir        <- file.path(root_dir, cancer_type)

if (!dir.exists(cdir)) stop("目录不存在: ", cdir)

cat("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
cat("\n开始解剖 —— 癌种:", cancer_type)
cat("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

# ---------- 通用解剖 ----------
inspect_seurat <- function(obj, title) {
  cat("\n\n--------------------------------------------------")
  cat("\n", title)
  cat("\n--------------------------------------------------")
  cat("\n类        :", class(obj)[1])
  cat("\n细胞×基因 :", ncol(obj), "×", nrow(obj))
  cat("\nAssays    :", paste(names(obj@assays), collapse = ", "))
  cat("\nReductions:", paste(names(obj@reductions), collapse = ", "))

  # assays
  for (a in names(obj@assays)) {
    mat <- GetAssayData(obj, assay = a, layer = "counts")
    cat("\n  Assay", a, "- counts 类:", class(mat), "维度:", dim(mat)[1], "×", dim(mat)[2])
    if (inherits(mat, "sparseMatrix"))
      cat(" 密度:", round(summary(mat)$nnz / prod(dim(mat)) * 100, 3), "%")
    cat("\n    前 5×5 真实值:")
    print(as.matrix(mat[1:5, 1:5]))
  }

  # meta.data
  meta <- obj@meta.data
  cat("\n  Meta.data 列类型 & NA 比例:")
  print(sapply(meta, class) %>% enframe(name = "col", value = "type") %>%
        mutate(NA_pct = round(colMeans(is.na(meta)) * 100, 2)))
  cat("\n  前 5 行:")
  print(meta[1:5, ])

  if ("fine_cell_type" %in% colnames(meta)) {
    cat("\n  fine_cell_type 前 10 频数:")
    print(table(meta$fine_cell_type) %>% sort(decreasing = TRUE) %>% head(10))
  }
}

inspect_csv <- function(csv, title) {
  cat("\n\n--------------------------------------------------")
  cat("\n", title)
  cat("\n--------------------------------------------------")
  if (!file.exists(csv)) {
    cat("文件不存在:", basename(csv))
    return()
  }
  df <- read.csv(csv, stringsAsFactors = FALSE)
  cat("维度:", dim(df), "\n列类型:")
  print(sapply(df, class) %>% enframe(name = "col", value = "type"))
  cat("前 6 行×前 5 列:")
  print(df[1:6, 1:min(5, ncol(df))])
}

# ===================== 1. 主对象 =====================
main_file <- file.path(cdir, paste0(cancer_type, "_fine_annotated.rds"))
if (file.exists(main_file)) {
  inspect_seurat(readRDS(main_file), paste0(cancer_type, "_fine_annotated.rds"))
} else {
  cat("\n⚠️ 缺失主对象:", main_file)
}

# ===================== 2. CSV 表 =====================
inspect_csv(file.path(cdir, paste0(cancer_type, "_cell_type_stats.csv")),
            paste0(cancer_type, "_cell_type_stats"))
inspect_csv(file.path(cdir, paste0(cancer_type, "_summary_stats.csv")),
            paste0(cancer_type, "_summary_stats"))

# ===================== 3. major 子对象 =====================
detail_files <- list.files(cdir, pattern = "_detailed\\.rds$", full.names = TRUE)
if (length(detail_files)) {
  for (f in detail_files) {
    inspect_seurat(readRDS(f), basename(f))
  }
} else {
  cat("\n⚠️ 未找到任何 *_detailed.rds 文件")
}

cat("\n\n==== 全部解剖完成 ====\n")