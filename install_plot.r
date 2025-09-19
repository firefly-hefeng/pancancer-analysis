## 强制第一条库路径
.libPaths("/mnt/public5/pancancercol/miniconda3/envs/hfa-plot/lib/R/library")
## ---------- 0. 准备 ----------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

## ---------- 1. 包列表 ----------
pkgs <- c(
  "Seurat", "ggplot2", "dplyr", "tidyverse", "ComplexHeatmap",
  "clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot", "msigdbr",
  "GSVA", "pheatmap", "corrplot", "VennDiagram", "UpSetR",
  "SingleR", "celldex", "scRNAseq", "scater", "gridExtra",
  "RColorBrewer", "viridis", "scales", "patchwork"
)

## 区分来源
cran_pkgs <- c("ggplot2", "dplyr", "tidyverse", "pheatmap", "corrplot",
               "VennDiagram", "UpSetR", "RColorBrewer", "viridis",
               "scales", "patchwork", "gridExtra")
bioc_pkgs <- setdiff(pkgs, cran_pkgs)

## 记录失败
fail_log <- character()

## ---------- 2. 安装函数 ----------
safe_install <- function(pkg, source = c("CRAN", "BIOC")) {
  source <- match.arg(source)
  msg <- tryCatch({
    if (source == "CRAN") {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    } else {
      BiocManager::install(pkg, ask = FALSE)
    }
    NULL
  }, error = function(e) {
    return(paste0("[", source, "] ", pkg, " -> ", conditionMessage(e)))
  })
  if (!is.null(msg)) fail_log <<- c(fail_log, msg)
}

## ---------- 3. 逐包安装 ----------
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) safe_install(p, "CRAN")
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) safe_install(p, "BIOC")
}

## ---------- 4. 加载并记录加载失败 ----------
load_fail <- character()
for (p in pkgs) {
  tryCatch(
    suppressPackageStartupMessages(library(p, character.only = TRUE)),
    error = function(e) load_fail <<- c(load_fail,
                                        paste0(p, " (load) -> ", conditionMessage(e)))
  )
}

## ---------- 5. 结果报告 ----------
if (length(fail_log)) {
  message("\n--- 安装失败的包 ---")
  writeLines(fail_log)
} else {
  message("\n--- 所有包均已安装成功！---")
}

if (length(load_fail)) {
  message("\n--- 加载失败的包 ---")
  writeLines(load_fail)
} else {
  message("\n--- 所有包均已加载成功！---")
}