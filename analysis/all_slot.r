#!/usr/bin/env Rscript
library(Seurat)
library(Matrix)
library(dplyr)
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
cat("\n【全槽解剖模式】癌种:", cancer_type)
cat("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

# ================= 修复版递归探针 =================
probe <- function(obj, prefix = "") {
  cls <- class(obj)[1]
  cat(prefix, "├─ 类:", cls, "\n")

  # Assay 专用槽展开
  if (inherits(obj, "Assay")) {
    for (sn in slotNames(obj)) {
      cat(prefix, "├─ @", sn, "\n")
      probe(slot(obj, sn), paste0(prefix, "│  "))
    }
    return(invisible())
  }

  # S4 通用槽展开
  if (isS4(obj)) {
    for (sn in slotNames(obj)) {
      cat(prefix, "├─ @", sn, "\n")
      probe(slot(obj, sn), paste0(prefix, "│  "))
    }
    return(invisible())
  }

  # data.frame
  if (is.data.frame(obj)) {
    cat(prefix, "├─ 维度:", nrow(obj), "×", ncol(obj), "\n")
    cat(prefix, "├─ 列:", paste(names(obj), collapse = ", "), "\n")
    for (cn in names(obj)) {
      cat(prefix, "   ├─ $", cn, "\n")
      probe(obj[[cn]], paste0(prefix, "   │  "))
    }
    return(invisible())
  }

  # Matrix / matrix
  if (is.matrix(obj) || inherits(obj, "Matrix")) {
    cat(prefix, "├─ 维度:", nrow(obj), "×", ncol(obj), "\n")
    if (inherits(obj, "sparseMatrix")) {
      nnz <- summary(obj)$nnz
      cat(prefix, "├─ 非零:", nnz,
          "| 密度:", round(nnz / prod(dim(obj)) * 100, 4), "%\n")
    }
    return(invisible())
  }

  # list / NULL
  if (is.list(obj)) {
    nms <- names(obj)
    if (is.null(nms)) nms <- paste0("[[", seq_along(obj), "]]")
    for (i in seq_along(obj)) {
      cat(prefix, "├─ $", nms[i], "\n")
      probe(obj[[i]], paste0(prefix, "│  "))
    }
    return(invisible())
  }

  # 原子向量 / 因子 / 其他
  len <- length(obj)
  cat(prefix, "├─ 长度:", len)
  if (len > 0) {
    if (is.numeric(obj)) cat(" | 范围:", range(obj, na.rm = TRUE))
    if (is.factor(obj))  cat(" | 水平:", nlevels(obj))
  }
  cat("\n")
}

# ================= 主对象 =================
main_file <- file.path(cdir, paste0(cancer_type, "_fine_annotated.rds"))
if (file.exists(main_file)) {
  seu <- readRDS(main_file)
  cat("\n【主对象】", basename(main_file), "\n")
  probe(seu)
} else {
  cat("\n⚠️ 缺失主对象:", main_file)
}

# ================= 子对象 =================
detail_files <- list.files(cdir, pattern = "_detailed\\.rds$", full.names = TRUE)
if (length(detail_files)) {
  for (f in detail_files) {
    cat("\n\n【子对象】", basename(f), "\n")
    probe(readRDS(f))
  }
} else {
  cat("\n⚠️ 未找到任何 *_detailed.rds 文件")
}

cat("\n\n==== 全槽解剖完成 ====\n")