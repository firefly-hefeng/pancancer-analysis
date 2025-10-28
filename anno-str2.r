# 安装和加载Seurat包
if (!requireNamespace("Seurat", quietly = TRUE)) {
    install.packages("Seurat")
}
library(Seurat)

# 设置工作目录到包含RDS文件的文件夹
setwd("/mnt/public7/pancancercol/seurat/GSE151530/rds")

# 加载RDS文件
files <- list.files(pattern = "\\.rds$")

# 遍历每个文件并查看其结构
for (file in files) {
  cat("Loading", file, "\n")
  data <- readRDS(file)
  
  # 检查是否为Seurat对象
  if (inherits(data, "Seurat")) {
    cat("Seurat object:", file, "\n")
    
    # 输出所有assay的名称
    assays <- names(assays(data))
    cat("Assays in", file, ":", paste(assays, collapse = ", "), "\n")
    
    # 遍历每个assay并查看其组成
    for (assay_name in assays) {
      cat("\nAssay:", assay_name, "\n")
      assay_data <- getAssay(data, assay_name)
      
      # 输出assay的组成
      cat("Slots in assay:", assay_name, ":", paste(names(assay_data), collapse = ", "), "\n")
      
      # 可以选择输出每个slot的内容，例如表达矩阵
      if ("counts" %in% names(assay_data)) {
        cat("Expression matrix (counts):\n")
        print(head(assay_data[["counts"]]))
      }
      if ("data" %in% names(assay_data)) {
        cat("Expression matrix (normalized data):\n")
        print(head(assay_data[["data"]]))
      }
      if ("scale.data" %in% names(assay_data)) {
        cat("Expression matrix (scaled data):\n")
        print(head(assay_data[["scale.data"]]))
      }
    }
  } else {
    cat("Not a Seurat object:", file, "\n")
  }
  
  cat("\n")
}