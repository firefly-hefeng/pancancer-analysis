library(Matrix)
library(SummarizedExperiment) 
library(BiocParallel)
register(MulticoreParam(workers = 16)) 
library(SingleR)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(stringr)

.libPaths('/bioapps/Rlibs/4.0.0')
library(celldex)
library(infercnv)
library(org.Hs.eg.db)

# 参数设置
gene_order_file <- "inferCNV_ref/gencode_v47_gene_pos-Copy1.infercnv.txt"
output_base_dir <- "cnv_output/"
base_data_dir <- '~/public7/seurat/'  # 样本数据的根目录
dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

# 读取样本信息表格
sample_info <- read.csv("inferCNV_ref/sample_info.csv", stringsAsFactors = FALSE)

# 准备SingleR参考数据
ref <- celldex::HumanPrimaryCellAtlasData()
gene_symbols <- rownames(ref)

ensembl_ids <- mapIds(
    org.Hs.eg.db,
    keys = gene_symbols,
    keytype = "SYMBOL",
    column = "ENSEMBL",
    multiVals = "first"
)

keep_genes <- !is.na(ensembl_ids)
ensembl_ids <- ensembl_ids[keep_genes]
ref <- ref[keep_genes, ]
rownames(ref) <- ensembl_ids

# 获取所有样本文件夹
sample_dirs <- list.dirs(path = base_data_dir, full.names = TRUE, recursive = FALSE)

# 主循环处理每个样本
for (i in 1:nrow(sample_info)) {
  sample_row <- sample_info[i, ]
  sample_id <- sample_row$ID
  file_prefix <- sample_row$file_name
  cancer_type <- sample_row$Cancer_general
  
  message("\n==== 处理样本: ", sample_id, " | 文件前缀: ", file_prefix, " | 癌症类型: ", cancer_type, " ====")
  
  # 创建癌症类型特定的输出目录
  cancer_output_dir <- file.path(output_base_dir, cancer_type)
  dir.create(cancer_output_dir, showWarnings = FALSE, recursive = TRUE)
    
    
  seurat_dir <- file.path(cancer_output_dir, "seurat_objects")
  annotated_rds_path <- file.path(seurat_dir, paste0(sample_id, "_", file_prefix, "_annotated.rds"))
  if (file.exists(annotated_rds_path)) {
    message(" 样本 ", sample_id, " 的annotated.rds已存在，跳过处理")
    next
  }
    
  
  # 查找匹配的样本文件夹
  matched_dir <- NULL
  for (dir in sample_dirs) {
    dir_name <- basename(dir)
    if (grepl(sample_id, dir_name, ignore.case = TRUE)) {
      matched_dir <- dir
      break
    }
  }
  
  if (is.null(matched_dir)) {
    warning("未找到匹配的文件夹: ", sample_id)
    next
  }
  
  # 构建Seurat对象文件路径 - 修改后的部分
  rds_file <- file.path(matched_dir, "rds", paste0(file_prefix, "_seurat.rds"))
  
  # 检查文件是否存在
  if (!file.exists(rds_file)) {
    warning("Seurat对象文件不存在: ", rds_file)
    next
  }
  
  # 读取Seurat对象
  seurat_obj <- readRDS(rds_file)
  
  # 添加样本级元数据
  seurat_obj$Cancer_general <- cancer_type
  seurat_obj$Sample_ID <- sample_id
  seurat_obj$TumorNormal <- sample_row$TumorNormal
  seurat_obj$Metastasis <- sample_row$Metastasis
  seurat_obj$File_prefix <- file_prefix
  
  
  # 确保数据已归一化
  if (is.null(LayerData(seurat_obj, assay = "RNA", layer = "data"))) {
    seurat_obj <- NormalizeData(seurat_obj)
  }
  
  # SingleR注释
  counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  
  singler_results <- SingleR(
    test = counts_mat,
    ref = ref,
    labels = ref$label.main,
    BPPARAM = bpparam()  
  )
  
  # 保存SingleR结果到metadata
  seurat_obj$singleR_labels <- singler_results$labels
  
  # 提取上皮和髓系细胞
  epi_cells <- colnames(seurat_obj)[
    grepl("Epithelial_cells|Keratinocytes|Hepatocytes|Neuroepithelial_cell", 
          seurat_obj$singleR_labels, ignore.case = TRUE)
  ]
  
  myel_cells <- colnames(seurat_obj)[
    grepl("Monocyte|Neutrophils|Macrophage|DC|GMP|Myelocyte|Pro-Myelocyte", 
          seurat_obj$singleR_labels, ignore.case = TRUE)
  ]
  
  if (length(epi_cells) < 10 || length(myel_cells) < 10) {
    warning("样本 ", sample_id, " 的Epithelial (", length(epi_cells), 
            ") 或 Myeloid (", length(myel_cells), ") 细胞过少，跳过")
    next
  }
  
  # 准备inferCNV输入
  tumor_counts <- counts_mat[, epi_cells, drop = FALSE]
  normal_counts <- counts_mat[, myel_cells, drop = FALSE]
  
  common_genes <- intersect(rownames(tumor_counts), rownames(normal_counts))
  combined_counts <- cbind(
    tumor_counts[common_genes, ], 
    normal_counts[common_genes, ]
  )
  
  # 创建细胞注释
  cell_anno <- data.frame(
    cell_id = colnames(combined_counts),
    group = c(
      rep("Epithelial", ncol(tumor_counts)),
      rep("Myeloid", ncol(normal_counts))
    ),
    stringsAsFactors = FALSE
  )
  
  # 运行inferCNV
  cnv_dir <- file.path(cancer_output_dir, paste0(sample_id, "_", file_prefix, "_inferCNV"))
  dir.create(cnv_dir, showWarnings = FALSE, recursive = TRUE)
  
  anno_file <- file.path(cnv_dir, "cell_annotations.txt")
  write.table(cell_anno, anno_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  infercnv_obj <- tryCatch(
    {
    CreateInfercnvObject(
      raw_counts_matrix = combined_counts,
      annotations_file = anno_file,
      gene_order_file = gene_order_file,
      ref_group_names = "Myeloid",
      chr_exclude = c("chrM", "chrY")
     )
    }, 
    error = function(e) {
    message("⚠️ 样本 ", sample_id, " CreateInfercnvObject 出错: ", conditionMessage(e))
    return(NULL)
   })

  if (is.null(infercnv_obj)){
    next
    }
  
    
  infercnv_obj <- tryCatch(
    {
      infercnv::run(
        infercnv_obj,
        cutoff = 0.1,
        out_dir = cnv_dir,
        cluster_by_groups = TRUE,
        plot_steps = FALSE,
        denoise = TRUE,
        HMM = TRUE,
        num_threads = 10,
        output_format = 'pdf'
      )
    },
    error = function(e) {
      message("⚠️ 样本 ", sample_id, " 运行 inferCNV 出错: ", conditionMessage(e))
      return(NULL) # 返回 NULL 表示跳过
    }
   )

  # 可选：如果 NULL 就直接跳过后续处理
  if (is.null(infercnv_obj)) {
    next
   }
    
  
  # 保存结果
  saveRDS(infercnv_obj, file.path(cnv_dir, "inferCNV_obj.rds"))
  
  # 保存完整的Seurat对象（包含所有注释）
  seurat_dir <- file.path(cancer_output_dir, "seurat_objects")
  dir.create(seurat_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(seurat_obj, file.path(seurat_dir, paste0(sample_id, "_", file_prefix, "_annotated.rds")))
  
  message("成功处理样本: ", sample_id, " | 结果保存至: ", cnv_dir)
}

message("\n==== 所有样本处理完成 ====")