#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
  library(dplyr)
})

# ============ 用户可配置区域 ============
# 整合方法: "harmony", "liger", "scvi", 或 "none"（若已在code2阶段做了强力批次校正，且只需合并和聚类）
integration_method <- "harmony"

# 归一化方法: "SCT" 或 "LogNormalize"
normalize_method <- "SCT"

# 是否进行子群细化分析
do_subclustering <- TRUE

# 物种与参考库（SingleR自动注释）
# 人类参考：HumanPrimaryCellAtlasData 或 MonacoImmuneData / BlueprintEncodeData 等
use_singler <- TRUE
singler_reference <- "HumanPrimaryCellAtlasData"  # 可选: "MonacoImmuneData","BlueprintEncodeData","DatabaseImmuneCellExpressionData"

# Azimuth（可选，互联网访问必要）
use_azimuth <- FALSE
azimuth_ref <- "Human_PBMC" # 示例参考；根据Azimuth官网选择合适癌种或组织参考

# UMAP维度、邻域大小、分辨率等超参
npcs <- 50
umap_neighbors <- 30
resolution_main <- 0.5

# 输出目录
dir.create("plots", showWarnings = FALSE)
dir.create("rds", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

# ============ 函数区 ============

read_all_seurat <- function(path = "rds", pattern = "_seurat\\.rds$"){
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if(length(files) == 0){
    stop("未在rds/下找到 *_seurat.rds 文件，请先运行code2.R")
  }
  message(sprintf("检测到 %d 个Seurat对象：", length(files)))
  print(basename(files))
  objs <- lapply(files, readRDS)
  # 确保orig.ident存在
  for(i in seq_along(objs)){
    if(is.null(objs[[i]]$orig.ident)){
      objs[[i]]$orig.ident <- objs[[i]]@project.name
    }
    # 加入sample字段
    if(is.null(objs[[i]]$sample)){
      objs[[i]]$sample <- objs[[i]]$orig.ident
    }
  }
  names(objs) <- gsub("_seurat\\.rds$", "", basename(files))
  objs
}

basic_qc_plots <- function(obj, prefix="merged"){
  p1 <- VlnPlot(obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"), group.by = "orig.ident",
                pt.size = 0.1, ncol=3) + NoLegend()
  ggsave(sprintf("plots/%s_qc_violin_by_sample.pdf", prefix), p1, width=12, height=4)
  
  p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    ggtitle("nCount vs nFeature")
  ggsave(sprintf("plots/%s_qc_scatter.pdf", prefix), p2, width=5, height=4)
}

perform_normalization <- function(obj, method="SCT"){
  if(method == "SCT"){
    suppressPackageStartupMessages(library(sctransform))
    # 建议按样本分块SCT，减少批次偏差
    obj <- SCTransform(obj, vars.to.regress = "percent.mt", method = "glmGamPoi", verbose = FALSE)
  }else{
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
    obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  }
  obj
}

run_pca <- function(obj, npcs=50){
  if("SCT" %in% names(obj@assays)){
    obj <- RunPCA(obj, npcs = npcs, assay = "SCT", verbose = FALSE)
  } else {
    obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  }
  ElbowPlot(obj, ndims = npcs)
  ggsave("plots/elbowplot.pdf", width=6, height=4)
  obj
}

integrate_harmony <- function(obj, npcs=50){
  suppressPackageStartupMessages(library(harmony))
  # 以样本作为批次
  obj <- RunHarmony(obj, group.by.vars = "orig.ident", reduction = "pca", dims.use = 1:npcs, verbose = TRUE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:npcs, n.neighbors = umap_neighbors, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:npcs)
  obj <- FindClusters(obj, resolution = resolution_main)
  obj
}

integrate_liger <- function(obj, npcs=50){
  # LIGER工作流（SeuratWrappers）
  suppressPackageStartupMessages(library(SeuratWrappers))
  suppressPackageStartupMessages(library(rliger))
  # 将对象按样本拆分再整合
  obj.list <- SplitObject(obj, split.by = "orig.ident")
  if("SCT" %in% names(obj@assays)){
    obj.list <- lapply(obj.list, function(x) {
      DefaultAssay(x) <- "SCT"
      x
    })
  }
  for(i in seq_along(obj.list)){
    obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE)
    obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], verbose = FALSE)
  }
  obj.liger <- RunOptimizeALS(obj.list, k = 30, lambda = 5)
  obj.liger <- RunQuantileNorm(obj.liger)
  obj <- as.Seurat(obj.liger)
  obj <- RunUMAP(obj, reduction = "iNMF", dims = 1:30, n.neighbors = umap_neighbors, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "iNMF", dims = 1:30)
  obj <- FindClusters(obj, resolution = resolution_main)
  obj
}

integrate_scvi <- function(obj, npcs=50){
  # 提示：最佳实践是在python里用scvi-tools完成整合。
  # 这里提供reticulate桥接的占位流程，若本机未配置，请改为在python完成后导回整合坐标。
  suppressPackageStartupMessages(library(reticulate))
  use_python(Sys.getenv("PYTHON_BIN", unset = "python"), required = FALSE)
  sc <- import("scanpy", convert = FALSE)
  scvi <- import("scvi", convert = FALSE)
  
  # 导出到Anndata
  tmp_h5ad <- "rds/tmp_for_scvi.h5ad"
  SeuratDisk::SaveH5Seurat(obj, filename = "rds/tmp_for_scvi.h5seurat", overwrite = TRUE)
  SeuratDisk::Convert("rds/tmp_for_scvi.h5seurat", dest = "h5ad", overwrite = TRUE)
  
  # 在python侧：构建AnnData -> scVI整合 -> 获取latent
  # 这里只是调用，需确保python环境可用
  py_run_string(sprintf("
import scanpy as sc
import scvi
adata = sc.read_h5ad('%s')
if 'orig.ident' not in adata.obs:
    adata.obs['orig.ident'] = 'sample'
scvi.model.SCVI.setup_anndata(adata, batch_key='orig.ident')
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
model.train(max_epochs=100, plan_kwargs={'lr':1e-3})
adata.obsm['X_scVI'] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep='X_scVI', n_neighbors=%d)
sc.tl.umap(adata)
adata.write_h5ad('%s')
", tmp_h5ad, umap_neighbors, "rds/tmp_scvi_out.h5ad"))
  
  # 读回
  ad <- SeuratDisk::Connect("rds/tmp_scvi_out.h5ad")
  # 将坐标嵌入回Seurat（示意）
  # 实际可通过zellkonverter/SeuratDisk读取obsm，并放入 obj@reductions$new
  # 简化：从h5ad提取UMAP坐标表，再赋值到obj@reductions[["umap"]]
  # 这里保留为示意，生产建议在python内完成后另存csv坐标再merge回R。
  message("注意：scVI桥接示意完成。生产中请在python内完成整合并返回坐标。")
  
  # 后续邻近图/聚类
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:npcs)
  obj <- FindClusters(obj, resolution = resolution_main)
  obj <- RunUMAP(obj, dims = 1:npcs, n.neighbors = umap_neighbors)
  obj
}

cluster_and_embed <- function(obj, npcs=50, method="harmony"){
  # PCA
  obj <- run_pca(obj, npcs = npcs)
  # 整合
  if(method == "harmony"){
    obj <- integrate_harmony(obj, npcs = npcs)
  } else if(method == "liger"){
    obj <- integrate_liger(obj, npcs = npcs)
  } else if(method == "scvi"){
    obj <- integrate_scvi(obj, npcs = npcs)
  } else {
    # 不整合，直接用PCA
    obj <- RunUMAP(obj, dims = 1:npcs, n.neighbors = umap_neighbors, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:npcs)
    obj <- FindClusters(obj, resolution = resolution_main)
  }
  obj
}

auto_annotation_singler <- function(obj){
  suppressPackageStartupMessages(library(SingleR))
  suppressPackageStartupMessages(library(celldex))
  # 选参考
  ref <- switch(singler_reference,
                "HumanPrimaryCellAtlasData" = celldex::HumanPrimaryCellAtlasData(),
                "MonacoImmuneData" = celldex::MonacoImmuneData(),
                "BlueprintEncodeData" = celldex::BlueprintEncodeData(),
                "DatabaseImmuneCellExpressionData" = celldex::DatabaseImmuneCellExpressionData(),
                celldex::HumanPrimaryCellAtlasData())
  DefaultAssay(obj) <- if("SCT" %in% names(obj@assays)) "SCT" else "RNA"
  mat <- GetAssayData(obj, slot = "data")
  pred <- SingleR(test = as.matrix(mat), ref = ref, labels = ref$label.main)
  obj$SingleR_label <- pred$labels
  write.csv(data.frame(cell=colnames(obj), label=obj$SingleR_label),
            file = "tables/SingleR_labels.csv", row.names = FALSE)
  obj
}

azimuth_annotation <- function(obj){
  # 需联网与正确参考库名称，示意
  suppressPackageStartupMessages(library(Azimuth))
  # 请根据Azimuth选择正确reference，如：reference = "pbmc_multimodal"
  obj <- Azimuth::RunAzimuth(obj, reference = azimuth_ref)
  obj$Azimuth_label <- obj$predicted.celltype.l2
  obj
}

plot_embeddings <- function(obj, prefix="merged"){
  p1 <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("Clusters")
  p2 <- DimPlot(obj, reduction = "umap", group.by = "orig.ident") + ggtitle("Samples")
  ggsave(sprintf("plots/%s_umap_clusters.pdf", prefix), p1, width=6, height=5)
  ggsave(sprintf("plots/%s_umap_by_sample.pdf", prefix), p2, width=6, height=5)
  if("SingleR_label" %in% colnames(obj@meta.data)){
    p3 <- DimPlot(obj, reduction = "umap", group.by = "SingleR_label", label = TRUE, repel = TRUE) + ggtitle("SingleR")
    ggsave(sprintf("plots/%s_umap_singler.pdf", prefix), p3, width=7, height=5)
  }
  if("Azimuth_label" %in% colnames(obj@meta.data)){
    p4 <- DimPlot(obj, reduction = "umap", group.by = "Azimuth_label", label = TRUE, repel = TRUE) + ggtitle("Azimuth")
    ggsave(sprintf("plots/%s_umap_azimuth.pdf", prefix), p4, width=7, height=5)
  }
}

export_results <- function(obj, prefix="merged"){
  saveRDS(obj, sprintf("rds/%s_cluster_annotated_seurat.rds", prefix))
  write.csv(obj@meta.data, sprintf("tables/%s_metadata.csv", prefix))
  # 写出平均表达与marker
  DefaultAssay(obj) <- if("SCT" %in% names(obj@assays)) "SCT" else "RNA"
  Idents(obj) <- obj$seurat_clusters
  markers <- FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1, return.thresh = 0.05)
  fwrite(markers, sprintf("tables/%s_markers_all_clusters.csv", prefix))
}

# 常用marker列表（用于人工校对，可根据癌种扩展）
canonical_markers <- list(
  Epithelial_Tumor = c("EPCAM","KRT8","KRT18","KRT19","KRT17"),
  T_CD8 = c("CD3D","CD3E","CD2","TRAC","CD8A","CD8B","GZMB","PRF1"),
  T_CD4 = c("CD3D","IL7R","CXCR4","CCR7"),
  T_Exhausted = c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","ENTPD1"),
  B_cell = c("MS4A1","CD79A","CD79B","MZB1","IGHG1"),
  Plasma = c("MZB1","SDC1","XBP1","IGHG1"),
  Myeloid_Macro = c("LYZ","CD68","MSR1","MRC1"),
  Myeloid_Mono = c("LYZ","S100A8","S100A9","FCGR3A"),
  DC = c("ITGAX","CD1C","CLEC9A","LAMP3"),
  Fibroblast = c("COL1A1","COL1A2","DCN","LUM","PDGFRB","COL3A1"),
  Endothelial = c("PECAM1","VWF","KDR","CLDN5","ESAM"),
  Mast = c("TPSAB1","CPA3","KIT"),
  NK = c("NKG7","KLRD1","KLRB1","GNLY","XCL1")
)

plot_canonical_markers <- function(obj, markers=canonical_markers, prefix="merged"){
  DefaultAssay(obj) <- if("SCT" %in% names(obj@assays)) "SCT" else "RNA"
  feats <- unique(unlist(markers))
  feats <- feats[feats %in% rownames(obj)]
  if(length(feats) > 0){
    p1 <- DotPlot(obj, features = feats, group.by = "seurat_clusters") + RotatedAxis()
    ggsave(sprintf("plots/%s_marker_dotplot.pdf", prefix), p1, width=10, height=6)
  }
  # 各大类挑代表基因画FeaturePlot
  for(cat in names(markers)){
    f <- intersect(markers[[cat]], rownames(obj))
    if(length(f)>0){
      p <- FeaturePlot(obj, features = head(f, 6), ncol = 3, order = TRUE)
      ggsave(sprintf("plots/%s_feature_%s.pdf", prefix, cat), p, width=8, height=6)
    }
  }
}

subcluster_lineage <- function(obj, lineage, features=NULL, resolution=0.4, prefix="merged"){
  cells <- WhichCells(obj, expression = celltype_major == lineage)
  if(length(cells) < 200) {
    message(sprintf("子群 '%s' 细胞过少(%d)，跳过子聚类。", lineage, length(cells)))
    return(NULL)
  }
  sub <- subset(obj, cells = cells)
  DefaultAssay(sub) <- if("SCT" %in% names(sub@assays)) "SCT" else "RNA"
  # 重新归一化与降维
  sub <- perform_normalization(sub, method = if("SCT" %in% names(obj@assays)) "SCT" else "LogNormalize")
  sub <- RunPCA(sub, npcs = 30, verbose = FALSE)
  sub <- RunUMAP(sub, dims = 1:30, n.neighbors = 30, verbose = FALSE)
  sub <- FindNeighbors(sub, dims = 1:30)
  sub <- FindClusters(sub, resolution = resolution)
  # markers
  Idents(sub) <- sub$seurat_clusters
  markers <- FindAllMarkers(sub, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
  fwrite(markers, sprintf("tables/%s_%s_sub_markers.csv", prefix, lineage))
  # 可视化
  p1 <- DimPlot(sub, reduction = "umap", label = TRUE) + ggtitle(sprintf("%s subclusters", lineage))
  ggsave(sprintf("plots/%s_%s_sub_umap.pdf", prefix, lineage), p1, width=6, height=5)
  # 返回
  saveRDS(sub, sprintf("rds/%s_%s_sub.rds", prefix, lineage))
  sub
}

assign_major_celltypes <- function(obj){
  # 基于SingleR/Azimuth与marker的规则融合示例
  major <- rep(NA_character_, ncol(obj))
  if("SingleR_label" %in% colnames(obj@meta.data)){
    major <- obj$SingleR_label
  }
  # 简单规则覆盖：若表达明显肿瘤上皮标志则设为Epithelial_Tumor
  DefaultAssay(obj) <- if("SCT" %in% names(obj@assays)) "SCT" else "RNA"
  expr <- GetAssayData(obj, slot = "data")
  epi_score <- Matrix::colMeans(as.matrix(expr[intersect(canonical_markers$Epithelial_Tumor, rownames(expr)), , drop=FALSE])) 
  t_score <- Matrix::colMeans(as.matrix(expr[intersect(c("CD3D","CD3E","TRAC"), rownames(expr)), , drop=FALSE]))
  b_score <- Matrix::colMeans(as.matrix(expr[intersect(c("MS4A1","CD79A"), rownames(expr)), , drop=FALSE]))
  myeloid_score <- Matrix::colMeans(as.matrix(expr[intersect(c("LYZ","S100A8","S100A9"), rownames(expr)), , drop=FALSE]))
  fibro_score <- Matrix::colMeans(as.matrix(expr[intersect(c("COL1A1","COL1A2","DCN"), rownames(expr)), , drop=FALSE]))
  endo_score <- Matrix::colMeans(as.matrix(expr[intersect(c("PECAM1","VWF","KDR"), rownames(expr)), , drop=FALSE]))
  
  # 阈值策略（可调）
  # 用相对最大得分决定
  scores <- rbind(epi=epi_score, Tcell=t_score, Bcell=b_score, Myeloid=myeloid_score, Fibro=fibro_score, Endo=endo_score)
  max_idx <- apply(scores, 2, which.max)
  label_map <- c("Epithelial_Tumor","T_cell","B_cell","Myeloid","Fibroblast","Endothelial")
  inferred <- label_map[max_idx]
  
  # 融合：优先规则打标签，否则保留SingleR
  major[is.na(major)] <- inferred[is.na(major)]
  obj$celltype_major <- major
  
  # 基于PDCD1/CTLA4等注释耗竭T
  if(all(c("PDCD1","LAG3","HAVCR2") %in% rownames(obj))){
    ex_score <- Matrix::colMeans(as.matrix(expr[c("PDCD1","LAG3","HAVCR2"), , drop=FALSE]))
    obj$T_exhausted <- ex_score
  }
  obj
}

# ============ 主流程 ============

set.seed(1)

objs <- read_all_seurat("rds", pattern = "_seurat\\.rds$")
# 合并
if(length(objs) == 1){
  seu <- objs[[1]]
} else {
  # 使用Seurat合并，不在此处做整合（整合在后续步骤进行）
  seu <- Reduce(function(a,b) merge(a, y=b, project = "PanCancer", add.cell.ids = c(a@project.name, b@project.name)), objs)
  # 统一orig.ident
  if(is.null(seu$orig.ident)) seu$orig.ident <- seu$sample
}

# 基本QC可视化
basic_qc_plots(seu, prefix = "merged_initial")

# 归一化与变量基因
seu <- perform_normalization(seu, method = normalize_method)

# 降维与整合+聚类
seu <- cluster_and_embed(seu, npcs = npcs, method = integration_method)

# 自动注释（可选）
if(use_singler){
  seu <- auto_annotation_singler(seu)
}
if(use_azimuth){
  seu <- azimuth_annotation(seu)
}

# 结合marker与规则，得到大类注释
seu <- assign_major_celltypes(seu)

# 可视化与导出
plot_embeddings(seu, prefix = "merged")
plot_canonical_markers(seu, canonical_markers, prefix = "merged")
export_results(seu, prefix = "merged")

# ============ 子群细化（可选）===========
if(do_subclustering){
  # 针对主要大类分别做子聚类，如T细胞、髓系
  t_sub <- subcluster_lineage(seu, "T_cell", resolution = 0.6, prefix = "merged")
  my_sub <- subcluster_lineage(seu, "Myeloid", resolution = 0.6, prefix = "merged")
  b_sub <- subcluster_lineage(seu, "B_cell", resolution = 0.6, prefix = "merged")
  epi_sub <- subcluster_lineage(seu, "Epithelial_Tumor", resolution = 0.6, prefix = "merged")
}

message("完成：聚类与细胞类型注释。输出已保存到 rds/, plots/, tables/。")