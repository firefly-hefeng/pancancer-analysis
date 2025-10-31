# quick_fix_annotations.R
# 快速修复明显的注释错误

library(Seurat)
library(dplyr)

# 快速修复函数
quick_fix_cell_types <- function(fine_cell_types) {
  # 定义明显错误的模式
  error_patterns <- list(
    # B细胞不应该有的亚型
    "B_cell_(Classical_monocytes|Intermediate_monocytes|Non_classical_monocytes|Myeloid_dendritic_cells|Plasmacytoid_dendritic_cells|Natural_killer_cells|Th1_cells|Progenitor_cells)" = "B_cell_Unassigned",
    
    # 内皮细胞不应该有的亚型  
    "Endothelial_cells_(Classical_monocytes|Intermediate_monocytes|Non_classical_monocytes|Myeloid_dendritic_cells|Plasmacytoid_dendritic_cells|Natural_killer_cells|Th1_cells|Progenitor_cells)" = "Endothelial_cells_Unassigned",
    
    # 上皮细胞不应该有的亚型
    "Epithelial_cells_(Classical_monocytes|Intermediate_monocytes|Non_classical_monocytes|Myeloid_dendritic_cells|Plasmacytoid_dendritic_cells|Natural_killer_cells|Th1_cells|Progenitor_cells|Plasmablasts)" = "Epithelial_cells_Unassigned",
    
    # 成纤维细胞不应该有的亚型
    "Fibroblasts_(Classical_monocytes|Intermediate_monocytes|Non_classical_monocytes|Myeloid_dendritic_cells|Plasmacytoid_dendritic_cells|Exhausted_B_cells|Th1_cells|Progenitor_cells)" = "Fibroblasts_Unassigned",
    
    # NK细胞不应该有的亚型
    "NK_cell_(Exhausted_B_cells|Intermediate_monocytes|Myeloid_dendritic_cells|Plasmablasts)" = "NK_cell_Unassigned"
  )
  
  fixed_types <- fine_cell_types
  
  for (pattern in names(error_patterns)) {
    replacement <- error_patterns[[pattern]]
    fixed_types <- gsub(pattern, replacement, fixed_types, perl = TRUE)
  }
  
  return(fixed_types)
}

# 批量修复所有癌种
quick_fix_all_cancers <- function() {
  output_base_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation"
  analysis_summary <- read.csv(file.path(output_base_dir, "analysis_summary.csv"), stringsAsFactors = FALSE)
  
  for (i in 1:nrow(analysis_summary)) {
    cancer_type <- analysis_summary$Cancer_Type[i]
    
    main_file <- file.path(output_base_dir, cancer_type, paste0(cancer_type, "_fine_annotated.rds"))
    
    if (file.exists(main_file)) {
      seurat_obj <- readRDS(main_file)
      
      # 备份原始注释
      seurat_obj$fine_cell_type_original <- seurat_obj$fine_cell_type
      
      # 快速修复
      seurat_obj$fine_cell_type <- quick_fix_cell_types(seurat_obj$fine_cell_type)
      
      # 保存修复版本
      fixed_file <- file.path(output_base_dir, cancer_type, paste0(cancer_type, "_fine_annotated_quickfix.rds"))
      saveRDS(seurat_obj, fixed_file)
      
      cat("Fixed:", cancer_type, "\n")
    }
  }
}

# 运行快速修复
quick_fix_all_cancers()