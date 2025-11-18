# ============================================================================
# KEGG数据库下载脚本
# 用途: 预先下载KEGG数据，避免在线访问问题
# ============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# 设置输出目录（与主脚本的output_dir一致）
output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/annotation_evaluation_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 设置下载参数
options(timeout = 600)  # 10分钟超时
options(download.file.method = "libcurl")

message("=" , rep("=", 70))
message("KEGG数据库下载工具")
message("=" , rep("=", 70))

# ============================================================================
# 方法1: 使用clusterProfiler下载KEGG数据
# ============================================================================
download_kegg_clusterProfiler <- function() {
  message("\n方法1: 使用clusterProfiler下载KEGG数据...")
  
  tryCatch({
    # 搜索人类organism
    kegg_org <- search_kegg_organism('Homo sapiens', by='scientific_name')
    message("✓ 找到KEGG物种代码: ", kegg_org$kegg_code[1])
    
    # 下载KEGG pathway数据
    message("正在下载KEGG pathway数据 (可能需要几分钟)...")
    
    kegg_data <- download_KEGG(species = "hsa", keggType = "KEGG", keyType = "kegg")
    
    # 保存数据
    kegg_file <- file.path(output_dir, "kegg_hsa_database.RData")
    save(kegg_data, file = kegg_file)
    
    message("✓ KEGG数据已保存到: ", kegg_file)
    message("  - KEGG pathways: ", length(unique(kegg_data$from)))
    
    return(TRUE)
    
  }, error = function(e) {
    message("✗ 方法1失败: ", e$message)
    return(FALSE)
  })
}

# ============================================================================
# 方法2: 使用KEGGREST下载
# ============================================================================
download_kegg_rest <- function() {
  message("\n方法2: 使用KEGGREST直接下载...")
  
  # 安装KEGGREST（如果未安装）
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    message("安装KEGGREST包...")
    BiocManager::install("KEGGREST")
  }
  
  library(KEGGREST)
  
  tryCatch({
    message("下载KEGG pathway列表...")
    kegg_pathways <- keggList("pathway", "hsa")
    message("✓ 获得 ", length(kegg_pathways), " 个pathways")
    
    message("下载pathway-gene关联...")
    kegg_genes <- keggLink("hsa", "pathway")
    message("✓ 获得 ", length(kegg_genes), " 个基因关联")
    
    # 整理数据格式
    kegg_data_rest <- list(
      pathways = kegg_pathways,
      genes = kegg_genes,
      download_date = Sys.Date()
    )
    
    # 保存
    kegg_rest_file <- file.path(output_dir, "kegg_rest_database.RData")
    save(kegg_data_rest, file = kegg_rest_file)
    
    message("✓ KEGGREST数据已保存到: ", kegg_rest_file)
    
    return(TRUE)
    
  }, error = function(e) {
    message("✗ 方法2失败: ", e$message)
    return(FALSE)
  })
}

# ============================================================================
# 方法3: 使用MSigDB获取KEGG基因集（推荐，最稳定）
# ============================================================================
download_msigdb_kegg <- function() {
  message("\n方法3: 从MSigDB下载KEGG基因集 (推荐)...")
  
  # 安装msigdbr（如果未安装）
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    message("安装msigdbr包...")
    install.packages("msigdbr")
  }
  
  library(msigdbr)
  library(dplyr)
  
  tryCatch({
    message("下载KEGG基因集...")
    msigdb_kegg <- msigdbr(
      species = "Homo sapiens", 
      category = "C2", 
      subcategory = "CP:KEGG"
    )
    
    message("✓ 获得 ", length(unique(msigdb_kegg$gs_name)), " 个KEGG pathways")
    message("✓ 包含 ", length(unique(msigdb_kegg$entrez_gene)), " 个基因")
    
    # 转换为clusterProfiler格式
    kegg_term2gene <- msigdb_kegg %>%
      dplyr::select(gs_name, entrez_gene) %>%
      as.data.frame()
    
    kegg_term2name <- msigdb_kegg %>%
      dplyr::select(gs_name, gs_description) %>%
      distinct() %>%
      as.data.frame()
    
    msigdb_data <- list(
      term2gene = kegg_term2gene,
      term2name = kegg_term2name,
      full_data = msigdb_kegg,
      download_date = Sys.Date()
    )
    
    # 保存
    msigdb_file <- file.path(output_dir, "msigdb_kegg_genesets.RData")
    save(msigdb_data, file = msigdb_file)
    
    message("✓ MSigDB KEGG数据已保存到: ", msigdb_file)
    
    # 显示前几个pathway
    message("\n示例pathway:")
    print(head(kegg_term2name, 3))
    
    return(TRUE)
    
  }, error = function(e) {
    message("✗ 方法3失败: ", e$message)
    return(FALSE)
  })
}

# ============================================================================
# 方法4: 下载Reactome数据（作为KEGG的替代）
# ============================================================================
download_reactome <- function() {
  message("\n方法4: 下载Reactome数据库...")
  
  # 安装ReactomePA（如果未安装）
  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    message("安装ReactomePA包...")
    BiocManager::install("ReactomePA")
  }
  
  library(ReactomePA)
  library(msigdbr)
  
  tryCatch({
    # 从MSigDB获取Reactome数据
    message("从MSigDB下载Reactome基因集...")
    msigdb_reactome <- msigdbr(
      species = "Homo sapiens", 
      category = "C2", 
      subcategory = "CP:REACTOME"
    )
    
    message("✓ 获得 ", length(unique(msigdb_reactome$gs_name)), " 个Reactome pathways")
    
    reactome_term2gene <- msigdb_reactome %>%
      dplyr::select(gs_name, entrez_gene) %>%
      as.data.frame()
    
    reactome_term2name <- msigdb_reactome %>%
      dplyr::select(gs_name, gs_description) %>%
      distinct() %>%
      as.data.frame()
    
    reactome_data <- list(
      term2gene = reactome_term2gene,
      term2name = reactome_term2name,
      full_data = msigdb_reactome,
      download_date = Sys.Date()
    )
    
    # 保存
    reactome_file <- file.path(output_dir, "reactome_genesets.RData")
    save(reactome_data, file = reactome_file)
    
    message("✓ Reactome数据已保存到: ", reactome_file)
    
    return(TRUE)
    
  }, error = function(e) {
    message("✗ 方法4失败: ", e$message)
    return(FALSE)
  })
}

# ============================================================================
# 主执行函数
# ============================================================================
main <- function() {
  message("\n开始下载数据库...")
  message("输出目录: ", output_dir)
  
  success_count <- 0
  
  # 尝试所有方法
  if (download_msigdb_kegg()) success_count <- success_count + 1
  
  Sys.sleep(2)  # 避免API限制
  
  if (download_reactome()) success_count <- success_count + 1
  
  Sys.sleep(2)
  
  if (download_kegg_rest()) success_count <- success_count + 1
  
  Sys.sleep(2)
  
  if (download_kegg_clusterProfiler()) success_count <- success_count + 1
  
  # 总结
  message("\n", rep("=", 70))
  message("下载完成!")
  message("成功下载: ", success_count, " / 4 个数据库")
  message(rep("=", 70))
  
  # 列出已下载的文件
  message("\n已下载文件:")
  files <- list.files(output_dir, pattern = "\\.RData$", full.names = TRUE)
  for (f in files) {
    size <- file.info(f)$size / 1024 / 1024  # MB
    message("  - ", basename(f), " (", round(size, 2), " MB)")
  }
  
  if (success_count > 0) {
    message("\n✓ 数据库下载成功！可以运行主分析脚本了。")
  } else {
    message("\n✗ 所有下载方法都失败了。请检查网络连接。")
  }
}

# 执行下载
main()