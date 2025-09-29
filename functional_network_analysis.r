# functional_network_analysis.R
# 轨迹相关基因的功能富集、转录因子调控网络和代谢分析

library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)
library(pathview)
library(enrichplot)
library(ggupset)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(ggraph)
library(tidygraph)

# 设置路径
trajectory_output_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation/trajectory_analysis_2"
functional_output_dir <- file.path(trajectory_output_dir, "functional_network_analysis")
dir.create(functional_output_dir, showWarnings = FALSE, recursive = TRUE)

# 创建子目录
sub_dirs <- c("enrichment_analysis", "tf_networks", "metabolic_analysis", 
              "integrated_networks", "pathway_plots", "network_plots")
for(dir in sub_dirs) {
  dir.create(file.path(functional_output_dir, dir), showWarnings = FALSE)
}

# ==================== 基因富集分析 ====================

# 综合基因富集分析
comprehensive_gene_enrichment <- function(gene_list, 
                                        background_genes = NULL,
                                        organism = "hsa",
                                        p_cutoff = 0.05,
                                        q_cutoff = 0.1) {
  
  message("开始综合基因富集分析...")
  
  if(length(gene_list) < 5) {
    stop("基因数量太少，无法进行有效的富集分析")
  }
  
  # 基因ID转换
  gene_df <- tryCatch({
    bitr(gene_list, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), 
         OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    # 如果是ENSEMBL ID，尝试转换
    bitr(gene_list, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), 
         OrgDb = org.Hs.eg.db)
  })
  
  if(is.null(gene_df) || nrow(gene_df) == 0) {
    stop("基因ID转换失败")
  }
  
  entrez_ids <- unique(gene_df$ENTREZID)
  gene_symbols <- unique(gene_df$SYMBOL)
  
  message(paste("成功转换", length(entrez_ids), "个基因ID"))
  
  enrichment_results <- list()
  
  # 1. GO富集分析
  message("执行GO富集分析...")
  
  # GO生物过程
  tryCatch({
    go_bp <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = p_cutoff,
                     qvalueCutoff = q_cutoff,
                     readable = TRUE)
    
    if(!is.null(go_bp) && nrow(go_bp@result) > 0) {
      enrichment_results$GO_BP <- go_bp
      message(paste("GO BP富集到", nrow(go_bp@result), "个条目"))
    }
  }, error = function(e) {
    message("GO BP富集失败:", e$message)
  })
  
  # GO分子功能
  tryCatch({
    go_mf <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = p_cutoff,
                     qvalueCutoff = q_cutoff,
                     readable = TRUE)
    
    if(!is.null(go_mf) && nrow(go_mf@result) > 0) {
      enrichment_results$GO_MF <- go_mf
    }
  }, error = function(e) {
    message("GO MF富集失败:", e$message)
  })
  
  # GO细胞组分
  tryCatch({
    go_cc <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = p_cutoff,
                     qvalueCutoff = q_cutoff,
                     readable = TRUE)
    
    if(!is.null(go_cc) && nrow(go_cc@result) > 0) {
      enrichment_results$GO_CC <- go_cc
    }
  }, error = function(e) {
    message("GO CC富集失败:", e$message)
  })
  
  # 2. KEGG通路富集
  message("执行KEGG通路富集...")
  tryCatch({
    kegg <- enrichKEGG(gene = entrez_ids,
                      organism = organism,
                      pvalueCutoff = p_cutoff,
                      qvalueCutoff = q_cutoff)
    
    if(!is.null(kegg) && nrow(kegg@result) > 0) {
      enrichment_results$KEGG <- kegg
      message(paste("KEGG富集到", nrow(kegg@result), "个通路"))
    }
  }, error = function(e) {
    message("KEGG富集失败:", e$message)
  })
  
  # 3. Reactome通路富集
  message("执行Reactome通路富集...")
  tryCatch({
    library(ReactomePA)
    reactome <- enrichPathway(gene = entrez_ids,
                            pvalueCutoff = p_cutoff,
                            qvalueCutoff = q_cutoff,
                            readable = TRUE)
    
    if(!is.null(reactome) && nrow(reactome@result) > 0) {
      enrichment_results$Reactome <- reactome
      message(paste("Reactome富集到", nrow(reactome@result), "个通路"))
    }
  }, error = function(e) {
    message("Reactome富集失败:", e$message)
  })
  
  # 4. WikiPathways富集
  message("执行WikiPathways富集...")
  tryCatch({
    library(clusterProfiler)
    wp <- enrichWP(gene = entrez_ids,
                  organism = "Homo sapiens",
                  pvalueCutoff = p_cutoff,
                  qvalueCutoff = q_cutoff)
    
    if(!is.null(wp) && nrow(wp@result) > 0) {
      enrichment_results$WikiPathways <- wp
    }
  }, error = function(e) {
    message("WikiPathways富集失败:", e$message)
  })
  
  # 5. MSigDB免疫签名富集
  message("执行免疫相关基因集富集...")
  tryCatch({
    library(msigdbr)
    
    # 获取免疫相关基因集
    immune_signatures <- msigdbr(species = "Homo sapiens", category = "C7")
    immune_t2g <- immune_signatures %>% 
      select(gs_name, entrez_gene) %>% 
      filter(grepl("T_CELL|TCELL|CD8|CD4|EXHAUSTION|ACTIVATION", gs_name, ignore.case = TRUE))
    
    if(nrow(immune_t2g) > 0) {
      immune_gsea <- enricher(gene = entrez_ids,
                             TERM2GENE = immune_t2g,
                             pvalueCutoff = p_cutoff,
                             qvalueCutoff = q_cutoff)
      
      if(!is.null(immune_gsea) && nrow(immune_gsea@result) > 0) {
        enrichment_results$Immune_Signatures <- immune_gsea
      }
    }
  }, error = function(e) {
    message("免疫签名富集失败:", e$message)
  })
  
  # 6. 疾病本体富集（DO）
  message("执行疾病本体富集...")
  tryCatch({
    library(DOSE)
    do_result <- enrichDO(gene = entrez_ids,
                         ont = "DO",
                         pvalueCutoff = p_cutoff,
                         qvalueCutoff = q_cutoff,
                         readable = TRUE)
    
    if(!is.null(do_result) && nrow(do_result@result) > 0) {
      enrichment_results$Disease_Ontology <- do_result
    }
  }, error = function(e) {
    message("DO富集失败:", e$message)
  })
  
  return(enrichment_results)
}

# 基于轨迹方向的分层富集分析
trajectory_directional_enrichment <- function(trajectory_genes, 
                                             correlation_threshold = 0.2) {
  
  message("执行轨迹方向性富集分析...")
  
  # 分离正相关和负相关基因
  positive_genes <- trajectory_genes$gene_short_name[
    trajectory_genes$correlation > correlation_threshold & 
    trajectory_genes$q_value < 0.05
  ]
  
  negative_genes <- trajectory_genes$gene_short_name[
    trajectory_genes$correlation < -correlation_threshold & 
    trajectory_genes$q_value < 0.05
  ]
  
  # 转换基因名
  positive_symbols <- convert_gene_ids(positive_genes)
  negative_symbols <- convert_gene_ids(negative_genes)
  
  results <- list()
  
  # 正相关基因富集（轨迹终点富集基因）
  if(length(positive_symbols) >= 5) {
    message(paste("分析", length(positive_symbols), "个正相关基因"))
    results$positive_correlation <- comprehensive_gene_enrichment(positive_symbols)
    results$positive_genes <- positive_symbols
  }
  
  # 负相关基因富集（轨迹起点富集基因）
  if(length(negative_symbols) >= 5) {
    message(paste("分析", length(negative_symbols), "个负相关基因"))
    results$negative_correlation <- comprehensive_gene_enrichment(negative_symbols)
    results$negative_genes <- negative_symbols
  }
  
  return(results)
}

# ==================== 转录因子调控网络分析 ====================

# 转录因子靶基因富集分析
tf_target_enrichment_analysis <- function(trajectory_genes, 
                                        correlation_threshold = 0.2) {
  
  message("分析转录因子-靶基因调控网络...")
  
  # 加载转录因子数据库
  library(dorothea)
  
  # 获取人类转录因子活性数据
  tf_regulons <- dorothea_hs
  
  # 筛选显著轨迹基因
  significant_genes <- trajectory_genes[
    abs(trajectory_genes$correlation) > correlation_threshold & 
    trajectory_genes$q_value < 0.05, 
  ]
  
  if(nrow(significant_genes) == 0) {
    message("没有显著的轨迹基因")
    return(NULL)
  }
  
  # 转换基因名
  gene_symbols <- convert_gene_ids(significant_genes$gene_short_name)
  
  if(length(gene_symbols) == 0) {
    message("基因名转换失败")
    return(NULL)
  }
  
  # TF富集分析
  tf_enrichment_results <- list()
  
  # 使用dorothea进行TF活性预测
  tryCatch({
    # 准备TF-target基因集
    tf_sets <- split(tf_regulons$target, tf_regulons$tf)
    
    # 对每个TF进行超几何检验
    tf_pvalues <- sapply(names(tf_sets), function(tf) {
      tf_targets <- tf_sets[[tf]]
      overlap <- intersect(gene_symbols, tf_targets)
      
      if(length(overlap) >= 2) {
        # 超几何检验
        phyper(length(overlap) - 1, 
               length(tf_targets), 
               20000 - length(tf_targets), 
               length(gene_symbols), 
               lower.tail = FALSE)
      } else {
        1
      }
    })
    
    # 多重检验校正
    tf_qvalues <- p.adjust(tf_pvalues, method = "BH")
    
    # 整理结果
    tf_results <- data.frame(
      TF = names(tf_pvalues),
      p_value = tf_pvalues,
      q_value = tf_qvalues,
      n_targets_in_trajectory = sapply(names(tf_sets), function(tf) {
        length(intersect(gene_symbols, tf_sets[[tf]]))
      }),
      total_targets = sapply(tf_sets, length),
      stringsAsFactors = FALSE
    )
    
    # 过滤显著的TF
    significant_tfs <- tf_results[tf_results$q_value < 0.05 & 
                                 tf_results$n_targets_in_trajectory >= 2, ]
    
    if(nrow(significant_tfs) > 0) {
      significant_tfs <- significant_tfs[order(significant_tfs$q_value), ]
      tf_enrichment_results$tf_enrichment <- significant_tfs
      
      message(paste("发现", nrow(significant_tfs), "个显著调控TF"))
    }
    
  }, error = function(e) {
    message("TF富集分析失败:", e$message)
  })
  
  # 构建TF调控网络
  if(length(tf_enrichment_results) > 0) {
    tf_network <- build_tf_regulatory_network(
      significant_tfs$TF, 
      gene_symbols, 
      tf_regulons
    )
    tf_enrichment_results$regulatory_network <- tf_network
  }
  
  return(tf_enrichment_results)
}

# 构建转录因子调控网络
build_tf_regulatory_network <- function(significant_tfs, target_genes, tf_regulons) {
  
  message("构建TF调控网络...")
  
  # 筛选相关的TF-target关系
  network_edges <- tf_regulons[
    tf_regulons$tf %in% significant_tfs & 
    tf_regulons$target %in% target_genes,
  ]
  
  if(nrow(network_edges) == 0) {
    message("没有足够的TF-target关系构建网络")
    return(NULL)
  }
  
  # 构建igraph网络
  library(igraph)
  
  g <- graph_from_data_frame(
    d = network_edges[, c("tf", "target")],
    directed = TRUE
  )
  
  # 计算网络属性
  network_stats <- list(
    n_nodes = vcount(g),
    n_edges = ecount(g),
    n_tfs = length(unique(network_edges$tf)),
    n_targets = length(unique(network_edges$target)),
    density = edge_density(g),
    avg_degree = mean(degree(g))
  )
  
  # 计算中心性指标
  V(g)$betweenness <- betweenness(g)
  V(g)$degree <- degree(g)
  V(g)$pagerank <- page_rank(g)$vector
  
  # 识别节点类型
  V(g)$node_type <- ifelse(V(g)$name %in% significant_tfs, "TF", "Target")
  
  return(list(
    graph = g,
    network_stats = network_stats,
    edges = network_edges
  ))
}

# ==================== 代谢轨迹分析 ====================

# 代谢通路活性分析
metabolic_trajectory_analysis <- function(seurat_obj, 
                                        pseudotime_values,
                                        metabolism_method = "AUCell") {
  
  message("分析代谢轨迹...")
  
  # 检查包的可用性
  if(!requireNamespace("scMetabolism", quietly = TRUE)) {
    message("scMetabolism包不可用，使用简化的代谢分析")
    return(simplified_metabolic_analysis(seurat_obj, pseudotime_values))
  }
  
  library(scMetabolism)
  
  # 使用scMetabolism计算代谢通路活性
  tryCatch({
    metabolism_scores <- sc.metabolism.Seurat(
      obj = seurat_obj, 
      method = metabolism_method,
      imputation = FALSE,
      ncores = 2
    )
    
    # 提取代谢评分
    metabolic_pathways <- rownames(metabolism_scores@assays$METABOLISM@data)
    
    # 分析代谢评分与pseudotime的相关性
    common_cells <- intersect(names(pseudotime_values), colnames(metabolism_scores))
    
    if(length(common_cells) < 20) {
      message("共同细胞数量不足")
      return(NULL)
    }
    
    metabolic_correlations <- data.frame()
    
    for(pathway in metabolic_pathways) {
      pathway_scores <- GetAssayData(metabolism_scores, assay = "METABOLISM")[pathway, common_cells]
      pt_values <- pseudotime_values[common_cells]
      
      # 计算相关性
      if(sum(!is.na(pathway_scores)) > 10 && sum(!is.na(pt_values)) > 10) {
        cor_test <- cor.test(pathway_scores, pt_values, method = "spearman")
        
        metabolic_correlations <- rbind(metabolic_correlations, data.frame(
          pathway = pathway,
          correlation = cor_test$estimate,
          p_value = cor_test$p.value,
          n_cells = length(common_cells),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # 多重检验校正
    if(nrow(metabolic_correlations) > 0) {
      metabolic_correlations$q_value <- p.adjust(metabolic_correlations$p_value, method = "BH")
      metabolic_correlations <- metabolic_correlations[order(abs(metabolic_correlations$correlation), 
                                                           decreasing = TRUE), ]
      
      # 识别显著的代谢通路
      significant_pathways <- metabolic_correlations[
        metabolic_correlations$q_value < 0.05 & 
        abs(metabolic_correlations$correlation) > 0.2,
      ]
      
      message(paste("发现", nrow(significant_pathways), "个显著的代谢轨迹通路"))
      
      return(list(
        all_correlations = metabolic_correlations,
        significant_pathways = significant_pathways,
        metabolism_scores = metabolism_scores
      ))
    }
    
  }, error = function(e) {
    message("scMetabolism分析失败:", e$message)
    return(simplified_metabolic_analysis(seurat_obj, pseudotime_values))
  })
}

# 简化的代谢基因分析
simplified_metabolic_analysis <- function(seurat_obj, pseudotime_values) {
  
  message("执行简化的代谢基因分析...")
  
  # 定义关键代谢通路基因
  metabolic_gene_sets <- list(
    Glycolysis = c("HK1", "HK2", "GPI", "PFKM", "ALDOA", "TPI1", "GAPDH", 
                   "PGK1", "PGAM1", "ENO1", "PKM", "LDHA", "LDHB"),
    
    Oxidative_Phosphorylation = c("NDUFB1", "NDUFB5", "SDHB", "UQCRC1", 
                                 "COX4I1", "ATP5A1", "ATP5B", "ATP5G1"),
    
    Fatty_Acid_Synthesis = c("ACACA", "FASN", "SCD", "ELOVL6", "FADS1", "FADS2"),
    
    Fatty_Acid_Oxidation = c("CPT1A", "CPT1B", "ACADM", "HADHA", "ACOX1"),
    
    Pentose_Phosphate = c("G6PD", "PGD", "RPIA", "RPE", "TALDO1", "TKT"),
    
    Glutaminolysis = c("GLS", "GLS2", "GLUD1", "GLUD2", "GOT1", "GOT2"),
    
    One_Carbon = c("SHMT1", "SHMT2", "MTHFD1", "MTHFD2", "TYMS", "DHFR"),
    
    TCA_Cycle = c("CS", "ACO2", "IDH1", "IDH2", "OGDH", "SUCLA2", "SDHA", "FH", "MDH2")
  )
  
  # 计算每个通路的活性评分
  pathway_scores <- list()
  pathway_correlations <- data.frame()
  
  common_cells <- intersect(names(pseudotime_values), colnames(seurat_obj))
  
  for(pathway_name in names(metabolic_gene_sets)) {
    genes <- metabolic_gene_sets[[pathway_name]]
    available_genes <- intersect(genes, rownames(seurat_obj))
    
    if(length(available_genes) >= 3) {
      # 计算通路评分（平均表达）
      expr_data <- GetAssayData(seurat_obj, layer = "data")[available_genes, common_cells, drop = FALSE]
      pathway_score <- colMeans(expr_data)
      
      # 与pseudotime相关性
      cor_test <- cor.test(pathway_score, pseudotime_values[common_cells], method = "spearman")
      
      pathway_correlations <- rbind(pathway_correlations, data.frame(
        pathway = pathway_name,
        correlation = cor_test$estimate,
        p_value = cor_test$p.value,
        n_genes = length(available_genes),
        n_cells = length(common_cells),
        stringsAsFactors = FALSE
      ))
      
      pathway_scores[[pathway_name]] <- pathway_score
    }
  }
  
  # 多重检验校正
  if(nrow(pathway_correlations) > 0) {
    pathway_correlations$q_value <- p.adjust(pathway_correlations$p_value, method = "BH")
    
    return(list(
      pathway_correlations = pathway_correlations,
      pathway_scores = pathway_scores,
      metabolic_gene_sets = metabolic_gene_sets
    ))
  } else {
    return(NULL)
  }
}

# ==================== 整合网络分析 ====================

# 构建基因-功能-TF整合网络
build_integrated_network <- function(enrichment_results, 
                                   tf_results, 
                                   metabolic_results,
                                   trajectory_genes) {
  
  message("构建整合调控网络...")
  
  # 准备节点和边的数据框
  nodes <- data.frame()
  edges <- data.frame()
  
  # 添加基因节点
  gene_symbols <- convert_gene_ids(trajectory_genes$gene_short_name)
  if(length(gene_symbols) > 0) {
    gene_nodes <- data.frame(
      id = gene_symbols,
      type = "Gene",
      correlation = trajectory_genes$correlation[match(names(gene_symbols), 
                                                      trajectory_genes$gene_short_name)],
      stringsAsFactors = FALSE
    )
    nodes <- rbind(nodes, gene_nodes)
  }
  
  # 添加TF节点和边
  if(!is.null(tf_results) && "tf_enrichment" %in% names(tf_results)) {
    tf_data <- tf_results$tf_enrichment
    
    # TF节点
    tf_nodes <- data.frame(
      id = tf_data$TF,
      type = "TF",
      significance = -log10(tf_data$q_value),
      stringsAsFactors = FALSE
    )
    
    nodes <- rbind(nodes, tf_nodes)
    
    # TF-Gene边
    if("regulatory_network" %in% names(tf_results) && 
       !is.null(tf_results$regulatory_network)) {
      tf_edges <- tf_results$regulatory_network$edges
      tf_edges_df <- data.frame(
        from = tf_edges$tf,
        to = tf_edges$target,
        type = "TF_regulation",
        stringsAsFactors = FALSE
      )
      edges <- rbind(edges, tf_edges_df)
    }
  }
  
  # 添加功能模块节点和边
  if(!is.null(enrichment_results) && length(enrichment_results) > 0) {
    for(db_name in names(enrichment_results)) {
      enrich_result <- enrichment_results[[db_name]]
      
      if(!is.null(enrich_result) && nrow(enrich_result@result) > 0) {
        # 选择top富集条目
        top_terms <- head(enrich_result@result[enrich_result@result$p.adjust < 0.05, ], 10)
        
        if(nrow(top_terms) > 0) {
          # 功能节点
          func_nodes <- data.frame(
            id = top_terms$ID,
            type = paste0("Function_", db_name),
            description = top_terms$Description,
            significance = -log10(top_terms$p.adjust),
            stringsAsFactors = FALSE
          )
          
          nodes <- rbind(nodes, func_nodes)
          
          # 基因-功能边
          for(i in 1:nrow(top_terms)) {
            term_genes <- strsplit(top_terms$geneID[i], "/")[[1]]
            gene_func_edges <- data.frame(
              from = term_genes,
              to = top_terms$ID[i],
              type = "Gene_function",
              stringsAsFactors = FALSE
            )
            edges <- rbind(edges, gene_func_edges)
          }
        }
      }
    }
  }
  
  # 添加代谢通路节点和边
  if(!is.null(metabolic_results)) {
    if("significant_pathways" %in% names(metabolic_results)) {
      metab_pathways <- metabolic_results$significant_pathways
      
      if(nrow(metab_pathways) > 0) {
        # 代谢通路节点
        metab_nodes <- data.frame(
          id = metab_pathways$pathway,
          type = "Metabolic_pathway",
          correlation = metab_pathways$correlation,
          significance = -log10(metab_pathways$q_value),
          stringsAsFactors = FALSE
        )
        
        nodes <- rbind(nodes, metab_nodes)
      }
    }
  }
  
  # 移除重复节点
  nodes <- nodes[!duplicated(nodes$id), ]
  
  # 构建igraph网络
  if(nrow(nodes) > 0 && nrow(edges) > 0) {
    library(igraph)
    
    # 过滤有效边
    valid_edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id, ]
    
    if(nrow(valid_edges) > 0) {
      g <- graph_from_data_frame(d = valid_edges, vertices = nodes, directed = TRUE)
      
      return(list(
        graph = g,
        nodes = nodes,
        edges = valid_edges,
        network_stats = list(
          n_nodes = vcount(g),
          n_edges = ecount(g),
          density = edge_density(g)
        )
      ))
    }
  }
  
  return(NULL)
}

# ==================== 可视化函数 ====================

# 生成富集分析可视化
generate_enrichment_plots <- function(enrichment_results, output_prefix) {
  
  library(enrichplot)
  library(ggplot2)
  
  plots <- list()
  
  for(db_name in names(enrichment_results)) {
    enrich_obj <- enrichment_results[[db_name]]
    
    if(!is.null(enrich_obj) && nrow(enrich_obj@result) > 0) {
      
      # 条形图
      tryCatch({
        p1 <- barplot(enrich_obj, showCategory = 15) +
          ggtitle(paste(db_name, "Enrichment"))
        plots[[paste0(db_name, "_barplot")]] <- p1
      }, error = function(e) {
        message(paste("生成", db_name, "条形图失败"))
      })
      
      # 点图
      tryCatch({
        p2 <- dotplot(enrich_obj, showCategory = 15) +
          ggtitle(paste(db_name, "Enrichment"))
        plots[[paste0(db_name, "_dotplot")]] <- p2
      }, error = function(e) {
        message(paste("生成", db_name, "点图失败"))
      })
      
      # 网络图（如果条目不太多）
      if(nrow(enrich_obj@result) <= 30) {
        tryCatch({
          p3 <- emapplot(enrich_obj, showCategory = 20)
          plots[[paste0(db_name, "_network")]] <- p3
        }, error = function(e) {
          message(paste("生成", db_name, "网络图失败"))
        })
      }
    }
  }
  
  # 保存图形
  if(length(plots) > 0) {
    pdf(paste0(output_prefix, "_enrichment_plots.pdf"), width = 12, height = 8)
    for(plot_name in names(plots)) {
      tryCatch({
        print(plots[[plot_name]])
      }, error = function(e) {
        message(paste("打印图形失败:", plot_name))
      })
    }
    dev.off()
  }
  
  return(plots)
}

# 生成TF调控网络可视化
generate_tf_network_plot <- function(tf_network, output_file) {
  
  if(is.null(tf_network) || is.null(tf_network$graph)) {
    return(NULL)
  }
  
  library(ggraph)
  library(tidygraph)
  
  g <- tf_network$graph
  
  # 设置节点颜色
  V(g)$color <- ifelse(V(g)$node_type == "TF", "red", "lightblue")
  V(g)$size <- ifelse(V(g)$node_type == "TF", 
                     V(g)$degree * 2 + 5, 
                     V(g)$degree + 3)
  
  # 创建ggraph图形
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(arrow = arrow(length = unit(2, 'mm')), 
                   alpha = 0.6, color = "gray50") +
    geom_node_point(aes(color = node_type, size = degree)) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    scale_color_manual(values = c("TF" = "red", "Target" = "lightblue")) +
    theme_graph() +
    labs(title = "Transcription Factor Regulatory Network",
         subtitle = paste("Nodes:", vcount(g), "Edges:", ecount(g))) +
    theme(legend.position = "bottom")
  
  # 保存图形
  ggsave(output_file, p, width = 12, height = 10, dpi = 300)
  
  return(p)
}

# ==================== 主分析函数 ====================

# 单个癌种-细胞类型的功能网络分析
analyze_single_cancer_cell_function <- function(cancer_type, cell_type) {
  
  message(paste("开始功能网络分析:", cancer_type, "-", cell_type))
  
  # 1. 加载轨迹分析结果
  results_file <- file.path(trajectory_output_dir, "monocle3_results",
                           paste0(cancer_type, "_", cell_type, "_trajectory_results.rds"))
  
  if(!file.exists(results_file)) {
    stop(paste("结果文件不存在:", results_file))
  }
  
  trajectory_results <- readRDS(results_file)
  
  # 2. 加载轨迹相关基因
  gene_file <- file.path(trajectory_output_dir, "gene_analysis",
                        paste0(cancer_type, "_", cell_type, "_trajectory_genes.csv"))
  
  if(!file.exists(gene_file)) {
    stop(paste("基因文件不存在:", gene_file))
  }
  
  trajectory_genes <- read.csv(gene_file, stringsAsFactors = FALSE)
  
  # 3. 提取pseudotime
  if(is.null(trajectory_results$cds) || 
     !"pseudotime" %in% colnames(colData(trajectory_results$cds))) {
    stop("Pseudotime数据不可用")
  }
  
  pseudotime_values <- colData(trajectory_results$cds)$pseudotime
  names(pseudotime_values) <- colnames(trajectory_results$cds)
  
  # 4. 基因富集分析
  message("执行基因富集分析...")
  
  # 提取显著基因
  significant_genes <- trajectory_genes[trajectory_genes$q_value < 0.05, ]
  gene_symbols <- convert_gene_ids(significant_genes$gene_short_name)
  
  enrichment_results <- NULL
  if(length(gene_symbols) >= 5) {
    enrichment_results <- comprehensive_gene_enrichment(gene_symbols)
    
    # 方向性富集分析
    directional_enrichment <- trajectory_directional_enrichment(trajectory_genes)
  } else {
    message("显著基因数量不足，跳过富集分析")
    directional_enrichment <- NULL
  }
  
  # 5. 转录因子调控网络分析
  message("执行TF调控网络分析...")
  tf_results <- tf_target_enrichment_analysis(trajectory_genes)
  
  # 6. 代谢轨迹分析
  message("执行代谢轨迹分析...")
  metabolic_results <- metabolic_trajectory_analysis(
    trajectory_results$seurat_obj, 
    pseudotime_values
  )
  
  # 7. 构建整合网络
  message("构建整合网络...")
  integrated_network <- build_integrated_network(
    enrichment_results, 
    tf_results, 
    metabolic_results,
    trajectory_genes
  )
  
  # 8. 生成可视化
  message("生成可视化...")
  
  output_prefix <- file.path(functional_output_dir, "enrichment_analysis",
                            paste0(cancer_type, "_", cell_type))
  
  # 富集分析图形
  if(!is.null(enrichment_results)) {
    enrichment_plots <- generate_enrichment_plots(enrichment_results, output_prefix)
  }
  
  # TF网络图形
  if(!is.null(tf_results) && "regulatory_network" %in% names(tf_results)) {
    tf_plot_file <- file.path(functional_output_dir, "network_plots",
                             paste0(cancer_type, "_", cell_type, "_tf_network.pdf"))
    tf_network_plot <- generate_tf_network_plot(tf_results$regulatory_network, tf_plot_file)
  }
  
  # 9. 整合结果
  final_results <- list(
    cancer_type = cancer_type,
    cell_type = cell_type,
    trajectory_genes = trajectory_genes,
    enrichment_analysis = enrichment_results,
    directional_enrichment = directional_enrichment,
    tf_analysis = tf_results,
    metabolic_analysis = metabolic_results,
    integrated_network = integrated_network,
    analysis_summary = list(
      n_trajectory_genes = nrow(trajectory_genes),
      n_significant_genes = sum(trajectory_genes$q_value < 0.05),
      n_enriched_terms = if(!is.null(enrichment_results)) {
        sum(sapply(enrichment_results, function(x) nrow(x@result)))
      } else 0,
      n_regulatory_tfs = if(!is.null(tf_results) && "tf_enrichment" %in% names(tf_results)) {
        nrow(tf_results$tf_enrichment)
      } else 0
    )
  )
  
  # 10. 保存结果
  results_file <- file.path(functional_output_dir, 
                           paste0(cancer_type, "_", cell_type, "_functional_analysis.rds"))
  saveRDS(final_results, results_file)
  
  # 生成报告
  generate_functional_analysis_report(final_results)
  
  message(paste("功能网络分析完成:", cancer_type, "-", cell_type))
  
  return(final_results)
}

# 批量功能网络分析
batch_functional_analysis <- function() {
  
  message("开始批量功能网络分析...")
  
  # 检测可用的轨迹分析结果
  results_dir <- file.path(trajectory_output_dir, "monocle3_results")
  result_files <- list.files(results_dir, pattern = "_trajectory_results\\.rds$")
  
  if(length(result_files) == 0) {
    stop("未找到轨迹分析结果文件")
  }
  
  # 解析文件名
  all_results <- list()
  
  for(file in result_files) {
    result_name <- gsub("_trajectory_results\\.rds$", "", file)
    parts <- strsplit(result_name, "_")[[1]]
    
    if(length(parts) >= 2) {
      cancer_type <- parts[1]
      cell_type <- paste(parts[-1], collapse = "_")
      
      message(paste("\n处理:", cancer_type, "-", cell_type))
      
      tryCatch({
        result <- analyze_single_cancer_cell_function(cancer_type, cell_type)
        all_results[[result_name]] <- result
        
      }, error = function(e) {
        message(paste("ERROR:", cancer_type, "-", cell_type, ":", e$message))
      })
    }
  }
  
  # 生成汇总分析
  if(length(all_results) > 0) {
    generate_batch_summary(all_results)
  }
  
  message("批量功能网络分析完成!")
  return(all_results)
}

# 生成功能分析报告
generate_functional_analysis_report <- function(results) {
  
  report_file <- file.path(functional_output_dir,
                          paste0(results$cancer_type, "_", results$cell_type, 
                                "_functional_report.txt"))
  
  sink(report_file)
  
  cat("=== 功能网络分析报告 ===\n")
  cat("癌种:", results$cancer_type, "\n")
  cat("细胞类型:", results$cell_type, "\n")
  cat("分析时间:", as.character(Sys.time()), "\n\n")
  
  # 基本统计
  cat("1. 基本统计\n")
  cat("===========\n")
  summary <- results$analysis_summary
  cat("轨迹相关基因总数:", summary$n_trajectory_genes, "\n")
  cat("显著基因数:", summary$n_significant_genes, "\n")
  cat("富集条目数:", summary$n_enriched_terms, "\n")
  cat("调控转录因子数:", summary$n_regulatory_tfs, "\n")
  
  # 富集分析结果
  if(!is.null(results$enrichment_analysis)) {
    cat("\n2. 功能富集分析\n")
    cat("================\n")
    
    for(db_name in names(results$enrichment_analysis)) {
      enrich_obj <- results$enrichment_analysis[[db_name]]
      if(!is.null(enrich_obj) && nrow(enrich_obj@result) > 0) {
        significant_terms <- enrich_obj@result[enrich_obj@result$p.adjust < 0.05, ]
        cat(sprintf("%s: %d个显著条目\n", db_name, nrow(significant_terms)))
        
        # 显示top 5条目
        if(nrow(significant_terms) > 0) {
          top_terms <- head(significant_terms, 5)
          for(i in 1:nrow(top_terms)) {
            cat(sprintf("  %d. %s (p=%.2e)\n", 
                        i, top_terms$Description[i], top_terms$p.adjust[i]))
          }
        }
        cat("\n")
      }
    }
  }
  
  # TF调控分析
  if(!is.null(results$tf_analysis) && "tf_enrichment" %in% names(results$tf_analysis)) {
    cat("3. 转录因子调控分析\n")
    cat("===================\n")
    tf_data <- results$tf_analysis$tf_enrichment
    cat("显著调控TF数:", nrow(tf_data), "\n")
    
    if(nrow(tf_data) > 0) {
      cat("Top 10 调控TF:\n")
      top_tfs <- head(tf_data, 10)
      for(i in 1:nrow(top_tfs)) {
        cat(sprintf("  %d. %s (靶基因数: %d, p=%.2e)\n",
                    i, top_tfs$TF[i], top_tfs$n_targets_in_trajectory[i], 
                    top_tfs$q_value[i]))
      }
    }
  }
  
  # 代谢分析
  if(!is.null(results$metabolic_analysis)) {
    cat("\n4. 代谢轨迹分析\n")
    cat("================\n")
    
    if("significant_pathways" %in% names(results$metabolic_analysis)) {
      metab_data <- results$metabolic_analysis$significant_pathways
      cat("显著代谢通路数:", nrow(metab_data), "\n")
      
      if(nrow(metab_data) > 0) {
        cat("显著代谢通路:\n")
        for(i in 1:nrow(metab_data)) {
          cat(sprintf("  %d. %s (相关性: %.3f, p=%.2e)\n",
                      i, metab_data$pathway[i], metab_data$correlation[i], 
                      metab_data$q_value[i]))
        }
      }
    } else if("pathway_correlations" %in% names(results$metabolic_analysis)) {
      metab_data <- results$metabolic_analysis$pathway_correlations
      significant_metab <- metab_data[metab_data$q_value < 0.05, ]
      cat("显著代谢通路数:", nrow(significant_metab), "\n")
      
      if(nrow(significant_metab) > 0) {
        for(i in 1:nrow(significant_metab)) {
          cat(sprintf("  %d. %s (相关性: %.3f, p=%.2e)\n",
                      i, significant_metab$pathway[i], 
                      significant_metab$correlation[i], 
                      significant_metab$q_value[i]))
        }
      }
    }
  }
  
  cat("\n=== 报告结束 ===\n")
  sink()
  
  message(paste("功能分析报告已保存:", report_file))
}

# 辅助函数：基因ID转换
convert_gene_ids <- function(gene_ids) {
  library(org.Hs.eg.db)
  
  # 首先尝试作为ENSEMBL ID转换
  tryCatch({
    symbols <- mapIds(org.Hs.eg.db, 
                     keys = gene_ids,
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
    
    valid_symbols <- symbols[!is.na(symbols)]
    return(valid_symbols)
    
  }, error = function(e) {
    # 如果失败，假设已经是SYMBOL
    return(gene_ids[gene_ids != "" & !is.na(gene_ids)])
  })
}

# ==================== 执行分析 ====================

if(!interactive()) {
  message("开始执行功能网络分析...")
  
  # 执行批量分析
  functional_results <- batch_functional_analysis()
  
  message("功能网络分析完成!")
  
} else {
  message("功能网络分析代码加载完成!")
  message("使用方法:")
  message("1. 单个分析: analyze_single_cancer_cell_function('BC', 'T_cells')")
  message("2. 批量分析: batch_functional_analysis()")
}

message("轨迹基因功能网络分析代码加载完成!")