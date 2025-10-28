library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(scales)
library(patchwork)

# Fine-grained Cell Type Annotation Quality Assessment and Visualization Analysis
# =============================================================================

# Set paths
analysis_dir <- "/mnt/public7/pancancercol/hefeng/cluster-annotation"
output_viz_dir <- file.path(analysis_dir, "visualization_analysis")
dir.create(output_viz_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Data loading function (supports both file types)
load_annotation_results <- function(analysis_dir) {
  message("Loading annotation results...")
  
  # Get all cancer type directories
  cancer_dirs <- list.dirs(analysis_dir, full.names = FALSE, recursive = FALSE)
  cancer_dirs <- cancer_dirs[cancer_dirs != "visualization_analysis"]
  
  results_list <- list()
  detailed_list <- list()
  summary_stats <- list()
  
  for (cancer_type in cancer_dirs) {
    cancer_dir <- file.path(analysis_dir, cancer_type)
    
    # Main file
    main_file <- file.path(cancer_dir, paste0(cancer_type, "_fine_annotated.rds"))
    # Detailed files
    detailed_files <- list.files(cancer_dir, pattern = paste0(cancer_type, "_.*_detailed\\.rds$"), full.names = TRUE)
    
    if (file.exists(main_file)) {
      message("  Loading ", cancer_type, " main file...")
      
      tryCatch({
        seurat_obj <- readRDS(main_file)
        results_list[[cancer_type]] <- seurat_obj
        
        # Collect basic statistics
        summary_stats[[cancer_type]] <- list(
          n_cells = ncol(seurat_obj),
          n_samples = length(unique(seurat_obj$sample_id)),
          n_fine_types = length(unique(seurat_obj$fine_cell_type)),
          n_major_types = length(unique(seurat_obj$major_cell_type)),
          cancer_general = unique(seurat_obj$Cancer_general)[1]
        )
        
        # Load detailed files (if exist)
        if (length(detailed_files) > 0) {
          message("    Found ", length(detailed_files), " detailed files")
          detailed_list[[cancer_type]] <- list()
          
          for (detailed_file in detailed_files) {
            file_name <- basename(detailed_file)
            # Extract file type identifier
            file_id <- gsub(paste0(cancer_type, "_(.*)_detailed\\.rds"), "\\1", file_name)
            
            tryCatch({
              detailed_obj <- readRDS(detailed_file)
              detailed_list[[cancer_type]][[file_id]] <- detailed_obj
              message("      Loaded detailed file: ", file_id)
            }, error = function(e) {
              warning("Failed to load detailed file ", detailed_file, ": ", e$message)
            })
          }
        }
        
      }, error = function(e) {
        warning("Failed to load ", cancer_type, ": ", e$message)
      })
    }
  }
  
  message("Successfully loaded ", length(results_list), " cancer types")
  message("Among them, ", length(detailed_list), " cancer types have detailed analysis files")
  
  return(list(
    main_data = results_list, 
    detailed_data = detailed_list,
    summary = summary_stats
  ))
}

# 2. Quality assessment based on main files
assess_main_annotation_quality <- function(seurat_obj, cancer_type) {
  message("Assessing ", cancer_type, " main annotation quality...")
  
  quality_metrics <- list()
  
  # 1. Cell type distribution assessment
  cell_type_counts <- table(seurat_obj$fine_cell_type)
  quality_metrics$cell_type_distribution <- list(
    n_types = length(cell_type_counts),
    gini_coefficient = calculate_gini_coefficient(cell_type_counts),
    min_cells_per_type = min(cell_type_counts),
    max_cells_per_type = max(cell_type_counts),
    median_cells_per_type = median(cell_type_counts)
  )
  
  # 2. Batch effect assessment (based on samples)
  if ("sample_id" %in% colnames(seurat_obj@meta.data)) {
    quality_metrics$batch_effect <- assess_batch_effect(seurat_obj)
  }
  
  # 3. Major vs fine type consistency
  quality_metrics$type_consistency <- assess_type_consistency(seurat_obj)
  
  # 4. Sample coverage assessment
  quality_metrics$sample_coverage <- assess_sample_coverage(seurat_obj)
  
  # 5. SingleR label consistency (if exists)
  if ("singleR_labels" %in% colnames(seurat_obj@meta.data)) {
    quality_metrics$singleR_consistency <- assess_singleR_consistency(seurat_obj)
  }
  
  return(quality_metrics)
}

# Helper assessment functions
calculate_gini_coefficient <- function(counts) {
  n <- length(counts)
  if (n <= 1) return(0)
  
  sorted_counts <- sort(counts)
  index <- 1:n
  gini <- (2 * sum(index * sorted_counts)) / (n * sum(sorted_counts)) - (n + 1) / n
  return(gini)
}

assess_batch_effect <- function(seurat_obj) {
  sample_celltype_table <- table(seurat_obj$sample_id, seurat_obj$fine_cell_type)
  
  if (min(dim(sample_celltype_table)) > 1 && sum(sample_celltype_table) > 0) {
    tryCatch({
      chisq_test <- chisq.test(sample_celltype_table)
      batch_effect_pvalue <- chisq_test$p.value
    }, error = function(e) {
      batch_effect_pvalue <- NA
    })
  } else {
    batch_effect_pvalue <- NA
  }
  
  celltype_sample_distribution <- apply(sample_celltype_table, 2, function(x) {
    if (sum(x) == 0) return(0)
    props <- x / sum(x)
    expected_prop <- 1 / length(x)
    variance <- sum((props - expected_prop)^2)
    return(variance)
  })
  
  return(list(
    batch_effect_pvalue = batch_effect_pvalue,
    mean_distribution_variance = mean(celltype_sample_distribution, na.rm = TRUE),
    sample_celltype_table = sample_celltype_table
  ))
}

assess_type_consistency <- function(seurat_obj) {
  major_fine_mapping <- table(seurat_obj$major_cell_type, seurat_obj$fine_cell_type)
  
  major_type_purity <- apply(major_fine_mapping, 1, function(x) {
    if (sum(x) == 0) return(0)
    max(x) / sum(x)
  })
  
  fine_type_consistency <- apply(major_fine_mapping, 2, function(x) {
    if (sum(x) == 0) return(0)
    max(x) / sum(x)
  })
  
  return(list(
    major_type_purity = major_type_purity,
    fine_type_consistency = fine_type_consistency,
    mean_major_purity = mean(major_type_purity, na.rm = TRUE),
    mean_fine_consistency = mean(fine_type_consistency, na.rm = TRUE),
    mapping_table = major_fine_mapping
  ))
}

assess_sample_coverage <- function(seurat_obj) {
  sample_stats <- seurat_obj@meta.data %>%
    group_by(sample_id) %>%
    summarise(
      n_cells = n(),
      n_major_types = length(unique(major_cell_type)),
      n_fine_types = length(unique(fine_cell_type)),
      .groups = 'drop'
    )
  
  return(list(
    n_samples = nrow(sample_stats),
    mean_cells_per_sample = mean(sample_stats$n_cells),
    mean_major_types_per_sample = mean(sample_stats$n_major_types),
    mean_fine_types_per_sample = mean(sample_stats$n_fine_types),
    sample_stats = sample_stats
  ))
}

assess_singleR_consistency <- function(seurat_obj) {
  # Assess consistency between singleR labels and final annotations
  consistency_table <- table(seurat_obj$singleR_labels, seurat_obj$major_cell_type)
  
  # Calculate consistency score
  total_cells <- sum(consistency_table)
  max_consistency <- sum(apply(consistency_table, 1, max))
  consistency_score <- max_consistency / total_cells
  
  return(list(
    consistency_score = consistency_score,
    consistency_table = consistency_table
  ))
}

# 3. Comprehensive statistics analysis function (simplified version)
generate_comprehensive_statistics <- function(results_list) {
  message("Generating comprehensive statistics...")
  
  # Statistics based on main files
  all_fine_types <- c()
  all_major_types <- c()
  cancer_celltype_matrix <- list()
  cancer_major_matrix <- list()
  
  for (cancer_type in names(results_list$main_data)) {
    seurat_obj <- results_list$main_data[[cancer_type]]
    
    # Fine type statistics
    fine_types <- table(seurat_obj$fine_cell_type)
    all_fine_types <- union(all_fine_types, names(fine_types))
    cancer_celltype_matrix[[cancer_type]] <- fine_types
    
    # Major type statistics
    major_types <- table(seurat_obj$major_cell_type)
    all_major_types <- union(all_major_types, names(major_types))
    cancer_major_matrix[[cancer_type]] <- major_types
  }
  
  # Create matrices
  fine_celltype_matrix <- matrix(0, nrow = length(names(results_list$main_data)), 
                                ncol = length(all_fine_types))
  rownames(fine_celltype_matrix) <- names(results_list$main_data)
  colnames(fine_celltype_matrix) <- all_fine_types
  
  major_celltype_matrix <- matrix(0, nrow = length(names(results_list$main_data)), 
                                 ncol = length(all_major_types))
  rownames(major_celltype_matrix) <- names(results_list$main_data)
  colnames(major_celltype_matrix) <- all_major_types
  
  # Fill matrices
  for (cancer_type in names(cancer_celltype_matrix)) {
    for (cell_type in names(cancer_celltype_matrix[[cancer_type]])) {
      fine_celltype_matrix[cancer_type, cell_type] <- cancer_celltype_matrix[[cancer_type]][cell_type]
    }
  }
  
  for (cancer_type in names(cancer_major_matrix)) {
    for (cell_type in names(cancer_major_matrix[[cancer_type]])) {
      major_celltype_matrix[cancer_type, cell_type] <- cancer_major_matrix[[cancer_type]][cell_type]
    }
  }
  
  # Quality assessment summary
  main_quality_summary <- list()
  
  for (cancer_type in names(results_list$main_data)) {
    # Main file quality assessment
    main_quality_summary[[cancer_type]] <- assess_main_annotation_quality(
      results_list$main_data[[cancer_type]], cancer_type
    )
  }
  
  return(list(
    fine_celltype_matrix = fine_celltype_matrix,
    major_celltype_matrix = major_celltype_matrix,
    main_quality_summary = main_quality_summary,
    all_fine_types = all_fine_types,
    all_major_types = all_major_types
  ))
}

# 4. Modified overview plots function
create_overview_plots <- function(results_list, stats, output_dir) {
  message("Creating overview plots...")
  
  summary_df <- do.call(rbind, lapply(names(results_list$summary), function(cancer) {
    data.frame(
      Cancer_Type = cancer,
      N_Cells = results_list$summary[[cancer]]$n_cells,
      N_Samples = results_list$summary[[cancer]]$n_samples,
      N_Fine_Types = results_list$summary[[cancer]]$n_fine_types,
      N_Major_Types = results_list$summary[[cancer]]$n_major_types,
      Has_Detailed = cancer %in% names(results_list$detailed_data)
    )
  }))
  
  # Cell count plot
  p1 <- ggplot(summary_df, aes(x = reorder(Cancer_Type, N_Cells), y = N_Cells, fill = Has_Detailed)) +
    geom_col(alpha = 0.7) +
    geom_text(aes(label = scales::comma(N_Cells)), hjust = -0.1, size = 3) +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "lightgray", "TRUE" = "steelblue"),
                     name = "Has Detailed Analysis") +
    labs(title = "Number of Cells per Cancer Type", x = "Cancer Type", y = "Number of Cells") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # Cell type count plot
  p2 <- ggplot(summary_df, aes(x = reorder(Cancer_Type, N_Fine_Types), y = N_Fine_Types)) +
    geom_col(fill = "coral", alpha = 0.7) +
    geom_text(aes(label = N_Fine_Types), hjust = -0.1, size = 3) +
    coord_flip() +
    labs(title = "Number of Fine Cell Types per Cancer Type", x = "Cancer Type", y = "Number of Fine Cell Types") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # Analysis method coverage plot (using base R instead of pivot_longer)
  main_only_count <- sum(!summary_df$Has_Detailed)
  with_detailed_count <- sum(summary_df$Has_Detailed)
  
  method_coverage <- data.frame(
    Analysis_Type = c("Main Only", "With Detailed"),
    Count = c(main_only_count, with_detailed_count)
  )
  
  p3 <- ggplot(method_coverage, aes(x = Analysis_Type, y = Count, fill = Analysis_Type)) +
    geom_col(alpha = 0.7) +
    geom_text(aes(label = Count), vjust = -0.3, size = 4) +
    scale_fill_manual(values = c("Main Only" = "lightcoral", "With Detailed" = "darkgreen")) +
    labs(title = "Analysis Method Coverage", x = "Analysis Type", y = "Number of Cancer Types") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "none")
  
  # Combined plot
  overview_plot <- (p1 / p2 / p3) + 
    plot_annotation(title = "Cancer Type Annotation Results Overview", 
                   theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  ggsave(file.path(output_dir, "01_overview_statistics.png"), 
         overview_plot, width = 12, height = 15, dpi = 300)
  
  return(summary_df)
}

# Simplified heatmap function
create_celltype_heatmap <- function(celltype_matrix, output_dir, matrix_type = "fine") {
  message("Creating ", matrix_type, " cell type distribution heatmap...")
  
  # Convert to proportions
  celltype_prop_matrix <- celltype_matrix / rowSums(celltype_matrix)
  celltype_prop_matrix[is.na(celltype_prop_matrix)] <- 0
  
  # Filtering strategy
  if (matrix_type == "fine") {
    celltype_filter <- apply(celltype_matrix, 2, function(x) {
      sum(x > 0) >= 2 | sum(x) > 50
    })
  } else {
    celltype_filter <- apply(celltype_matrix, 2, function(x) {
      sum(x > 0) >= 1
    })
  }
  
  filtered_matrix <- celltype_prop_matrix[, celltype_filter, drop = FALSE]
  
  if (ncol(filtered_matrix) == 0) {
    message("No cell types meet the criteria")
    return(NULL)
  }
  
  # Use basic heatmap
  tryCatch({
    png(file.path(output_dir, paste0("02_", matrix_type, "_celltype_heatmap.png")), 
        width = max(1200, ncol(filtered_matrix) * 50), 
        height = max(800, nrow(filtered_matrix) * 60),
        res = 100)
    
    par(mar = c(10, 6, 4, 2))
    
    heatmap(as.matrix(filtered_matrix),
            col = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
            main = paste0(toupper(substring(matrix_type, 1, 1)), substring(matrix_type, 2), 
                         " Cell Type Distribution Heatmap"),
            xlab = "Cell Types",
            ylab = "Cancer Types",
            cexRow = 0.8,
            cexCol = 0.6,
            margins = c(10, 6))
    
    dev.off()
    message(matrix_type, " heatmap created successfully")
    
  }, error = function(e) {
    dev.off()
    message(matrix_type, " heatmap creation failed: ", e$message)
  })
  
  return(filtered_matrix)
}

# Quality assessment plots (based on main files)
create_quality_assessment_plots <- function(main_quality_summary, output_dir) {
  message("Creating quality assessment plots...")
  
  quality_df <- do.call(rbind, lapply(names(main_quality_summary), function(cancer) {
    qs <- main_quality_summary[[cancer]]
    
    data.frame(
      Cancer_Type = cancer,
      Gini_Coefficient = qs$cell_type_distribution$gini_coefficient,
      N_Cell_Types = qs$cell_type_distribution$n_types,
      Mean_Major_Purity = ifelse(is.null(qs$type_consistency), NA, 
                                qs$type_consistency$mean_major_purity),
      N_Samples = ifelse(is.null(qs$sample_coverage), NA,
                        qs$sample_coverage$n_samples),
      SingleR_Consistency = ifelse(is.null(qs$singleR_consistency), NA,
                                  qs$singleR_consistency$consistency_score)
    )
  }))
  
  # Gini coefficient plot
  p1 <- ggplot(quality_df, aes(x = reorder(Cancer_Type, Gini_Coefficient), y = Gini_Coefficient)) +
    geom_col(fill = "orange", alpha = 0.7) +
    geom_text(aes(label = round(Gini_Coefficient, 3)), hjust = -0.1, size = 3) +
    coord_flip() +
    labs(title = "Cell Type Distribution Uniformity (Gini Coefficient)", 
         subtitle = "Lower values indicate more uniform distribution",
         x = "Cancer Type", y = "Gini Coefficient") +
    theme_minimal()
  
  # Type purity plot
  if (any(!is.na(quality_df$Mean_Major_Purity))) {
    p2 <- ggplot(quality_df[!is.na(quality_df$Mean_Major_Purity), ], 
                aes(x = reorder(Cancer_Type, Mean_Major_Purity), y = Mean_Major_Purity)) +
      geom_col(fill = "purple", alpha = 0.7) +
      geom_text(aes(label = round(Mean_Major_Purity, 3)), hjust = -0.1, size = 3) +
      coord_flip() +
      labs(title = "Major Cell Type Purity", 
           subtitle = "Higher values indicate better consistency within major types",
           x = "Cancer Type", y = "Average Purity") +
      theme_minimal()
  } else {
    p2 <- ggplot() + ggtitle("Type Purity Data Not Available") + theme_void()
  }
  
  # SingleR consistency plot
  if (any(!is.na(quality_df$SingleR_Consistency))) {
    p3 <- ggplot(quality_df[!is.na(quality_df$SingleR_Consistency), ], 
                aes(x = reorder(Cancer_Type, SingleR_Consistency), y = SingleR_Consistency)) +
      geom_col(fill = "darkgreen", alpha = 0.7) +
      geom_text(aes(label = round(SingleR_Consistency, 3)), hjust = -0.1, size = 3) +
      coord_flip() +
      labs(title = "SingleR Label Consistency", 
           subtitle = "Higher values indicate better agreement with SingleR annotations",
           x = "Cancer Type", y = "Consistency Score") +
      theme_minimal()
  } else {
    p3 <- ggplot() + ggtitle("SingleR Consistency Data Not Available") + theme_void()
  }
  
  quality_plot <- (p1 / p2 / p3) + 
    plot_annotation(title = "Annotation Quality Assessment", 
                   theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  ggsave(file.path(output_dir, "03_quality_assessment.png"), 
         quality_plot, width = 14, height = 18, dpi = 300)
  
  return(quality_df)
}

# Individual cancer type plots
create_individual_cancer_plots <- function(results_list, output_dir) {
  message("Creating detailed analysis plots for each cancer type...")
  
  individual_dir <- file.path(output_dir, "individual_cancers")
  dir.create(individual_dir, showWarnings = FALSE)
  
  for (cancer_type in names(results_list$main_data)) {
    message("  Creating detailed plots for ", cancer_type, "...")
    
    seurat_obj <- results_list$main_data[[cancer_type]]
    cancer_output_dir <- file.path(individual_dir, cancer_type)
    dir.create(cancer_output_dir, showWarnings = FALSE)
    
    # 1. Cell type proportion plots
    create_celltype_proportion_plots(seurat_obj, cancer_type, cancer_output_dir)
    
    # 2. Sample comparison plots
    if (length(unique(seurat_obj$sample_id)) > 1) {
      create_sample_comparison_plots(seurat_obj, cancer_type, cancer_output_dir)
    }
    
    # 3. Major vs fine type relationship plot
    create_major_fine_relationship_plot(seurat_obj, cancer_type, cancer_output_dir)
  }
}

create_celltype_proportion_plots <- function(seurat_obj, cancer_type, output_dir) {
  # Cell type counts
  celltype_counts <- table(seurat_obj$fine_cell_type)
  celltype_df <- data.frame(
    CellType = names(celltype_counts),
    Count = as.numeric(celltype_counts),
    Proportion = as.numeric(celltype_counts) / ncol(seurat_obj)
  )
  celltype_df <- celltype_df[order(celltype_df$Count, decreasing = TRUE), ]
  
  # Bar plot
  p1 <- ggplot(celltype_df, aes(x = reorder(CellType, Count), y = Count)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = Count), hjust = -0.1, size = 3) +
    coord_flip() +
    labs(title = paste0(cancer_type, " - Cell Type Counts"), x = "Cell Type", y = "Number of Cells") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  # Major type distribution
  major_counts <- table(seurat_obj$major_cell_type)
  major_df <- data.frame(
    CellType = names(major_counts),
    Count = as.numeric(major_counts),
    Proportion = as.numeric(major_counts) / ncol(seurat_obj)
  )
  
  p2 <- ggplot(major_df, aes(x = "", y = Count, fill = CellType)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    ggtitle(paste0(cancer_type, " - Major Cell Type Proportions")) +
    theme(legend.position = "right")
  
  # Combined plot
  proportion_combined <- p1 | p2
  
  ggsave(file.path(output_dir, paste0(cancer_type, "_celltype_proportions.png")), 
         proportion_combined, width = 16, height = 8, dpi = 300)
}

create_sample_comparison_plots <- function(seurat_obj, cancer_type, output_dir) {
  # Sample-cell type cross table
  sample_celltype_table <- table(seurat_obj$sample_id, seurat_obj$fine_cell_type)
  sample_celltype_prop <- prop.table(sample_celltype_table, margin = 1)
  
  # Convert to long format
  sample_celltype_df <- as.data.frame(sample_celltype_prop)
  colnames(sample_celltype_df) <- c("Sample", "CellType", "Proportion")
  
  # Stacked bar plot
  p1 <- ggplot(sample_celltype_df, aes(x = Sample, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity") +
    labs(title = paste0(cancer_type, " - Cell Type Proportions Across Samples"), 
         x = "Sample", y = "Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")  # Remove legend to save space
  
  ggsave(file.path(output_dir, paste0(cancer_type, "_sample_comparison.png")), 
         p1, width = 12, height = 8, dpi = 300)
}

create_major_fine_relationship_plot <- function(seurat_obj, cancer_type, output_dir) {
  # Major type vs fine type relationship
  relationship_table <- table(seurat_obj$major_cell_type, seurat_obj$fine_cell_type)
  
  # Convert to long format
  relationship_df <- as.data.frame(relationship_table)
  colnames(relationship_df) <- c("Major_Type", "Fine_Type", "Count")
  relationship_df <- relationship_df[relationship_df$Count > 0, ]
  
  p1 <- ggplot(relationship_df, aes(x = Major_Type, y = Count, fill = Fine_Type)) +
    geom_bar(stat = "identity") +
    labs(title = paste0(cancer_type, " - Major vs Fine Cell Type Relationship"), 
         x = "Major Cell Type", y = "Number of Cells") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          legend.text = element_text(size = 6))
  
  ggsave(file.path(output_dir, paste0(cancer_type, "_major_fine_relationship.png")), 
         p1, width = 14, height = 8, dpi = 300)
}

# Report generation function
generate_analysis_report <- function(results_list, stats, quality_df, output_dir) {
  message("Generating analysis report...")
  
  # Count cancer types with detailed analysis
  n_detailed <- length(results_list$detailed_data)
  
  report_content <- c(
    "# Fine-grained Cell Type Annotation Analysis Report",
    "",
    paste0("Generated on: ", Sys.time()),
    "",
    "## 1. Data Overview",
    "",
    paste0("- Number of cancer types analyzed: ", length(results_list$main_data)),
    paste0("- Number of cancer types with detailed analysis: ", n_detailed),
    paste0("- Total number of cells: ", sum(sapply(results_list$summary, function(x) x$n_cells))),
    paste0("- Total number of samples: ", sum(sapply(results_list$summary, function(x) x$n_samples))),
    paste0("- Total number of fine cell types: ", length(stats$all_fine_types)),
    paste0("- Total number of major cell types: ", length(stats$all_major_types)),
    "",
    "## 2. Analysis Methods",
    "",
    "- Main file analysis: Based on merged annotation results",
    "- Detailed file analysis: Contains Monaco, Blueprint, UMAP, Harmony detailed results",
    "",
    "## 3. Output Files Description",
    "",
    "- `01_overview_statistics.png`: Overview statistics plots",
    "- `02_fine_celltype_heatmap.png`: Fine cell type distribution heatmap", 
    "- `02_major_celltype_heatmap.png`: Major cell type distribution heatmap",
    "- `03_quality_assessment.png`: Quality assessment plots",
    "- `individual_cancers/`: Detailed analysis plots for each cancer type",
    "- `summary_statistics.csv`: Summary statistics table",
    "- `quality_metrics.csv`: Quality metrics table",
    "- `fine_celltype_matrix.csv`: Fine cell type count matrix",
    "- `major_celltype_matrix.csv`: Major cell type count matrix",
    ""
  )
  
  # Add detailed information for each cancer type
  for (cancer in names(results_list$summary)) {
    s <- results_list$summary[[cancer]]
    report_content <- c(report_content,
      paste0("### ", cancer),
      paste0("- Number of cells: ", s$n_cells),
      paste0("- Number of samples: ", s$n_samples),
      paste0("- Number of fine cell types: ", s$n_fine_types),
      paste0("- Number of major cell types: ", s$n_major_types),
      ""
    )
  }
  
  # Save report
  writeLines(report_content, file.path(output_dir, "ANALYSIS_REPORT.md"))
  
  # Save data tables
  write.csv(do.call(rbind, lapply(names(results_list$summary), function(cancer) {
    s <- results_list$summary[[cancer]]
    data.frame(Cancer_Type = cancer, 
               N_Cells = s$n_cells,
               N_Samples = s$n_samples, 
               N_Fine_Types = s$n_fine_types,
               N_Major_Types = s$n_major_types,
               Has_Detailed = cancer %in% names(results_list$detailed_data))
  })), file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
  
  write.csv(quality_df, file.path(output_dir, "quality_metrics.csv"), row.names = FALSE)
  write.csv(stats$fine_celltype_matrix, file.path(output_dir, "fine_celltype_matrix.csv"))
  write.csv(stats$major_celltype_matrix, file.path(output_dir, "major_celltype_matrix.csv"))
}

# Main execution function (simplified version)
main_analysis <- function() {
  message("Starting comprehensive analysis and visualization...")
  
  # 1. Load data (main files)
  results_list <- load_annotation_results(analysis_dir)
  
  if (length(results_list$main_data) == 0) {
    stop("No valid annotation results found")
  }
  
  # 2. Generate comprehensive statistics
  stats <- generate_comprehensive_statistics(results_list)
  
  # 3. Create visualization plots
  message("Creating visualization plots...")
  
  # Overview plots
  summary_df <- create_overview_plots(results_list, stats, output_viz_dir)
  
  # Cell type heatmaps
  filtered_fine_matrix <- create_celltype_heatmap(stats$fine_celltype_matrix, output_viz_dir, "fine")
  filtered_major_matrix <- create_celltype_heatmap(stats$major_celltype_matrix, output_viz_dir, "major")
  
  # Quality assessment plots (based on main files)
  quality_df <- create_quality_assessment_plots(stats$main_quality_summary, output_viz_dir)
  
  # Individual cancer type plots
  create_individual_cancer_plots(results_list, output_viz_dir)
  
  # 4. Generate report
  generate_analysis_report(results_list, stats, quality_df, output_viz_dir)
  
  message("Analysis completed! Results saved in: ", output_viz_dir)
  
  return(list(
    results = results_list,
    statistics = stats,
    quality_metrics = quality_df,
    summary = summary_df
  ))
}

# Run analysis
if (!interactive()) {
  analysis_results <- main_analysis()
}