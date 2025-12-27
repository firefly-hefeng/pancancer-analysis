# Single-cell Sequencing Data Analysis Script (R Version)
# ========================================================

# Install and load required packages
required_packages <- c("readxl", "ggplot2", "dplyr", "tidyr", "RColorBrewer", 
                       "pheatmap", "gridExtra", "scales", "writexl")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Set working directory and output folder
output_dir <- "data_distribution"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Set plotting parameters
theme_set(theme_bw(base_size = 12))
options(scipen = 999)  # Disable scientific notation

cat("================================================================================\n")
cat("Single-cell Sequencing Data Analysis Report\n")
cat("================================================================================\n\n")

# Read data
file_path <- '/mnt/public7/pancancercol/hefeng/cluster-annotation/all_data.xlsx'
df <- read_excel(file_path)

# Cancer to tissue mapping dictionary
cancer_to_tissue <- list(
  # Breast
  'BC' = 'Breast',
  'Breast Cancer' = 'Breast',
  
  # Liver
  'LIHC' = 'Liver',
  
  # Pancreas
  'PAAD' = 'Pancreas',
  'PNT' = 'Pancreas',
  'PC' = 'Pancreas',
  
  # Skin
  'Melan' = 'Skin',
  'BCC' = 'Skin',
  'SCC、BCC' = 'Skin',
  'cSCC' = 'Skin',
  'MCC' = 'Skin',
  
  # Prostate
  'PRAD' = 'Prostate',
  'CRPC' = 'Prostate',
  'PCA' = 'Prostate',
  
  # Lung
  'LUSC' = 'Lung',
  'LUAD' = 'Lung',
  'LCLC' = 'Lung',
  'SCLC' = 'Lung',
  'NSCLC' = 'Lung',
  'LC' = 'Lung',
  'LSCC' = 'Lung',
  
  # Ovary
  'OV' = 'Ovary',
  
  # Colorectum
  'COADREAD' = 'Colorectum',
  'CRC' = 'Colorectum',
  'COAD' = 'Colorectum',
  'READ' = 'Colorectum',
  
  # Head and Neck
  'HNSC' = 'Head and Neck',
  'HNSCC' = 'Head and Neck',
  'OCC' = 'Head and Neck',
  
  # Kidney
  'KIRC' = 'Kidney',
  'RCC' = 'Kidney',
  'KIRP' = 'Kidney',
  'WT' = 'Kidney',
  
  # Blood/Bone Marrow
  'MM' = 'Blood/Bone Marrow',
  'RRMM' = 'Blood/Bone Marrow',
  'BCL' = 'Blood/Bone Marrow',
  'NHL' = 'Blood/Bone Marrow',
  
  # Lymphoid
  'CTCL' = 'Lymphoid',
  
  # Bladder
  'BLCA' = 'Bladder',
  'BCa' = 'Bladder',
  
  # Brain
  'GBM' = 'Brain',
  'other_glioma' = 'Brain',
  'GLI' = 'Brain',
  'HGG' = 'Brain',
  'Glioma' = 'Brain',
  'MB' = 'Brain',
  'DMG' = 'Brain',
  'EPN' = 'Brain',
  
  # Neuroblast
  'NB' = 'Neuroblast',
  
  # Bile Duct
  'CHOL' = 'Bile Duct',
  'hilar CCA' = 'Bile Duct',
  'CCA' = 'Bile Duct',
  'BTC' = 'Bile Duct',
  
  # Endometrium/Cervix
  'UCEC' = 'Endometrium',
  'CESC' = 'Cervix',
  
  # Neuroendocrine
  'NEN' = 'Neuroendocrine',
  'NETs' = 'Neuroendocrine',
  
  # Esophagus
  'ESCA' = 'Esophagus',
  'ESCC' = 'Esophagus',
  
  # Thyroid
  'THCA' = 'Thyroid',
  'FTC' = 'Thyroid',
  
  # Sarcoma
  'SARC' = 'Sarcoma',
  'DSRCT' = 'Sarcoma',
  'GIST' = 'Sarcoma',
  'MS' = 'Sarcoma',
  
  # Stomach
  'STAD' = 'Stomach',
  'GAC' = 'Stomach',
  'GC' = 'Stomach',
  'PSC' = 'Stomach',
  
  # Liver Metastasis
  'LM' = 'Liver Metastasis',
  
  # Epithelial
  'EP' = 'Epithelial',
  
  # Retina
  'RB' = 'Retina',
  
  # Testis
  'TGCT' = 'Testis',
  
  # Mesothelium
  'MESO' = 'Mesothelium',
  
  # Giant Cell Tumor of Bone
  'GCTB' = 'Bone',
  
  # Vestibular Schwannoma
  'VS' = 'Nerve',
  
  # Appendix
  'Appendiceal Tumors' = 'Appendix',
  
  # Multiple Chondromatosis
  'MCA' = 'Cartilage',
  
  # Other
  'other' = 'Other',
  'Cancer_general' = 'Other'
)

# Create tissue mapping function
map_cancer_to_tissue <- function(cancer_name) {
  if (is.na(cancer_name)) {
    return('Unknown')
  }
  
  cancer_name <- trimws(as.character(cancer_name))
  
  if (cancer_name %in% names(cancer_to_tissue)) {
    return(cancer_to_tissue[[cancer_name]])
  }
  
  return('Other')
}

# Apply mapping
df$Tissue_Type <- sapply(df$Cancer_general, map_cancer_to_tissue)

# Check unmapped cancers
unmapped_cancers <- unique(df$Cancer_general[df$Tissue_Type == 'Other'])
unmapped_cancers <- unmapped_cancers[!is.na(unmapped_cancers)]

if (length(unmapped_cancers) > 0) {
  cat("\nWarning: The following cancer types were not mapped to tissue types:\n")
  for (cancer in unmapped_cancers) {
    cat(paste0("  - ", cancer, "\n"))
  }
  cat("\n")
}

# ============ Basic Statistics ============
cat(sprintf("\nTotal Entries (Sample Count): %d\n", nrow(df)))
cat(sprintf("Total Cells: %s\n", format(sum(df$num_cells, na.rm = TRUE), big.mark = ",")))
cat(sprintf("Average Cells per Sample: %.0f\n", mean(df$num_cells, na.rm = TRUE)))
cat(sprintf("Median Cell Count: %.0f\n", median(df$num_cells, na.rm = TRUE)))

# Cancer type statistics
cancer_counts <- as.data.frame(table(df$Cancer_general))
colnames(cancer_counts) <- c("Cancer_Type", "Count")
cancer_counts <- cancer_counts %>% arrange(desc(Count))

cat(sprintf("\nNumber of Cancer Types: %d\n", nrow(cancer_counts)))
cat("\nSample Count by Cancer Type:\n")
print(cancer_counts)

# Tissue type statistics
tissue_counts <- as.data.frame(table(df$Tissue_Type))
colnames(tissue_counts) <- c("Tissue_Type", "Count")
tissue_counts <- tissue_counts %>% arrange(desc(Count))

cat(sprintf("\n\nNumber of Tissue Types: %d\n", nrow(tissue_counts)))
cat("\nSample Count by Tissue Type:\n")
print(tissue_counts)

# Tissue type cell statistics
tissue_cells <- df %>%
  group_by(Tissue_Type) %>%
  summarise(
    Total_Cells = sum(num_cells, na.rm = TRUE),
    Average_Cells = mean(num_cells, na.rm = TRUE),
    Sample_Count = n()
  ) %>%
  arrange(desc(Total_Cells))

cat("\nCell Statistics by Tissue Type:\n")
print(tissue_cells)

# Cancer to tissue mapping table
cancer_tissue_mapping <- df %>%
  group_by(Cancer_general, Tissue_Type) %>%
  summarise(Sample_Count = n(), .groups = 'drop') %>%
  arrange(Tissue_Type, desc(Sample_Count))

cat("\nCancer to Tissue Type Mapping:\n")
print(cancer_tissue_mapping)

# Cancer cell statistics
cancer_cells <- df %>%
  group_by(Cancer_general) %>%
  summarise(
    Total_Cells = sum(num_cells, na.rm = TRUE),
    Average_Cells = mean(num_cells, na.rm = TRUE),
    Sample_Count = n()
  ) %>%
  arrange(desc(Total_Cells))

# ============ Save Summary Data ============
# Basic statistics
basic_stats <- data.frame(
  Metric = c('Total Entries (Samples)', 'Total Cells', 'Average Cells', 
             'Median Cells', 'Cancer Types', 'Tissue Types'),
  Value = c(nrow(df), sum(df$num_cells, na.rm = TRUE), mean(df$num_cells, na.rm = TRUE),
           median(df$num_cells, na.rm = TRUE), nrow(cancer_counts), nrow(tissue_counts))
)
write.csv(basic_stats, file.path(output_dir, 'summary_basic_statistics.csv'), 
          row.names = FALSE, fileEncoding = "UTF-8")

# Cancer statistics
write.csv(cancer_counts, file.path(output_dir, 'summary_cancer_sample_counts.csv'), 
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(cancer_cells, file.path(output_dir, 'summary_cancer_cell_statistics.csv'), 
          row.names = FALSE, fileEncoding = "UTF-8")

# Tissue statistics
write.csv(tissue_counts, file.path(output_dir, 'summary_tissue_sample_counts.csv'), 
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tissue_cells, file.path(output_dir, 'summary_tissue_cell_statistics.csv'), 
          row.names = FALSE, fileEncoding = "UTF-8")

# Cancer-tissue mapping
write.csv(cancer_tissue_mapping, file.path(output_dir, 'summary_cancer_to_tissue_mapping.csv'), 
          row.names = FALSE, fileEncoding = "UTF-8")

# Save unmapped cancer list
if (length(unmapped_cancers) > 0) {
  write.csv(data.frame(Unmapped_Cancer = unmapped_cancers), 
            file.path(output_dir, 'unmapped_cancers.csv'), 
            row.names = FALSE, fileEncoding = "UTF-8")
}

cat(sprintf("\nData summary has been saved to the %s folder\n", output_dir))

# ============ Visualization Section ============

# Figure 1: All Cancer Types Sample Distribution (Horizontal Bar Chart)
cancer_counts_sorted <- cancer_counts %>% arrange(Count)
cancer_counts_sorted$Cancer_Type <- factor(cancer_counts_sorted$Cancer_Type, 
                                           levels = cancer_counts_sorted$Cancer_Type)

p1 <- ggplot(cancer_counts_sorted, aes(x = Cancer_Type, y = Count, fill = Cancer_Type)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = Count), hjust = -0.2, size = 3.5, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(nrow(cancer_counts_sorted))) +
  labs(
    title = "Sample Distribution Across Cancer Types",
    x = "Cancer Type",
    y = "Number of Samples"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "Fig1_all_cancer_sample_distribution.png"), 
       plot = p1, width = 10, height = max(8, nrow(cancer_counts_sorted) * 0.4), dpi = 300)
cat(sprintf("\nSaved: %s/Fig1_all_cancer_sample_distribution.png\n", output_dir))

# Figure 2: Top 10 Cancer Types Sample Distribution (Vertical Bar Chart)
top10_cancers <- cancer_counts %>% arrange(desc(Count)) %>% head(10)
top10_cancers$Cancer_Type <- factor(top10_cancers$Cancer_Type, 
                                     levels = top10_cancers$Cancer_Type)

p2 <- ggplot(top10_cancers, aes(x = Cancer_Type, y = Count, fill = Cancer_Type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, show.legend = FALSE) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Top 10 Cancer Types by Sample Count",
    x = "Cancer Type",
    y = "Number of Samples"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "Fig2_top10_cancer_sample_distribution.png"), 
       plot = p2, width = 12, height = 6, dpi = 300)
cat(sprintf("Saved: %s/Fig2_top10_cancer_sample_distribution.png\n", output_dir))

# Figure 3: Tissue Type Sample Distribution (Horizontal Bar Chart)
tissue_counts_sorted <- tissue_counts %>% arrange(Count)
tissue_counts_sorted$Tissue_Type <- factor(tissue_counts_sorted$Tissue_Type, 
                                           levels = tissue_counts_sorted$Tissue_Type)

p3 <- ggplot(tissue_counts_sorted, aes(x = Tissue_Type, y = Count, fill = Tissue_Type)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = Count), hjust = -0.2, size = 3.5, fontface = "bold") +
  coord_flip() +
  scale_fill_viridis_d() +
  labs(
    title = "Sample Distribution Across Tissue Types",
    x = "Tissue Type",
    y = "Number of Samples"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "Fig3_tissue_sample_distribution.png"), 
       plot = p3, width = 10, height = max(6, nrow(tissue_counts_sorted) * 0.5), dpi = 300)
cat(sprintf("Saved: %s/Fig3_tissue_sample_distribution.png\n", output_dir))

# Figure 4: Tissue Type Cell Distribution
tissue_total_cells <- df %>%
  group_by(Tissue_Type) %>%
  summarise(Total_Cells = sum(num_cells, na.rm = TRUE)) %>%
  arrange(Total_Cells)

tissue_total_cells$Tissue_Type <- factor(tissue_total_cells$Tissue_Type, 
                                         levels = tissue_total_cells$Tissue_Type)

p4 <- ggplot(tissue_total_cells, aes(x = Tissue_Type, y = Total_Cells, fill = Tissue_Type)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = format(Total_Cells, big.mark = ",")), 
            hjust = -0.1, size = 3, fontface = "bold") +
  coord_flip() +
  scale_fill_viridis_d(option = "plasma") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Total Cell Count Distribution Across Tissue Types",
    x = "Tissue Type",
    y = "Total Number of Cells"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "Fig4_tissue_cell_distribution.png"), 
       plot = p4, width = 10, height = max(6, nrow(tissue_total_cells) * 0.5), dpi = 300)
cat(sprintf("Saved: %s/Fig4_tissue_cell_distribution.png\n", output_dir))

# Figure 5: Tissue Type Proportion Pie Chart
tissue_counts_pie <- tissue_counts %>% 
  mutate(Percentage = Count / sum(Count) * 100,
         Label = paste0(Tissue_Type, "\n", round(Percentage, 1), "%"))

p5 <- ggplot(tissue_counts_pie, aes(x = "", y = Count, fill = Tissue_Type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set3"))(nrow(tissue_counts_pie))) +
  labs(
    title = "Tissue Type Sample Proportion",
    fill = "Tissue Type"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 20)),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "Fig5_tissue_pie_chart.png"), 
       plot = p5, width = 12, height = 8, dpi = 300)
cat(sprintf("Saved: %s/Fig5_tissue_pie_chart.png\n", output_dir))

# Figure 6: Top 10 Cancer Types Cell Distribution
top10_cancer_names <- (cancer_counts %>% arrange(desc(Count)) %>% head(10))$Cancer_Type
top10_cells <- df %>%
  filter(Cancer_general %in% top10_cancer_names) %>%
  group_by(Cancer_general) %>%
  summarise(Total_Cells = sum(num_cells, na.rm = TRUE))

# Arrange by original top 10 order
top10_cells$Cancer_general <- factor(top10_cells$Cancer_general, levels = top10_cancer_names)
top10_cells <- top10_cells %>% arrange(Cancer_general)

p6 <- ggplot(top10_cells, aes(x = Cancer_general, y = Total_Cells, fill = Cancer_general)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, show.legend = FALSE) +
  geom_text(aes(label = format(Total_Cells, big.mark = ",")), 
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Top 10 Cancer Types by Cell Count",
    x = "Cancer Type",
    y = "Number of Cells"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "Fig6_top10_cancer_cell_distribution.png"), 
       plot = p6, width = 12, height = 6, dpi = 300)
cat(sprintf("Saved: %s/Fig6_top10_cancer_cell_distribution.png\n", output_dir))

# Figure 7: Cancer-Tissue Heatmap
cancer_tissue_matrix <- df %>%
  group_by(Tissue_Type, Cancer_general) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Cancer_general, values_from = Count, values_fill = 0)

# Convert to matrix format
matrix_data <- as.matrix(cancer_tissue_matrix[, -1])
rownames(matrix_data) <- cancer_tissue_matrix$Tissue_Type

png(file.path(output_dir, "Fig7_cancer_tissue_heatmap.png"), 
    width = max(1600, ncol(matrix_data) * 35), 
    height = max(800, nrow(matrix_data) * 50), res = 300)

pheatmap(matrix_data,
         color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 9,
         fontsize_number = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Cancer-Tissue Distribution Heatmap",
         angle_col = 45,
         border_color = "grey60")

dev.off()
cat(sprintf("Saved: %s/Fig7_cancer_tissue_heatmap.png\n", output_dir))

# Figure 8: Tissue Type Sample Count vs Cell Count Comparison (Dual Y-axis)
tissue_stats <- df %>%
  group_by(Tissue_Type) %>%
  summarise(
    Sample_Count = n(),
    Total_Cells = sum(num_cells, na.rm = TRUE)
  ) %>%
  arrange(desc(Sample_Count))

tissue_stats$Tissue_Type <- factor(tissue_stats$Tissue_Type, 
                                   levels = tissue_stats$Tissue_Type)

# Calculate scaling factor for dual Y-axis display
scale_factor <- max(tissue_stats$Total_Cells) / max(tissue_stats$Sample_Count)

p8 <- ggplot(tissue_stats, aes(x = Tissue_Type)) +
  geom_bar(aes(y = Sample_Count, fill = "Sample Count"), stat = "identity", 
           position = position_nudge(x = -0.2), width = 0.4, alpha = 0.8) +
  geom_bar(aes(y = Total_Cells / scale_factor, fill = "Total Cells"), stat = "identity", 
           position = position_nudge(x = 0.2), width = 0.4, alpha = 0.8) +
  scale_y_continuous(
    name = "Number of Samples",
    sec.axis = sec_axis(~.*scale_factor, name = "Total Number of Cells",
                        labels = scales::comma)
  ) +
  scale_fill_manual(values = c("Sample Count" = "steelblue", "Total Cells" = "coral"),
                    name = "") +
  labs(
    title = "Sample Count vs Cell Count by Tissue Type",
    x = "Tissue Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "Fig8_tissue_sample_cell_comparison.png"), 
       plot = p8, width = 14, height = 6, dpi = 300)
cat(sprintf("Saved: %s/Fig8_tissue_sample_cell_comparison.png\n", output_dir))

# Completion message
cat("\n================================================================================\n")
cat("Analysis Complete! All charts and data summaries have been saved.\n")
cat("================================================================================\n")
cat(sprintf("\nAll files have been saved to folder: %s/\n", output_dir))
cat("\nGenerated file list:\n")
cat("\n  CSV Files:\n")
cat("    - summary_basic_statistics.csv (Basic statistics)\n")
cat("    - summary_cancer_sample_counts.csv (Sample counts by cancer type)\n")
cat("    - summary_cancer_cell_statistics.csv (Cell statistics by cancer type)\n")
cat("    - summary_tissue_sample_counts.csv (Sample counts by tissue type)\n")
cat("    - summary_tissue_cell_statistics.csv (Cell statistics by tissue type)\n")
cat("    - summary_cancer_to_tissue_mapping.csv (Cancer to tissue mapping table)\n")
if (length(unmapped_cancers) > 0) {
  cat("    - unmapped_cancers.csv (List of unmapped cancer types)\n")
}
cat("\n  Chart Files:\n")
cat("    - Fig1_all_cancer_sample_distribution.png (All cancer types sample distribution)\n")
cat("    - Fig2_top10_cancer_sample_distribution.png (Top 10 cancer types sample distribution)\n")
cat("    - Fig3_tissue_sample_distribution.png (Tissue type sample distribution)\n")
cat("    - Fig4_tissue_cell_distribution.png (Tissue type cell count distribution)\n")
cat("    - Fig5_tissue_pie_chart.png (Tissue type proportion pie chart)\n")
cat("    - Fig6_top10_cancer_cell_distribution.png (Top 10 cancer types cell distribution)\n")
cat("    - Fig7_cancer_tissue_heatmap.png (Cancer-tissue distribution heatmap)\n")
cat("    - Fig8_tissue_sample_cell_comparison.png (Tissue sample count vs cell count comparison)\n")
cat("================================================================================\n") 