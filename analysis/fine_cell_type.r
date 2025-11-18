# 将输出文件保存在fine_annotation目录内
extract_fine_cell_types_from_stats <- function(
  base_dir = "/mnt/public7/pancancercol/hefeng/cluster-annotation/fine_annotation",
  output_file = "fine_cell_types_list.txt"
) {
  
  message("开始从cell_type_stats.csv文件提取fine_cell_type...")
  message("扫描目录: ", base_dir)
  
  all_fine_types <- c()
  
  # 获取所有癌种目录
  cancer_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  
  message("找到 ", length(cancer_dirs), " 个癌种目录\n")
  
  for (cancer_dir in cancer_dirs) {
    cancer_type <- basename(cancer_dir)
    
    # 构建stats文件路径
    stats_file <- file.path(cancer_dir, paste0(cancer_type, "_cell_type_stats.csv"))
    
    if (file.exists(stats_file)) {
      message("读取: ", cancer_type)
      
      tryCatch({
        # 读取CSV文件
        stats_data <- read.csv(stats_file, row.names = 1, check.names = FALSE)
        
        # 提取行名（即fine_cell_type）
        fine_types <- rownames(stats_data)
        
        # 过滤掉Other和Unknown
        fine_types <- fine_types[!fine_types %in% c("Other", "Unknown")]
        
        all_fine_types <- c(all_fine_types, fine_types)
        
        message("  发现 ", length(fine_types), " 种细胞类型")
        
      }, error = function(e) {
        warning("读取 ", stats_file, " 失败: ", e$message)
      })
    } else {
      message("文件不存在: ", stats_file)
    }
  }
  
  # 去重并排序
  all_fine_types <- unique(all_fine_types)
  all_fine_types <- sort(all_fine_types)
  
  message("\n总共发现 ", length(all_fine_types), " 种不同的fine_cell_type")
  
  # 保存为txt格式（R向量格式）
  # 保存在fine_annotation目录内
  output_path <- file.path(base_dir, output_file)
  
  # 写入文件头
  cat('cell_types_list <- c(\n', file = output_path)
  
  # 按主要细胞类型分组写入
  current_major <- ""
  
  for (i in 1:length(all_fine_types)) {
    cell_type <- all_fine_types[i]
    
    # 提取主要类型（下划线前的部分）
    major_type <- strsplit(cell_type, "_")[[1]][1]
    
    # 如果是新的主要类型，添加注释
    if (major_type != current_major) {
      if (i > 1) {
        cat('\n', file = output_path, append = TRUE)
      }
      cat(paste0('  # ', major_type, '\n'), file = output_path, append = TRUE)
      current_major <- major_type
    }
    
    # 写入细胞类型
    if (i < length(all_fine_types)) {
      cat(paste0('  "', cell_type, '",\n'), file = output_path, append = TRUE)
    } else {
      cat(paste0('  "', cell_type, '"\n'), file = output_path, append = TRUE)
    }
  }
  
  cat(')\n', file = output_path, append = TRUE)
  
  message("\n✓ 结果已保存到: ", output_path)
  
  # 打印统计信息
  message("\n=== 统计信息 ===")
  major_type_counts <- table(sapply(strsplit(all_fine_types, "_"), `[`, 1))
  for (mt in names(sort(major_type_counts, decreasing = TRUE))) {
    message("  ", mt, ": ", major_type_counts[mt], " 种亚型")
  }
  
  # 返回向量
  return(all_fine_types)
}

# 使用方法
fine_types <- extract_fine_cell_types_from_stats(
  base_dir = "/mnt/public7/pancancercol/hefeng/cluster-annotation/fine_annotation",
  output_file = "fine_cell_types_list.txt"
)

# 打印前30个作为预览
message("\n=== 前30个细胞类型 ===")
print(head(fine_types, 30))

# 打印总数
message("\n总计: ", length(fine_types), " 种细胞类型")