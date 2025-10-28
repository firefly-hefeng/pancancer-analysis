#!/usr/bin/env Rscript
## 目的：遍历所有 <cancer_type>_cell_type_stats.csv
##       提取所有出现过的 fine_cell_type 名字（去重+排序）
## 用法：Rscript list_fine_cell_types.R

library(tidyverse)

## 1. 找到当前目录下所有 *cell_type_stats.csv 文件
## 递归查找当前目录及所有子目录
csv_files <- list.files(path = ".",
                        pattern = "*cell_type_stats.csv",
                        recursive = TRUE,
                        full.names = TRUE)
if (length(csv_files) == 0) {
  stop("未找到任何 *cell_type_stats.csv 文件，请确认运行目录正确！")
}

## 2. 批量读取并合并
all_df <- map_dfr(csv_files, function(f) {
  # 文件名格式：<cancer>_cell_type_stats.csv
  cancer <- tools::file_path_sans_ext(basename(f)) %>% str_remove("_cell_type_stats")
  read_csv(f, col_types = cols()) %>% 
    select(fine_cell_type = 1) %>%   # 第一列永远是 fine_cell_type
    mutate(cancer = cancer)
})

## 3. 去重+排序
fine_types <- sort(unique(all_df$fine_cell_type))

## 4. 屏幕打印
cat(">>> 所有出现过的 fine_cell_type 名称（去重后）：\n")
writeLines(fine_types)

## 5. 可选：保存为文本 / Excel
write_lines(fine_types, "all_fine_cell_types.txt")
readr::write_csv(tibble(fine_cell_type = fine_types),
                 "all_fine_cell_types.csv")

cat("\n结果已写入：all_fine_cell_types.txt / .csv\n")