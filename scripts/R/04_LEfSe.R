#!/usr/bin/env Rscript

# ==============================================================================
# 04_LEfSe.R — LEfSe LDA条形图 (基于input.res结果)
# ==============================================================================
# 读取LEfSe的 input.res 结果文件，用ggplot2重新绘制
# 正文版: LDA > 4 (精简，约40个feature)
# 补充版: LDA > 2 (全部，约91个feature)
# 只展示属和种水平（去掉门纲目科株等冗余层级）
#
# 用法:
#   Rscript 04_LEfSe.R \
#     --input temp/input.res \
#     --output ~/thesis_project/figures/paper_en/Fig2-7 \
#     --stat_dir ~/thesis_project/data/processed \
#     --style paper
# ==============================================================================

options(warn = -1)

site <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only=TRUE, quietly=TRUE)))) {
  install.packages("optparse", repos=site); require("optparse", character.only=TRUE)
}

option_list <- list(
  make_option(c("-i", "--input"),  type="character", default="temp/input.res"),
  make_option(c("-o", "--output"), type="character", default="Fig2-7"),
  make_option(c("--stat_dir"),    type="character", default="."),
  make_option(c("-w", "--width"),  type="numeric",   default=180),
  make_option(c("-e", "--height"), type="numeric",   default=200),
  make_option(c("-s", "--style"),  type="character",  default="paper")
)
opts <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages({
  library(ggplot2); library(dplyr); library(stringr)
}))

# ---- 字体 ----
sys_fonts <- tryCatch(unique(systemfonts::system_fonts()$family), error = function(e) character(0))
if (opts$style == "thesis") {
  FONT <- ifelse("SimSun" %in% sys_fonts, "SimSun",
          ifelse("宋体" %in% sys_fonts, "宋体", "serif"))
} else {
  FONT <- ifelse("Arial" %in% sys_fonts, "Arial", "sans")
}
cat("Style:", opts$style, "Font:", FONT, "\n")

# ---- 配色 (方案B) ----
GROUP_COLORS <- c(CK = "#1c9cbd", A = "#e3943d", B = "#92181e")

# ==============================================================================
# 1. 读取 LEfSe 结果
# ==============================================================================
res <- read.table(opts$input, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                  fill=TRUE, quote="")
colnames(res) <- c("Feature", "LogMaxMean", "Group", "LDA", "pvalue")

# 只保留有显著差异的
res <- res %>% filter(Group != "" & Group != "-" & !is.na(LDA) & LDA != "")
res$LDA <- as.numeric(res$LDA)
res$pvalue <- as.numeric(res$pvalue)
res <- res %>% filter(!is.na(LDA))

cat("总显著特征:", nrow(res), "\n")

# ---- 提取分类层级和简短名称 ----
res$n_levels <- str_count(res$Feature, "\\.") + 1
res$tax_level <- case_when(
  res$n_levels == 1 ~ "Kingdom",
  res$n_levels == 2 ~ "Phylum",
  res$n_levels == 3 ~ "Class",
  res$n_levels == 4 ~ "Order",
  res$n_levels == 5 ~ "Family",
  res$n_levels == 6 ~ "Genus",
  res$n_levels == 7 ~ "Species",
  res$n_levels == 8 ~ "Strain",
  TRUE ~ "Other"
)

# 提取最后一级名称
res$short_name <- sapply(strsplit(res$Feature, "\\."), function(x) tail(x, 1))
# 去掉前缀 k__ p__ c__ 等
res$short_name <- gsub("^[kpcofgst]__", "", res$short_name)

res$Group <- factor(res$Group, levels = c("CK", "A", "B"))

cat("各组: ", table(res$Group), "\n")
cat("各层级: ", table(res$tax_level), "\n")

# ==============================================================================
# 2. 保存统计结果
# ==============================================================================
dir.create(opts$stat_dir, showWarnings=FALSE, recursive=TRUE)

# 保存完整结果表
stat_file <- file.path(opts$stat_dir, "LEfSe_significant_features.txt")
res_out <- res %>%
  select(Feature, short_name, tax_level, Group, LDA, pvalue) %>%
  arrange(Group, desc(LDA))
write.table(res_out, file=stat_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("统计已保存:", stat_file, "\n")

# ==============================================================================
# 3. LDA 条形图绘制函数
# ==============================================================================

draw_lda_bar <- function(data, title_suffix = "", letter = "a") {

  # B组LDA取负值（左侧），CK和A取正值（右侧）
  # 按惯例：enriched in不同组的feature分左右
  # 这里按组着色，LDA值都向右画，按组分段排列

  # 按组和LDA排序
  data <- data %>%
    arrange(Group, LDA) %>%
    mutate(short_name = factor(short_name, levels = unique(short_name)))

  p <- ggplot(data, aes(x = short_name, y = LDA, fill = Group)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = GROUP_COLORS) +
    labs(
      x = NULL,
      y = "LDA Score (log 10)",
      title = paste0("(", letter, ") ", title_suffix),
      fill = "Enriched in"
    ) +
    theme_classic(base_family = FONT) +
    theme(
      text = element_text(family = FONT),
      axis.text.y = element_text(family = FONT, size = 7, color = "black"),
      axis.text.x = element_text(family = FONT, size = 8, color = "black"),
      axis.title = element_text(family = FONT, size = 9),
      legend.text = element_text(family = FONT, size = 8),
      legend.title = element_text(family = FONT, size = 9, face = "bold"),
      legend.position = "top",
      plot.title = element_text(family = FONT, size = 10, face = "bold"),
      plot.margin = margin(5, 10, 5, 5, "mm")
    )

  return(p)
}

# ==============================================================================
# 4. 正文版: 只展示属+种水平, LDA > 4
# ==============================================================================
cat("\n--- 正文版 (Genus + Species, LDA > 4) ---\n")

main_data <- res %>%
  filter(tax_level %in% c("Genus", "Species") & LDA > 4)
cat("  特征数:", nrow(main_data), "\n")

p_main <- draw_lda_bar(main_data, "Genus & Species level, LDA > 4", "a")

# 高度根据特征数动态调整
main_h <- max(80, nrow(main_data) * 5 + 30)

ggsave(paste0(opts$output, "_LDA_main.pdf"), p_main,
       width = opts$width, height = main_h, units = "mm")
ggsave(paste0(opts$output, "_LDA_main.png"), p_main,
       width = opts$width, height = main_h, units = "mm", dpi = 300)
cat("  ->", paste0(opts$output, "_LDA_main.pdf/png\n"))

# ==============================================================================
# 5. 完整版: 所有分类水平, LDA > 2
# ==============================================================================
cat("\n--- 完整版 (All levels, LDA > 2) ---\n")

full_data <- res %>% filter(LDA > 2)
cat("  特征数:", nrow(full_data), "\n")

p_full <- draw_lda_bar(full_data, "All levels, LDA > 2", "a")

full_h <- max(120, nrow(full_data) * 4.5 + 30)

ggsave(paste0(opts$output, "_LDA_full.pdf"), p_full,
       width = opts$width, height = full_h, units = "mm")
ggsave(paste0(opts$output, "_LDA_full.png"), p_full,
       width = opts$width, height = full_h, units = "mm", dpi = 300)
cat("  ->", paste0(opts$output, "_LDA_full.pdf/png\n"))

# ==============================================================================
# 6. 属+种水平, LDA > 3 (推荐正文使用)
# ==============================================================================
cat("\n--- 推荐版 (Genus + Species, LDA > 3) ---\n")

rec_data <- res %>%
  filter(tax_level %in% c("Genus", "Species") & LDA > 3)
cat("  特征数:", nrow(rec_data), "\n")

p_rec <- draw_lda_bar(rec_data, "Genus & Species level, LDA > 3", "a")

rec_h <- max(100, nrow(rec_data) * 5 + 30)

ggsave(paste0(opts$output, "_LDA_genus_species.pdf"), p_rec,
       width = opts$width, height = rec_h, units = "mm")
ggsave(paste0(opts$output, "_LDA_genus_species.png"), p_rec,
       width = opts$width, height = rec_h, units = "mm", dpi = 300)
cat("  ->", paste0(opts$output, "_LDA_genus_species.pdf/png\n"))

cat("\nDone!\n")
cat("\n建议:\n")
cat("  正文使用: Fig2-7_LDA_genus_species (属+种, LDA>3, ", nrow(rec_data), "个feature)\n")
cat("  补充材料: Fig2-7_LDA_full (全层级, LDA>2, ", nrow(full_data), "个feature)\n")
