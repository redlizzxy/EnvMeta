#!/usr/bin/env Rscript

# ==============================================================================
# 02_alpha_diversity.R — α多样性箱线图 (Kraken2 + MetaPhlAn4 双数据源)
# ==============================================================================
# 上排 (Kraken2): Richness / Chao1 / ACE — 丰富度估计
# 下排 (MetaPhlAn4): Observed Species / Shannon / Simpson — 多样性指数
# 箱体浅色填充 + 原色散点 + 样品标签 + Wilcoxon p值
#
# 用法:
#   Rscript 02_alpha_diversity.R \
#     --kraken2 ~/meta2/result/kraken2/tax_count.alpha \
#     --metaphlan4 ~/meta2/result/metaphlan4/alpha.txt \
#     --design ~/meta2/result/metadata.txt \
#     --output ~/thesis_project/figures/paper_en/Fig2-4 \
#     --stat_dir ~/thesis_project/data/processed \
#     --style paper
# ==============================================================================

options(warn = -1)

site <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only=TRUE, quietly=TRUE)))) {
  install.packages("optparse", repos=site); require("optparse", character.only=TRUE)
}

option_list <- list(
  make_option(c("--kraken2"),     type="character", default="kraken2/tax_count.alpha"),
  make_option(c("--metaphlan4"),  type="character", default="metaphlan4/alpha.txt"),
  make_option(c("-d", "--design"),type="character", default="metadata.txt"),
  make_option(c("-n", "--group"), type="character", default="Group"),
  make_option(c("-o", "--output"),type="character", default="Fig2-4"),
  make_option(c("--stat_dir"),   type="character", default="."),
  make_option(c("-w", "--width"), type="numeric",   default=89),
  make_option(c("-e", "--height"),type="numeric",   default=75),
  make_option(c("-s", "--style"), type="character",  default="paper")
)
opts <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages({
  library(ggplot2); library(dplyr); library(ggpubr); library(ggrepel); library(patchwork)
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

# ---- 配色 ----
GROUP_COLORS <- c(CK = "#1c9cbd", A = "#e3943d", B = "#92181e")
GROUP_FILLS  <- c(CK = "#d4e5ef", A = "#fbd8b5", B = "#d0916e")

# ---- 主题 ----
theme_alpha <- function() {
  theme_classic(base_family = FONT) +
    theme(
      text = element_text(family = FONT),
      axis.text = element_text(color = "black", family = FONT, size = 9),
      axis.title = element_text(family = FONT, size = 10),
      legend.text = element_text(family = FONT, size = 8),
      legend.title = element_text(family = FONT, size = 9, face = "bold"),
      plot.title = element_text(family = FONT, size = 10, face = "bold", hjust = 0),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
}

comparisons <- list(c("CK", "A"), c("CK", "B"), c("A", "B"))

# ==============================================================================
# 1. 读取数据
# ==============================================================================
metadata <- read.table(opts$design, header=TRUE, row.names=1, sep="\t",
                       comment.char="", stringsAsFactors=FALSE)

# 样品名映射
design_full <- read.table(opts$design, header=TRUE, sep="\t", stringsAsFactors=FALSE)
name_map <- setNames(
  paste0(design_full[[opts$group]], "_", design_full$Replicate),
  design_full$SampleID
)

# ---- Kraken2 ----
kr2 <- read.table(opts$kraken2, header=TRUE, row.names=1, sep="\t", comment.char="")
common_kr2 <- intersect(rownames(kr2), rownames(metadata))
kr2 <- kr2[common_kr2, , drop=FALSE]
kr2_meta <- metadata[common_kr2, , drop=FALSE]
df_kr2 <- cbind(Sample = rownames(kr2), kr2, Group = kr2_meta[[opts$group]])
df_kr2$Group <- factor(df_kr2$Group, levels = c("CK", "A", "B"))
df_kr2$Label <- ifelse(df_kr2$Sample %in% names(name_map), name_map[df_kr2$Sample], df_kr2$Sample)
cat("Kraken2 样本:", nrow(df_kr2), " 指数:", colnames(kr2), "\n")

# ---- MetaPhlAn4 ----
mp4 <- read.table(opts$metaphlan4, header=TRUE, row.names=1, sep="\t", comment.char="")
common_mp4 <- intersect(rownames(mp4), rownames(metadata))
mp4 <- mp4[common_mp4, , drop=FALSE]
mp4_meta <- metadata[common_mp4, , drop=FALSE]
df_mp4 <- cbind(Sample = rownames(mp4), mp4, Group = mp4_meta[[opts$group]])
df_mp4$Group <- factor(df_mp4$Group, levels = c("CK", "A", "B"))
df_mp4$Label <- ifelse(df_mp4$Sample %in% names(name_map), name_map[df_mp4$Sample], df_mp4$Sample)
cat("MetaPhlAn4 样本:", nrow(df_mp4), " 指数:", colnames(mp4), "\n")

# ==============================================================================
# 2. 统计检验 (两个数据源都做)
# ==============================================================================
dir.create(opts$stat_dir, showWarnings=FALSE, recursive=TRUE)
stat_file <- file.path(opts$stat_dir, "alpha_wilcoxon_results.txt")

do_stats <- function(df, indices, source_name) {
  cat("\n========== ", source_name, " ==========\n")
  for (idx in indices) {
    if (!idx %in% colnames(df)) next
    cat("--- ", idx, " ---\n")
    for (g in c("CK", "A", "B")) {
      vals <- df[df$Group == g, idx]
      cat(sprintf("  %s: %.4f +/- %.4f (n=%d)\n", g, mean(vals, na.rm=TRUE), sd(vals, na.rm=TRUE), length(vals)))
    }
    cat("  Pairwise Wilcoxon:\n")
    for (comp in comparisons) {
      v1 <- df[df$Group == comp[1], idx]
      v2 <- df[df$Group == comp[2], idx]
      wt <- tryCatch(wilcox.test(v1, v2, exact=FALSE), error = function(e) list(p.value=NA))
      cat(sprintf("    %s vs %s: p = %.4f\n", comp[1], comp[2], wt$p.value))
    }
    kt <- kruskal.test(as.formula(paste0("`", idx, "` ~ Group")), data = df)
    cat(sprintf("  Kruskal-Wallis: chi2 = %.3f, p = %.4f\n\n", kt$statistic, kt$p.value))
  }
}

sink(stat_file)
cat("=== Alpha多样性 Wilcoxon秩和检验结果 ===\n")
cat("生成时间:", format(Sys.time()), "\n")
do_stats(df_kr2, c("richness", "chao1", "ACE"), "Kraken2")
do_stats(df_mp4, c("observed_species", "shannon", "simpson"), "MetaPhlAn4")
sink()
cat("统计已保存:", stat_file, "\n")

# ==============================================================================
# 3. 绘图函数
# ==============================================================================

draw_alpha <- function(df, idx_name, display_name, letter) {

  p <- ggplot(df, aes(x = Group, y = .data[[idx_name]])) +
    # 箱体浅色填充
    geom_boxplot(aes(fill = Group), width = 0.5, linewidth = 0.4,
                 outlier.shape = NA, color = "black", alpha = 0.9) +
    # 散点原色
    geom_jitter(aes(color = Group), width = 0.12, size = 2.2, alpha = 0.85) +
    # 样品名标签
    geom_text_repel(aes(label = Label, color = Group),
                    size = 2.0, family = FONT, show.legend = FALSE,
                    max.overlaps = 20, seed = 42,
                    segment.color = "grey70", segment.size = 0.2,
                    nudge_x = 0.25, direction = "y") +
    # Wilcoxon p值
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.format",
      size = 2.3, family = FONT,
      step.increase = 0.08,
      tip.length = 0.01,
      bracket.size = 0.3
    ) +
    scale_fill_manual(values = GROUP_FILLS) +
    scale_color_manual(values = GROUP_COLORS) +
    labs(x = NULL, y = display_name, title = paste0("(", letter, ")")) +
    theme_alpha() +
    theme(legend.position = "none")

  return(p)
}

# ==============================================================================
# 4. 单图
# ==============================================================================
cat("\n--- 单图 ---\n")
W <- opts$width; H <- opts$height

# 上排: Kraken2
kr2_indices <- list(richness = "Richness", chao1 = "Chao1", ACE = "ACE")
kr2_indices <- kr2_indices[names(kr2_indices) %in% colnames(df_kr2)]

# 下排: MetaPhlAn4
mp4_indices <- list(observed_species = "Observed Species", shannon = "Shannon", simpson = "Simpson")
mp4_indices <- mp4_indices[names(mp4_indices) %in% colnames(df_mp4)]

all_single <- list()
ltr_i <- 1

for (idx_name in names(kr2_indices)) {
  ltr <- letters[ltr_i]
  p <- draw_alpha(df_kr2, idx_name, kr2_indices[[idx_name]], ltr)
  all_single[[paste0("kr2_", idx_name)]] <- p
  ggsave(paste0(opts$output, "_", ltr, "_", idx_name, "_kraken2.pdf"), p, width=W, height=H, units="mm")
  ggsave(paste0(opts$output, "_", ltr, "_", idx_name, "_kraken2.png"), p, width=W, height=H, units="mm", dpi=300)
  cat("  ->", ltr, idx_name, "(Kraken2)\n")
  ltr_i <- ltr_i + 1
}

for (idx_name in names(mp4_indices)) {
  ltr <- letters[ltr_i]
  p <- draw_alpha(df_mp4, idx_name, mp4_indices[[idx_name]], ltr)
  all_single[[paste0("mp4_", idx_name)]] <- p
  ggsave(paste0(opts$output, "_", ltr, "_", idx_name, "_metaphlan4.pdf"), p, width=W, height=H, units="mm")
  ggsave(paste0(opts$output, "_", ltr, "_", idx_name, "_metaphlan4.png"), p, width=W, height=H, units="mm", dpi=300)
  cat("  ->", ltr, idx_name, "(MetaPhlAn4)\n")
  ltr_i <- ltr_i + 1
}

# ==============================================================================
# 5. 合并图 (2行3列)
# ==============================================================================
cat("\n--- 合并图 ---\n")

# 上排标题
title_kr2 <- ggplot() + theme_void() +
  annotate("text", x=0.5, y=0.5, label="Kraken2-based (k-mer alignment)",
           size=3.5, family=FONT, fontface="bold.italic", color="grey30") +
  theme(plot.margin = margin(2, 0, 0, 0, "mm"))

title_mp4 <- ggplot() + theme_void() +
  annotate("text", x=0.5, y=0.5, label="MetaPhlAn4-based (marker gene)",
           size=3.5, family=FONT, fontface="bold.italic", color="grey30") +
  theme(plot.margin = margin(2, 0, 0, 0, "mm"))

# 组装上排 (Kraken2)
kr2_plots <- list()
for (i in seq_along(kr2_indices)) {
  kr2_plots[[i]] <- draw_alpha(df_kr2, names(kr2_indices)[i], kr2_indices[[i]], letters[i])
}
row_kr2 <- wrap_plots(kr2_plots, nrow = 1)

# 组装下排 (MetaPhlAn4)
mp4_plots <- list()
offset <- length(kr2_indices)
for (i in seq_along(mp4_indices)) {
  mp4_plots[[i]] <- draw_alpha(df_mp4, names(mp4_indices)[i], mp4_indices[[i]], letters[offset + i])
}
row_mp4 <- wrap_plots(mp4_plots, nrow = 1)

# 总合并: 标题 + 行 + 标题 + 行
p_combined <- title_kr2 / row_kr2 / title_mp4 / row_mp4 +
  plot_layout(heights = c(0.06, 1, 0.06, 1))

comb_w <- W * 3.2
comb_h <- H * 2.5

ggsave(paste0(opts$output, "_combined.pdf"), p_combined,
       width = comb_w, height = comb_h, units = "mm")
ggsave(paste0(opts$output, "_combined.png"), p_combined,
       width = comb_w, height = comb_h, units = "mm", dpi = 300)
cat("  ->", paste0(opts$output, "_combined.pdf/png\n"))

cat("\nDone!\n")
