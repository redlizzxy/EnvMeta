#!/usr/bin/env Rscript

# ==============================================================================
# 01_tax_stackplot.R — 基于 tax_stackplot_fixed.R (amplicon::tax_stackplot)
# 改动: 配色→STACK_12 | 样品名→A_1 | 去unclassified | 字体 | +combined拼合
# 调用: Rscript 01_tax_stackplot.R --input X.txt --design metadata.txt \
#         --group Group --output prefix --legend 10 --width 120 --height 90
# ==============================================================================

options(warn = -1)

site <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only=TRUE, quietly=TRUE)))) {
  install.packages("optparse", repos=site); require("optparse", character.only=TRUE)
}

option_list <- list(
  make_option(c("-i", "--input"),   type="character", default="metaphlan4/Phylum.txt"),
  make_option(c("-d", "--design"),  type="character", default="metadata.txt"),
  make_option(c("-n", "--group"),   type="character", default="Group"),
  make_option(c("-o", "--output"),  type="character", default=""),
  make_option(c("-l", "--legend"),  type="numeric",   default=10),
  make_option(c("-c", "--color"),   type="character", default="custom"),
  make_option(c("-w", "--width"),   type="numeric",   default=120),
  make_option(c("-e", "--height"),  type="numeric",   default=90),
  make_option(c("-s", "--style"),   type="character", default="paper"),
  make_option(c("-p", "--percent_dir"), type="character", default="")
)
opts <- parse_args(OptionParser(option_list=option_list))
if (opts$output == "") opts$output <- opts$input
if (opts$percent_dir == "") opts$percent_dir <- dirname(opts$output)

suppressWarnings(suppressMessages(library(amplicon)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(RColorBrewer)))

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
STACK_12 <- c(
  "#bbc4e4", "#d7e7af", "#bf9d6d", "#a4d6c1", "#bce2e8", "#b5b5b6",
  "#0084c2", "#3ab483", "#00978c", "#21825e", "#6aa3d2", "#3d62ad"
)

font_theme <- theme(
  text = element_text(family = FONT),
  axis.text = element_text(family = FONT, color = "black"),
  axis.title = element_text(family = FONT),
  legend.text = element_text(family = FONT),
  legend.title = element_text(family = FONT)
)

# ==============================================================================
# 读取数据
# ==============================================================================
metadata <- read.table(opts$design, header=TRUE, row.names=1, sep="\t",
                       comment.char="", stringsAsFactors=FALSE)
taxonomy <- read.table(opts$input, header=TRUE, row.names=1, sep="\t",
                       comment.char="", quote="", check.names=FALSE)

cat("原始维度:", dim(taxonomy), "\n")

# ---- 去 unclassified ----
keep <- !grepl("(?i)(unclassif|unassign|unknown)", rownames(taxonomy), perl=TRUE)
cat("去除unclassified:", sum(!keep), "行\n")
taxonomy <- taxonomy[keep, , drop=FALSE]

# ---- 样品名重命名: 8_1 → A_1 ----
design_full <- read.table(opts$design, header=TRUE, sep="\t", stringsAsFactors=FALSE)
name_map <- setNames(
  paste0(design_full[[opts$group]], "_", design_full$Replicate),
  design_full$SampleID
)
colnames(taxonomy) <- ifelse(colnames(taxonomy) %in% names(name_map),
                             name_map[colnames(taxonomy)], colnames(taxonomy))
rownames(metadata) <- ifelse(rownames(metadata) %in% names(name_map),
                             name_map[rownames(metadata)], rownames(metadata))

cat("样品:", colnames(taxonomy), "\n")
cat("维度:", dim(taxonomy), "\n")

# ==============================================================================
# 绘图（仅 >=2 行）
# ==============================================================================
if (nrow(taxonomy) < 2) { cat("跳过: 不足2行\n"); quit(status=0) }

dir.create(opts$percent_dir, showWarnings=FALSE, recursive=TRUE)
pct_prefix <- file.path(opts$percent_dir, basename(opts$output))

W <- opts$width; H <- opts$height

# ==============================================================================
# Sample 版
# ==============================================================================
cat("\n--- Sample ---\n")
p <- tax_stackplot(taxonomy, metadata, topN=opts$legend,
                   groupID=opts$group, style="sample", sorted="abundance")
p <- p + scale_fill_manual(values=STACK_12) + font_theme

ggsave(paste0(opts$output, ".sample.pdf"), p, width=W, height=H, units="mm")
ggsave(paste0(opts$output, ".sample.png"), p, width=W, height=H, units="mm", dpi=300)
cat("  →", paste0(opts$output, ".sample.pdf/png\n"))

# 提取百分比（原脚本逻辑）
plot_data_sample <- ggplot2::ggplot_build(p)$data[[1]]
sample_species <- unique(plot_data_sample$fill)
sample_names   <- unique(as.character(plot_data_sample$x))
pct_s <- matrix(NA, nrow=length(sample_species), ncol=length(sample_names))
rownames(pct_s) <- sample_species; colnames(pct_s) <- sample_names
for (i in 1:nrow(plot_data_sample)) {
  pct_s[plot_data_sample$fill[i], as.character(plot_data_sample$x[i])] <-
    (plot_data_sample$ymax[i] - plot_data_sample$ymin[i]) * 100
}
write.table(data.frame(Species=rownames(pct_s), pct_s, check.names=FALSE),
            file=paste0(pct_prefix, ".species_percentage.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ==============================================================================
# Group 版
# ==============================================================================
cat("\n--- Group ---\n")
p_group <- tax_stackplot(taxonomy, metadata, topN=opts$legend,
                         groupID=opts$group, style="group", sorted="abundance")
p_group <- p_group + scale_fill_manual(values=STACK_12) + font_theme

ggsave(paste0(opts$output, ".group.pdf"), p_group, width=W, height=H, units="mm")
ggsave(paste0(opts$output, ".group.png"), p_group, width=W, height=H, units="mm", dpi=300)
cat("  →", paste0(opts$output, ".group.pdf/png\n"))

# 提取组百分比
plot_data_g <- ggplot2::ggplot_build(p_group)$data[[1]]
gp_species <- unique(plot_data_g$fill)
gp_names   <- unique(as.character(plot_data_g$x))
pct_g <- matrix(NA, nrow=length(gp_species), ncol=length(gp_names))
rownames(pct_g) <- gp_species; colnames(pct_g) <- gp_names
for (i in 1:nrow(plot_data_g)) {
  pct_g[plot_data_g$fill[i], as.character(plot_data_g$x[i])] <-
    (plot_data_g$ymax[i] - plot_data_g$ymin[i]) * 100
}
write.table(data.frame(Species=rownames(pct_g), pct_g, check.names=FALSE),
            file=paste0(pct_prefix, ".group_mean_percentage.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ==============================================================================
# Combined 版：用 cowplot 把 sample.png 和 group.png 拼合
# 两图各自用原始尺寸生成，拼合时 sample 占 60% group 占 40%
# ==============================================================================
cat("\n--- Combined ---\n")

if (!requireNamespace("cowplot", quietly=TRUE)) install.packages("cowplot", repos=site)
library(cowplot)

p_s <- p + theme(legend.position = "none")
p_g <- p_group

p_combined <- plot_grid(p_s, p_g, nrow=1, rel_widths=c(0.58, 0.42))

# combined 总宽 = sample宽度 / 0.58（让sample部分保持原始大小）
W_comb <- round(W / 0.58)
ggsave(paste0(opts$output, ".combined.pdf"), p_combined,
       width=W_comb, height=H, units="mm")
ggsave(paste0(opts$output, ".combined.png"), p_combined,
       width=W_comb, height=H, units="mm", dpi=300)
cat("  →", paste0(opts$output, ".combined.pdf/png\n"))

cat("\n百分比文件:\n")
cat("  ", paste0(pct_prefix, ".species_percentage.txt\n"))
cat("  ", paste0(pct_prefix, ".group_mean_percentage.txt\n"))
cat("\nDone!\n")
