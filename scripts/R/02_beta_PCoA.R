#!/usr/bin/env Rscript

# ==============================================================================
# 02_beta_PCoA.R — PCoA + PERMANOVA (无椭圆版)
# ==============================================================================

options(warn = -1)

site <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only=TRUE, quietly=TRUE)))) {
  install.packages("optparse", repos=site); require("optparse", character.only=TRUE)
}

option_list <- list(
  make_option(c("--dist"),      type="character", default="metaphlan4/beta_bray.txt"),
  make_option(c("-d", "--design"), type="character", default="metadata.txt"),
  make_option(c("-n", "--group"),  type="character", default="Group"),
  make_option(c("-o", "--output"), type="character", default="Fig2-5_PCoA"),
  make_option(c("--stat_dir"),    type="character", default="."),
  make_option(c("-w", "--width"),  type="numeric", default=120),
  make_option(c("-e", "--height"), type="numeric", default=90),
  make_option(c("-s", "--style"),  type="character", default="paper")
)
opts <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages({
  library(vegan); library(ggplot2); library(ggrepel)
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
# 1. 读取
# ==============================================================================
dist_matrix <- as.matrix(read.table(opts$dist, header=TRUE, row.names=1, check.names=FALSE))
metadata <- read.table(opts$design, header=TRUE, row.names=1, sep="\t",
                       comment.char="", stringsAsFactors=FALSE)

common <- intersect(rownames(dist_matrix), rownames(metadata))
dist_matrix <- dist_matrix[common, common]
metadata <- metadata[common, , drop=FALSE]
dist_obj <- as.dist(dist_matrix)

# ---- 样品名映射: 8_1 → A_1 ----
design_full <- read.table(opts$design, header=TRUE, sep="\t", stringsAsFactors=FALSE)
name_map <- setNames(
  paste0(design_full[[opts$group]], "_", design_full$Replicate),
  design_full$SampleID
)

cat("样本数:", length(common), "\n")

# ==============================================================================
# 2. PCoA
# ==============================================================================
pcoa_res <- cmdscale(dist_obj, k=2, eig=TRUE)
eigenvalues <- pcoa_res$eig
pos_eig <- eigenvalues[eigenvalues > 0]
var_exp <- pos_eig / sum(pos_eig)

pcoa_df <- data.frame(
  PCo1 = pcoa_res$points[, 1],
  PCo2 = pcoa_res$points[, 2],
  Sample = rownames(pcoa_res$points),
  stringsAsFactors = FALSE
)
pcoa_df <- merge(pcoa_df, cbind(Sample=rownames(metadata), metadata), by="Sample")
pcoa_df$Group <- factor(pcoa_df[[opts$group]], levels = c("CK", "A", "B"))

# 样品标签
pcoa_df$Label <- ifelse(pcoa_df$Sample %in% names(name_map),
                        name_map[pcoa_df$Sample], pcoa_df$Sample)

cat("PCo1:", round(var_exp[1]*100, 2), "%\n")
cat("PCo2:", round(var_exp[2]*100, 2), "%\n")

# ==============================================================================
# 3. PERMANOVA
# ==============================================================================
cat("\n--- PERMANOVA ---\n")
dir.create(opts$stat_dir, showWarnings=FALSE, recursive=TRUE)
stat_file <- file.path(opts$stat_dir, "PERMANOVA_results.txt")

sink(stat_file)
cat("=== PERMANOVA分析结果 (组间差异检验) ===\n\n")

set.seed(123)
adonis_res <- adonis2(dist_obj ~ metadata[[opts$group]], permutations=999)
print(adonis_res)

overall_R2 <- adonis_res$R2[1]
overall_P  <- adonis_res$`Pr(>F)`[1]

cat("\n=== 两两比较结果 ===\n")
groups <- c("CK", "A", "B")
for (i in 1:(length(groups)-1)) {
  for (j in (i+1):length(groups)) {
    g1 <- groups[i]; g2 <- groups[j]
    idx <- metadata[[opts$group]] %in% c(g1, g2)
    sub_dist <- as.dist(dist_matrix[idx, idx])
    sub_meta <- metadata[idx, , drop=FALSE]
    pair_res <- adonis2(sub_dist ~ sub_meta[[opts$group]], permutations=999)
    cat(sprintf("%s vs %s: R²=%.3f, p=%.4f\n", g1, g2, pair_res$R2[1], pair_res$`Pr(>F)`[1]))
  }
}

cat("\n=== PCoA解释度 ===\n")
for (k in 1:min(5, length(var_exp))) {
  cat(sprintf("PCo%d: %.2f%%\n", k, var_exp[k]*100))
}
cat(sprintf("前2轴累计: %.2f%%\n", sum(var_exp[1:2])*100))
sink()

cat("PERMANOVA已保存:", stat_file, "\n")
cat(readLines(stat_file), sep="\n")

# ==============================================================================
# 4. 绘图
# ==============================================================================
cat("\n--- 绘图 ---\n")

annot_text <- sprintf("PERMANOVA\nR² = %.3f, P = %.3f", overall_R2, overall_P)

W <- opts$width; H <- opts$height

# ---- 带标签版 ----
p <- ggplot(pcoa_df, aes(x=PCo1, y=PCo2, color=Group)) +
  geom_point(size=3.5, alpha=0.9) +
  geom_text_repel(aes(label=Label), size=2.5, max.overlaps=20,
                  family=FONT, show.legend=FALSE, seed=42,
                  segment.color="grey70", segment.size=0.3) +
  scale_color_manual(values=GROUP_COLORS) +
  labs(
    x = paste0("PCoA1 (", round(var_exp[1]*100, 2), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2]*100, 2), "%)"),
    color = "Group"
  ) +
  # PERMANOVA 右上角
  annotate("text", x=Inf, y=Inf, label=annot_text,
           hjust=1.05, vjust=1.3, size=2.8, family=FONT,
           fontface="italic", color="grey30") +
  theme_classic(base_family=FONT) +
  theme(
    text = element_text(family=FONT),
    axis.text = element_text(color="black", family=FONT),
    axis.title = element_text(family=FONT),
    legend.text = element_text(family=FONT),
    legend.title = element_text(family=FONT, face="bold"),
    legend.position = "right",
    plot.margin = margin(8, 5, 5, 5, "mm")
  )

ggsave(paste0(opts$output, ".pdf"), p, width=W, height=H, units="mm")
ggsave(paste0(opts$output, ".png"), p, width=W, height=H, units="mm", dpi=300)
cat("  →", paste0(opts$output, ".pdf/png\n"))

# ---- 无标签版 (备用，适合小论文) ----
p_clean <- ggplot(pcoa_df, aes(x=PCo1, y=PCo2, color=Group)) +
  geom_point(size=3.5, alpha=0.9) +
  scale_color_manual(values=GROUP_COLORS) +
  labs(
    x = paste0("PCoA1 (", round(var_exp[1]*100, 2), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2]*100, 2), "%)"),
    color = "Group"
  ) +
  annotate("text", x=Inf, y=Inf, label=annot_text,
           hjust=1.05, vjust=1.3, size=2.8, family=FONT,
           fontface="italic", color="grey30") +
  theme_classic(base_family=FONT) +
  theme(
    text = element_text(family=FONT),
    axis.text = element_text(color="black", family=FONT),
    axis.title = element_text(family=FONT),
    legend.text = element_text(family=FONT),
    legend.title = element_text(family=FONT, face="bold"),
    legend.position = "right",
    plot.margin = margin(8, 5, 5, 5, "mm")
  )

ggsave(paste0(opts$output, "_clean.pdf"), p_clean, width=W, height=H, units="mm")
ggsave(paste0(opts$output, "_clean.png"), p_clean, width=W, height=H, units="mm", dpi=300)
cat("  →", paste0(opts$output, "_clean.pdf/png\n"))

cat("\nDone!\n")
