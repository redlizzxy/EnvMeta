#!/usr/bin/env Rscript

# ==============================================================================
# 03_RDA.R — RDA排序图 + 环境因子显著性检验
# ==============================================================================
# 输入: Species丰度表 + 环境因子表
# 输出: RDA排序图 (箭头加粗, F/P值标注) + 统计结果txt
#
# 用法:
#   Rscript 03_RDA.R \
#     --species ~/meta2/result/metaphlan4/Species.txt \
#     --env ~/meta2/result/elementdata/env_factors.txt \
#     --output ~/thesis_project/figures/paper_en/Fig2-6_RDA \
#     --stat_dir ~/thesis_project/data/processed \
#     --style paper
# ==============================================================================

options(warn = -1)

site <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only=TRUE, quietly=TRUE)))) {
  install.packages("optparse", repos=site); require("optparse", character.only=TRUE)
}

option_list <- list(
  make_option(c("--species"),    type="character", default="metaphlan4/Species.txt"),
  make_option(c("--env"),        type="character", default="elementdata/env_factors.txt"),
  make_option(c("-o", "--output"),type="character", default="Fig2-6_RDA"),
  make_option(c("--stat_dir"),   type="character", default="."),
  make_option(c("-w", "--width"), type="numeric",   default=130),
  make_option(c("-e", "--height"),type="numeric",   default=100),
  make_option(c("-s", "--style"), type="character",  default="paper")
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

# ---- 配色 ----
GROUP_COLORS <- c(CK = "#1c9cbd", A = "#e3943d", B = "#92181e")

# ==============================================================================
# 1. 读取数据
# ==============================================================================

# 物种矩阵 (列=样品, 行=物种)
sp_raw <- read.table(opts$species, header=TRUE, sep="\t", row.names=1,
                     comment.char="", quote="", check.names=FALSE)
cat("物种矩阵:", dim(sp_raw), "\n")

# 去 unclassified
keep <- !grepl("(?i)(unclassif|unassign|unknown)", rownames(sp_raw), perl=TRUE)
sp_raw <- sp_raw[keep, , drop=FALSE]

# 转置: 行=样品, 列=物种
sp <- as.data.frame(t(sp_raw))
cat("转置后:", dim(sp), "(样品 x 物种)\n")

# 环境因子
env_raw <- read.table(opts$env, header=TRUE, sep="\t", stringsAsFactors=FALSE)
cat("环境因子:", dim(env_raw), "\n")

# ---- 样品名映射: Species.txt用2_1格式, env用CK_1格式 ----
# 建立 2_1 → CK_1 的映射
id_map <- c(
  "2_1" = "CK_1", "2_5" = "CK_5", "2_7" = "CK_7",
  "8_1" = "A_1",  "8_2" = "A_2",  "8_3" = "A_3",  "8_4" = "A_4",
  "9_1" = "B_1",  "9_2" = "B_2",  "9_3" = "B_3"
)

# 重命名 sp 行名
rownames(sp) <- ifelse(rownames(sp) %in% names(id_map),
                       id_map[rownames(sp)], rownames(sp))

# env 以 SampleID 为行名
rownames(env_raw) <- env_raw$SampleID

# 对齐
common <- intersect(rownames(sp), rownames(env_raw))
sp <- sp[common, ]
env <- env_raw[common, c("pH", "Eh", "TOC", "Total_As")]
group_vec <- env_raw[common, "Group"]
group_vec <- factor(group_vec, levels = c("CK", "A", "B"))

cat("对齐样本:", length(common), "\n")

# 标准化环境因子
env_scaled <- as.data.frame(scale(env))

# ==============================================================================
# 2. RDA 分析
# ==============================================================================
cat("\n--- RDA ---\n")

# Hellinger 变换物种数据
sp_hell <- decostand(sp, method = "hellinger")

# RDA
rda_res <- rda(sp_hell ~ ., data = env_scaled)

# 解释度
rda_sum <- summary(rda_res)
rda1_pct <- round(rda_sum$cont$importance[2, 1] * 100, 2)
rda2_pct <- round(rda_sum$cont$importance[2, 2] * 100, 2)
cat("RDA1:", rda1_pct, "%  RDA2:", rda2_pct, "%\n")

# 总体显著性
set.seed(123)
anova_overall <- anova.cca(rda_res, permutations = 999)
cat("Overall: F =", round(anova_overall$F[1], 3), " P =", anova_overall$`Pr(>F)`[1], "\n")

# 各环境因子显著性 (逐项检验)
set.seed(123)
anova_terms <- anova.cca(rda_res, by = "terms", permutations = 999)
cat("\n各因子检验:\n")
print(anova_terms)

# ==============================================================================
# 3. 统计结果保存
# ==============================================================================
dir.create(opts$stat_dir, showWarnings=FALSE, recursive=TRUE)
stat_file <- file.path(opts$stat_dir, "RDA_results.txt")

sink(stat_file)
cat("=== RDA分析结果 ===\n")
cat("生成时间:", format(Sys.time()), "\n\n")
cat("解释度: RDA1 =", rda1_pct, "%  RDA2 =", rda2_pct, "%\n\n")
cat("=== 总体显著性 (ANOVA) ===\n")
print(anova_overall)
cat("\n=== 各环境因子显著性 (by terms) ===\n")
print(anova_terms)

# Mantel 检验
cat("\n=== Mantel检验 ===\n")
sp_dist <- vegdist(sp_hell, method = "bray")
for (env_name in colnames(env)) {
  env_dist <- dist(scale(env[, env_name, drop=FALSE]))
  mt <- mantel(sp_dist, env_dist, permutations = 999)
  cat(sprintf("%s: r = %.3f, p = %.4f\n", env_name, mt$statistic, mt$signif))
}
sink()

cat("统计已保存:", stat_file, "\n")

# ==============================================================================
# 4. 绘图
# ==============================================================================
cat("\n--- 绘图 ---\n")

# 提取坐标
site_scores <- as.data.frame(scores(rda_res, display = "sites", choices = 1:2))
site_scores$Sample <- rownames(site_scores)
site_scores$Group <- group_vec

# 环境因子箭头坐标
env_scores <- as.data.frame(scores(rda_res, display = "bp", choices = 1:2))
env_scores$Factor <- rownames(env_scores)

# 提取各因子的F值和P值（用于标注）
f_vals <- anova_terms$F[1:nrow(env_scores)]
p_vals <- anova_terms$`Pr(>F)`[1:nrow(env_scores)]

# 构建标注文本
env_scores$label <- paste0(
  env_scores$Factor,
  "\n(F=", round(f_vals, 2),
  ", P=", ifelse(p_vals < 0.001, "<0.001",
          ifelse(p_vals < 0.01, sprintf("%.3f", p_vals),
                 sprintf("%.3f", p_vals))), ")"
)

# 标注中加星号
env_scores$sig <- ifelse(p_vals < 0.001, "***",
                  ifelse(p_vals < 0.01, "**",
                  ifelse(p_vals < 0.05, "*", "")))

env_scores$label_short <- paste0(env_scores$Factor, " ", env_scores$sig)

# 箭头缩放因子（让箭头长度适配散点分布）
arrow_scale <- max(abs(site_scores[, 1:2])) / max(abs(env_scores[, 1:2])) * 0.85

# 绘图
p <- ggplot() +
  # 样品散点
  geom_point(data = site_scores, aes(x = RDA1, y = RDA2, color = Group),
             size = 3.5, alpha = 0.9) +
  # 样品标签
  geom_text_repel(data = site_scores, aes(x = RDA1, y = RDA2, label = Sample, color = Group),
                  size = 2.2, family = FONT, show.legend = FALSE,
                  max.overlaps = 20, seed = 42,
                  segment.color = "grey70", segment.size = 0.2) +
  # 环境因子箭头（加粗）
  geom_segment(data = env_scores,
               aes(x = 0, y = 0, xend = RDA1 * arrow_scale, yend = RDA2 * arrow_scale),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               linewidth = 0.9, color = "#333333") +
  # 环境因子标签（带 F/P 值）
  geom_text_repel(data = env_scores,
                  aes(x = RDA1 * arrow_scale * 1.12,
                      y = RDA2 * arrow_scale * 1.12,
                      label = label_short),
                  size = 3, family = FONT, fontface = "bold",
                  color = "#333333", seed = 42,
                  box.padding = 0.5, point.padding = 0.3,
                  segment.color = NA) +
  # 配色
  scale_color_manual(values = GROUP_COLORS) +
  # 轴标签
  labs(
    x = paste0("RDA1 (", rda1_pct, "%)"),
    y = paste0("RDA2 (", rda2_pct, "%)"),
    color = "Group"
  ) +
  # 参考线
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.3) +
  # 主题
  theme_classic(base_family = FONT) +
  theme(
    text = element_text(family = FONT),
    axis.text = element_text(color = "black", family = FONT),
    axis.title = element_text(family = FONT),
    legend.text = element_text(family = FONT),
    legend.title = element_text(family = FONT, face = "bold"),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5, "mm")
  )

W <- opts$width; H <- opts$height

ggsave(paste0(opts$output, ".pdf"), p, width=W, height=H, units="mm")
ggsave(paste0(opts$output, ".png"), p, width=W, height=H, units="mm", dpi=300)
cat("  ->", paste0(opts$output, ".pdf/png\n"))

# 带详细F/P标注的版本
p_detail <- p
# 替换简短标签为详细标签
p_detail$layers[[4]] <- NULL  # 去掉短标签层
p_detail <- p_detail +
  geom_label_repel(data = env_scores,
                   aes(x = RDA1 * arrow_scale * 1.12,
                       y = RDA2 * arrow_scale * 1.12,
                       label = label),
                   size = 2.3, family = FONT, fontface = "bold",
                   color = "#333333", fill = "white", alpha = 0.85,
                   label.size = 0.2, seed = 42,
                   box.padding = 0.5, segment.color = NA)

ggsave(paste0(opts$output, "_detail.pdf"), p_detail, width=W*1.15, height=H, units="mm")
ggsave(paste0(opts$output, "_detail.png"), p_detail, width=W*1.15, height=H, units="mm", dpi=300)
cat("  ->", paste0(opts$output, "_detail.pdf/png\n"))

cat("\nDone!\n")
