# ==============================================================================
# 论文绘图统一R主题配置 (theme_thesis.R) — 最终版 2026-03-12
# source("config/theme_thesis.R")
# ==============================================================================

library(ggplot2)
library(RColorBrewer)

GROUPS <- list(CK = c("2_1","2_5","2_7"), A = c("8_1","8_2","8_3","8_4"), B = c("9_1","9_2","9_3"))
GROUP_ORDER <- c("CK", "A", "B")

# 三组核心配色 (方案B Earth-Warm, 用于柱状图/PCoA/元素循环)
GROUP_COLORS <- c(CK = "#1c9cbd", A = "#e3943d", B = "#92181e")
GROUP_COLORS_LIGHT <- c(CK = "#d4e5ef", A = "#fbd8b5", B = "#d0916e")

# 元素循环配色
ELEMENT_COLORS <- c(arsenic = "#92181e", sulfur = "#e3943d", iron = "#a35626", nitrogen = "#1c9cbd")

# 堆叠图12色 (方案A Ocean-Forest, 重排: 高丰度→浅色, 低丰度→深色)
STACK_12 <- c(
  "#bbc4e4", "#d7e7af", "#bf9d6d", "#a4d6c1", "#bce2e8", "#b5b5b6",
  "#0084c2", "#3ab483", "#00978c", "#21825e", "#6aa3d2", "#3d62ad"
)
PALETTE_12 <- STACK_12

# 热图色阶 (方案D)
# 用法: scale_fill_distiller(palette = "RdYlBu", direction = -1)
HEATMAP_PALETTE <- "RdYlBu"

# ---- 毕业论文主题: 宋体(中文) + TNR(英文/数字) ----
theme_thesis <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = "serif") + theme(
    text = element_text(family = "serif"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black", family = "serif"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.background = element_blank(), legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(5, 5, 5, 5, "mm"))
}

# ---- SCI小论文主题: Arial (sans-serif) ----
theme_paper <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "sans") + theme(
    text = element_text(family = "sans"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black", size = 7, family = "sans"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.background = element_blank(), legend.key = element_blank(),
    legend.text = element_text(size = 6), legend.title = element_text(size = 7),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(5, 5, 5, 5, "mm"))
}

save_fig <- function(p, fig_id, desc_cn = "", desc_en = "",
                     thesis_w = 180, thesis_h = 120,
                     paper_w = 90, paper_h = 70, fig_dir = ".") {
  thesis_path <- file.path(fig_dir, "figures/thesis_cn", paste0(fig_id, "_", desc_cn, ".pdf"))
  ggsave(thesis_path, p, width = thesis_w, height = thesis_h, units = "mm", dpi = 300)
  paper_path <- file.path(fig_dir, "figures/paper_en", paste0(fig_id, "_", desc_en, ".pdf"))
  ggsave(paper_path, p, width = paper_w, height = paper_h, units = "mm", dpi = 600)
}

sig_label <- function(p) ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns")))

cat("theme_thesis.R loaded (final 2026-03-12)\n")
