# R 包安装脚本
# Rscript scripts/R/install_packages.R

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))

pkgs <- c("ggplot2", "vegan", "pheatmap", "RColorBrewer", "reshape2",
           "tidyverse", "ape", "igraph", "ggraph")

bioc_pkgs <- c("ggtree", "ggtreeExtra")

# CRAN 包
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
    cat("Installed:", p, "\n")
  } else {
    cat("Already installed:", p, "\n")
  }
}

# Bioconductor 包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE)
    cat("Installed:", p, "\n")
  } else {
    cat("Already installed:", p, "\n")
  }
}

cat("\n✓ 所有R包安装完成\n")
