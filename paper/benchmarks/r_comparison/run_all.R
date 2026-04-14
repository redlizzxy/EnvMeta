# ==============================================================================
# run_all.R — 跑原 R 脚本生成 EnvMeta 对照所需的参考 PDF
# ==============================================================================
# 在 RStudio 中：
#   setwd("d:/workdata/envmeta")
#   source("paper/benchmarks/r_comparison/run_all.R")
# 或者 Terminal：
#   cd d:/workdata/envmeta
#   Rscript paper/benchmarks/r_comparison/run_all.R
# ==============================================================================

setwd("d:/workdata/envmeta")

OUT <- "paper/benchmarks/r_comparison"
DATA <- "data/raw"
META <- file.path(DATA, "metadata.txt")

# ---- 1. 堆叠图（Genus 和 Species 两版）----
cat("\n>>> 1. 堆叠图 Genus\n")
system2("Rscript", c(
  "scripts/R/01_tax_stackplot.R",
  "--input",  file.path(DATA, "Genus.txt"),
  "--design", META,
  "--output", file.path(OUT, "stackplot/Fig2-3_Genus"),
  "--legend", "10",
  "--style",  "paper"
))

cat("\n>>> 1. 堆叠图 Species\n")
system2("Rscript", c(
  "scripts/R/01_tax_stackplot.R",
  "--input",  file.path(DATA, "Species.txt"),
  "--design", META,
  "--output", file.path(OUT, "stackplot/Fig2-3_Species"),
  "--legend", "10",
  "--style",  "paper"
))

# ---- 2. α 多样性 ----
cat("\n>>> 2. α 多样性\n")
system2("Rscript", c(
  "scripts/R/02_alpha_diversity.R",
  "--kraken2",    file.path(DATA, "tax_count.alpha"),
  "--metaphlan4", file.path(DATA, "alpha.txt"),
  "--design",     META,
  "--output",     file.path(OUT, "alpha/Fig2-4"),
  "--stat_dir",   file.path(OUT, "alpha"),
  "--style",      "paper"
))

# ---- 3. β 多样性 PCoA ----
cat("\n>>> 3. PCoA\n")
system2("Rscript", c(
  "scripts/R/02_beta_PCoA.R",
  "--dist",     file.path(DATA, "beta_bray.txt"),
  "--design",   META,
  "--output",   file.path(OUT, "pcoa/Fig2-5_PCoA"),
  "--stat_dir", file.path(OUT, "pcoa"),
  "--style",    "paper"
))

# ---- 4. RDA ----
cat("\n>>> 4. RDA\n")
system2("Rscript", c(
  "scripts/R/03_RDA.R",
  "--species",  file.path(DATA, "Species.txt"),
  "--env",      file.path(DATA, "env_factors.txt"),
  "--output",   file.path(OUT, "rda/Fig2-6_RDA"),
  "--stat_dir", file.path(OUT, "rda"),
  "--style",    "paper"
))

# ---- 5. LEfSe（从 input.res 绘图）----
cat("\n>>> 5. LEfSe\n")
system2("Rscript", c(
  "scripts/R/04_LEfSe.R",
  "--input",    file.path(DATA, "input.res"),
  "--output",   file.path(OUT, "lefse/Fig2-7"),
  "--stat_dir", file.path(OUT, "lefse"),
  "--style",    "paper"
))

cat("\n=== 全部完成 ===\n")
cat("PDF 位置：", normalizePath(OUT), "\n")
