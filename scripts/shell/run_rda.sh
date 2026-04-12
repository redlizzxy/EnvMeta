#!/bin/bash
# ==============================================================================
# run_rda.sh — RDA排序图 + 环境因子检验
# cd ~/thesis_project && bash scripts/shell/run_rda.sh
# ==============================================================================

set -euo pipefail

GREEN='\033[0;32m'; BLUE='\033[0;34m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

META2="${1:-$HOME/meta2/result}"
SPECIES="$META2/metaphlan4/Species.txt"
ENV="$META2/elementdata/env_factors.txt"

OUT_PAPER="$PROJECT_ROOT/figures/paper_en"
OUT_THESIS="$PROJECT_ROOT/figures/thesis_cn"
STAT_DIR="$PROJECT_ROOT/data/processed"
RSCRIPT="$PROJECT_ROOT/scripts/R/03_RDA.R"

mkdir -p "$OUT_PAPER" "$OUT_THESIS" "$STAT_DIR"

[[ ! -f "$SPECIES" ]] && echo "ERROR: $SPECIES not found" && exit 1
[[ ! -f "$ENV" ]]     && echo "ERROR: $ENV not found" && exit 1

info "=========================================="
info "RDA排序图"
info "=========================================="

# SCI版
info "Fig2-6 — paper..."
Rscript "$RSCRIPT" \
    --species "$SPECIES" --env "$ENV" \
    --output "$OUT_PAPER/Fig2-6_RDA" \
    --stat_dir "$STAT_DIR" \
    --width 130 --height 100 --style paper
ok "  paper"

# 毕业论文版
info "Fig2-6 — thesis..."
Rscript "$RSCRIPT" \
    --species "$SPECIES" --env "$ENV" \
    --output "$OUT_THESIS/Fig2-6_RDA" \
    --stat_dir "$STAT_DIR" \
    --width 130 --height 100 --style thesis
ok "  thesis"

echo ""
ok "完成！"
echo ""
echo "=== 图表 ==="
ls "$OUT_PAPER"/Fig2-6*.pdf "$OUT_THESIS"/Fig2-6*.pdf 2>/dev/null | while read f; do echo "  $f"; done
echo ""
echo "=== 统计 ==="
[[ -f "$STAT_DIR/RDA_results.txt" ]] && cat "$STAT_DIR/RDA_results.txt"
