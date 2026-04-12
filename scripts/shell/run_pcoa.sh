#!/bin/bash
# ==============================================================================
# run_pcoa.sh — PCoA + PERMANOVA
# cd ~/thesis_project && bash scripts/shell/run_pcoa.sh
# ==============================================================================

set -euo pipefail

GREEN='\033[0;32m'; BLUE='\033[0;34m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

META2_RESULT="${1:-$HOME/meta2/result}"
DIST="$META2_RESULT/metaphlan4/beta_bray.txt"
DESIGN="$META2_RESULT/metadata.txt"

OUT_PAPER="$PROJECT_ROOT/figures/paper_en"
OUT_THESIS="$PROJECT_ROOT/figures/thesis_cn"
STAT_DIR="$PROJECT_ROOT/data/processed"
RSCRIPT="$PROJECT_ROOT/scripts/R/02_beta_PCoA.R"

mkdir -p "$OUT_PAPER" "$OUT_THESIS" "$STAT_DIR"

[[ ! -f "$DIST" ]] && echo "ERROR: $DIST not found" && exit 1
[[ ! -f "$DESIGN" ]] && echo "ERROR: $DESIGN not found" && exit 1

info "=========================================="
info "PCoA + PERMANOVA"
info "=========================================="

# SCI版 (Arial)
info "Fig2-5 — paper..."
Rscript "$RSCRIPT" \
    --dist "$DIST" --design "$DESIGN" --group Group \
    --output "$OUT_PAPER/Fig2-5_PCoA" \
    --stat_dir "$STAT_DIR" \
    --width 120 --height 90 --style paper
ok "  paper"

# 毕业论文版 (SimSun+TNR)
info "Fig2-5 — thesis..."
Rscript "$RSCRIPT" \
    --dist "$DIST" --design "$DESIGN" --group Group \
    --output "$OUT_THESIS/Fig2-5_PCoA" \
    --stat_dir "$STAT_DIR" \
    --width 120 --height 90 --style thesis
ok "  thesis"

echo ""
ok "完成！"
echo ""
echo "=== 图表 ==="
ls "$OUT_PAPER"/Fig2-5*.pdf "$OUT_THESIS"/Fig2-5*.pdf 2>/dev/null | while read f; do echo "  $f"; done
echo ""
echo "=== 统计结果 ==="
ls "$STAT_DIR"/PERMANOVA*.txt 2>/dev/null | while read f; do echo "  $f"; done
