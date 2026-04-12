#!/bin/bash
# ==============================================================================
# run_alpha.sh — α多样性 (Kraken2 + MetaPhlAn4)
# cd ~/thesis_project && bash scripts/shell/run_alpha.sh
# ==============================================================================

set -euo pipefail

GREEN='\033[0;32m'; BLUE='\033[0;34m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

META2="${1:-$HOME/meta2/result}"
KR2_ALPHA="$META2/kraken2/tax_count.alpha"
MP4_ALPHA="$META2/metaphlan4/alpha.txt"
DESIGN="$META2/metadata.txt"

OUT_PAPER="$PROJECT_ROOT/figures/paper_en"
OUT_THESIS="$PROJECT_ROOT/figures/thesis_cn"
STAT_DIR="$PROJECT_ROOT/data/processed"
RSCRIPT="$PROJECT_ROOT/scripts/R/02_alpha_diversity.R"

mkdir -p "$OUT_PAPER" "$OUT_THESIS" "$STAT_DIR"

[[ ! -f "$KR2_ALPHA" ]] && echo "ERROR: $KR2_ALPHA not found" && exit 1
[[ ! -f "$MP4_ALPHA" ]] && echo "ERROR: $MP4_ALPHA not found" && exit 1
[[ ! -f "$DESIGN" ]]    && echo "ERROR: $DESIGN not found" && exit 1

info "=========================================="
info "α多样性 (Kraken2 + MetaPhlAn4)"
info "=========================================="

# SCI版 (Arial)
info "Fig2-4 — paper..."
Rscript "$RSCRIPT" \
    --kraken2 "$KR2_ALPHA" --metaphlan4 "$MP4_ALPHA" \
    --design "$DESIGN" --group Group \
    --output "$OUT_PAPER/Fig2-4" \
    --stat_dir "$STAT_DIR" \
    --width 89 --height 75 --style paper
ok "  paper"

# 毕业论文版 (SimSun+TNR)
info "Fig2-4 — thesis..."
Rscript "$RSCRIPT" \
    --kraken2 "$KR2_ALPHA" --metaphlan4 "$MP4_ALPHA" \
    --design "$DESIGN" --group Group \
    --output "$OUT_THESIS/Fig2-4" \
    --stat_dir "$STAT_DIR" \
    --width 89 --height 75 --style thesis
ok "  thesis"

echo ""
ok "完成！"
echo ""
echo "=== 图表 ==="
for dir in "$OUT_PAPER" "$OUT_THESIS"; do
    ls "$dir"/Fig2-4*.pdf 2>/dev/null | while read f; do echo "  $(basename $f)"; done
done
echo ""
echo "=== 统计 ==="
[[ -f "$STAT_DIR/alpha_wilcoxon_results.txt" ]] && head -30 "$STAT_DIR/alpha_wilcoxon_results.txt"
