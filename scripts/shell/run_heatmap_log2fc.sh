#!/bin/bash
# ==============================================================================
# run_heatmap_log2fc.sh — Fig2-8 热图 + Fig2-9 log2FC
# cd ~/thesis_project && bash scripts/shell/run_heatmap_log2fc.sh
# ==============================================================================

set -euo pipefail

GREEN='\033[0;32m'; BLUE='\033[0;34m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

META2="${1:-$HOME/meta2/result}"
SPF="$META2/eggnog/eggnog.KEGG_ko.TPM.spf"
DESIGN="$META2/metadata.txt"

OUT_PAPER="$PROJECT_ROOT/figures/paper_en"
OUT_THESIS="$PROJECT_ROOT/figures/thesis_cn"
STAT_DIR="$PROJECT_ROOT/data/processed"
PYSCRIPT="$PROJECT_ROOT/scripts/python/05_gene_heatmap_log2fc.py"

mkdir -p "$OUT_PAPER" "$OUT_THESIS" "$STAT_DIR"

[[ ! -f "$SPF" ]] && echo "ERROR: $SPF not found" && exit 1
[[ ! -f "$DESIGN" ]] && echo "ERROR: $DESIGN not found" && exit 1

info "=========================================="
info "Fig2-8 热图 + Fig2-9 log2FC"
info "=========================================="

# SCI版
info "paper..."
python3 "$PYSCRIPT" \
    -i "$SPF" --metadata "$DESIGN" \
    -o "$OUT_PAPER" --stat_dir "$STAT_DIR" \
    --style paper
ok "  paper"

# 毕业论文版
info "thesis..."
python3 "$PYSCRIPT" \
    -i "$SPF" --metadata "$DESIGN" \
    -o "$OUT_THESIS" --stat_dir "$STAT_DIR" \
    --style thesis
ok "  thesis"

echo ""
ok "完成！"
echo ""
echo "=== 图表 ==="
for dir in "$OUT_PAPER" "$OUT_THESIS"; do
    tag=$(basename "$dir")
    ls "$dir"/Fig2-8*.pdf "$dir"/Fig2-9*.pdf "$dir"/FigS*.pdf 2>/dev/null | \
        while read f; do echo "  [$tag] $(basename $f)"; done
done
echo ""
echo "=== 统计 ==="
ls "$STAT_DIR"/gene_log2fc*.txt 2>/dev/null | while read f; do echo "  $(basename $f)"; done
