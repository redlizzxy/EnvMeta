#!/bin/bash
# ==============================================================================
# run_stackplots.sh
# cd ~/thesis_project && bash scripts/shell/run_stackplots.sh
# ==============================================================================

set -euo pipefail

GREEN='\033[0;32m'; BLUE='\033[0;34m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

META2_RESULT="${1:-$HOME/meta2/result}"
INPUT_DIR="$META2_RESULT/metaphlan4"
DESIGN="$META2_RESULT/metadata.txt"

OUT_THESIS="$PROJECT_ROOT/figures/thesis_cn"
OUT_PAPER="$PROJECT_ROOT/figures/paper_en"
PCT_DIR="$PROJECT_ROOT/data/processed"
RSCRIPT="$PROJECT_ROOT/scripts/R/01_tax_stackplot.R"

mkdir -p "$OUT_THESIS" "$OUT_PAPER" "$PCT_DIR"

[[ ! -f "$DESIGN" ]] && echo "ERROR: $DESIGN not found" && exit 1
[[ ! -d "$INPUT_DIR" ]] && echo "ERROR: $INPUT_DIR not found" && exit 1

info "=========================================="
info "堆叠图 (Phylum / Genus / Species)"
info "=========================================="

declare -A FIGNUM
FIGNUM[Phylum]="Fig2-1"; FIGNUM[Genus]="Fig2-2"; FIGNUM[Species]="Fig2-3"

for tax in Phylum Genus Species; do
    INPUT="$INPUT_DIR/${tax}.txt"
    [[ ! -f "$INPUT" ]] && echo "  WARN: $INPUT not found" && continue
    fignum="${FIGNUM[$tax]}"

    # SCI版 (Arial) —— 用原始参数 --width 120 --height 90
    info "${fignum} ${tax} — paper..."
    Rscript "$RSCRIPT" \
        --input "$INPUT" --design "$DESIGN" --group Group \
        --output "$OUT_PAPER/${fignum}_${tax}" \
        --legend 10 --width 120 --height 90 --style paper \
        --percent_dir "$PCT_DIR"
    ok "  paper"

    # 毕业论文版 (SimSun+TNR)
    info "${fignum} ${tax} — thesis..."
    Rscript "$RSCRIPT" \
        --input "$INPUT" --design "$DESIGN" --group Group \
        --output "$OUT_THESIS/${fignum}_${tax}" \
        --legend 10 --width 120 --height 90 --style thesis \
        --percent_dir "$PCT_DIR"
    ok "  thesis"
    echo ""
done

ok "全部完成！"
echo ""
echo "=== 图表 ==="
for dir in "$OUT_PAPER" "$OUT_THESIS"; do
    tag=$(basename "$dir")
    ls "$dir"/Fig2-*.pdf 2>/dev/null | while read f; do echo "  [$tag] $(basename $f)"; done
done
echo ""
echo "=== 百分比数据 ==="
ls "$PCT_DIR"/*percentage*.txt 2>/dev/null | while read f; do echo "  $(basename $f)"; done
