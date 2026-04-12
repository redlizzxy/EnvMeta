#!/bin/bash
# ==============================================================================
# run_lefse.sh — LEfSe LDA条形图 (基于已有的input.res)
# cd ~/thesis_project && bash scripts/shell/run_lefse.sh
# ==============================================================================

set -euo pipefail

GREEN='\033[0;32m'; BLUE='\033[0;34m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

META2="${1:-$HOME/meta2/result}"
INPUT_RES="$META2/temp/input.res"

# 如果temp/input.res不存在，尝试其他路径
if [[ ! -f "$INPUT_RES" ]]; then
    INPUT_RES="$META2/metaphlan4/input.res"
fi
if [[ ! -f "$INPUT_RES" ]]; then
    # 也检查result上级目录
    INPUT_RES="$HOME/meta2/temp/input.res"
fi

OUT_PAPER="$PROJECT_ROOT/figures/paper_en"
OUT_THESIS="$PROJECT_ROOT/figures/thesis_cn"
STAT_DIR="$PROJECT_ROOT/data/processed"
RSCRIPT="$PROJECT_ROOT/scripts/R/04_LEfSe.R"

mkdir -p "$OUT_PAPER" "$OUT_THESIS" "$STAT_DIR"

[[ ! -f "$INPUT_RES" ]] && echo "ERROR: input.res not found. 请指定路径: bash run_lefse.sh /path/to/meta2/result" && exit 1

info "=========================================="
info "LEfSe LDA条形图"
info "输入: $INPUT_RES"
info "=========================================="

# SCI版
info "Fig2-7 — paper..."
Rscript "$RSCRIPT" \
    --input "$INPUT_RES" \
    --output "$OUT_PAPER/Fig2-7" \
    --stat_dir "$STAT_DIR" \
    --width 180 --height 200 --style paper
ok "  paper"

# 毕业论文版
info "Fig2-7 — thesis..."
Rscript "$RSCRIPT" \
    --input "$INPUT_RES" \
    --output "$OUT_THESIS/Fig2-7" \
    --stat_dir "$STAT_DIR" \
    --width 180 --height 200 --style thesis
ok "  thesis"

# 复制原始Cladogram PDF（LEfSe原图直接用）
CLADO="$META2/metaphlan4/lefse_cladogram.pdf"
if [[ -f "$CLADO" ]]; then
    cp "$CLADO" "$OUT_PAPER/Fig2-7b_cladogram.pdf"
    cp "$CLADO" "$OUT_THESIS/Fig2-7b_cladogram.pdf"
    ok "Cladogram copied"
else
    echo "  WARN: lefse_cladogram.pdf not found at $CLADO"
fi

echo ""
ok "完成！"
echo ""
echo "=== 图表 ==="
ls "$OUT_PAPER"/Fig2-7*.pdf "$OUT_THESIS"/Fig2-7*.pdf 2>/dev/null | while read f; do echo "  $f"; done
echo ""
echo "=== 统计 ==="
[[ -f "$STAT_DIR/LEfSe_significant_features.txt" ]] && \
    echo "  $(wc -l < "$STAT_DIR/LEfSe_significant_features.txt") features saved"
