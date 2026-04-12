#!/bin/bash
# run_MAG_quality.sh — Fig3-1 MAG质量评估散点图
# cd ~/thesis_project && bash scripts/shell/run_MAG_quality.sh

set -euo pipefail
GREEN='\033[0;32m'; BLUE='\033[0;34m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PYSCRIPT="$PROJECT_ROOT/scripts/python/06_MAG_quality.py"

# 数据路径
CHECKM2="$HOME/meta2/binresult/result/checkm2/quality_report.tsv"
BAC120="$HOME/meta2/binresult/temp/gtdb_all/classify/tax.bac120.summary.tsv"
AR53="$HOME/meta2/binresult/temp/gtdb_all/classify/tax.ar53.summary.tsv"
KEYSTONE="$HOME/meta2/binresult/result/keystone_phylogeny/keystone_species.txt"

OUT_PAPER="$PROJECT_ROOT/figures/paper_en"
OUT_THESIS="$PROJECT_ROOT/figures/thesis_cn"
STAT_DIR="$PROJECT_ROOT/data/processed"

mkdir -p "$OUT_PAPER" "$OUT_THESIS" "$STAT_DIR"

for f in "$CHECKM2" "$BAC120" "$AR53" "$KEYSTONE"; do
    [[ ! -f "$f" ]] && echo "ERROR: $f not found" && exit 1
done

info "Fig3-1 MAG Quality Assessment"

info "paper..."
python3 "$PYSCRIPT" \
    --checkm2 "$CHECKM2" --bac120 "$BAC120" --ar53 "$AR53" \
    --keystone "$KEYSTONE" \
    -o "$OUT_PAPER" --stat_dir "$STAT_DIR" --style paper
ok "paper"

info "thesis..."
python3 "$PYSCRIPT" \
    --checkm2 "$CHECKM2" --bac120 "$BAC120" --ar53 "$AR53" \
    --keystone "$KEYSTONE" \
    -o "$OUT_THESIS" --stat_dir "$STAT_DIR" --style thesis
ok "thesis"

echo ""
ok "完成！"
ls "$OUT_PAPER"/Fig3-1*.pdf "$OUT_THESIS"/Fig3-1*.pdf 2>/dev/null | while read f; do echo "  $(basename $(dirname $f))/$(basename $f)"; done
ls "$STAT_DIR"/MAG_quality*.txt 2>/dev/null | while read f; do echo "  $(basename $f)"; done
