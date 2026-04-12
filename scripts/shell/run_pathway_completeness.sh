#!/bin/bash
# run_pathway_completeness.sh — Fig3-5 通路完整度热图+气泡图
# cd ~/thesis_project && bash scripts/shell/run_pathway_completeness.sh

set -euo pipefail
GREEN='\033[0;32m'; BLUE='\033[0;34m'; RED='\033[0;31m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }
fail() { echo -e "${RED}[FAIL]${NC} $1"; exit 1; }

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PYSCRIPT="$PROJECT_ROOT/scripts/python/08_pathway_completeness.py"

# ============================================================
# 数据路径
# ============================================================
KEGG_TARGET="$HOME/meta2/binresult/result/eggnog/kegg_target_only.tsv"
BAC120="$HOME/meta2/binresult/temp/gtdb_all/classify/tax.bac120.summary.tsv"
AR53="$HOME/meta2/binresult/temp/gtdb_all/classify/tax.ar53.summary.tsv"
ABUNDANCE="$HOME/meta2/binresult/result/coverm/abundance.tsv"
KEYSTONE="$HOME/meta2/binresult/result/keystone_phylogeny/keystone_species.txt"

OUT_DIR="$PROJECT_ROOT/figures"
STAT_DIR="$PROJECT_ROOT/data/processed"

# ============================================================
# 检查文件
# ============================================================
info "检查输入文件..."
for f in "$KEGG_TARGET" "$BAC120" "$AR53" "$ABUNDANCE" "$KEYSTONE"; do
    if [[ ! -f "$f" ]]; then
        fail "文件不存在: $f"
    fi
done
ok "所有输入文件就绪"

# ============================================================
# 创建输出目录
# ============================================================
mkdir -p "$OUT_DIR/paper_en" "$OUT_DIR/thesis_cn" "$STAT_DIR"

# ============================================================
# 运行
# ============================================================
info "Fig3-5 Pathway Completeness"

PYTHONPATH="$PROJECT_ROOT" python3 "$PYSCRIPT" \
    --ko     "$KEGG_TARGET" \
    --bac    "$BAC120" \
    --ar     "$AR53" \
    --abund  "$ABUNDANCE" \
    --ks     "$KEYSTONE" \
    --outdir "$OUT_DIR" \
    --statdir "$STAT_DIR"

ok "完成！"

echo ""
echo "输出文件:"
find "$OUT_DIR" -name "Fig3-5*" -type f | sort | while read f; do
    echo "  $(basename $(dirname $f))/$(basename $f)"
done
ls "$STAT_DIR"/Fig3-5*.* 2>/dev/null | while read f; do
    echo "  data/processed/$(basename $f)"
done
