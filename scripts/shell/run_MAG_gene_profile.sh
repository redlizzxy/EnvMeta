#!/bin/bash
# run_MAG_gene_profile.sh — Fig3-3 MAG元素循环基因谱热图
# cd ~/thesis_project && bash scripts/shell/run_MAG_gene_profile.sh

set -euo pipefail
GREEN='\033[0;32m'; BLUE='\033[0;34m'; RED='\033[0;31m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }
fail() { echo -e "${RED}[FAIL]${NC} $1"; exit 1; }

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PYSCRIPT="$PROJECT_ROOT/scripts/python/06_MAG_gene_profile.py"

# ============================================================
# 数据路径（按WSL实际路径配置）
# ============================================================
ITOL_HEATMAP="$HOME/meta2/binresult/result/itol/07_element_KO_heatmap.txt"
KEGG_TARGET="$HOME/meta2/binresult/result/eggnog/kegg_target_only.tsv"
BAC120="$HOME/meta2/binresult/temp/gtdb_all/classify/tax.bac120.summary.tsv"
AR53="$HOME/meta2/binresult/temp/gtdb_all/classify/tax.ar53.summary.tsv"
ABUNDANCE="$HOME/meta2/binresult/result/coverm/abundance.tsv"
KEYSTONE="$HOME/meta2/binresult/result/keystone_phylogeny/keystone_species.txt"

OUT_DIR="$PROJECT_ROOT/figures"

# ============================================================
# 检查文件
# ============================================================
info "检查输入文件..."
for f in "$ITOL_HEATMAP" "$KEGG_TARGET" "$BAC120" "$AR53" "$ABUNDANCE" "$KEYSTONE"; do
    if [[ ! -f "$f" ]]; then
        fail "文件不存在: $f"
    fi
done
ok "所有输入文件就绪"

# ============================================================
# 创建输出目录
# ============================================================
mkdir -p "$OUT_DIR/paper_en" "$OUT_DIR/thesis_cn"

# ============================================================
# 运行
# ============================================================
info "Fig3-3 MAG Element Cycle Gene Profiles"

info "生成 EN + CN 版本..."
PYTHONPATH="$PROJECT_ROOT" python3 "$PYSCRIPT" \
    --itol   "$ITOL_HEATMAP" \
    --ko     "$KEGG_TARGET" \
    --bac    "$BAC120" \
    --ar     "$AR53" \
    --abund  "$ABUNDANCE" \
    --ks     "$KEYSTONE" \
    --outdir "$OUT_DIR"

ok "完成！"

echo ""
echo "输出文件:"
ls "$OUT_DIR"/paper_en/Fig3-3*.pdf "$OUT_DIR"/paper_en/Fig3-3*.tiff \
   "$OUT_DIR"/thesis_cn/Fig3-3*.pdf 2>/dev/null | while read f; do
    echo "  $(basename $(dirname $f))/$(basename $f)"
done
