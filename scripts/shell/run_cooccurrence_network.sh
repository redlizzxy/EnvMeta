#!/bin/bash
# run_cooccurrence_network.sh — Fig3-7 共现网络图
# cd ~/thesis_project && bash scripts/shell/run_cooccurrence_network.sh

set -euo pipefail
GREEN='\033[0;32m'; BLUE='\033[0;34m'; RED='\033[0;31m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }
fail() { echo -e "${RED}[FAIL]${NC} $1"; exit 1; }

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PYSCRIPT="$PROJECT_ROOT/scripts/python/09_cooccurrence_network.py"

# ============================================================
# 数据路径
# ============================================================
ABUNDANCE="$HOME/meta2/binresult/result/coverm/abundance.tsv"
BAC120="$HOME/meta2/binresult/temp/gtdb_all/classify/tax.bac120.summary.tsv"
AR53="$HOME/meta2/binresult/temp/gtdb_all/classify/tax.ar53.summary.tsv"
KEYSTONE="$HOME/meta2/binresult/result/keystone_phylogeny/keystone_species.txt"
NODES="$HOME/meta2/binresult/result/downstream_analysis/network_node_stats.txt"
MODULES="$HOME/meta2/binresult/result/downstream_analysis/network_module_info.txt"

OUT_DIR="$PROJECT_ROOT/figures"
STAT_DIR="$PROJECT_ROOT/data/processed"

# ============================================================
# 检查文件
# ============================================================
info "检查输入文件..."
for f in "$ABUNDANCE" "$BAC120" "$AR53" "$KEYSTONE"; do
    if [[ ! -f "$f" ]]; then
        fail "文件不存在: $f"
    fi
done

# nodes和modules可能在不同路径，尝试查找
for f in "$NODES" "$MODULES"; do
    if [[ ! -f "$f" ]]; then
        # 尝试thesis_project下
        alt="$PROJECT_ROOT/data/raw/$(basename $f)"
        if [[ -f "$alt" ]]; then
            info "使用替代路径: $alt"
        else
            fail "文件不存在: $f (也不在 $alt)"
        fi
    fi
done
ok "所有输入文件就绪"

# ============================================================
# 安装依赖
# ============================================================
pip install networkx --break-system-packages -q 2>/dev/null || true

# ============================================================
# 运行
# ============================================================
mkdir -p "$OUT_DIR/paper_en" "$OUT_DIR/thesis_cn" "$STAT_DIR"

info "Fig3-7 Co-occurrence Network"

PYTHONPATH="$PROJECT_ROOT" python3 "$PYSCRIPT" \
    --abund  "$ABUNDANCE" \
    --bac    "$BAC120" \
    --ar     "$AR53" \
    --ks     "$KEYSTONE" \
    --nodes  "$NODES" \
    --modules "$MODULES" \
    --outdir "$OUT_DIR" \
    --statdir "$STAT_DIR"

ok "完成！"

echo ""
echo "输出文件:"
find "$OUT_DIR" -name "Fig3-7*" -type f | sort | while read f; do
    echo "  $(basename $(dirname $f))/$(basename $f)"
done
ls "$STAT_DIR"/Fig3-7*.* 2>/dev/null | while read f; do
    echo "  data/processed/$(basename $f)"
done
