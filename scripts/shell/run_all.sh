#!/bin/bash
# ==============================================================================
# run_all.sh — 一键运行分析 + 同步到 D:\workdata
# ==============================================================================
# 用法:
#   bash scripts/shell/run_all.sh                    # 全部运行+同步
#   bash scripts/shell/run_all.sh --skip-analysis    # 只同步
#   bash scripts/shell/run_all.sh --chapter 2        # 只跑第二章
# ==============================================================================

set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
fail() { echo -e "${RED}[FAIL]${NC} $1"; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

LOG_DIR="$PROJECT_ROOT/logs"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="$LOG_DIR/run_${TIMESTAMP}.log"

SKIP_ANALYSIS=false
TARGET_CHAPTER="all"

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-analysis) SKIP_ANALYSIS=true; shift ;;
        --chapter)       TARGET_CHAPTER="$2"; shift 2 ;;
        *) shift ;;
    esac
done

info "=========================================="
info "论文分析 + 同步系统"
info "=========================================="
info "项目: $PROJECT_ROOT"
info "日志: $LOG_FILE"
echo ""

# 环境检查
python3 -c "import matplotlib, numpy, pandas, scipy" 2>/dev/null && ok "Python packages" || \
    { warn "安装缺失包..."; pip install matplotlib numpy pandas scipy seaborn --break-system-packages -q; }

run_script() {
    local script="$1" desc="$2" ext="${1##*.}"
    info "  → $desc ($script)"
    if [[ "$ext" == "py" ]]; then
        PYTHONPATH="$PROJECT_ROOT" python3 "$script" >> "$LOG_FILE" 2>&1 && ok "  $desc" || fail "  $desc (see $LOG_FILE)"
    elif [[ "$ext" == "R" ]]; then
        Rscript "$script" >> "$LOG_FILE" 2>&1 && ok "  $desc" || fail "  $desc (see $LOG_FILE)"
    fi
}

if [[ "$SKIP_ANALYSIS" == false ]]; then
    info "开始运行分析..."

    if [[ "$TARGET_CHAPTER" == "all" || "$TARGET_CHAPTER" == "1" ]]; then
        info "--- 第一章 ---"
        [[ -f scripts/python/01_physicochemical.py ]] && run_script scripts/python/01_physicochemical.py "理化指标" || warn "跳过"
    fi

    if [[ "$TARGET_CHAPTER" == "all" || "$TARGET_CHAPTER" == "2" ]]; then
        info "--- 第二章 ---"
        for s in scripts/R/01_tax_stackplot.R scripts/python/02_alpha_diversity.py \
                 scripts/R/02_beta_PCoA.R scripts/R/03_RDA_Mantel.R \
                 scripts/python/03_gene_heatmap.py scripts/python/04_gene_log2fc.py; do
            [[ -f "$s" ]] && run_script "$s" "$(basename $s)" || warn "$(basename $s) 不存在"
        done
    fi

    if [[ "$TARGET_CHAPTER" == "all" || "$TARGET_CHAPTER" == "3" ]]; then
        info "--- 第三章 ---"
        for s in scripts/python/05_MAG_quality.py scripts/R/05_phylogenetic_tree.R \
                 scripts/python/06_MAG_heatmap.py scripts/python/07_pathway_completeness.py \
                 scripts/python/08_species_contribution.py; do
            [[ -f "$s" ]] && run_script "$s" "$(basename $s)" || warn "$(basename $s) 不存在"
        done
    fi

    ok "分析完成"
    echo ""
fi

# 同步到 D:\workdata
info "同步到 Windows..."
bash "$SCRIPT_DIR/sync_to_win.sh"

ok "全部完成！"
