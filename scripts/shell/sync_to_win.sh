#!/bin/bash
# ==============================================================================
# sync_to_win.sh — 同步论文产出到 D:\workdata
# ==============================================================================
#
# 目标结构 (D:\workdata):
#   figures/          ← 全部图表（thesis_cn + paper_en + supplementary）
#   scripts/          ← 全部代码（Python + R + Shell + config）
#   data/             ← 处理后数据 + 表格（不含原始大文件）
#   docs/             ← 文档（README、图表清单、软著材料）
#
# 用法:
#   bash scripts/shell/sync_to_win.sh                  # 同步到 D:\workdata
#   bash scripts/shell/sync_to_win.sh /mnt/e/backup    # 同步到自定义路径
#   bash scripts/shell/sync_to_win.sh --dry-run         # 预览不执行
# ==============================================================================

set -euo pipefail

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# 参数解析
DRY_RUN=""
WIN_TARGET="${THESIS_WIN_TARGET:-/mnt/d/workdata}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run) DRY_RUN="--dry-run"; shift ;;
        *)         WIN_TARGET="$1"; shift ;;
    esac
done

[[ -n "$DRY_RUN" ]] && warn "预览模式，不实际复制"

info "=========================================="
info "同步到 Windows: $WIN_TARGET"
info "源目录: $PROJECT_ROOT"
info "=========================================="

# 检查 D 盘挂载
if [[ "$WIN_TARGET" == /mnt/d* && ! -d "/mnt/d" ]]; then
    echo -e "\033[0;31m[FAIL]\033[0m /mnt/d 不可访问"
    echo "  尝试: sudo mount -t drvfs D: /mnt/d"
    exit 1
fi

# 创建目标结构
mkdir -p "$WIN_TARGET"/{figures/{thesis_cn,paper_en,supplementary},scripts/{python,R,shell,config},data/{raw,processed,tables},docs}

RSYNC_OPTS="-av --delete --exclude='__pycache__' --exclude='*.pyc' --exclude='.git' $DRY_RUN"

sync_dir() {
    local src="$1" dst="$2" desc="$3"
    if [[ -d "$src" ]] && [[ -n "$(ls -A "$src" 2>/dev/null)" ]]; then
        eval rsync $RSYNC_OPTS "$src/" "$dst/" > /dev/null 2>&1
        local count=$(find "$dst" -type f 2>/dev/null | wc -l)
        ok "$desc ($count files)"
    else
        warn "$desc: empty ($src)"
    fi
}

# ---- figures/ ----
info "--- 图表 ---"
sync_dir "$PROJECT_ROOT/figures/thesis_cn"     "$WIN_TARGET/figures/thesis_cn"     "毕业论文图"
sync_dir "$PROJECT_ROOT/figures/paper_en"      "$WIN_TARGET/figures/paper_en"      "小论文图"
sync_dir "$PROJECT_ROOT/figures/supplementary" "$WIN_TARGET/figures/supplementary" "补充图"

# ---- scripts/ ----
info "--- 代码 ---"
sync_dir "$PROJECT_ROOT/scripts/python"  "$WIN_TARGET/scripts/python"  "Python脚本"
sync_dir "$PROJECT_ROOT/scripts/R"       "$WIN_TARGET/scripts/R"       "R脚本"
sync_dir "$PROJECT_ROOT/scripts/shell"   "$WIN_TARGET/scripts/shell"   "Shell脚本"

# config 单独复制到 scripts/config/（让D盘的代码包自洽）
if [[ -d "$PROJECT_ROOT/config" ]]; then
    eval rsync $RSYNC_OPTS "$PROJECT_ROOT/config/" "$WIN_TARGET/scripts/config/" > /dev/null 2>&1
    ok "配置文件 → scripts/config/"
fi

# ---- data/ ----
info "--- 数据 ---"
# 只同步 processed（脚本产出），不同步 raw（大文件，通过软链接指向meta2）
sync_dir "$PROJECT_ROOT/data/processed"  "$WIN_TARGET/data/processed"  "处理后数据"
sync_dir "$PROJECT_ROOT/tables"          "$WIN_TARGET/data/tables"     "表格"

# raw 中只复制小文件（<5MB），跳过大文件和软链接目标
if [[ -d "$PROJECT_ROOT/data/raw" ]]; then
    mkdir -p "$WIN_TARGET/data/raw"
    find "$PROJECT_ROOT/data/raw" -maxdepth 1 -type f -size -5M | while read -r f; do
        cp "$f" "$WIN_TARGET/data/raw/" 2>/dev/null
    done
    # 也复制软链接指向的小文件
    find "$PROJECT_ROOT/data/raw" -maxdepth 1 -type l | while read -r lnk; do
        target=$(readlink -f "$lnk")
        if [[ -f "$target" ]] && [[ $(stat -c%s "$target") -lt 5242880 ]]; then
            cp "$target" "$WIN_TARGET/data/raw/$(basename "$lnk")" 2>/dev/null
        fi
    done
    ok "原始数据(小文件<5MB)"
fi

# ---- docs/ ----
info "--- 文档 ---"
sync_dir "$PROJECT_ROOT/docs" "$WIN_TARGET/docs" "文档"

# 根目录文件
for f in README.md requirements.txt LICENSE; do
    [[ -f "$PROJECT_ROOT/$f" ]] && cp "$PROJECT_ROOT/$f" "$WIN_TARGET/"
done

echo ""
info "=========================================="
info "同步完成！"
info "=========================================="
echo ""

# 统计
if [[ -z "$DRY_RUN" ]]; then
    info "D:\\workdata\\ 内容统计:"
    echo "  figures/   $(find "$WIN_TARGET/figures" -type f 2>/dev/null | wc -l) files"
    echo "  scripts/   $(find "$WIN_TARGET/scripts" -type f 2>/dev/null | wc -l) files"
    echo "  data/      $(find "$WIN_TARGET/data" -type f 2>/dev/null | wc -l) files"
    echo "  docs/      $(find "$WIN_TARGET/docs" -type f 2>/dev/null | wc -l) files"
    echo "  total:     $(du -sh "$WIN_TARGET" | cut -f1)"
    echo ""
    echo "  Windows路径: D:\\workdata\\"
fi
