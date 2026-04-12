#!/bin/bash
# ==============================================================================
# init_project.sh — 初始化 thesis_project
# ==============================================================================
#
# 功能:
#   1. 创建完整目录结构
#   2. 建立 data/raw/ → ~/meta2/ 子目录的软链接
#   3. 检查Python/R环境
#   4. 验证D盘可写
#
# 用法:
#   cd ~/thesis_project
#   bash scripts/shell/init_project.sh
#
# 软链接说明:
#   脚本会读取下方 LINKS 数组，为每个条目在 data/raw/ 下创建软链接。
#   格式: "链接名|源路径"
#   源路径支持相对于 ~/meta2/ 的路径，也支持绝对路径。
#
#   ★ 你需要根据自己的 meta2 目录结构修改下方 LINKS 数组 ★
#
# ==============================================================================

set -euo pipefail

GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
fail() { echo -e "${RED}[FAIL]${NC} $1"; }

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
META2="$HOME/meta2"
WIN_TARGET="/mnt/d/workdata"

echo ""
echo "========================================================"
echo "  thesis_project 初始化"
echo "========================================================"
echo "  项目目录:  $PROJECT_ROOT"
echo "  数据源:    $META2"
echo "  输出目标:  $WIN_TARGET (D:\\workdata)"
echo "========================================================"
echo ""

# ============================================================
# 1. 创建目录结构
# ============================================================
info "创建目录结构..."

dirs=(
    config
    scripts/python scripts/R scripts/shell
    data/raw data/processed
    figures/thesis_cn figures/paper_en figures/supplementary
    tables docs docs/copyright logs
)

for d in "${dirs[@]}"; do
    mkdir -p "$PROJECT_ROOT/$d"
done
ok "目录结构已创建"

# ============================================================
# 2. 建立软链接 (data/raw/ → meta2/)
# ============================================================
info "建立数据软链接..."

if [[ ! -d "$META2" ]]; then
    fail "$META2 不存在！请确认路径"
    exit 1
fi

# ================================================================
# ★★★ 在这里配置你的软链接 ★★★
#
# 格式: "链接名|meta2下的相对路径"
#
# 例如:
#   "kraken2|kraken2"
#     → 创建 data/raw/kraken2 → ~/meta2/kraken2
#
#   "metadata.txt|result/metadata.txt"
#     → 创建 data/raw/metadata.txt → ~/meta2/result/metadata.txt
#
# 等你确认 meta2 的目录结构后，在下方添加条目即可。
# ================================================================

LINKS=(
    # "链接名|meta2下的相对路径"
    "elementdata|result/elementdata"
    # ---------- 后续添加更多条目 ----------
    # "kraken2|kraken2"
    # "binning_result|binning/result"
    # "metadata.txt|result/metadata.txt"
)

cd "$PROJECT_ROOT/data/raw"

link_count=0
for entry in "${LINKS[@]}"; do
    # 跳过注释行和空行
    [[ "$entry" =~ ^#.*$ || -z "$entry" ]] && continue

    link_name="${entry%%|*}"
    source_rel="${entry##*|}"
    source_abs="$META2/$source_rel"

    if [[ -e "$source_abs" ]]; then
        # 如果软链接已存在且指向正确，跳过
        if [[ -L "$link_name" ]] && [[ "$(readlink -f "$link_name")" == "$(readlink -f "$source_abs")" ]]; then
            ok "  $link_name → $source_rel (已存在)"
        else
            ln -sfn "$source_abs" "$link_name"
            ok "  $link_name → $source_rel"
        fi
        link_count=$((link_count + 1))
    else
        warn "  $link_name → $source_rel (源不存在: $source_abs)"
    fi
done

cd "$PROJECT_ROOT"

if [[ $link_count -eq 0 ]]; then
    warn "未配置任何软链接"
    echo ""
    echo "  请编辑 scripts/shell/init_project.sh 中的 LINKS 数组，"
    echo "  添加你需要链接的 meta2 子目录/文件。"
    echo ""
    echo "  或者手动创建:"
    echo "    cd ~/thesis_project/data/raw"
    echo "    ln -s ~/meta2/你的子目录 链接名"
    echo ""
else
    ok "$link_count 个软链接已建立"
fi

# ============================================================
# 3. 检查环境
# ============================================================
echo ""
info "检查环境..."

# Python
if command -v python3 &>/dev/null; then
    py_ver=$(python3 --version 2>&1)
    ok "Python: $py_ver"

    # 检查核心包
    missing=""
    for pkg in matplotlib numpy pandas scipy seaborn; do
        python3 -c "import $pkg" 2>/dev/null || missing="$missing $pkg"
    done
    if [[ -z "$missing" ]]; then
        ok "Python核心包已安装"
    else
        warn "缺少Python包:$missing"
        echo "  安装: pip install$missing --break-system-packages"
    fi
else
    fail "Python3 未安装"
fi

# R
if command -v Rscript &>/dev/null; then
    r_ver=$(Rscript --version 2>&1 | head -1)
    ok "R: $r_ver"
else
    warn "R 未安装 (堆叠图/PCoA/RDA/ggtree 需要)"
fi

# 字体
if [[ -f "/mnt/c/Windows/Fonts/arial.ttf" ]]; then
    ok "Windows字体可访问 (Arial/SimSun/TNR)"
else
    warn "Windows字体不可访问，将使用Linux替代字体"
    echo "  Liberation Sans ≡ Arial, Liberation Serif ≡ TNR"
fi

# ============================================================
# 4. 验证 D 盘
# ============================================================
echo ""
info "检查 D 盘..."

if [[ -d "/mnt/d" ]]; then
    # 尝试写入
    test_file="/mnt/d/.thesis_write_test_$$"
    if touch "$test_file" 2>/dev/null; then
        rm -f "$test_file"
        ok "D盘可写"
        mkdir -p "$WIN_TARGET"
        ok "D:\\workdata\\ 已创建"
    else
        warn "D盘只读，同步可能失败"
        echo "  尝试: sudo mount -t drvfs D: /mnt/d -o metadata"
    fi
else
    warn "/mnt/d 不存在"
    echo "  如果D盘已挂载到其他位置，设置环境变量:"
    echo "  export THESIS_WIN_TARGET=/mnt/你的路径/workdata"
fi

# ============================================================
# 5. 完成
# ============================================================
echo ""
echo "========================================================"
ok "初始化完成！"
echo "========================================================"
echo ""
echo "  接下来:"
echo "  1. 编辑 scripts/shell/init_project.sh 中的 LINKS 数组"
echo "     配置 meta2 → data/raw 的软链接"
echo "  2. 重新运行: bash scripts/shell/init_project.sh"
echo "  3. 开始写绘图脚本，每张图完成后:"
echo "     bash scripts/shell/sync_to_win.sh"
echo ""
