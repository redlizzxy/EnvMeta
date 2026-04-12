#!/bin/bash
# ==============================================================================
# pack_copyright.sh — 软件著作权申请专用打包
# ==============================================================================
# 软著要求：提交前30页 + 后30页源代码（每页50行，共3000行）
# 本脚本自动：
#   1. 合并所有 .py 和 .R 脚本为单个文件
#   2. 添加统一的文件头（软件名称、版本、作者等）
#   3. 截取前30页和后30页
#   4. 输出 Word 可用的纯文本
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

OUT_DIR="$PROJECT_ROOT/docs/copyright"
mkdir -p "$OUT_DIR"

# ---- 软件信息（根据实际修改） ----
SW_NAME="砷选冶渣微生物群落宏基因组分析系统"
SW_NAME_EN="Arsenic Smelting Slag Metagenomics Analysis System"
SW_VERSION="V1.0"
SW_AUTHOR="[你的姓名]"
SW_DATE="2026年"
LINES_PER_PAGE=50

echo "=========================================="
echo "软件著作权代码打包"
echo "=========================================="

# ---- 合并所有源代码 ----
MERGED="$OUT_DIR/merged_source.txt"

cat > "$MERGED" << EOF
$SW_NAME $SW_VERSION
$SW_NAME_EN
版权所有 (C) $SW_DATE $SW_AUTHOR

================================================================
源代码清单
================================================================

EOF

# 依次合并 Python 和 R 文件
file_count=0
for ext in py R; do
    find "$PROJECT_ROOT/scripts" "$PROJECT_ROOT/config" -name "*.$ext" -type f | sort | while read -r f; do
        rel_path="${f#$PROJECT_ROOT/}"
        echo "" >> "$MERGED"
        echo "# ================================================================" >> "$MERGED"
        echo "# 文件: $rel_path" >> "$MERGED"
        echo "# ================================================================" >> "$MERGED"
        echo "" >> "$MERGED"
        cat "$f" >> "$MERGED"
        echo "" >> "$MERGED"
        file_count=$((file_count + 1))
    done
done

TOTAL_LINES=$(wc -l < "$MERGED")
echo "合并完成: $TOTAL_LINES 行"

# ---- 截取前30页和后30页 ----
FRONT_LINES=$((30 * LINES_PER_PAGE))  # 1500行
BACK_LINES=$((30 * LINES_PER_PAGE))   # 1500行

FRONT_FILE="$OUT_DIR/前30页源代码.txt"
BACK_FILE="$OUT_DIR/后30页源代码.txt"
FULL_FILE="$OUT_DIR/完整源代码.txt"

# 前30页
head -n $FRONT_LINES "$MERGED" > "$FRONT_FILE"
echo "前30页: $(wc -l < "$FRONT_FILE") 行 → $FRONT_FILE"

# 后30页
if [[ $TOTAL_LINES -gt $BACK_LINES ]]; then
    tail -n $BACK_LINES "$MERGED" > "$BACK_FILE"
else
    cp "$MERGED" "$BACK_FILE"
fi
echo "后30页: $(wc -l < "$BACK_FILE") 行 → $BACK_FILE"

# 完整源代码
cp "$MERGED" "$FULL_FILE"
echo "完整版: $TOTAL_LINES 行 → $FULL_FILE"

# ---- 代码统计 ----
STATS_FILE="$OUT_DIR/代码统计.txt"
cat > "$STATS_FILE" << EOF
$SW_NAME 代码统计报告
==============================
生成时间: $(date '+%Y-%m-%d %H:%M')

文件统计:
  Python脚本: $(find scripts config -name "*.py" -type f 2>/dev/null | wc -l) 个
  R脚本:      $(find scripts config -name "*.R" -type f 2>/dev/null | wc -l) 个
  Shell脚本:  $(find scripts -name "*.sh" -type f 2>/dev/null | wc -l) 个

代码行数:
  总行数:     $TOTAL_LINES 行
  前30页:     $FRONT_LINES 行
  后30页:     $BACK_LINES 行

各文件行数:
EOF

find scripts config -name "*.py" -o -name "*.R" -o -name "*.sh" | sort | while read -r f; do
    lines=$(wc -l < "$f")
    printf "  %-50s %5d 行\n" "$f" "$lines" >> "$STATS_FILE"
done

echo ""
echo "=========================================="
echo "软著打包完成！"
echo "=========================================="
echo "输出目录: $OUT_DIR/"
ls -la "$OUT_DIR/"
echo ""
echo "提交材料："
echo "  1. ${FRONT_FILE##*/} （前30页）"
echo "  2. ${BACK_FILE##*/} （后30页）"
echo "  3. ${STATS_FILE##*/} （代码统计，备查）"
