#!/bin/bash
# collect_and_sync.sh — 收集所有脚本依赖的输入文件到thesis_project，然后同步到D盘
# cd ~/thesis_project && bash scripts/shell/collect_and_sync.sh

set -euo pipefail
GREEN='\033[0;32m'; BLUE='\033[0;34m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; NC='\033[0m'
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[  OK]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
fail() { echo -e "${RED}[FAIL]${NC} $1"; }

PROJECT_ROOT="$HOME/thesis_project"
RAW_DIR="$PROJECT_ROOT/data/raw"
META2="$HOME/meta2/result"
META2_BIN="$HOME/meta2/binresult"
mkdir -p "$RAW_DIR"

# ============================================================
# 第一步：收集所有输入文件
# ============================================================
info "========== 收集输入文件到 $RAW_DIR =========="
echo ""

copied=0; skipped=0; missing=0

copy_file() {
    local src="$1"
    local dst_name="$2"
    local dst="$RAW_DIR/$dst_name"
    if [[ -f "$src" ]]; then
        if [[ ! -f "$dst" ]] || [[ "$src" -nt "$dst" ]]; then
            cp "$src" "$dst"
            ok "  $dst_name"
            copied=$((copied + 1))
        else
            echo "  [SKIP] $dst_name (已是最新)"
            skipped=$((skipped + 1))
        fi
    else
        fail "  $dst_name — 源不存在: $src"
        missing=$((missing + 1))
    fi
}

# ---- 第二章：reads-based分析 ----
info "第二章 reads-based 输入文件"

# run_stackplots.sh: Phylum/Genus/Species.txt + metadata.txt
copy_file "$META2/metaphlan4/Phylum.txt"    "Phylum.txt"
copy_file "$META2/metaphlan4/Genus.txt"     "Genus.txt"
copy_file "$META2/metaphlan4/Species.txt"   "Species.txt"
copy_file "$META2/metadata.txt"             "metadata.txt"

# run_alpha.sh: alpha多样性
copy_file "$META2/kraken2/tax_count.alpha"  "tax_count.alpha"
copy_file "$META2/metaphlan4/alpha.txt"     "alpha.txt"

# run_pcoa.sh: beta多样性
copy_file "$META2/metaphlan4/beta_bray.txt" "beta_bray.txt"

# run_rda.sh: 物种+环境因子
copy_file "$META2/elementdata/env_factors.txt" "env_factors.txt"

# run_heatmap_log2fc.sh: KO TPM数据
copy_file "$META2/eggnog/eggnog.KEGG_ko.TPM.spf" "eggnog.KEGG_ko.TPM.spf"

# run_lefse.sh: LEfSe输入（多个候选路径）
lefse_found=false
for try_path in "$META2/temp/input.res" "$META2/metaphlan4/input.res" "$HOME/meta2/temp/input.res"; do
    if [[ -f "$try_path" ]]; then
        copy_file "$try_path" "input.res"
        lefse_found=true
        break
    fi
done
[[ "$lefse_found" == false ]] && { fail "  input.res — 所有候选路径均不存在"; missing=$((missing + 1)); }

echo ""

# ---- 第三章：MAG分析 ----
info "第三章 MAG 输入文件"

# GTDB分类（多个sh共用）
copy_file "$META2_BIN/temp/gtdb_all/classify/tax.bac120.summary.tsv" "tax_bac120_summary.tsv"
copy_file "$META2_BIN/temp/gtdb_all/classify/tax.ar53.summary.tsv"   "tax_ar53_summary.tsv"

# CoverM丰度（多个sh共用）
copy_file "$META2_BIN/result/coverm/abundance.tsv" "abundance.tsv"

# Keystone物种列表（多个sh共用）
copy_file "$META2_BIN/result/keystone_phylogeny/keystone_species.txt" "keystone_species.txt"

# run_MAG_quality.sh: CheckM2质量报告
copy_file "$META2_BIN/result/checkm2/quality_report.tsv" "quality_report.tsv"

# run_MAG_gene_profile.sh: iTOL热图 + KEGG目标KO
# iTOL文件可能在多个位置
itol_found=false
for try_path in \
    "$META2_BIN/result/itol/07_element_KO_heatmap.txt" \
    "$META2_BIN/result/element_cycles/07_element_KO_heatmap.txt" \
    "$PROJECT_ROOT/data/raw/07_element_KO_heatmap.txt"; do
    if [[ -f "$try_path" ]]; then
        copy_file "$try_path" "07_element_KO_heatmap.txt"
        itol_found=true
        break
    fi
done
[[ "$itol_found" == false ]] && { fail "  07_element_KO_heatmap.txt — 未找到"; missing=$((missing + 1)); }

# kegg_target_only.tsv 可能在多个位置
kegg_found=false
for try_path in \
    "$META2_BIN/result/eggnog/kegg_target_only.tsv" \
    "$PROJECT_ROOT/kegg_target_only.tsv" \
    "$PROJECT_ROOT/data/raw/kegg_target_only.tsv"; do
    if [[ -f "$try_path" ]]; then
        copy_file "$try_path" "kegg_target_only.tsv"
        kegg_found=true
        break
    fi
done
[[ "$kegg_found" == false ]] && { fail "  kegg_target_only.tsv — 未找到"; missing=$((missing + 1)); }

# run_cooccurrence_network.sh: 网络分析文件
copy_file "$META2_BIN/result/downstream_analysis/network_node_stats.txt"  "network_node_stats.txt"
copy_file "$META2_BIN/result/downstream_analysis/network_module_info.txt" "network_module_info.txt"

echo ""
info "收集完成: 复制=$copied, 跳过=$skipped, 缺失=$missing"

# ============================================================
# 第二步：验证完整性
# ============================================================
echo ""
info "========== 验证文件完整性 =========="

EXPECTED_FILES=(
    # 第二章 (10个)
    "Phylum.txt" "Genus.txt" "Species.txt" "metadata.txt"
    "tax_count.alpha" "alpha.txt" "beta_bray.txt"
    "env_factors.txt" "eggnog.KEGG_ko.TPM.spf" "input.res"
    # 第三章 (9个)
    "tax_bac120_summary.tsv" "tax_ar53_summary.tsv"
    "abundance.tsv" "keystone_species.txt" "quality_report.tsv"
    "07_element_KO_heatmap.txt" "kegg_target_only.tsv"
    "network_node_stats.txt" "network_module_info.txt"
)

present=0; absent=0
for f in "${EXPECTED_FILES[@]}"; do
    if [[ -f "$RAW_DIR/$f" ]]; then
        present=$((present + 1))
    else
        fail "  缺失: $f"
        absent=$((absent + 1))
    fi
done
echo ""
if [[ $absent -eq 0 ]]; then
    ok "全部 ${#EXPECTED_FILES[@]} 个输入文件就绪！"
else
    warn "$present/${#EXPECTED_FILES[@]} 个文件就绪, $absent 个缺失"
fi

# ============================================================
# 第三步：检查脚本和输出
# ============================================================
echo ""
info "========== 检查脚本和输出 =========="

py_count=$(find "$PROJECT_ROOT/scripts/python" -name "*.py" 2>/dev/null | wc -l)
r_count=$(find "$PROJECT_ROOT/scripts/R" -name "*.R" 2>/dev/null | wc -l)
sh_count=$(find "$PROJECT_ROOT/scripts/shell" -name "*.sh" 2>/dev/null | wc -l)
info "脚本: Python=$py_count, R=$r_count, Shell=$sh_count"

mkdir -p "$PROJECT_ROOT/figures/paper_en" "$PROJECT_ROOT/figures/thesis_cn" "$PROJECT_ROOT/data/processed"
n_en=$(find "$PROJECT_ROOT/figures/paper_en" -name "Fig*" -type f 2>/dev/null | wc -l)
n_cn=$(find "$PROJECT_ROOT/figures/thesis_cn" -name "Fig*" -type f 2>/dev/null | wc -l)
n_stat=$(find "$PROJECT_ROOT/data/processed" -name "Fig*" -type f 2>/dev/null | wc -l)
info "图表: EN=$n_en, CN=$n_cn, 统计=$n_stat"

# ============================================================
# 第四步：同步到Windows
# ============================================================
echo ""
WIN_DIR="/mnt/d/workdata/thesis_project"
info "========== 同步到 $WIN_DIR =========="

rsync -av --delete \
    --exclude='__pycache__' \
    --exclude='*.pyc' \
    --exclude='.ipynb_checkpoints' \
    "$PROJECT_ROOT/" "$WIN_DIR/"

ok "同步完成！"

# ============================================================
# 汇总
# ============================================================
echo ""
echo "============================================"
echo "  输入文件: $RAW_DIR/ (${#EXPECTED_FILES[@]}个)"
echo "  WSL项目:  $PROJECT_ROOT"
echo "  Windows:  D:\\workdata\\thesis_project"
echo "============================================"
echo ""
echo "data/raw/ 文件列表:"
ls -lh "$RAW_DIR/" | tail -n +2
