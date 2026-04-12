#关于R修改sh后运行不了的问题——是因为R（或RStudio）保存时把换行符从Unix的LF改成了Windows的CRLF，bash不认\r。以后遇到这个问题，运行前执行一下：
bashsed -i 's/\r$//' scripts/shell/run_cooccurrence_network.sh
#或者用VSCode/Notepad++编辑sh文件，不要用R/RStudio。
# cn、en论文相关字体安装方法1：从Windows复制字体到Linux（推荐）
mkdir -p ~/.local/share/fonts
cp /mnt/c/Windows/Fonts/arial.ttf ~/.local/share/fonts/
cp /mnt/c/Windows/Fonts/arialbd.ttf ~/.local/share/fonts/
cp /mnt/c/Windows/Fonts/times.ttf ~/.local/share/fonts/
cp /mnt/c/Windows/Fonts/timesbd.ttf ~/.local/share/fonts/
cp /mnt/c/Windows/Fonts/simsun.ttc ~/.local/share/fonts/

# 刷新字体缓存
fc-cache -fv ~/.local/share/fonts

# 清除matplotlib字体缓存（关键！）
rm -rf ~/.cache/matplotlib

# 验证
python3 -c "
import matplotlib.font_manager as fm
fm._load_fontmanager(try_read_cache=False)
avail = sorted({f.name for f in fm.fontManager.ttflist})
for f in avail:
    if any(k in f.lower() for k in ['arial','times','simsun']):
        print(f)
"
#前置
cd ~/thesis_project
bash scripts/shell/init_project.sh
#理化指标
pip install openpyxl
pip install -r requirements.txt
cd ~/thesis_project
python3 scripts/python/01_physicochemical.py ~/meta2/result/elementdata/ordata.xlsx
#物种组成堆叠图重绘
cd ~/thesis_project
bash scripts/shell/run_stackplots.sh
#PCoA
cd ~/thesis_project
bash scripts/shell/run_pcoa.sh
#alpha
cd ~/thesis_project
bash scripts/shell/run_alpha.sh
#RDA
cd ~/thesis_project
bash scripts/shell/run_rda.sh
#lefse
cd ~/thesis_project
bash scripts/shell/run_lefse.sh
#heatmap gene
cd ~/thesis_project
bash scripts/shell/run_heatmap_log2fc.sh
#第三章输入文件查找
# CheckM2结果
ls ~/meta2/result/binning/result/checkm2/
# 或
find ~/meta2/result -name "quality_report*" -o -name "checkm2*" 2>/dev/null | head

# GTDB-Tk分类
find ~/meta2/result -name "*gtdbtk*" -o -name "*classify*" -o -name "*taxonomy*" 2>/dev/null | head

# MAG统计
ls ~/meta2/result/binning/result/
# 或
find ~/meta2/result -name "*bin*stat*" -o -name "*mag*stat*" 2>/dev/null | head

# Keystone列表
ls ~/meta2/result/binning/keystone_phylogeny/
find ~/meta2/binresult/temp/gtdb_classify/ -name "*.tsv" | head
ls ~/meta2/binresult/temp/gtdb_classify/
ls ~/meta2/binresult/temp/drep95/
find ~/meta2/binresult -name "*keystone*" -o -name "*key*" 2>/dev/null | head
find ~/meta2/binresult -name "*.summary.tsv" 2>/dev/null | head
find ~/meta2/binresult/temp/gtdb_classify -type f 2>/dev/null | head -20
find ~/meta2/binresult/temp/gtdb_all -type f 2>/dev/null | head -20

cat ~/meta2/binresult/result/keystone_phylogeny/keystone_taxonomy.txt | head -10
cat ~/meta2/binresult/result/keystone_phylogeny/keystone_species.txt | head -10
find ~/meta2/binresult/temp/gtdb_classify -type f 2>/dev/null | head -20
find ~/meta2/binresult/temp/gtdb_all -type f 2>/dev/null | head -20
#MAG质量
cd ~/thesis_project
bash scripts/shell/run_MAG_quality.sh
#前置文件查找
# 1. GTDB-Tk生成的树文件
find ~/meta2/binresult/temp/gtdb_all -name "*.tree" | head -10
find ~/meta2/binresult/temp/gtdb_classify -name "*.tree" | head -10

# 2. 看看keystone_phylogeny目录里有什么
ls ~/meta2/binresult/result/keystone_phylogeny/

# 3. itol目录（你可能已经准备过iTOL文件）
ls ~/meta2/binresult/result/itol/ 2>/dev/null

# 4. coverm丰度数据
ls ~/meta2/binresult/result/coverm/ 2>/dev/null
head -3 ~/meta2/binresult/result/coverm/*.tsv 2>/dev/null | head -20

# 5. 元素循环基因数据（MAG级别）
ls ~/meta2/binresult/result/element_cycles/ 2>/dev/null
head -3 ~/meta2/binresult/result/element_cycles/*.txt 2>/dev/null | head -20
# 检查eggnog目录里的MAG功能注释
ls ~/meta2/binresult/result/eggnog/
head -5 ~/meta2/binresult/result/eggnog/keystone_gene_profiles.txt
head -5 ~/meta2/binresult/result/eggnog/keystone_species_functions.txt

# 检查element_cycles里有没有其他文件
find ~/meta2/binresult/result/element_cycles/ -type f
cat ~/meta2/binresult/result/element_cycles/target_KOs.txt

# 检查是否有MAG级别的KO注释（eggnog-mapper输出）
find ~/meta2/binresult -name "*annotations*" -o -name "*KEGG*" -o -name "*.emapper*" 2>/dev/null | head -10

head -3 ~/meta2/binresult/result/eggnog/KEGG_ko_list.tsv | cat
# 看看列是：MAG名 KO列表？还是已经是矩阵格式？
wc -l ~/meta2/binresult/result/eggnog/KEGG_ko_list.tsv
head -3 ~/meta2/binresult/result/element_cycles/target_KOs.txt
head -3 ~/meta2/binresult/result/eggnog/KEGG_ko_list.tsv | cut -f1,3
# 用cat -A看有没有特殊字符/空格
cat -A ~/meta2/binresult/result/element_cycles/target_KOs.txt | head -5

# 看列数
awk '{print NF, $0}' ~/meta2/binresult/result/element_cycles/target_KOs.txt | head -5

# 第1步：提取干净的KO号列表
grep -oP 'K\d{5}' \
    ~/meta2/binresult/result/element_cycles/target_KOs.txt \
    > /tmp/ko_ids_clean.txt

wc -l /tmp/ko_ids_clean.txt   # 预期：57

# 第2步：过滤KEGG注释文件
awk 'NR==FNR{ko[$1]=1; next}
     FNR==1{print; next}
     {
       split($3,a,",");
       for(i in a){
         sub(/^ko:/,"",a[i]);
         if(a[i] in ko){print; break}
       }
     }' \
    /tmp/ko_ids_clean.txt \
    ~/meta2/binresult/result/eggnog/KEGG_ko_list.tsv \
    > kegg_target_only.tsv

wc -l kegg_target_only.tsv   # 预期：几千行

# 提取MAG ID和分类信息（取属名或种名作为标签）
cat ~/meta2/binresult/temp/gtdb_all/classify/tax.bac120.summary.tsv \
    ~/meta2/binresult/temp/gtdb_all/classify/tax.ar53.summary.tsv | \
    grep "^Mx_All" | cut -f1,2 > ~/thesis_project/mag_taxonomy_labels.tsv

# 同步到Windows
cp ~/thesis_project/mag_taxonomy_labels.tsv /mnt/d/workdata/

cp ~/meta2/binresult/temp/gtdb_all/classify/tax.bac120.summary.tsv ~/thesis_project/
cp ~/meta2/binresult/temp/gtdb_all/classify/tax.ar53.summary.tsv ~/thesis_project/
#3-3物种热图和子图
cd ~/thesis_project
bash scripts/shell/run_MAG_gene_profile.sh
#3-4top30mag
cd ~/thesis_project
bash scripts/shell/run_MAG_abundance_heatmap.sh
#3-5通路完整图和子图
cd ~/thesis_project
bash scripts/shell/run_pathway_completeness.sh
#网络图
cd ~/thesis_project
pip install networkx  # 如果没装的话
bash scripts/shell/run_cooccurrence_network.sh

#同步输入输出文件，在终端中运行
cd ~/thesis_project
bash scripts/shell/collect_and_sync.sh
