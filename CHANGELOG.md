# 更新日志

本项目遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/) 格式，
版本号采用 [Semantic Versioning](https://semver.org/lang/zh-CN/) 约定（主.次.补丁）。

> 新手用户：在 README 看到提示有新版本时，本地执行 `git pull origin master` 即可更新。
> 完整流程见 [README "更新到新版本"](README.md#-更新到新版本内测阶段) 一节。

---

## [0.8.2] — 2026-05-08

> R/Python 11 图侧侧对照工作完成 + RDA 数值对齐修复。

### 🐛 修复

- **RDA 数值与 R vegan 不一致**：`skbio.stats.ordination.rda()` 与
  `vegan::rda()` 的 eigval 归一化方式不同，导致 inertia 数值差 16-20×、
  环境因子 ANOVA F 差 3.5×、p 值方向反转（R: pH p=0.002 ↔ EnvMeta: 0.216）。
  修法：保留 `skbio.rda()` 输出的 ordination axes（site / biplot scores），
  但 inertia + ANOVA F 改用 SS-based 公式（`SS(Y) / (n-1)`，对标 vegan）。
  `explained_ref` 默认从 `"constrained"` 改为 `"total"`（对标 R vegan
  `summary()$cont$importance[2,]` 默认）。修复后所有 F / r / 解释度数值
  与 R vegan 4 位精度对齐。
  （[863e886](https://github.com/redlizzxy/EnvMeta/commit/863e886)）

### 📚 文档 / 验证

- **R/Python 侧侧对照 11 图全部完成**（`paper/benchmarks/validation/`）：
  - 数值精确一致（5 图）：alpha / pcoa / mag_quality / stackplot / RDA（修复后）
  - 算法等价（6 图）：lefse / gene_heatmap / log2fc / gene_profile / pathway / mag_heatmap
  - 5 个 R 脚本副本（强制 sans 字体）+ 11 个 README + 论文引用模板
  - 撤回 4 个误判"伪 bug"标记（log2fc 方向 / mag_heatmap 标尺 / gene_profile KO 过滤 / pathway KB 升级，均确认非 bug）
  （[ebe88b5](https://github.com/redlizzxy/EnvMeta/commit/ebe88b5)）

- **Paper 3 投稿前清单**（`paper/manuscript/paper3_pre_submission_checklist.md`）：
  跟踪 4 大任务进度（R 对照 ✅ / 第二数据集 ⬜ / 英文 README ⬜ / Methods 起草 ⬜）

### 🧪 测试

- pytest **293/293 全绿**（修复未破坏现有功能）；test_rda 4 个测试全部通过

---

## [0.8.1] — 2026-04-21

> Mac 端首批内测反馈修复版。聚焦 macOS 安装 / 假说评分 / 交互 HTML 三处实测 bug。

### 🐛 修复

- **HTML 交互导出**：切换分组（A / B / CK）后化学物-通路连线失效（拖拽不跟随、
  hover 不高亮）。根因为 `renderCycle` 清理清单漏了 `#em-chem-links` 容器，
  D3 data join key 复用导致 `.enter()` 返回空集合，DOM 残留旧 `<line>` 绑着
  过期 chemical 引用。修补一行清理 + 加回归测试断言五个兄弟容器都被清理。
  （[6667ae3](https://github.com/redlizzxy/EnvMeta/commit/6667ae3)）

- **假说评分（macOS）**：上传 YAML 文件评分时报 `[Errno 63] File name too long: '#`。
  根因为 `load_hypothesis()` 对字符串入参先 `Path(src).exists()` 试探是否文件路径，
  YAML 正文常达数 KB → macOS `stat()` 抛 `ENAMETOOLONG`，Python 3.11 的
  `Path._IGNORED_ERRNOS` 白名单不含该 errno → 异常上浮到 UI。修法：启发式
  （含换行 / ≥1024 字符 → 直判 YAML 正文）+ try/except OSError 双重保护。
  （[8085d14](https://github.com/redlizzxy/EnvMeta/commit/8085d14)）

### 📚 文档

- **小白安装指南**新增 3 个 Mac 安装常见问题：
  - Q1.1：`CondaToSNonInteractiveError` —— Anaconda 2024 ToS 接受门槛
  - Q1.2：`ResolutionImpossible: protobuf` —— Apple Silicon pip 旧 resolver
    选到无 arm64 wheel 的 protobuf；解法是先 `pip install --upgrade pip`
  - Q1.3：首次启动 `Welcome to Streamlit! ... Email:` —— 留空回车跳过
  （[4742b6c](https://github.com/redlizzxy/EnvMeta/commit/4742b6c)）

- **README 新增"更新到新版本"章节**：本地用户内测期 `git pull origin master`
  完整工作流 + 3 个常见场景（Already up to date / 本地 stash 冲突 / 新依赖）。
  （[977ef83](https://github.com/redlizzxy/EnvMeta/commit/977ef83)）

### 🧹 维护

- 包版本 `__version__` 由历史遗留的 `0.1.0` 同步到 `0.8.1`，与 README 长期使用
  的 "v0.8" 内测版叙事对齐。HTML 导出顶部的 "EnvMeta v…" 徽章会显示新版本。
- 测试：291 → **293 全绿**（新增 2 个回归测试）

---

## [0.8.0] — 2026-04-19

> v0.8 内测版（"Sunday Sprint"）—— 投稿前的功能完整版。

### ✨ 新增

- **新手落地包（S8-ux）**：数据准备指南（11 上游工具 → EnvMeta 输入映射）+
  图表选择向导（按研究问题反查推荐分析）+ 14 图"如何解读"expander +
  样例数据一键加载
- **导出中心（T1）**：4 tab 统一入口（图表 / Bundle / 脚本 / 文档）+ 批量 ZIP
- **HTML 交互导出（T2）**：400-550 KB 独立 HTML（D3 v7 inline 嵌入）+
  per-group 切换 + 化学物精确锚点 + 化学式上下标格式化（v1.3 精修）
- **在线部署**：Streamlit Cloud 自动部署 + 本地 Windows / Linux 行为一致性修复
- **内测素材**：腾讯问卷 + 海报 + 部署指南 + 小白安装指南

### 📚 文档

- CLAUDE.md 拆分：核心规范保留在 [CLAUDE.md](CLAUDE.md)，历史日志移到
  [DEBUG_NOTES.md](DEBUG_NOTES.md)（精简前 1969 行 → 264 行）

---

## [0.7.0] — 2026-04-13 至 2026-04-18

> Phase 3 核心：循环图 + 假说评分 + Fork Bundle + KEGG-driven KB。

### ✨ 新增

- **生物地球化学循环图 v2**：Mockup 10 合并细胞布局 + 跨元素化学物耦合 +
  keystone ★ 标注 + 4 元素 (As/N/S/Fe) × 18 通路 × 57 KO 自动推断
- **统计可信度三件套（S1/S2/S3.5）**：去偏 + 999 次置换 null_p + Saltelli 权重敏感度
- **机制假说 YAML 评分器（S3）**：5 类 claim（pathway_active / coupling_possible /
  env_correlation / keystone_in_pathway / group_contrast）+ Bradford Hill required veto
- **Fork Bundle（S4）**：打包 KB + YAML + config + KEGG 快照 → zip
- **KEGG-driven KB（S5）**：`envmeta kb-build` CLI 从 KEGG snapshot 重建知识库

---

## [0.5.0] — 2026-04-08 前后

> Phase 1 + Phase 2：14 图分析引擎完整版。

### ✨ 新增

- **Phase 1**：7 张 Reads-based 图（堆叠图 / α / PCoA / RDA / LEfSe / 基因热图 /
  log2FC）+ 代码生成器
- **Phase 2**：5 张 MAG-based 图（MAG 质量 / 丰度热图 / 通路完整度 / 基因谱 /
  共现网络 Gephi-prep）+ 4 图共享统一 4 层参数面板

---

## [0.1.0] — 项目骨架

> Phase 0：环境 + 知识库 + 11 文件类型识别 + Streamlit 基础。

---

## 版本号约定

- **主版本号 (X.0.0)**：架构性变化（如 Phase 4 插件框架上线时升 1.0.0）
- **次版本号 (0.X.0)**：新功能 / 新分析图 / 新模块
- **补丁号 (0.0.X)**：bug 修复 / 文档 / 跨平台兼容修复

发表论文后会标记 v1.0.0 + 申请 Zenodo DOI 永久归档。
