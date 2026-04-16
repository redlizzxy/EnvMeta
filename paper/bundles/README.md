# Fork Bundle — 论文-EnvMeta 绑定发布协议

## 这是什么

Fork Bundle 是一个 `.zip` 文件，打包一篇论文的 EnvMeta 分析所需的全部**非数据**
产物：

```
paper_bundle.zip
├── manifest.yaml           # 论文 + bundle 元数据
├── kb/
│   ├── elements.json       # 这篇论文的专属知识库
│   └── kegg_snapshot.json  # 可选，KEGG 快照（溯源）
├── hypotheses/
│   └── *.yaml              # 论文的机制假说清单
├── config/
│   └── cycle_params.yaml   # 分析参数（复现图表用）
└── README.md               # 如何使用本 bundle
```

**不包含原始数据**（体积 / 敏感度），数据指针写在 `manifest.paper_doi` 或
manifest 里的自定义字段。

## 为什么需要 Bundle

CLAUDE.md 的"可定制框架 + 论文绑定发布"原则要求：**每篇论文都有一份自己的
EnvMeta 配置**，读者 fork 即可复现。不依赖中央仓库 / 云服务。

**审稿场景**：审稿人拿到 bundle → 启动 EnvMeta → 一键加载 → 跑一遍 → 验证图
表可复现 → 不用问作者任何配置细节。

## 如何制作

### A. 在 EnvMeta 界面里（推荐）

1. 启动 `streamlit run app.py`
2. 循环图页完成参数调节后，展开 **"📦 Fork Bundle"**
3. 右列「导出当前状态为 Bundle」：填名字 / 作者 / DOI → 点"创建 Bundle"
4. 下载 `.zip`

界面会自动把：
- 当前的 `elements.json`
- 已上传的假说 YAML
- 当前所有循环图参数

打进 zip。

### B. 命令行（CI / 自动化）

```bash
python -m envmeta bundle-create \
  --output my_paper.zip \
  --kb envmeta/geocycle/knowledge_base/elements.json \
  --hypothesis paper/hypotheses/arsenic_steel_slag.yaml \
  --kegg-snapshot envmeta/geocycle/kegg_snapshot.json \
  --config my_params.yaml \
  --name "My Paper 2026" \
  --author "..." \
  --paper-doi "10.xxx/my"
```

### C. Python API

```python
from envmeta.tools.bundle import create_bundle

create_bundle(
    "my_paper.zip",
    kb_path="envmeta/geocycle/knowledge_base/elements.json",
    hypothesis_paths=["paper/hypotheses/arsenic_steel_slag.yaml"],
    config={"completeness_threshold": 50, ...},
    name="My Paper 2026",
    author="...",
    paper_doi="10.xxx/my",
)
```

## 如何使用别人的 Bundle

### 界面加载

1. 启动 EnvMeta
2. 循环图页 → **"📦 Fork Bundle"** → 左列"加载 Bundle" → 上传 zip
3. 看 summary（作者 / DOI / 版本兼容性）
4. 随后正常上传你自己的 6 个数据文件 → 按 bundle 的参数生成图表

### 命令行检查

```bash
python -m envmeta bundle-inspect my_paper.zip
```

输出示例：
```
Bundle: my_paper.zip
  name:         My Paper 2026
  author:       ...
  paper_doi:    10.xxx/my
  envmeta ver:  bundle=0.1.0 runtime=0.1.0
  files:        8
    kb:         yes
    hypotheses: 2
```

## 版本兼容性

`manifest.envmeta_version` 记录制作 bundle 时的 EnvMeta 版本。加载时与运行时
对比：
- **一致**：完全兼容
- **不一致**：UI 显示 ⚠️ warning 但仍允许加载；未知字段被忽略不报错。这是
  有意设计（**向后容忍**），让老 bundle 在新版本可用

## 示例

[arsenic_steel_slag_bundle.zip](arsenic_steel_slag_bundle.zip) 是"钢渣驱动
铁氧化-砷固定"研究的实测 bundle，可直接加载参考。

## 设计不做的事

- ❌ 加密 / 签名（社区可用 GPG 自己签；EnvMeta 不介入信任链）
- ❌ 在线注册 / 检索（与"fork-based 分布式"原则冲突）
- ❌ 自动 DOI 申请（去 Zenodo 走标准流程）
- ❌ 捆绑原始数据（体积 / 伦理 / 版权）

## 相关 CLI

| 命令 | 作用 |
|---|---|
| `envmeta bundle-create` | 打包 |
| `envmeta bundle-inspect` | 检查 |
| `envmeta kb-build` | 从 KEGG 快照生成 KB（kb/elements.json 来源）|
| `envmeta hypothesis-validate` | 校验假说 YAML |

典型工作流：
```
kb-build  →  hypothesis-validate  →  bundle-create  →  (发布)  →  bundle-inspect (审稿端)
```
