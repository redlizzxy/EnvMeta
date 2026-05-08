# Zenodo DOI 操作步骤（待办：投稿时执行）

> **状态**：待定，等到准备投稿 iMeta 时打 v1.0.0 release 一起做
> **创建日期**：2026-05-08
> **预计工时**：30 min（首次集成）+ 10 min（DOI 嵌入 README）

---

## 为什么留到投稿时做

1. **Zenodo 一旦归档版本，metadata 改不了**（标题、作者、关键词等）。投稿时
   论文标题 + 作者列表才会最终敲定。
2. **DOI 应该指向投稿对应的代码版本**（v1.0.0），而不是内测中间版（v0.8.x）。
3. iMeta 接收后，仓库可能再有 review 改动，那时再发 v1.0.1 / v1.1 二次归档。

---

## 5 步操作（投稿前 1 周做）

### Step 1 — 登录 Zenodo（GitHub OAuth）

访问：<https://zenodo.org/login/?next=/account/settings/github/>

点 "Log in with GitHub" → 授权 Zenodo 读取你的公开仓库列表。

### Step 2 — 打开 EnvMeta 仓库的 Zenodo 集成开关

进入 GitHub 集成设置页：<https://zenodo.org/account/settings/github/>

在仓库列表里找到 `redlizzxy/EnvMeta`，把右侧 toggle 从 **OFF** 拨到 **ON**。

> ⚠️ **关键**：这一步必须在创建 GitHub release **之前**做。Zenodo 只归档"开
> 关打开后"创建的 release。已存在的 release（如当前 v0.8.2）不会被自动归档。

### Step 3 — 在 GitHub 创建 v1.0.0 release

访问：<https://github.com/redlizzxy/EnvMeta/releases/new>

填写：

- **Choose a tag**：选 `v1.0.0`（或本地先 `git tag -a v1.0.0 -m "..."` 再 push）
- **Release title**：`EnvMeta v1.0.0 — Methodology paper submission release`
- **Description**：复制 [CHANGELOG.md](../../CHANGELOG.md) 对应版本段，加投稿信息：

  ```markdown
  EnvMeta v1.0.0 — Submission release for the methodology paper:
  "EnvMeta: A visualization platform for environmental metagenomics with
   automated biogeochemical-cycle inference and YAML hypothesis scoring"
  Submitted to iMeta on YYYY-MM-DD.

  ## Highlights

  - 14 publication-grade plots (7 Reads-based + 5 MAG-based + 1 cycle figure + 1 hypothesis scorer)
  - KEGG-driven cycle inference (4 elements × 18 pathways × 57 KOs)
  - YAML hypothesis scorer with 5 claim types + permutation null-p + weight sensitivity
  - Standalone interactive HTML SI (D3.js inlined, 400-550 KB, fully offline)
  - Fork Bundle for paper-tool binding
  - 293/293 unit tests green; R/Python side-by-side validation across 11 figures
  ```

- 不要勾 "This is a pre-release"
- 点 "Publish release"

### Step 4 — 等待 Zenodo 自动归档

通常 1-3 分钟。归档成功后，注册的 Zenodo 邮箱会收到通知邮件。

去 <https://zenodo.org/account/settings/github/> 查看仓库行右侧出现 DOI badge，
形如：

```
DOI: 10.5281/zenodo.XXXXXXX
```

### Step 5 — 完善 Zenodo metadata + 嵌入 README

进入对应记录页（点 DOI badge → "View on Zenodo"），点 "Edit"：

- **Title**：`EnvMeta v1.0.0` → 改为论文标题（投稿时确定）
- **Authors**：补全（默认只有 GitHub username；加 ORCID + 单位）
- **Description**：补论文摘要（200-300 字）
- **Keywords**：metagenomics, biogeochemistry, microbial ecology, KEGG, hypothesis testing, visualization
- **License**：MIT（已自动从 LICENSE 文件读取）
- **Related identifiers**（重要）：投稿后加论文 DOI 作为 "is supplemented by"
- 点 "Save" → "Publish"

得到最终 DOI 后，告诉 Claude，会自动加 badge 到 README.md + README_CN.md：

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

通常加在 README 顶部 logo / tagline 下方。

---

## 论文里如何引用

iMeta Methods 末尾或 Data Availability 段落：

> EnvMeta source code is available at <https://github.com/redlizzxy/EnvMeta>
> (release v1.0.0, archived at Zenodo: <https://doi.org/10.5281/zenodo.XXXXXXX>),
> licensed under MIT.

---

## 后续维护

- 论文 review 中如有代码修改：发 v1.0.1 release，Zenodo 会自动归档第二个 DOI
  （新版本 DOI 与 v1.0.0 不同，但二者通过 "concept DOI" 关联）
- 论文接收后：发 v1.1.0 → 第三个 DOI，metadata 里加论文 DOI 链接
- 重大架构升级（Phase 4 插件框架）：发 v2.0.0
