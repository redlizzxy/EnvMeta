# EnvMeta Publication Strategy — Path A 改良版

> **存档日期**：2026-05-09
> **决策依据**：基于用户硬约束（2026 内 ≥ 1 篇接收 + 2027-05 前所有文章见刊）+ 课题
> 论文当前状态（仅大纲 + As 形态待重测 + 机制争议未决）的可行性分析。
> **替代方案讨论**：Path A / B / C / D 详见对话上下文（2026-05-09 conversation）。

---

## 选定方案：**Path A 改良版（EnvMeta 先 + bioRxiv preprint + 课题论文并行起草）**

### 决策核心

- **课题论文起草要 3-6 月**（As 重测 + 机制争议）→ 起步就晚，无法在 2026 年内接收
- **EnvMeta paper 起草需 3-4 周**（修 5 大 Major Issues + 6 张图）→ 可在 2026 年内接收
- 因此 **EnvMeta 必须先发**才能满足"2026 内 ≥ 1 篇接收"的硬指标

### 关键加速器：bioRxiv preprint

- 投 bioRxiv 后立即拿 DOI（2-3 天上线）
- 课题论文起稿时可立即引用 EnvMeta（不需要等 iMeta 接收）
- iMeta 投稿不受影响（明确允许 preprint）
- 国内可用，无政策限制（但需偶尔科学上网或用 ChinaXiv 镜像）

---

## 时间表（按月）

```
2026-05      [Week 1-3]  修 Mock Review 5 大 Major Issues + 文字部分
                          课题论文 As 形态重测启动（并行）

2026-06      [Week 4]    F1/F3/F5/F6/F7/F10 6 张图绘制
             [Week 5]    §5.7.1 Vancouver → iMeta 引用格式转换
             [Week 6]    EnvMeta paper 上 bioRxiv ⭐（关键节点）

2026-06      [Week 7-8]  EnvMeta 投 iMeta（cover letter 含 bioRxiv DOI + 课题论文 in prep）

2026-07~09   课题论文起草（用 EnvMeta bioRxiv DOI 作为 Methods 引用）
             EnvMeta 第一轮 review

2026-10      EnvMeta first-round revision response
             课题论文初稿完成 + 导师审稿

2026-11~12   EnvMeta 接收 ⭐ (满足 Q2 "2026 接收" 硬指标)
             课题论文投环境微生物期刊（cite EnvMeta accepted）

2027-01~02   EnvMeta 见刊（iMeta 2027 早期）⭐
             课题论文 review

2027-03~05   课题论文 revision + 接收 + 见刊（满足 "所有文章 2027-05 前见刊" 截止）
```

---

## 课题论文目标期刊候选（2026 Q4 投稿时再定）

| 期刊 | IF (2024) | 平均周期 | 适配度 |
|---|---|---|---|
| **Environmental Science & Technology (ES&T)** | 11.4 | 4-6 月 | ⭐⭐⭐ 砷修复机制契合 |
| **ISME Journal** | 11.0 | 5-7 月 | ⭐⭐⭐ MAG-based mechanism 契合 |
| **Microbiome** | 13.8 | 5-7 月 | ⭐⭐ 注重微生物组 + 机制 |
| **FEMS Microbiology Ecology** | 4.2 | 3-4 月 | ⭐ 速度快但 IF 低 |
| **Front. Microbiol.** | 4.0 | 3-4 月 | ⭐ 速度快但 IF 低 |

> **建议**：投稿前 1 月再决定，看导师意见 + 课题论文写完后的实际亮点。

---

## bioRxiv 投稿操作 checklist（投稿前 1 周做）

- [ ] 准备 manuscript PDF（含 Title / Authors / Affiliations / Abstract / Highlights / 全文 / Figures / References）
- [ ] 注册 bioRxiv 账号（biorxiv.org）+ ORCID 关联
- [ ] 选择 subject area: "Bioinformatics" 或 "Microbiology"
- [ ] 上传 PDF + supplementary materials（Tables / 大图 PDF / SI HTML）
- [ ] 填 author info（中文姓名 + 拼音 + 单位 + 邮箱）
- [ ] **License 选择**：CC-BY 4.0（最宽松，符合 iMeta 要求）
- [ ] 选 "Submit to journal" → 选 iMeta（这步告诉 bioRxiv 你打算投 iMeta，bioRxiv 会自动 forward 给 iMeta — 如果选了 iMeta 的 transfer 服务）
  - 或选 "Will submit elsewhere"（自己 manual 投 iMeta）
- [ ] 提交 → 等 2-3 天审核 → 上线 + 拿 DOI
- [ ] 拿到 DOI 后立即更新 GitHub README + 简历 + 项目报告

---

## 风险评估 + 应急

| 风险 | 概率 | 应急方案 |
|---|---|---|
| EnvMeta iMeta 一审 reject | 30% | 改投 *Bioinformatics* / *Front Microbiol*；bioRxiv preprint 不受影响 |
| 课题论文 As 重测延迟 | 50% | 推迟课题论文投稿到 2027-Q1，仍能 2027-05 前见刊（紧但可行）|
| iMeta 审稿超期 | 20% | 投稿后 4 月没消息可发邮件催；preprint 状态保持 |
| bioRxiv 国内访问慢 | 70% | 用 ChinaXiv 中国镜像（imeta.science 也支持） |
| 课题论文被拒 | 30% | 不影响 EnvMeta；改投速度更快期刊 |

---

## 维护记录

| 日期 | 事项 |
|---|---|
| 2026-05-09 | 方案确定（基于 mock review v0.9.1 + 用户约束 Q1-Q3）|
