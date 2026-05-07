# v1 问卷阶段性发现存档（2026-05-07，n=2，校准版）

> **状态**：归档暂存，**已用原始文本标签 CSV 校正**，**不急着检查或改 app**。
> 等 v2 问卷拿到 n ≥ 8 后再综合分析。
>
> **本档历史**：
> - 2026-05-07 v1：基于编码 CSV 写的初版（部分误读 Likert 选项 vs ABC 选项混淆，已纠正）
> - 2026-05-07 v2：基于原始文本标签 CSV 重写，**这是当前正确版本**

---

## ⚠️ 我之前的解读错误（已纠正）

编码 CSV 里数字"1, 2, 3..."同时表示 ① Likert 1-5 分数 ② 选项 ABCDE 的索引（1=A）。
我之前把后者当成前者，导致虚报红色信号：

| 题号 | 我之前错说 | 真实情况 |
|---|---|---|
| Q15 上下文推荐 ??? 含义 | "都给 1（不懂）" ❌ | **都选 A（一眼能看出）** ✅ |
| Q44 数据准备文档完整度 | "都给 1（最低）" ❌ | **都给 5（最高）** |
| Q57 help tooltip 有用度 | "都给 1" ❌ | 陈源 5 / 匿名 4 |

那 3 条假阳性已删除。

---

## 数据来源

- 问卷：v1 设计（87 题，预期 30-40 min，实测 30+ min）
- 回收：n=2，共 2 份导出格式
  - 编码 CSV：`responses_2026-05-07_n2_encoded.csv`
  - 原始 CSV（文本标签）：`responses_2026-05-07_n2_raw.csv` ← 当前分析依据
- 时间窗：2026-04-21 ~ 2026-04-23
- **回收率：~20%（10 名目标 → 2 份）—— 这是最关键的元信号**

## 受访者画像

| ID | 来源 | 背景 | 任务完成 | 耗时 |
|---|---|---|---|---|
| #1 | Win Chrome / 课题组 | 陈源 / 环境工程硕士 2 / **零代码 + 零生信基础** | 7/8（跳过假说评分整页）| 30 min |
| #2 | iPhone 微信 / 课题组 | 匿名 / 硕士 2 / **生信基础（跑过 GTDB-Tk / CheckM）+ 中等代码** | 7/8（跳过假说评分整页）| 30 min |

---

## 🚨 用户提供的 3 个偏差校准（关键！）

读分析前必须先记住的"偏差因子"：

### 偏差 1：社交期望偏差（同门师弟师妹给面子）

**两位受访者都是用户的师弟师妹**，礼貌打高分是常态。原始数据的"高分"应当**整体下调 1-2 档**才接近真实评价：

| 指标 | 原始报告 | 校准估计（去社交偏差）| 真实含义 |
|---|---|---|---|
| SUS 均值 | 71.3（"Good" >68）| **~62-65**（接近边界）| 可用性中等偏下 |
| NPS | +50（1 promoter + 1 passive）| **~0**（弱 passive）| 推荐意愿不强 |
| Q66 卖点价值 | 4.5 | **~3.5**（受访方向匹配偏差）| 中性 |
| 各种 4-5 分 Likert | 看似很好 | 校准后 3-4 | 中等 |

**论文里不能直接报这些数字** —— 审稿人会看到 n=2 + 课题组内一秒识破。

### 偏差 2：图表加列建议是噪声，不是信号

Q63 / Q67 问"想加什么图"，用户列出的火山图、维恩图、桑基图、UpSet plot 等 ——
**用户未必真的知道这些图用于什么分析、自己研究里需不需要**。这些是"听过名字
的图"，不是"经验上需要的图"。

→ **这些功能扩展建议是 noise，不应该作为优先级输入到 backlog**

### 偏差 3：回收率 20% 本身是最重要的信号

- 内测群 ~10 人发了 v1，**8 人沉默**
- 沉默 ≠ 没意见
- 沉默 = "**没有动机抽 30-40 分钟来填**"
- 即使后续 v2 投到网络小范围用户，回收率仍可能低
- 反映的根本问题：
  1. **EnvMeta 是工具型软件**，必须在用户"正在面对那个具体痛点"时才有动机
  2. **测序数据下游可视化是小众且时段性需求**（一个课题跑一次，不是日常）
  3. **进入门槛高** —— 用户得先有 MAG 数据才能跑全套，符合条件的人少

---

## 🔴 校准后的真负面信号（基于实际数据 + 偏差校准）

以下 4 条是**真的需要关注的问题**：

### 1. HTML 离线部分功能失效（匿名实测）

- 数据：Q28 匿名选 B"部分功能失效"（陈源选 A"全可用"）
- 可信度：⭐⭐⭐⭐ 高 —— 唯一一个具体可复现的 product bug
- 推测：iPhone 微信内置浏览器（WKWebView + Safari iOS 18.7）可能对 D3 v7 的某些 API 兼容性差，或者文件下载到 iPhone 后用别的 app 打开时部分功能（比如 SVG 导出按钮）失效
- 后续：建议在 iOS Safari + iPhone 微信浏览器跑一次 HTML，确认是哪类失效（导出按钮？拖拽？hover？）
- ⚠️ 不要混淆：用户上一轮自测是**Mac/PC 桌面端**，那个 fix 已修复（commit 6667ae3 切组后 chem-link 失效）。本条是**移动端 / WebView**的另一个潜在问题

### 2. 循环图生成性能慢（两人独立确认）

- 数据：陈源 Q80 + 匿名 Q67 都提到"生成循环图等图的时间过长"
- 可信度：⭐⭐⭐⭐ 高 —— 唯一两人独立同向负面反馈
- 推测：循环图推断管线包括去偏 + 999 次置换 + 敏感度扫描，对 168 MAG × 57 KO 数据集算 30-60 秒
- 后续：profile `envmeta/geocycle/inference.py`，看哪一步占大头；可加 progress bar 改善感知（不一定要真的优化算法）

### 3. 零基础用户读不懂 PCoA 解读 expander

- 数据：陈源（环境工程，零生信，零代码）Q21 选 C "部分能，关键术语不懂"，Q22 不懂 Bray-Curtis
- 对比：匿名（生信背景）Q21 选 A "完全能独立解读"，并独立写出统计学正确的解读
- 含义：**"如何解读"expander 对生信用户有效，对真零基础用户失效**
- 后续：解读 expander 可加"零基础速查"分支（用日常语言讲 Bray-Curtis = "看两个样本物种组成有多不一样的距离"），但不是急事

### 4. 核心卖点对零基础用户没吸引力

- 数据：
  - Q66 循环图+假说价值：陈源选 B"独特但暂时用不上"
  - Q82 替代工作流：陈源选 D"暂时不替代任何东西"，匿名选 B"替代我自己写的 R/Python 脚本"
- 含义：**EnvMeta 真正能感受到价值的用户群只有"会写脚本的研究生"**
- 这不是 bug，是市场定位 —— 论文 Methods 里要诚实说明这个 user persona

---

## 🟢 校准后仍稳健的绿色信号

去除社交偏差后还能站住的：

1. **匿名（生信）独立写出 PCoA 标准解读**："解释度主要看 PC1，总解释了 70 多，PERMANOVA 群间显著分离"
   → "如何解读"expander **对生信用户**有效。论文可作为 case study evidence，但**必须明确标注用户类型**
2. **匿名 Q43 awk 命令答错**（选 Python，正确 awk）—— 即使生信背景用户也没在文档里找到 awk 命令；说明 docs/data_preparation_zh.md 里的 awk 命令需要更醒目（标题加 ⭐ 或加 example block）
3. **Q1.1 找到样例数据按钮**：陈源 5 秒 / 匿名 30 秒 —— 首屏一键加载达到设计目标
4. **匿名用了 EnvMeta 后选了 B"替代我自己写的 R/Python 脚本"** —— 对目标用户群的真实价值得到验证

---

## 📋 待办优先级（投稿前）

按 ROI 排序：

| 优先级 | 任务 | 工时 | 理由 |
|---|---|---|---|
| 🟥 P0 | **iOS Safari / 微信浏览器 HTML 兼容性自检** | 1-2h | 唯一具体可复现的 bug 报告 |
| 🟥 P0 | **第二数据集复现**（Oak Ridge 铀 / 阿拉斯加冻土）| 2-3 天 | iMeta 投稿硬指标 |
| 🟧 P1 | **循环图性能 profile + 加 progress bar** | 半天-1 天 | 双人独立确认的负面 |
| 🟧 P1 | **English README + LICENSE + Zenodo DOI** | 4-6h | iMeta 投稿硬指标 |
| 🟨 P2 | docs/data_preparation_zh.md 里 awk 命令加视觉强调 | 1h | 单人证据但明确可改 |
| 🟦 P3 | "如何解读"expander 加零基础速查分支 | 半天 | 单人证据，可缓 |

---

## 🎯 给论文 user study 部分的诚实建议

n=2 + 社交偏差不能写 SUS / NPS 数字（审稿人秒拒）。改写成 **case study 风格**：

```
We conducted preliminary usability evaluations with 2 graduate student
users with contrasting backgrounds (one with no programming or
bioinformatics experience, environmental engineering major; one with
intermediate scripting and CLI experience in metagenomics). Both
completed 7/8 tasks within 30 minutes. The non-programmer successfully
loaded sample data, generated PCoA plots and exported interactive HTML,
but reported difficulty interpreting statistical terms (e.g.,
Bray-Curtis distance). The bioinformatics-experienced user
independently reproduced the standard PCoA interpretation including
PC1 explanation rate (~70%) and PERMANOVA significance, suggesting
that in-app interpretation expanders effectively support users with
prior statistical exposure. Both users independently reported that
geochemical cycle figure rendering was perceived as slow (an explicit
performance optimization target for v0.9). The intermediate user
indicated that EnvMeta would replace her custom R/Python scripts for
this analysis, validating the tool's value proposition for the target
user group.

Limitations: The two evaluators were members of the authors' research
group, introducing social desirability bias. SUS/NPS scores were
collected but are not reported as the sample size (n=2) is below the
threshold for meaningful aggregate statistics. A larger external
recruitment is currently underway via online survey targeting
non-affiliated environmental microbiology researchers (preliminary
recruitment URL provided in repository README).
```

不报数字 + 讲两个用户的故事 + 承认 limitation + 指向后续招募 —— 这比硬撑 n=2 SUS 安全 100 倍。

---

## 后续动作

1. ✅ v1 反馈存档（含校正版分析）
2. ✅ 设计 v2 双轨问卷（[online_survey_design_v2.md](online_survey_design_v2.md)）
3. ✅ v2 投放小范围网络用户
4. ⏳ 等 v2 拿到 n ≥ 8 后再综合分析
5. ⏳ iOS / 微信浏览器 HTML 兼容性自检（独立任务）
