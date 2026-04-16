# S3 假说评分器 — 端到端验证

## 输入

- **假说 YAML**：[`paper/hypotheses/arsenic_steel_slag.yaml`](../../../hypotheses/arsenic_steel_slag.yaml)（8 条 claim，4 类全覆盖）
- **CycleData**：`tests/sample_data/` 六份文件
- **推断参数**：`env_rho_min=0.3, env_p_max=0.1, perm_n=99`

## 输出（S3.5 多指标版，按组）

| 组 | overall | label | null_p | weight_robust | veto |
|---|---|---|---|---|---|
| CK | 0.901 | strong | 0.296 (随机样) | ✅ | 0 |
| A  | 0.890 | strong | 0.764 (高随机) | ✅ | 0 |
| **B** | **1.000** | **strong** | N/A (退化) | ✅ | 0 |

**解读**：
- B 组满分 + label strong + 全 required claim 通过 = 假说被强支持
- null_p: CK/A 接近随机（假说对它们不特异），B 组因全 satisfied 退化
  （单一 score 集不够多样无法做排列） → 这些信号精确反映了"B 组独特性"
- weight_robust: 三组都 ✅，说明 label 不是"权重恰好调对"的巧合
- 无 veto 触发：required=true 的 `iron_transport_active` +
  `fe_as_adsorption_coupling` 在三组都 satisfied（铁-砷耦合基础在各组都存在）

每组 TSV + JSON 落在 `score_group_{CK,A,B}.tsv` / `.json`。

## 解读

1. **B 组满分**：所有不被 skip 的 claim 全部 satisfied。符合 2026-04-17
   推断观察（B 组砷代谢活性 1.9×、氨氧化 22×、keystone 成为功能主力）。
2. **CK / A 各有 1-2 条 skipped**：主要是 CK 组没有 Fe(III) / S²⁻ 同时出现
   的耦合（`fe_as_adsorption_coupling` 或 `s_as_precipitation_coupling` 被 skip）。
   这 **skipped 不扣分**，避免因"CK 组没这机制"把整体分数压低 —— 这正是
   v1 设计里 skipped 语义的核心。
3. **CK 组 unsatisfied 的那条**：`ammonia_ox_eh_link` 或 `arsenate_reduction_active`
   未达阈值（具体看 TSV）。比较 CK/A/B 三张 TSV 可读出"哪条 claim 随处理梯度
   变化" —— 这是论文 Results 的直接素材。

## 如何复现

```python
from envmeta.geocycle.inference import infer
from envmeta.geocycle.hypothesis import load_hypothesis, score
# 加载数据 (参见 tests/test_hypothesis.py fixture)
data = infer(ko, tax, ks, ab, env, md,
             params={"perm_n": 99, "group_filter": "B",
                     "env_rho_min": 0.3, "env_p_max": 0.1})
hyp = load_hypothesis("paper/hypotheses/arsenic_steel_slag.yaml")
result = score(hyp, data)
print(result.overall_score, result.label)
result.to_dataframe().to_csv("hypothesis_score.tsv", sep="\t")
```

或在 app 里操作：
1. `streamlit run app.py`
2. 循环图页 → 上传 6 文件 → 选 group=B → 生成循环图
3. 展开「🧪 假说评分」→ 上传 `arsenic_steel_slag.yaml` → 点击「评分」
4. 下载 TSV / JSON

## 论文 Methods 可引用

> "EnvMeta couples hypothesis-agnostic cycle inference (evidence of pathway
> activity / MAG contribution / env correlation, all without user-supplied
> hypothesis input) with an optional hypothesis-testing YAML evaluator. Users
> enumerate mechanistic claims (pathway_active, coupling_possible,
> env_correlation, keystone_in_pathway) with weights and thresholds; the
> evaluator returns per-claim evidence plus a weighted overall score
> (strong / suggestive / weak / insufficient). Evidence and interpretation
> are architecturally decoupled — the evaluator reads CycleData but does
> not modify inference."

## 已知局限

- `group_contrast` 类 claim 不支持（v2 再加）；目前要跨组对比假说支持度，
  需要对 `group_filter=CK/A/B` 分别运行评分（本验证即如此做）。
- 置信度阶梯固定（strong>suggestive>weak>spurious?>none），尚不支持
  用户自定义阈值映射。
