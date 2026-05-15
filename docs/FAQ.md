# 常见问题（FAQ）

> 在线版 / 本地版常见问题与解决方法。安装相关问题见
> [install_for_beginners.md](install_for_beginners.md)；数据格式问题见
> [data_preparation_zh.md](data_preparation_zh.md)。

---

## 🌐 在线版

### Q1：第一次打开网页很慢，要等多久？

**答**：Streamlit Cloud 免费版有"休眠 → 唤醒"机制：

- 如果应用 7 天内没人访问 → 进入休眠
- 第一次访问会看到「App is sleeping…」+「Yes, get this app back up!」按钮
- 点按钮后**等 30-90 秒**就好，期间 Cloud 在容器里重新启动 Python 进程
- 唤醒后 1 小时内访问都很快

**如果等了 2 分钟还没动静**：刷新一次（Ctrl+F5），或者直接换本地装（见 Q3）。

---

### Q2：网页显示「😟 Oh no. Error running app」怎么办？

**答**：这是 Streamlit Cloud 应用崩溃页面，**不是浏览器问题**。可能原因：

1. **依赖冲突**：Streamlit Cloud 默认 Python 版本升级后，某些库 wheel 不可用
   导致 build 失败。我们已通过 `runtime.txt` 锁定 Python 3.11；如果还遇到，
   说明又出新破坏，请到 [GitHub Issues](https://github.com/redlizzxy/EnvMeta/issues) 反馈。
2. **OOM（内存不足）**：免费版只有 1 GB RAM，多人同时跑大数据集会爆。见 Q4。
3. **临时网络抖动**：等 5-10 分钟自动恢复。

**快速验证应用是否在跑**：访问 `https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/healthz`，
返回 `{"status":"ok"}` 说明 Cloud 是好的，重新刷主页即可。

**如果一直崩**：直接装本地版（见 Q3），10 分钟搞定。

---

### Q3：本地装难吗？

**答**：不难。我们准备了 [小白安装指南](install_for_beginners.md)，
零基础 20 分钟装完：

1. 装 Miniconda（一键 installer）
2. `git clone https://github.com/redlizzxy/EnvMeta.git`
3. `conda create -n envmeta python=3.11 -y && conda activate envmeta`
4. `pip install -r requirements.txt`
5. `streamlit run app.py`

浏览器自动开 <http://localhost:8501>。**本地版没有内存限制 + 没有休眠 + 不会
被别人并发挤崩**，强烈推荐严肃使用者本地装。

---

### Q4：几个人能同时用在线版？会被挤崩吗？

**答**：Streamlit Cloud 免费版的硬约束是 **1 GB 内存 / 1 vCPU / 单容器多
用户共用**。EnvMeta 每用户内存占用（加载示例数据 + 跑几张图）大致：

| 同时在线人数 | 预计内存 | 状态 |
|---|---|---|
| 1-3 人 | 250-500 MB | ✅ 顺畅 |
| 5 人 | 700-900 MB | ⚠️ 紧 |
| 8+ 人 | 接近 / 超 1 GB | 🔴 大概率 OOM 崩溃 |

**为缓解并发崩溃**，v0.9.x 起在线版默认载入的是 **30 MAG 测试样例**（`tests/sample_data_demo/`），
而非完整 168 MAG 数据集。内存占用降到约 1/5，10 人**轻交互**也能撑住。

**师弟师妹批量测试场景建议**：

- **错峰**：2-3 人一组测，间隔 5-10 分钟。
- **本地装**：所有需要"跑完整 pipeline + 生成图保存"的人都装本地版。
- **在线版仅做** "瞄一眼 UI"、"看看长啥样" 的入口。

---

### Q5：在线版的示例数据是什么？能用它复现论文吗？

**答**：在线版载入的是 **30 MAG 子集**（`tests/sample_data_demo/`），从砷渣
修复研究完整数据集（168 MAG）抽样而来。**仅用于功能演示，不反映原研究的科学结论**。

要复现 EnvMeta 论文（Paper 3）的图：

- 本地装 + 用 `tests/sample_data/` 完整 168 MAG 数据集。
- 或下载论文 **Fork Bundle**（投稿时发布的 Zenodo DOI），含完整 KB / YAML /
  config / 数据，一键加载。

外部数据集（Liu 2023 / Grettenberger 2021 / Ayala 2020 等）数据 + 假说 YAML
随 Paper 3 Fork Bundle 发布（投稿后 Zenodo DOI）；本地完整 archive 见
`software/papers/paper3_envmeta/paper_archive/benchmarks/external/`（仅作者本机）。

---

### Q6：上传我自己的数据后，刷新一下就丢了？

**答**：是的。Streamlit 的 `session_state` 是按浏览器会话存的，刷新页面 = 新
会话 = 数据清空。**这是 Streamlit 框架的设计**，不是 EnvMeta bug。

**保存工作的办法**：

- 跑完分析后**导出 PDF / SVG / PNG / TSV**（每张图右下都有下载按钮）。
- 用「**Fork Bundle**」打包当前 KB + YAML + 上传的文件 → 下次加载 Bundle 恢复。
- 严肃使用走**本地版**，文件存本地不丢。

---

## 💻 本地版

### Q7：本地装好了，`streamlit run app.py` 报错 `ModuleNotFoundError: No module named 'xxx'`？

**答**：通常是依赖没装全。检查：

1. 确认在正确的 conda 环境：`conda activate envmeta`
2. 重装依赖：`pip install -r requirements.txt`
3. 还不行 → 删环境重来：`conda env remove -n envmeta && conda create -n envmeta python=3.11 -y && conda activate envmeta && pip install -r requirements.txt`

### Q8：Mac 装 scikit-bio / scipy 报编译错？

**答**：见 [小白安装指南 Q1-Q3](install_for_beginners.md#-常见问题)（conda ToS /
pip protobuf resolver / Streamlit 欢迎邮箱 三个坑）。

### Q9：matplotlib 中文字体显示方块？

**答**：Windows 默认装的 matplotlib 没含中文字体配置。两个方法：

1. **不用中文**：所有图改英文标签（最简）。
2. **装中文字体**：在 `~/.matplotlib/matplotlibrc` 加：

   ```
   font.sans-serif: SimHei, Arial Unicode MS, DejaVu Sans
   axes.unicode_minus: False
   ```

---

## 📊 分析

### Q10：「文件管理」里我的数据没被识别？

**答**：见 [数据准备指南](data_preparation_zh.md)。EnvMeta 通过表头规则识别
11 种文件类型，未识别通常是：

- 列名不在已知 alias 列表里（如把 `Sample_ID` 写成了 `样本ID`）→ 改列名。
- 文件编码异常（gb2312 / GBK）→ 另存为 UTF-8。
- 表头有合并行 / 多级表头 → 改成单行表头。

如果改不了 → 在「文件管理」面板**手动指定文件类型**也行（绕过自动识别）。

### Q11：循环图生成慢 / 卡死？

**答**：`cycle_diagram` 的 S2 置换检验跑 999 次 shuffle × Spearman × N_pathway ×
N_env，大数据集会慢：

| 数据规模 | 耗时（本地）|
|---|---|
| 30 MAG × 10 样本 × 4 env | < 1 秒 |
| 168 MAG × 10 样本 × 14 env | ~2 秒 |
| 1000 MAG × 87 样本 × 4 env | ~15 秒 |

如果 > 30 秒还没完，多半是浏览器不响应了，刷新重来一次即可。

### Q12：YAML 假说评分器返回 `INSUFFICIENT` 怎么办？

**答**：见 [假说写作指南](hypothesis_writing_guide.md)。常见原因：

- 某条 `required: true` 的 claim 不满足 → 触发 veto，整体被否决。
- 数据的 KEGG 注释覆盖不全（如只有 ROCker 14-gene 注释）→ 多条 claim 被
  skipped，分母减小但 required veto 仍触发。
- 阈值（`min_completeness`、`min_dominance_fraction`）设得过严。

每条 claim 的 `status` / `score` / `explanation` 字段会告诉你具体哪里没满足。

---

## 🐛 报 bug / 提需求

- GitHub Issues：<https://github.com/redlizzxy/EnvMeta/issues>
- 想要新功能 → 标 `enhancement` label
- 报 bug → 附上**操作步骤 + 错误截图 + 你用的什么数据集**

---

## 🔗 相关文档

- [README.md](../README.md) — 英文项目首页
- [README_CN.md](../README_CN.md) — 中文项目首页
- [install_for_beginners.md](install_for_beginners.md) — 小白安装指南
- [data_preparation_zh.md](data_preparation_zh.md) — 上游工具 → EnvMeta 输入映射
- [hypothesis_writing_guide.md](hypothesis_writing_guide.md) — YAML 假说写作教程
- [performance_zh.md](performance_zh.md) — 性能与硬件配置建议
