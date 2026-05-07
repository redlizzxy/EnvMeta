# 历史学位论文绘图配置（thesis legacy plotting configs）

> 来自原研究的 R / matplotlib 主题文件，**不被 EnvMeta app 使用**。
> EnvMeta 的绘图风格独立维护在 `envmeta/params/` 和 `envmeta/analysis/_mag_common.py`。

## 文件

- `__init__.py` — 模块占位
- `plot_config.py` — matplotlib 配色 / 字体配置（论文 EN 版）
- `theme_thesis.R` — ggplot2 主题（论文中文版）

## 用途

仅供论文复现 + 与 EnvMeta 渲染对照（论文 Methods 里"R 侧侧对比验证"用得到）。

## 整理时间

2026-05-07 从根目录 `config/` 迁移到 `paper/config_thesis/`。
