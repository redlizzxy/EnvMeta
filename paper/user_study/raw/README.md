# 用户问卷原始数据

> 不入 git 公共仓库（里面有受访者邮箱 / 姓名 / IP）。
>
> 此目录被 [.gitignore](../../../.gitignore) 排除（`paper/user_study/raw/` 已加，
> 见父目录 `.gitignore` 规则）。

## 现有文件

- `responses_2026-05-07_n2_encoded.csv` — v1 问卷前 2 份回收（编码版）
  - 来源：腾讯问卷后台 2026-05-07 导出
  - 编码：UTF-8 with BOM（开头 `EF BB BF`）
  - 格式：每行一受访者，闭卷题=数字编码，开放题=原文文本
  - 已分析：见 [../v1_archived_issues.md](../v1_archived_issues.md)（信号存档，未深度分析）

## 命名约定

```
responses_<YYYY-MM-DD>_n<N>_<encoded|raw>.csv
```

例：
- `responses_2026-05-15_n8_encoded.csv` — v2 第一轮回收（n=8）
- `responses_2026-05-25_n12_encoded.csv` — v2 第二轮（n=12，最终）

## 隐私处理

- 个人邮箱 / 姓名 / IP / 微信 OpenID 在论文 SI 里**全部脱敏**或不公开
- 论文公开的是聚合 SUS 分布 / NPS 分布 / 任务完成率，不公开个体逐行回答
- 致谢页姓名以参与者**自填**为准（默认匿名）
