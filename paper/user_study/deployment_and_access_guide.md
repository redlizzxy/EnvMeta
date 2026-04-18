# EnvMeta 部署 + 访问指南

> 作者用来把 EnvMeta 部署到公网给师弟师妹访问，+ 师弟师妹访问说明。
> 场景：10 人线上体验 + 问卷，周日晚发射。

---

## Part 1：作者部署（~30 min）— 推荐 Streamlit Cloud

### 为什么选 Streamlit Cloud

- ✅ **零维护**：推送代码自动部署
- ✅ **公网 HTTPS URL**：师弟师妹直接访问，无需安装
- ✅ **个人免费**：1 GB RAM / 无限带宽（够 10 人同时体验）
- ✅ **支持 private repo**：不用设为 public
- ✅ **Python 3.11/3.12**：EnvMeta 原生支持

### 前置检查清单

- [ ] GitHub 仓库 `redlizzxy/EnvMeta` 最新 commit 已推送
- [ ] 根目录有 `requirements.txt` ✅ 已有
- [ ] 根目录有 `app.py` ✅ 已有
- [ ] `tests/sample_data/` 文件齐全（师弟师妹测试要用）
- [ ] 整个 repo < 1 GB（现有大约 50 MB，远小于）

### 部署步骤

**1. 访问 https://share.streamlit.io**
- 点击右上「Sign up」→ 选「Continue with GitHub」
- 授权 Streamlit 访问你的 GitHub（**选择给它访问 EnvMeta repo 的权限**）

**2. 创建新应用**
- 回到首页点「New app」
- 表单填写：
  - **Repository**: `redlizzxy/EnvMeta`
  - **Branch**: `master`
  - **Main file path**: `app.py`
  - **App URL (可自定义)**: 建议填 `envmeta` 或 `envmeta-beta`，会得到 `https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/`
  - 「Advanced settings」→ Python version 选 **3.11**

**3. 点击「Deploy」**
- 第一次部署约 5-10 分钟（安装依赖）
- 看日志面板实时输出，如果 `pip install` 失败会显示错误

**4. 常见部署问题**

| 问题 | 解决 |
|---|---|
| `scikit-bio` 安装失败 | requirements.txt 加一行 `cython` 放 scikit-bio 之前 |
| 启动时 ModuleNotFoundError | 检查 requirements.txt 有没有漏掉本地开发时用的依赖 |
| RAM 不够 | Advanced settings 提高到 2.5 GB（免费层可选） |
| 启动慢 > 30 秒 | 正常，第一次 warm-up 要几分钟 |

**5. 测试公网访问**
- 部署成功后得到 URL：`https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/`
- 用手机 4G 访问（断开 WiFi 验证真的是公网可访问）
- 走一遍任务 1 验证样例数据加载能用

### 部署后维护

- **推新 commit 到 master 会自动重新部署**（约 2-3 min）
- 如需重启：Streamlit Cloud dashboard → Reboot
- 查看访问日志：dashboard → Logs

---

## Part 2：备选方案（若 Streamlit Cloud 不行）

### 方案 B：本地 + ngrok 隧道（作者电脑持续开机）

**优势**：不依赖任何第三方托管，保留 private repo
**劣势**：作者电脑必须一直开着 + 联网

```bash
# 1. 正常启动 streamlit
conda activate envmeta
cd D:/workdata/envmeta
streamlit run app.py  # 默认 http://localhost:8501

# 2. 另开终端装 ngrok 并暴露
# 下载 ngrok: https://ngrok.com/download（注册免费号）
ngrok authtoken YOUR_TOKEN
ngrok http 8501

# 得到类似 https://abcd-xx-xx-xx.ngrok-free.app 的公网 URL
```

把 ngrok URL 发给师弟师妹即可。**ngrok 免费版每次重启 URL 会变**，需要每天发一次新链接；或升级到 ngrok 付费版 $8/月 拿固定域名。

### 方案 C：实验室服务器 Docker 部署

若有实验室内网服务器（或校内公网 IP）：

```bash
# 服务器端
git clone https://github.com/redlizzxy/EnvMeta.git
cd EnvMeta
docker build -t envmeta .  # 需要先写 Dockerfile（未提供）
docker run -p 8501:8501 envmeta
# 访问 http://服务器IP:8501
```

**注意**：EnvMeta 现在没有 Dockerfile。临时写一个：
```dockerfile
FROM python:3.11-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
EXPOSE 8501
CMD ["streamlit", "run", "app.py", "--server.address=0.0.0.0"]
```

### 方案 D：本地安装（师弟师妹愿意装）

最不友好但最独立的方案。对生信有基础的师兄师姐可用。详见 Part 3「本地安装」章节。

---

## Part 3：师弟师妹访问指南（复制到问卷里）

### 标准流程（推荐所有人）

**1. 访问 EnvMeta**

用电脑浏览器（Chrome / Edge / Firefox）打开：

```
https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/
```

（首次访问若显示「App is sleeping」，点击「Wake up」等 30 秒即可）

**2. 加载样例数据**

在首页找到「🚀 快速体验」区域 → 点击「📦 加载砷渣修复示例数据」按钮 → 会自动跳转到「文件管理」页面并加载 18 个识别好的文件。

**3. 按问卷任务操作**

回到手机问卷，按 8 个任务描述在电脑浏览器里操作，并截图上传到问卷。

**4. 遇到问题怎么办**

- 页面转圈超过 60 秒 → 刷新浏览器
- 报错 / 截图异常 → 跳过该任务，在问卷备注里说明 + 截图
- 全程不懂 → 微信群@师兄 / 发邮件 [作者邮箱]

### 高级：本地安装（可选，仅生信基础用户）

**适用场景**：希望用自己的数据跑，或离线使用。

**1. 装 Miniconda**

Windows：https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
Mac：https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh

安装时**勾选「Add to PATH」**。

**2. 克隆 + 安装**

```bash
git clone https://github.com/redlizzxy/EnvMeta.git
cd EnvMeta
conda create -n envmeta python=3.11 -y
conda activate envmeta
pip install -r requirements.txt
```

**3. 启动**

```bash
conda activate envmeta
cd EnvMeta
streamlit run app.py
```

浏览器自动打开 `http://localhost:8501`。

**4. 用自己的数据**

读 `docs/data_preparation_zh.md`（在「导出中心 → 文档」tab 可下载），对齐 EnvMeta 输入格式。

---

## Part 4：维护清单（作者测试期间）

### 发布当天

- [ ] GitHub 最新 commit 推到 master
- [ ] Streamlit Cloud 部署成功，公网可访问
- [ ] 自己手机 + 电脑各走一遍完整 30 分钟流程
- [ ] 腾讯问卷二维码 + 短链保存
- [ ] 微信群发送「版本 A 海报文案 + 二维码」

### 测试期间（1 周）

- [ ] 每日看 Streamlit Cloud dashboard 确认 app 运行正常
- [ ] 每日看腾讯问卷后台回收数
- [ ] 第 3 天发催收（版本 A 催收模板）
- [ ] 第 6 天发最后 24h 提醒

### 收尾（1 周后）

- [ ] 导出腾讯问卷数据（Excel）
- [ ] 计算 SUS 得分 + NPS 分 + 关键 Likert 题均值
- [ ] 质性反馈按 4 象限分类（P0/P1/P2/P3）
- [ ] 写入 `paper/user_study/results.md`
- [ ] 论文 Methods/Results 段落初稿
- [ ] Streamlit Cloud app 保留运行（作 SI demo 用）

---

## Part 5：常见问答

**Q1: Streamlit Cloud 每天只能 1000 次访问？**
- 不，免费层个人用无访问次数硬限制，只有资源使用限制。10 人测试无压力。

**Q2: 师弟师妹截图上传到腾讯问卷能保存吗？**
- 腾讯问卷免费版附件题单题 ≤ 3 张 / 单张 ≤ 5 MB，够用。

**Q3: 作者中途改代码怎么办？**
- 建议**冻结 master**，师兄测试期间不推非必要改动。紧急 bug 直接推 master，Streamlit Cloud 会 3 min 内自动更新。

**Q4: Streamlit Cloud 首次访问慢怎么办？**
- 免费 tier 的 app 闲置 7 天会 sleep，首次访问需要 30-60 秒 wake up。测试期间可以每天访问一次保持活跃。

**Q5: 有没有办法让师弟师妹用自己的数据跑？**
- 测试阶段不建议（数据格式差异大，支持成本高）。用内置样例数据即可。正式版发布后指引他们按 `data_preparation_zh.md` 自备数据。

---

## 立即执行路径

**最快路径**（今晚部署 + 明天发射）：

```
18:30 - 19:00   Streamlit Cloud 登录 + 新建 app + 开始部署
19:00 - 19:30   等部署 + 自己测试访问
19:30 - 20:00   二维码生成 + 海报文案收尾
20:00           微信群发射问卷链接 + 工具 URL
```

祝周日汇报顺利！🎯
