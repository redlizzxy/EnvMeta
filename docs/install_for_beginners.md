# EnvMeta 本地安装指南（新手版）

> 👤 **适合谁**：没装过 conda / 没用过命令行 / 也想把 EnvMeta 装在自己电脑上的师弟师妹
> ⏱ **预计耗时**：20-30 分钟（其中等下载 10 分钟）
> 💾 **磁盘空间**：需要 **约 3 GB**

## 🔍 开始前 60 秒自检

### 你需要什么？
- 一台 Windows 10/11 或 Mac 的电脑（Linux 也可，本教程以 Windows 为例）
- 有管理员权限（能装软件）
- 网络能正常访问（建议校园网 / 家里 WiFi，不要用手机热点，有条件上外网）

### 不需要什么？
- ❌ 不需要会编程
- ❌ 不需要懂 Python
- ❌ 不需要 GitHub 账号

---

## 🚦 三步走

**第 1 步**：装 Miniconda（一个管理 Python 环境的工具）
**第 2 步**：下载 EnvMeta 代码 + 自动装好所有依赖（复制粘贴 3 行命令）
**第 3 步**：打开 EnvMeta（以后每次就这 2 行命令）

---

## 第 1 步：安装 Miniconda（~10 分钟）

### 1.1 下载

**Windows 用户**点这个链接下载：
https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe

**Mac 用户**点这个：
- 苹果 M 芯片（2020 年之后的 Mac）：https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.pkg
- Intel 芯片：https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.pkg

### 1.2 安装（Windows 版，Mac 用户看 1.3）

双击下载的 `.exe` 文件，一路点「Next」，**但遇到这个界面一定要打勾**：

```
☐ Add Miniconda3 to my PATH environment variable    ← 这条默认没勾，你要勾上！
☑ Register Miniconda3 as my default Python          ← 默认已勾
```

如果安装程序警告说"不推荐加 PATH"，**忽略警告**，继续勾选。这对新手非常重要。

点「Install」，等进度条跑完，点「Finish」。

### 1.3 Mac 安装

双击 `.pkg` 文件，一路点「继续」→「同意」→「安装」，输入 Mac 密码，等完成。

### 1.4 验证装好了

**Windows**：按 `Windows 键 + R` → 输入 `cmd` → 回车，黑色终端打开
**Mac**：按 `Cmd + 空格` → 输入 `terminal` → 回车

在终端里输入：

```
conda --version
```

回车。如果看到类似 `conda 24.x.x` 这样的输出，**装好了**。

如果提示「不是内部或外部命令」（Windows）或「command not found」（Mac）：
- 关掉终端，**重新打开一个新的**再试一次
- 若还不行：Windows 在开始菜单搜「Anaconda Prompt」打开那个黑色终端；Mac 重启终端

---

## 第 2 步：下载 EnvMeta + 安装依赖（~10 分钟）

**完全不需要理解代码。复制粘贴下面 5 行命令，一次一行，每行回车后等它跑完再敲下一行。**

> ⚠️ Windows 用户请用第 1.4 步打开的那个黑色终端（`cmd` 或 `Anaconda Prompt`）。
> Mac 用户继续用 Terminal。

**命令 1**：创建 EnvMeta 专用环境

```
conda create -n envmeta python=3.11 -y
```

看到「done」并且命令行显示「$」或「>」提示符就是成功。耗时约 1 分钟。

**命令 2**：激活环境

```
conda activate envmeta
```

成功后你会看到命令提示符前面多了 `(envmeta)`。

**命令 3**：下载 EnvMeta 代码到桌面

Windows 版：
```
cd %USERPROFILE%\Desktop
git clone https://github.com/redlizzxy/EnvMeta.git
cd EnvMeta
```

Mac 版：
```
cd ~/Desktop
git clone https://github.com/redlizzxy/EnvMeta.git
cd EnvMeta
```

如果提示 `git: command not found`（Mac 可能会），先装 git：
- Mac：命令行输入 `xcode-select --install`，弹窗点同意，装完 5 分钟
- Windows：Miniconda 正常都带 git；若没带，手动去 https://git-scm.com/download/win 下载装一下

**命令 4**：安装 EnvMeta 的所有依赖（**这一步最耗时，约 5-8 分钟**）

```
pip install -r requirements.txt
```

屏幕会滚动很多行，看到最后一行 `Successfully installed ...` 就是成功。

**命令 5**：验证装好了

```
python -c "import streamlit; print('EnvMeta 依赖 OK')"
```

看到 `EnvMeta 依赖 OK` 就是装好了。

---

## 第 3 步：启动 EnvMeta（每次用就这几行）

还在刚才那个终端里（黑色窗口没关过的话），输入：

```
streamlit run app.py
```

浏览器会**自动**打开 `http://localhost:8501`，看到 EnvMeta 首页 = 启动成功。

### 以后每次想用 EnvMeta：

**打开终端**，依次输入（每行回车）：

Windows：
```
conda activate envmeta
cd %USERPROFILE%\Desktop\EnvMeta
streamlit run app.py
```

Mac：
```
conda activate envmeta
cd ~/Desktop/EnvMeta
streamlit run app.py
```

### 关闭 EnvMeta

- 关掉浏览器
- 回到终端，按 `Ctrl + C`（Windows）或 `Control + C`（Mac）停止服务
- 关闭终端

---

## 🆘 常见问题

### Q1：`conda: command not found`

**原因**：Miniconda 没加到 PATH，或者没重启终端。

**解决**：
- 关掉所有终端窗口
- Windows 从开始菜单搜「Anaconda Prompt」打开那个终端
- Mac 彻底退出 Terminal（Cmd+Q）再打开

### Q2：`pip install` 中途卡住 / 下载慢

**原因**：国内访问 PyPI 慢。

**解决**：换清华源：

```
pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple
```

### Q3：浏览器没自动打开

**原因**：浏览器设置问题。

**解决**：手动复制终端显示的 URL（通常是 `http://localhost:8501`）粘贴到浏览器地址栏。

### Q4：`streamlit run app.py` 报错 "address already in use"

**原因**：8501 端口被占用。

**解决**：

```
streamlit run app.py --server.port 8502
```

然后浏览器访问 `http://localhost:8502`。

### Q5：装完但页面白屏 / 加载一直转圈

**原因**：第一次加载 Streamlit 需要几秒；或者有模块没装全。

**解决**：
- 等 30 秒再刷新
- 仍不行：回终端看报错，**最后一行**报错信息发给师兄

---

## 💬 实在搞不定怎么办

不要反复折腾浪费时间！两个选择：

1. **用在线版**（推荐）：直接访问 `https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/`，啥都不用装
2. **找我**：微信群截图报错，或发邮件 [18872605913@163.com]；把下面 3 样发给我就能远程定位：
   - 你的操作系统（Windows 10/11 还是 Mac）
   - `conda --version` 的输出
   - 报错的完整截图（**最后一行很重要**）

---

## 🎉 装好了之后干什么

1. 首页点「📦 加载砷渣修复示例数据」体验完整流程
2. 打开左侧「📚 图表选择向导」用研究问题找分析
3. 读「导出中心 → 文档」里的 `data_preparation_zh.md`，了解如何用自己的数据

祝使用愉快！
