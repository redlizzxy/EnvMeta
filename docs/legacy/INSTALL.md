# EnvMeta 环境安装指南

## 第一步：安装 Miniconda

1. 浏览器打开下载地址：
   https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe

2. 运行安装程序，按以下选项操作：
   - 安装类型：Just Me（仅当前用户）
   - 安装路径：默认 C:\Users\REDLIZZ\miniconda3（不用改）
   - ✅ 勾选 "Add Miniconda3 to my PATH environment variable"
   - ✅ 勾选 "Register Miniconda3 as my default Python"
   - 点击 Install，等待完成

3. 安装完成后，**关闭并重新打开终端**（让 PATH 生效）

## 第二步：创建项目 conda 环境

打开终端（PowerShell 或 CMD），依次执行：

```bash
# 创建 envmeta 环境（Python 3.11）
conda create -n envmeta python=3.11 -y

# 激活环境
conda activate envmeta

# 进入项目目录
cd D:\workdata\envmeta

# 安装所有依赖
pip install -r requirements.txt

# 验证安装
python -c "import streamlit; print(f'Streamlit {streamlit.__version__} OK')"
python -c "from envmeta import __version__; print(f'EnvMeta v{__version__} OK')"
```

## 第三步：启动应用

```bash
conda activate envmeta
cd D:\workdata\envmeta
streamlit run app.py
```

浏览器会自动打开 http://localhost:8501，看到 EnvMeta 主页面即成功。

## 日常使用

每次打开终端后，只需：
```bash
conda activate envmeta
cd D:\workdata\envmeta
streamlit run app.py
```
