此代码为《Erdös-Sos 猜想在二分图变体中的证明》论文中，为处理小图情况而设计的批量化处理程序
# 双部图 Turán 数计算器（批量模式）**Bipartite Turán Number Calculator — Batch Mode**

一个高效的 Python 工具，用于计算**二部图极值数** `ex(n, m; H)`，支持批量处理任意多个 (n, m) 组合，自动生成极值图可视化与结构分析报告。

##  主要特性
- 支持任意禁图 H（经典的矩阵表示输入）
- 使用 ILP（PuLP + HiGHS/CBC）精确求解
- 批量计算任意范围的 (n, m)
- 每组结果自动生成**高清可视化 PNG**（带颜色渐变度数、连通分量高亮、结构标注）
- 自动结构分析：连通分量数、是否为完全二部图 K_{a,b}、缺失边/密度
- 生成**汇总表格图片**（含验证状态、时间统计）
- 中文界面与注释，适合中文用户直接使用
---
## 📦 安装依赖
pip install pulp networkx matplotlib numpy
推荐使用 Python 3.9+。
---
## 🚀 使用方法
```
python bipartite_turan_batch.py
```
按照提示依次输入：
1. 禁图 H 的行数 p 和列数 q
2. H 的邻接矩阵（空格分隔）
3. (n, m) 范围（支持多种格式）
4. 输出文件夹地址
5. ILP 时间限制（秒）
---
## 📝 支持的 (n, m) 输入格式
| 输入示例         | 说明                                   |
|------------------|----------------------------------------|
| `4 6`            | 单个组合                               |
| `3-6 4-8`        | 所有 n∈[3,6]、m∈[4,8] 且 n≤m 的对     |
| `3,4,5 5,6,7`    | 笛卡尔积（自动保留 n≤m）               |
| `3-7 *`          | 每个 n，m 从 n 到 n+4（自动上三角）    
---
## 📊 输出内容
运行后会在指定目录生成：
- `ex_{n}x{m}.png` —— 每个案例的极值图（带结构标注）
- `_summary_table.png` —— 所有案例的汇总表格
---
## 使用案例
目标：研究路径图P7在n≤m≤6的turan数与极值图
第一步：将P7转换为H矩阵（有多种表示方式）
【1 0 1 0】
【0 1 0 1】
【0 0 1 1】
第二步：运行程序并按照提示输入
<img width="856" height="398" alt="Snipaste_2026-03-24_11-12-36" src="https://github.com/user-attachments/assets/56c2c1ed-9276-41a0-877b-947d5f55b835" />
第三步：回车运行程序，结果展示：
<img width="1165" height="1000" alt="image" src="https://github.com/user-attachments/assets/821c9eda-8f3a-49b2-b769-0d141ebacbfe" />

<img width="1376" height="803" alt="image" src="https://github.com/user-attachments/assets/4fc5926b-33d6-4fb5-b015-c4887ea02854" />

<img width="1871" height="915" alt="image" src="https://github.com/user-attachments/assets/8e471f94-ad4a-47d0-865b-b2037e2bbb0a" />
