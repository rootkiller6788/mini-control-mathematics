# Mini Control Mathematics（迷你控制数学）

一套**从零开始、零依赖的 C 语言实现**，涵盖控制理论所依赖的工程数学基础。每个子模块对应 MIT 及 Stanford（和其他顶尖大学）课程，覆盖复分析、常微分方程、差分方程、积分变换、线性代数、最优化、数值方法以及随机滤波——全部面向控制工程师构建。

## 子模块

| 子模块 | 主题 | 核心课程 |
|--------|------|----------|
| [mini-complex-analysis](mini-complex-analysis/) | 解析函数、Cauchy-Riemann、围道积分、留数演算、级数、Laplace/Fourier/Z 变换 | MIT 18.04, Stanford MATH 116 |
| [mini-difference-equations](mini-difference-equations/) | 离散时间动力学、差分算子、Lyapunov 稳定性、状态空间、Z 变换 | MIT 6.241, Stanford ENGR207B |
| [mini-laplace-z-transform](mini-laplace-z-transform/) | Laplace 变换、Tustin 双线性映射、Z 域数字滤波器设计（IIR/FIR） | MIT 6.003, Stanford EE 264 |
| [mini-linear-algebra-control](mini-linear-algebra-control/) | 能控性、能观性、矩阵分解、特征值方法、状态空间几何 | MIT 6.241, Stanford ENGR 205 |
| [mini-numerical-methods-control](mini-numerical-methods-control/) | 数值线性代数、ODE 积分、特征值计算、Riccati 求解器、非线性方法、插值 | MIT 18.335, Stanford CME 200 |
| [mini-ode-fundamentals](mini-ode-fundamentals/) | 初值问题理论、线性 ODE 系统、基本解矩阵、状态转移矩阵、控制应用 | MIT 18.03, Stanford CME 102 |
| [mini-optimization-basics](mini-optimization-basics/) | 凸分析、梯度方法、Newton 型方法（BFGS/DFP）、约束优化、单纯形法、内点法 | MIT 6.7201, Stanford EE 364a |
| [mini-probability-stochastic](mini-probability-stochastic/) | Kalman 滤波（线性、扩展、无迹）、蒙特卡洛方法、Markov 链、Poisson 过程 | MIT 6.241, Stanford AA 273 |

## 设计理念

- **零外部依赖** — 纯 C99/C11，仅使用标准库头文件
- **自包含子模块** — 每个子模块拥有独立的 `include/`、`src/`、`CMakeLists.txt` 和冒烟测试
- **理论到代码的映射** — 每个模块内联引用教材章节和讲义
- **控制工程导向** — 所有数学内容均以控制理论为动机（Nyquist 判据、极点配置、LQR、状态估计）

## 构建

每个子模块独立构建。使用 CMake 构建：

```bash
cd mini-complex-analysis
mkdir build && cd build
cmake ..
make
./smoke_test
```

需要 **C99 兼容编译器**和 **CMake ≥ 3.14**。

## 项目结构

```
0. mini-control-mathematics/
├── mini-complex-analysis/            # 解析函数、围道积分、留数、积分变换
├── mini-difference-equations/        # 离散时间系统、稳定性、状态空间、Z 变换
├── mini-laplace-z-transform/         # Laplace 变换、Tustin 映射、数字滤波器
├── mini-linear-algebra-control/      # 能控性、能观性、矩阵分解
├── mini-numerical-methods-control/   # ODE 积分、特征值求解器、Riccati、非线性方法
├── mini-ode-fundamentals/            # 初值问题理论、线性 ODE 系统、状态转移矩阵
├── mini-optimization-basics/         # 凸优化、梯度/Newton 方法、约束优化
├── mini-probability-stochastic/      # Kalman 滤波、蒙特卡洛、Markov 链、随机过程
├── .gitignore
├── README.md
└── README-CN.md
```

## 许可证

MIT
