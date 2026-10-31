# 线性代数速查手册（考前复习版）
# Linear Algebra Quick Reference Handbook

> 从 51 节双语讲稿提炼的**最小记忆单元**：每个概念一行，背完这张表 = 掌握全书骨架。
> 配套讲稿：`Linear Algebra and Its Applications(5th version)-双语.md`

---

## 第 1 章 线性方程组

| 概念 | 一句话 |
|---|---|
| 线性方程组解的三种可能 | 无解 / 唯一解 / 无穷多解（绝无"恰两解"） |
| 相容 ⟺ | 阶梯形无 $[0\cdots0\ |\ b\neq0]$ 行 |
| 自由变量 ⟺ | 无穷多解（在相容前提下） |
| 初等行变换 | 替换不变、交换变号、倍乘乘 k（对行列式的影响） |
| RREF | 唯一（附录 A 证明）；主元 1、上下全 0 |
| 线性组合 | $c_1\mathbf{v}_1 + \cdots + c_p\mathbf{v}_p$ |
| Span | 所有线性组合；问 b∈Span ⟺ 方程组有解 |
| $A\mathbf{x}$ | $A$ 各列以 $\mathbf{x}$ 分量为权的组合 |
| Theorem 4 (1.4) | $A\mathbf{x}=\mathbf{b}$ ∀b 有解 ⟺ 每行有主元 |
| 齐次解集 | 过原点（子空间）；非齐次 = 特解 + 齐次通解（平移） |

**关键流程**：写增广矩阵 → 行化简 → 判相容 → RREF → 基本变量用自由变量表示 → 参数向量形式

## 第 2 章 矩阵代数

| 概念 | 一句话 |
|---|---|
| 矩阵乘法 | $(AB)_{ij}$ = A 第 i 行 · B 第 j 列；源于复合变换 |
| 两大怪脾气 | 无交换律、无消去律 |
| 2×2 逆 | $A^{-1} = \frac{1}{ad-bc}\begin{bmatrix} d & -b \\ -c & a \end{bmatrix}$，$\det A = ad-bc$ |
| 求逆算法 | $[A\ |\ I] \sim [I\ |\ A^{-1}]$ |
| **IMT** | 可逆 ⟺ 12+ 等价条件（见文末完整版） |
| 初等矩阵 | 左乘 = 行变换；均可逆 |
| 分块乘法 | 块当数用；A 列分法 = B 行分法 |
| 列-行展开 | $AB = \sum_k \text{col}_k(A)\,\text{row}_k(B)$ |
| LU 分解 | $A = LU$（L 单位下三角存倍数）；存在 ⟺ 消元不换行；解方程：前代 + 回代 |
| 列昂惕夫 | $(I-C)\mathbf{x} = \mathbf{d}$；$(I-C)^{-1} \approx I + C + C^2 + \cdots$ |
| 子空间 | 含零向量、对加法/数乘封闭 |
| Col A | 各列张成；$\mathbf{b} \in \operatorname{Col} A$ ⟺ $A\mathbf{x}=\mathbf{b}$ 有解 |
| Nul A | $A\mathbf{x}=\mathbf{0}$ 解集 |
| 基 | 线性无关 + 张成 |
| Col A 的基 | **A 的主元列**（不是阶梯形的列！） |
| Nul A 的基 | 解 $A\mathbf{x}=\mathbf{0}$ 参数向量形式 |
| dim / rank | 基向量数 / $\dim\operatorname{Col} A$ = 主元列数 |
| **秩定理** | $\operatorname{rank} A + \dim\operatorname{Nul} A = n$（列数） |
| 基定理 | dim H = p ⟹ p 个无关向量自动是基 |
| 换基 | $[\mathbf{x}]_{\mathcal{C}} = P_{\mathcal{C}\leftarrow\mathcal{B}}[\mathbf{x}]_{\mathcal{B}}$；$(P_{\mathcal{C}\leftarrow\mathcal{B}})^{-1} = P_{\mathcal{B}\leftarrow\mathcal{C}}$ |

## 第 3 章 行列式

| 概念 | 一句话 |
|---|---|
| $\det$ 定义 | 余因子展开（任意行/列，选零多的） |
| 三角阵 | $\det$ = 对角元乘积 |
| 行变换影响 | 替换不变、交换变号、倍乘乘 k |
| **Theorem 4** | **可逆 ⟺ $\det \neq 0$** |
| $\det A^T$ | $= \det A$ |
| $\det(AB)$ | $= \det A \cdot \det B$；$\det A^{-1} = 1/\det A$ |
| ⚠️ 警示 | $\det(A+B) \neq \det A + \det B$ |
| 面积/体积 | 平行四边形/六面体 = $\|\det A\|$ |
| 变换伸缩 | $T(S)$ 面积 = $\|\det A\| \times$ S 面积 |
| 克拉默法则 | $x_i = \det A_i(\mathbf{b}) / \det A$（符号计算美，数值少用） |
| 伴随矩阵 | $A^{-1} = \operatorname{adj}A / \det A$ |

## 第 4 章 向量空间

| 概念 | 一句话 |
|---|---|
| 向量空间 | 集合 + 加法/数乘 + 10 公理（元素不必是数组） |
| 子空间 | 含零、封闭（继承全部公理） |
| Theorem 1 (4.1) | Span{任意向量组} 是子空间（最通用构造法） |
| ker T / range T | 核 = "压成零的输入"；值域 = "全部输出"（均为子空间） |
| 基 | 无关 + 张成；标准基 $\{\mathbf{e}_i\}$、$\{1,t,\dots,t^n\}$ |
| 张成集定理 | 张成集删冗余向量可收缩为基 |
| 坐标 $[\mathbf{x}]_{\mathcal{B}}$ | 基下的坐标；$p$ 维子空间 ↔ $\mathbb{R}^p$ 同构（Theorem 8） |
| 维数 | 所有基大小相同（Theorem 10） |
| 扩充定理 | 无关集可扩充为基；dim H ≤ dim V |
| 秩定理 (4.6) | 行秩 = 列秩 = rank A |
| 四子空间 | Row/Nul A 在 $\mathbb{R}^n$、Col/Nul Aᵀ 在 $\mathbb{R}^m$；A 在 Row A 上一一映到 Col A |
| 换基矩阵 | 列 = $[\mathbf{b}_j]_{\mathcal{C}}$ |

## 第 5 章 特征值

| 概念 | 一句话 |
|---|---|
| 特征向量 | $A\mathbf{x} = \lambda\mathbf{x}$（$\mathbf{x} \neq \mathbf{0}$）：方向不变的向量 |
| 特征值求法 | 解 $\det(A - \lambda I) = 0$（⚠️ 不能靠行化简求 λ） |
| 特征空间 | $\operatorname{Nul}(A - \lambda I)$ |
| 三角阵 | 特征值 = 对角元 |
| 0 是特征值 ⟺ | A 不可逆 |
| Theorem 2 (5.1) | 不同特征值的特征向量无关 |
| 代数重数 | 特征多项式根的重数 |
| 相似 | $B = P^{-1}AP$：同特征值（Theorem 4）；同特征值 ≠ 相似 |
| **对角化** | $A = PDP^{-1}$ ⟺ n 个无关特征向量 ⟺ Σ特征空间维数 = n（Theorem 5/7） |
| 对角化四步 | 特征值 → 特征向量 → 组 P → 组 D（对应！验算 AP=PD） |
| $A^k$ | $= PD^kP^{-1}$ |
| B-矩阵 | $[T]_{\mathcal{B}} = P^{-1}AP$：相似 = 同变换不同基（Theorem 8） |
| 复特征值 | $a \pm bi$ ⟺ 旋转 φ + 缩放 r；共轭成对 |
| Theorem 9 (5.5) | $A = PCP^{-1}$，$C = \begin{bmatrix} a & -b \\ b & a \end{bmatrix}$，$P = [\operatorname{Re}\mathbf{v}\ \operatorname{Im}\mathbf{v}]$ |
| 动力系统 | $\mathbf{x}_k = \sum c_j\lambda_j^k\mathbf{v}_j$；支配特征值定长期 |
| 吸引子/鞍点 | 两 \|λ\| < 1 / 一小一大；r=1 椭圆；复特征值 r<1 螺旋入 |
| 微分方程 | $\mathbf{x}' = A\mathbf{x}$ 通解 $= \sum c_j\mathbf{v}_je^{\lambda_j t}$ |
| 幂法 | 反复乘 A → 最大特征值方向（归一化） |

## 第 6 章 正交性

| 概念 | 一句话 |
|---|---|
| 内积 | $\mathbf{u}\cdot\mathbf{v} = \mathbf{u}^T\mathbf{v}$；正交 ⟺ 内积 0 |
| Theorem 3 (6.1) | $(\operatorname{Row} A)^\perp = \operatorname{Nul} A$、$(\operatorname{Col} A)^\perp = \operatorname{Nul} A^T$ |
| 勾股定理 | $\mathbf{u}\perp\mathbf{v}$ ⟺ $\|\mathbf{u}+\mathbf{v}\|^2 = \|\mathbf{u}\|^2 + \|\mathbf{v}\|^2$ |
| 正交基坐标 | $c_j = \mathbf{y}\cdot\mathbf{u}_j / \mathbf{u}_j\cdot\mathbf{u}_j$（一行公式！） |
| 投影公式 | $\operatorname{proj}_L\mathbf{y} = \frac{\mathbf{y}\cdot\mathbf{u}}{\mathbf{u}\cdot\mathbf{u}}\mathbf{u}$ |
| 正交矩阵 | 列标准正交 ⟺ $U^TU = I$ ⟺ $U^{-1} = U^T$；保长度保角度 |
| 分解定理 | $\mathbf{y} = \hat{\mathbf{y}} + \mathbf{z}$（$\hat{\mathbf{y}} \in W$、$\mathbf{z} \in W^\perp$）唯一 |
| **最佳逼近** | $\hat{\mathbf{y}} = \operatorname{proj}_W\mathbf{y}$ 是 W 中最近点（勾股定理证明） |
| Gram-Schmidt | $\mathbf{v}_{k+1} = \mathbf{x}_{k+1} - \operatorname{proj}_{W_k}\mathbf{x}_{k+1}$ |
| QR 分解 | $A = QR$（Q 正交列、R 上三角）；$R = Q^TA$ |
| **正规方程** | $A^TA\hat{\mathbf{x}} = A^T\mathbf{b}$（最小二乘 ⟺ 此解，Thm 13） |
| 唯一性 | ⟺ A 列无关 ⟺ $A^TA$ 可逆 |
| QR 解 LS | 解 $R\hat{\mathbf{x}} = Q^T\mathbf{b}$（回代，稳定） |
| 回归 | 设计矩阵 X；直线/抛物线/平面统一为 $X\boldsymbol\beta \approx \mathbf{y}$ |
| 积分内积 | $\langle f, g\rangle = \int_a^b fg$；傅里叶系数 = 内积/范数² |
| 柯西-施瓦茨 | $\|\langle\mathbf{u},\mathbf{v}\rangle\| \leq \|\mathbf{u}\|\|\mathbf{v}\|$ |

## 第 7 章 对称矩阵与二次型

| 概念 | 一句话 |
|---|---|
| Theorem 1 (7.1) | 对称 + 特征值不同 ⟹ 特征向量正交 |
| 谱定理 | 对称 ⟺ 可正交对角化（实特征值 + 正交特征向量基） |
| 谱分解 | $A = \sum\lambda_j\mathbf{u}_j\mathbf{u}_j^T$（秩 1 投影加权和） |
| 主轴定理 | $\mathbf{x} = P\mathbf{y}$ 消交叉项：$Q = \lambda_1y_1^2 + \cdots$ |
| 正定性 | 正定 ⟺ 特征值全正（Theorem 5）；顺序主子式判别更快 |
| 约束优化 | $\max_{\|\mathbf{x}\|=1}\mathbf{x}^TA\mathbf{x} = \lambda_{\max}$（Theorem 6） |
| 逐级约束 | 排除前 k−1 方向后最大 = $\lambda_k$（PCA 原理，Theorem 8） |
| **SVD** | $A = U\Sigma V^T$（任意矩阵！） |
| 奇异值 | $\sigma_i = \sqrt{\lambda_i(A^TA)}$；rank = 非零 σ 数 |
| 四子空间 | Row A = span{v₁..vᵣ}、Nul A = span{vᵣ₊₁..}、Col A = span{u₁..uᵣ}、Nul Aᵀ = span{uᵣ₊₁..} |
| 伪逆 | $A^+ = V_rD^{-1}U_r^T$；$\hat{\mathbf{x}} = A^+\mathbf{b}$ 最小范数 LS 解 |
| Eckart-Young | 截断 SVD = 最佳秩 k 近似（图像压缩/降噪/PCA） |

## 第 8 章 向量空间的几何

| 概念 | 一句话 |
|---|---|
| 仿射组合 | 权和 = 1（直线 = 两点仿射包） |
| Theorem 1 (8.1) | 仿射 ⟺ 减锚点后线性 |
| Theorem 4 (8.1) | 仿射组合 ⟺ 齐次形式 $[\mathbf{y};1]$ 的线性组合 |
| 仿射无关 | 平移后线性无关（Theorem 5 三判据） |
| 重心坐标 | 解齐次系统；全正 ⟺ 在三角形内（GPU 插值基础） |
| 凸组合 | 权和 1 + 全非负；凸包 = 最小凸集 |
| Carathéodory | 凸包每点 = 至多 n+1 个点的凸组合 |
| 超平面 | $\mathbf{n}\cdot\mathbf{x} = d$（法向量定方向） |
| 分离定理 | 紧致不相交 ⟹ 可严格分离（Theorem 12/13） |
| 极点 | n 个线性无关约束取等号（Theorem 15） |
| **LP 定理** | 线性泛函在多胞形有界 ⟺ 在极点取最大（Theorem 16） |
| Bézier | 二项分布权（凸组合）；过首末点、在凸包内 |
| G⁰/G¹/C¹ | 重合 / 切线同向 / 切向量全等 |
| de Casteljau | 逐次中点细分（数值稳定） |


---

## 附录：核心公式卡

**行列式**
- 2×2：$\det = ad - bc$；3×3：对角线法则（仅 3×3！）
- 化简求值：$\det A = (-1)^r \times$（主元积）

**特征/对角化**
- 特征方程：$\det(A - \lambda I) = 0$
- $A^k = PD^kP^{-1}$；$\det A = \prod\lambda_i$；$\operatorname{tr} A = \sum\lambda_i$

**正交/最小二乘**
- 投影：$\operatorname{proj}_W\mathbf{y} = UU^T\mathbf{y}$（标准正交基）
- 正规方程：$A^TA\hat{\mathbf{x}} = A^T\mathbf{b}$
- 误差：$\|\mathbf{b} - A\hat{\mathbf{x}}\|$ = b 到 Col A 的距离

**SVD**
- $\sigma_i = \sqrt{\lambda_i(A^TA)}$，$\kappa(A) = \sigma_1/\sigma_r$
- $A = \sum_{i=1}^r \sigma_i\mathbf{u}_i\mathbf{v}_i^T$

## IMT 完整版（24 条等价条件）

$A$ 可逆 ⟺ **(a)** $A^{-1}$ 存在 ⟺ **(b)** $A \sim I_n$ ⟺ **(c)** n 个主元 ⟺ **(d)** $A\mathbf{x}=\mathbf{0}$ 仅平凡解 ⟺ **(e)** 列线性无关 ⟺ **(f)** $\mathbf{x}\mapsto A\mathbf{x}$ 一对一 ⟺ **(g)** $A\mathbf{x}=\mathbf{b}$ ∀b 有解 ⟺ **(h)** 列张成 $\mathbb{R}^n$ ⟺ **(i)** 变换映上 ⟺ **(j,k)** 左/右逆存在 ⟺ **(l)** $A^T$ 可逆 ⟺ **(m)** 列构成基 ⟺ **(n)** $\operatorname{Col} A = \mathbb{R}^n$ ⟺ **(o)** $\dim\operatorname{Col} A = n$ ⟺ **(p)** $\operatorname{rank} A = n$ ⟺ **(q)** $\operatorname{Nul} A = \{\mathbf{0}\}$ ⟺ **(r)** $\dim\operatorname{Nul} A = 0$ ⟺ **(s,t)** 0 非特征值、$\det A \neq 0$ ⟺ **(u–x)** $(\operatorname{Col} A)^\perp = \{\mathbf{0}\}$、$(\operatorname{Nul} A)^\perp = \mathbb{R}^n$、$\operatorname{Row} A = \mathbb{R}^n$、n 个非零奇异值

## 考前 10 分钟流程

1. **看到"解方程组"**：增广矩阵 → 阶梯形 → 判相容 → 参数向量形式（先问自由变量）
2. **看到"求逆/可逆性"**：IMT 任挑一条验证（最常用：det、主元、特征值 0）
3. **看到"对角化"**：特征值 → 特征向量 → 数维数（够 n 个才能）→ 组 P、D
4. **看到"最优近似/回归"**：正规方程（或 QR）；误差 = 垂距
5. **看到"求最值方向"**：最大特征值方向（瑞利商）；主成分 = SVD 前 k 列
6. **看到"点在不在三角形内"**：重心坐标全正

---

*本手册从 51 节讲稿提炼 · 2026-10-31*
