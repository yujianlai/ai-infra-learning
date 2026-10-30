# Linear Algebra and Its Applications (5th Edition) 双语精读
# 《线性代数及其应用（第5版）》双语讲解 · 第1章

> 面向零基础读者。每个概念都先用中文直觉解释，再给英文原文精讲，公式与图片与教材完全一致。
> Based on Lay, Lay & McDonald, *Linear Algebra and Its Applications*, 5th ed., Pearson 2016.

---

## 第一章 线性代数中的线性方程
## Chapter 1 Linear Equations in Linear Algebra

**导入例子：经济学与工程学中的线性模型**
**INTRODUCTORY EXAMPLE: Linear Models in Economics and Engineering**

1949 年，哈佛大学的列昂惕夫（Wassily Leontief，1973 年诺贝尔经济学奖得主）把美国经济划分为 500 个"部门"，为每个部门写出一个线性方程，描述它的产出如何分配给其他部门。由于数据量巨大，他将其化简为 **42 个方程、42 个未知量** 的方程组，交给 Harvard Mark II 计算机求解。计算机嗡嗡作响 56 个小时后得出了答案——这是计算机最早的大规模数学建模应用之一。

> 核心思想：现实问题中变量太多，模型通常是**线性的**（linear），即可以用**线性方程组（systems of linear equations）**来描述。线性代数的重要性，与计算能力的增长同步提升。

---

## 1.1 线性方程组 Systems of Linear Equations

### 什么是线性方程？

**中文讲解**：含有变量 $x_1, \dots, x_n$ 的**线性方程**是可以写成下面形式的方程：

$$a_1x_1 + a_2x_2 + \cdots + a_nx_n = b$$

其中 $b$ 和系数 $a_1, \dots, a_n$ 是已知实数或复数。直观判断标准：**变量只以一次幂出现，没有相乘、没有根号**。

是线性的（可整理成标准形式）：

$$4x_1 - 5x_2 + 2 = x_1 \;\Rightarrow\; 3x_1 - 5x_2 = -2, \qquad x_2 = 2(\sqrt{6} - x_1) + x_3 \;\Rightarrow\; 2x_1 + x_2 - x_3 = 2\sqrt{6}$$

不是线性的：

$$4x_1 - 5x_2 = x_1x_2 \quad (\text{出现 } x_1x_2 \text{ 乘积}), \qquad x_2 = 2\sqrt{x_1} - 6 \quad (\text{出现根号})$$

**英文原文**：
> A **linear equation** in the variables $x_1, \dots, x_n$ is an equation that can be written in the form $a_1x_1 + a_2x_2 + \cdots + a_nx_n = b$, where $b$ and the **coefficients** $a_1, \dots, a_n$ are real or complex numbers, usually known in advance.

### 什么是线性方程组？

**中文讲解**：**线性方程组**（也叫线性系统 linear system）就是**同一组变量**满足的一个或多个线性方程的集合。例如：

$$2x_1 - x_2 + 1.5x_3 = 8 \qquad x_1 \quad\;\; - 4x_3 = -7$$

- **解（solution）**：一组数 $(s_1, s_2, \dots, s_n)$，代入后使每个方程都成立。例如 $(5, 6.5, 3)$ 就是上面方程组的一个解。
- **解集（solution set）**：所有解的集合。
- **等价（equivalent）**：两个方程组解集完全相同，就称它们等价。

### 两个方程、两个变量的几何意义

求两个二元线性方程的解，就是找**两条直线的交点**：

$$x_1 - 2x_2 = -1 \qquad -x_1 + 3x_2 = 3$$

解是唯一点 $(3, 2)$：

![FIGURE 1 唯一解](img/fig01_one_solution.png)

*FIGURE 1 Exactly one solution. 图1 唯一解：两直线相交于一点*

但两条直线也可能**平行**（无解）或**重合**（无穷多解）：

![FIGURE 2 无解与无穷多解](img/fig02_no_or_infinite.png)

*FIGURE 2 (a) No solution. (b) Infinitely many solutions. 图2 (a) 无解：两直线平行；(b) 无穷多解：两直线重合*

**重要事实（教材黄框，1.2 节将证明）**：
> A system of linear equations has **1.** no solution, or **2.** exactly one solution, or **3.** infinitely many solutions.
> 线性方程组只有三种可能：**无解**、**唯一解**、**无穷多解**。

- **相容（consistent）**：有唯一解或无穷多解
- **不相容（inconsistent）**：无解

### 矩阵记号 Matrix Notation

**中文讲解**：把方程组的核心信息压缩进一个矩形数组——**矩阵（matrix）**。对系统：

$$x_1 - 2x_2 + \;\;x_3 = \;0 \qquad 2x_2 - 8x_3 = \;8 \qquad 5x_1 \quad\;\; - 5x_3 = 10$$

系数矩阵与增广矩阵：

$$\begin{bmatrix} 1 & -2 & 1 \\ 0 & 2 & -8 \\ 5 & 0 & -5 \end{bmatrix} \qquad \left[\begin{array}{rrrr} 1 & -2 & 1 & 0 \\ 0 & 2 & -8 & 8 \\ 5 & 0 & -5 & 10 \end{array}\right]$$

- **系数矩阵（coefficient matrix）**：只装未知数的系数
- **增广矩阵（augmented matrix）**：系数矩阵右边再加一列常数项
- **规模（size）**：$m \times n$ 表示 $m$ 行 $n$ 列。上面的增广矩阵是 $3 \times 4$。

### 解线性方程组：消元思想

**基本策略**：把原系统换成一个**等价的、更容易解的系统**（解集不变）。

**三种基本运算**（对应增广矩阵的三种**初等行变换 elementary row operations**）：

1. **Replacement 替换**：把某一行替换为它自己与另一行的倍数之和
2. **Interchange 交换**：交换两行
3. **Scaling 倍乘**：把某一行所有元素乘以一个非零常数

> 关键性质：初等行变换都是**可逆的（reversible）**。
> 若两个线性系统的增广矩阵**行等价（row equivalent）**，则两系统解集相同。

**例1（Example 1）求解系统**：

$$\begin{aligned} x_1 - 2x_2 + \;\;x_3 &= 0 \\ 2x_2 - 8x_3 &= 8 \\ 5x_1 \quad\;\; - 5x_3 &= 10 \end{aligned}$$

**步骤（左边是方程，右边是对应增广矩阵，教材并排列出）**：

第1步：用第1方程消去其他方程的 $x_1$。第3行 += (−5)×第1行：

$$\left[\begin{array}{rrrr} 1 & -2 & 1 & 0 \\ 0 & 2 & -8 & 8 \\ 0 & 10 & -10 & 10 \end{array}\right]$$

第2步：第2行 × $\frac{1}{2}$，使 $x_2$ 系数变为 1：

$$\left[\begin{array}{rrrr} 1 & -2 & 1 & 0 \\ 0 & 1 & -4 & 4 \\ 0 & 10 & -10 & 10 \end{array}\right]$$

第3步：用第2行消去第3行的 $10x_2$。第3行 += (−10)×第2行，得 $30x_3 = -30$，再 × $\frac{1}{30}$：

$$\left[\begin{array}{rrrr} 1 & -2 & 1 & 0 \\ 0 & 1 & -4 & 4 \\ 0 & 0 & 1 & -1 \end{array}\right] \quad \text{（三角形式 triangular form）}$$

第4步（**回代 backward**，从下往上）：用 $x_3$ 消去上方方程中的 $x_3$ 项；再用 $x_2$ 消去第1方程中的 $-2x_2$。最终：

$$\left[\begin{array}{rrrr} 1 & 0 & 0 & 1 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & -1 \end{array}\right] \;\Rightarrow\; x_1 = 1,\quad x_2 = 0,\quad x_3 = -1$$

**验证**：$1(1) - 2(0) + 1(-1) = 0$；$2(0) - 8(-1) = 8$；$5(1) - 5(-1) = 10$

**几何意义**：每个方程在三维空间中确定一个平面，三个平面交于唯一点 $(1, 0, -1)$：

![三平面交于一点](img/fig03_planes_one_point.png)

### 存在性与唯一性问题 Existence and Uniqueness Questions

**TWO FUNDAMENTAL QUESTIONS ABOUT A LINEAR SYSTEM 关于线性系统的两个基本问题**：
1. 系统相容吗？即解**存在**吗？（Is the system consistent?）
2. 如果解存在，它是唯一的吗？（Is the solution unique?）

**例2（Example 2）**：例1 的系统已化为三角形式。由 $x_3 = -1$ 可唯一确定 $x_2$，再唯一确定 $x_1$——**解存在且唯一，系统相容**。

**例3（Example 3）判断系统是否相容**：

$$x_2 - \;4x_3 = 8 \qquad 2x_1 - 3x_2 + \;2x_3 = 1 \qquad 4x_1 - 8x_2 + 12x_3 = 1$$

增广矩阵行变换：交换1、2行 → 第3行 += (−2)×第1行 → 第3行 += 2×第2行：

$$\left[\begin{array}{rrrr} 2 & -3 & 2 & 1 \\ 0 & 1 & -4 & 8 \\ 0 & 0 & 0 & 15 \end{array}\right] \;\Rightarrow\; 0 = 15 \;\text{矛盾！}$$

最后一行 $0 = 15$ 是 $0x_1+0x_2+0x_3=15$ 的简写，永远不成立。**系统不相容（无解）**。

![三平面无公共点](img/fig04_planes_no_point.png)

> 记住：三角形式中若出现 $[\,0\;\cdots\;0\;\;b\,]$（$b\neq0$）这样的行，系统就不相容。

**数值注记（Numerical Note）**：现实中计算机用**浮点算术（floating point arithmetic）**，数字以 ±.d1...dp×10^r 存储，会产生舍入误差（roundoff error），但通常不会造成问题。

### 练习示例（1.1 Exercises 精选）

**判断 $(3,4,-2)$ 是否为下列系统的解**（Practice Problem 3）：

$$\begin{aligned} 5x_1 - \;x_2 + 2x_3 &= 7 \\ -2x_1 + 6x_2 + 9x_3 &= 0 \\ -7x_1 + 5x_2 - 3x_3 &= -7 \end{aligned}$$

代入 $x_1=3, x_2=4, x_3=-2$：
- $5(3)-(4)+2(-2) = 15-4-4 = 7$ 成立
- $-2(3)+6(4)+9(-2) = -6+24-18 = 0$ 成立
- $-7(3)+5(4)-3(-2) = -21+20+6 = 5 \neq -7$ 不成立

**不是解**。几何上：该点在前两个平面的交线上，但不在第三个平面上：

![点(3,4,-2)只满足前两个方程](img/fig05_planes_partial.png)

**$h, k$ 取何值时系统相容？**（Practice Problem 4）：
$$2x_1 - x_2 = h \qquad -6x_1 + 3x_2 = k$$
第2行 += 3×第1行，得 $0 = k + 3h$。**当且仅当 $k + 3h = 0$ 时系统相容**。

**练习题 1**：
$$x_1 + 5x_2 = 7 \qquad -2x_1 - 7x_2 = -5$$
第2行 += 2×第1行：$3x_2 = 9 \Rightarrow x_2 = 3$，回代 $x_1 = 7 - 15 = -8$。**解：$(-8, 3)$**

---

## 1.2 行化简与阶梯形 Row Reduction and Echelon Forms

**中文讲解**：上一节的消元法很"手工"，这一节把它升级成一套**标准化算法**，能分析任意线性方程组，并严格回答"解是否存在、是否唯一"两个问题。

### 阶梯形与简化阶梯形 Echelon Form

**DEFINITION 定义**——矩阵是**阶梯形（echelon form / 行阶梯形 row echelon form）**，当且仅当满足：

1. 所有非零行都在全零行之上
2. 每一行的**主元（leading entry）**所在列，严格位于上一行主元的右边
3. 主元下方的元素全为 0

若 additionally 满足下面两条，则为**简化阶梯形（reduced echelon form / RREF）**：

4. 每个非零行的主元是 1
5. 每个主元 1 是它所在列中唯一的非零元

**直觉**：主元呈"阶梯（staircase）"状向右下方排列。

- **阶梯矩阵 / 简化阶梯矩阵**：分别指处于这两种形态的矩阵
- **行化简（row reduced）**：用初等行变换把矩阵变形的过程

**THEOREM 1（Uniqueness of the Reduced Echelon Form 简化阶梯形的唯一性）**：
> Each matrix is row equivalent to **one and only one** reduced echelon matrix.
> 每个矩阵行等价于**唯一一个**简化阶梯矩阵。

### 主元位置 Pivot Positions

**DEFINITION 定义**：矩阵 $A$ 中的**主元位置（pivot position）**是对应其简化阶梯形中主元 1 的位置；包含主元位置的列叫**主元列（pivot column）**。

**例2（Example 2）**：把 $A$ 行化简到阶梯形并找出主元列：

$$A = \begin{bmatrix} 0 & -3 & -6 & 4 & 9 \\ -1 & -2 & -1 & 3 & 1 \\ -2 & -3 & 0 & 3 & -1 \\ 1 & 4 & 5 & -9 & -7 \end{bmatrix}$$

最左非零列是主元列，把非零主元换到顶部（交换1、4行）→ 主元下方消零 → 对子矩阵重复……最终：

$$\begin{bmatrix} 1 & 4 & 5 & -9 & -7 \\ 0 & 2 & 4 & -6 & -6 \\ 0 & 0 & 0 & -5 & 0 \\ 0 & 0 & 0 & 0 & 0 \end{bmatrix} \quad\text{一般形式：}\; \begin{bmatrix} \blacksquare & * & * & * & * \\ 0 & \blacksquare & * & * & * \\ 0 & 0 & 0 & \blacksquare & * \\ 0 & 0 & 0 & 0 & 0 \end{bmatrix}$$

$A$ 的主元列是**第 1、2、4 列**。

> 注意区分：**pivot（主元）** 是化简过程中用于消零的那个非零数（本例为 1, 2, −5），它不一定等于 $A$ 原位置上的元素。

### 行化简算法 The Row Reduction Algorithm（高斯消元）

**前向阶段（forward phase，步骤1–4）→ 阶梯形；后向阶段（backward phase，步骤5）→ 简化阶梯形**

- **STEP 1**：从最左非零列开始，这是主元列，主元位置在顶端
- **STEP 2**：在该列选一个非零元作主元，必要时换行使它位于主元位置
- **STEP 3**：用行替换把主元下方的位置全部化为 0
- **STEP 4**：盖住已处理的行，对剩余子矩阵重复步骤1–3，直到没有非零行可处理
- **STEP 5（后向阶段）**：从最右主元开始，从右往左、从下往上，把每个主元上方也化为 0；若主元不是 1，先用倍乘把它变成 1

**例3（Example 3）完整示范**：

$$\begin{bmatrix} 0 & 3 & -6 & 6 & 4 & -5 \\ 3 & -7 & 8 & -5 & 8 & 9 \\ 3 & -9 & 12 & -9 & 6 & 15 \end{bmatrix}$$

1. 最左非零列是第1列 → 交换1、3行使主元 3 到顶部
2. 第2行 += (−1)×第1行
3. 子矩阵中第2列是新主元列，主元 2；第3行 += (−3/2)×第2行：

$$\begin{bmatrix} 3 & -9 & 12 & -9 & 6 & 15 \\ 0 & 2 & -4 & 4 & 2 & -6 \\ 0 & 0 & 0 & 0 & 1 & 4 \end{bmatrix} \;\text{（阶梯形）}$$

4. 后向阶段：第1行 += (−6)×第3行，第2行 += (−2)×第3行；第2行 ×1/2；第1行 += 9×第2行；第1行 ×1/3：

$$\begin{bmatrix} 1 & 0 & -2 & 3 & 0 & -24 \\ 0 & 1 & -2 & 2 & 0 & -7 \\ 0 & 0 & 0 & 0 & 1 & 4 \end{bmatrix} \;\text{（简化阶梯形，唯一）}$$

> **数值注记**：计算机通常选列中绝对值最大的元素作主元，这叫**部分主元法（partial pivoting）**，可减小舍入误差。

### 线性系统的解 Solutions of Linear Systems

**基本变量与自由变量（核心概念！）**：
- **基本变量（basic variables）**：对应于主元列的变量
- **自由变量（free variables）**：其余变量，可以任意取值

例如增广矩阵化简为：

$$\left[\begin{array}{rrrr} 1 & 0 & -5 & 1 \\ 0 & 1 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{array}\right] \;\Rightarrow\; \begin{aligned} x_1 \quad\;\; - 5x_3 &= 1 \\ x_2 + \;\;x_3 &= 4 \\ 0 &= 0 \end{aligned}$$

$x_1, x_2$ 是基本变量，$x_3$ 是自由变量。把基本变量用自由变量表示：

$$\begin{cases} x_1 = 1 + 5x_3 \\ x_2 = 4 - x_3 \\ x_3 \text{ 自由（free）} \end{cases}$$

$x_3$ 每取一个值就得到一个解：$x_3=0 \Rightarrow (1,4,0)$；$x_3=1 \Rightarrow (6,3,1)$。**不同的 $x_3$ 选择给出不同的解，每个解都由 $x_3$ 的取值决定**——这就是"无穷多解"的具体样子。这种写法叫解集的**参数描述（parametric description）**。

**例4（Example 4）**：增广矩阵（已化到阶梯形）：

$$\left[\begin{array}{rrrrrr} 1 & 6 & 2 & -5 & -2 & -4 \\ 0 & 0 & 2 & -8 & -1 & 3 \\ 0 & 0 & 0 & 0 & 1 & 7 \end{array}\right]$$

继续化成简化阶梯形：

$$\sim \left[\begin{array}{rrrrrr} 1 & 6 & 0 & 3 & 0 & 0 \\ 0 & 0 & 1 & -4 & 0 & 5 \\ 0 & 0 & 0 & 0 & 1 & 7 \end{array}\right]$$

主元列是 1、3、5，基本变量 $x_1, x_3, x_5$；自由变量 $x_2, x_4$。通解：

$$\begin{cases} x_1 = -6x_2 - 3x_4 \\ x_2 \text{ 自由} \\ x_3 = 5 + 4x_4 \\ x_4 \text{ 自由} \\ x_5 = 7 \end{cases}$$

> 写通解时必须**每个变量都写出**，并明确标注哪些是自由变量。$\{x_1 = -3+8x_3,\; x_2 = -1+x_3,\; x_3 = 1+x_2\}$ 是**错误**的写法，因为它错误地暗示 $x_2$ 和 $x_3$ 都是自由变量。

**回代（Back-Substitution）**：计算机程序通常只化到（非简化）阶梯形，再用回代求解：从最下方的方程解出最后一个基本变量，逐层向上代入。手算时教材建议**统一使用简化阶梯形**，矩阵格式能大幅减少出错。

**数值注记**：化到阶梯形约需 $2n^3/3$ 次浮点运算（flops），而继续化到简化阶梯形最多只需 $n^2$ 次——前向阶段远比后向阶段耗时。

### 存在性与唯一性定理 THEOREM 2

**Existence and Uniqueness Theorem 存在唯一性定理**：
> 线性系统**相容**，当且仅当增广矩阵的最右列**不是**主元列——即其阶梯形中没有 $[\,0\;\cdots\;0\;\;b\,]$（$b\neq 0$）这样的行。
> 若系统相容，则解集要么 (i) **唯一**（无自由变量），要么 (ii) **无穷多**（至少有一个自由变量）。

**用行化简解线性系统的标准流程（USING ROW REDUCTION TO SOLVE A LINEAR SYSTEM）**：
1. 写出增广矩阵
2. 行化简到阶梯形，判断相容性；若无解，停止
3. 继续化到简化阶梯形
4. 写出对应的方程组
5. 把每个非零方程改写为：基本变量 = 自由变量的表达式

**几何图示**：两个平面相交于一条直线，通解就是这条交线的参数描述：

![通解=两平面交线](img/fig06_line_intersection.png)

*The general solution of the system of equations is the line of intersection of the two planes. 方程组的通解就是两平面的交线。*

### 练习示例（1.2 精选）

**练习 3**：行化简到简化阶梯形并圈出主元位置。

$$\begin{bmatrix} 1 & 2 & 3 & 4 \\ 4 & 5 & 6 & 7 \\ 6 & 7 & 8 & 9 \end{bmatrix} \Rightarrow \begin{bmatrix} 1 & 0 & -1 & -2 \\ 0 & 1 & 2 & 3 \\ 0 & 0 & 0 & 0 \end{bmatrix}$$

主元在第1、2列（原矩阵中对应位置圈出）。

**练习 19**：$x_1 + hx_2 = 2,\;\; 4x_1 + 8x_2 = k$
- (a) 无解：第2行是第1行的倍数但常数不符，如 $h = 4$ 且 $k \neq 8$
- (b) 唯一解：$h \neq 4$
- (c) 无穷多解：$h = 4$ 且 $k = 8$（第2行 = 4×第1行）

**真/假判断（T/F）精选**：
- "初等行变换不改变增广矩阵对应线性系统的解集" True
- "阶梯形是唯一的" False（简化阶梯形才唯一，见 Theorem 1）
- "有自由变量则必有无穷多解" 需先相容才成立；若不相容则解集为空，故为 False（需小心！）
- "化到阶梯形叫前向阶段（forward phase）" True

---

## 1.3 向量方程 Vector Equations

**中文讲解**：用**向量（vector）**的概念和记号来描述线性系统的性质，能把"方程组"与"几何"联系起来。本节开始，**向量 = 有序数列（ordered list of numbers）**。

### 二维空间中的向量 Vectors in $\mathbb{R}^2$

**英文原文**：
> A matrix with only one column is called a **column vector**, or simply a **vector**. The set of all vectors with two entries is denoted by $\mathbb{R}^2$ (read "r-two"). The $\mathbb{R}$ stands for the real numbers that appear as entries, and the exponent 2 indicates that each vector contains two entries.

- 只有一列的矩阵叫**列向量（column vector）**，简称**向量**：

$$\mathbf{u} = \begin{bmatrix} 3 \\ -1 \end{bmatrix}, \qquad \mathbf{v} = \begin{bmatrix} .2 \\ .3 \end{bmatrix}, \qquad \mathbf{w} = \begin{bmatrix} w_1 \\ w_2 \end{bmatrix}$$

- 所有"两个分量的向量"组成的集合记作 $\mathbb{R}^2$（r-two）
- **相等**：对应分量相等才算相等。$\begin{bmatrix} 4 \\ 7 \end{bmatrix} \neq \begin{bmatrix} 7 \\ 4 \end{bmatrix}$——向量是**有序**的！

（后续 1.3 节内容——向量的几何表示、线性组合、Span 等——将在下一讲继续。To be continued...）

---

## 本章小结 Chapter Summary

| 概念 Concept | 一句话理解 |
|---|---|
| 线性方程 linear equation | 变量的"一次式"：$a_1x_1+\cdots+a_nx_n=b$ |
| 增广矩阵 augmented matrix | 系数 + 常数项拼成的表格，行变换的舞台 |
| 初等行变换 row operations | 替换 / 交换 / 倍乘——三种不改变解集的变换 |
| 阶梯形 / 简化阶梯形 echelon / RREF | 主元呈阶梯排布；RREF 唯一 |
| 主元列 pivot column | 决定基本变量与自由变量的关键 |
| 自由变量 free variable | 可任意取值，若系统相容则有无穷多解 |
| 相容 consistent | 有解，当且仅当无 $[0\cdots0\;b]$（$b\neq0$）行 |

> **致初学者（A Note to Students）**：线性代数中**概念与计算同等重要**。数值练习帮你检验基本程序，但今后的职业生涯中，计算交给计算机，你要做的是**选择计算方法、解读结果、并向他人解释**。坚持每天学习——线性代数是一门"语言"，像学外语一样需要日积月累。
