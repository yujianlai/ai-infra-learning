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

$$
a_1x_1 + a_2x_2 + \cdots + a_nx_n = b
$$

其中 $b$ 和系数 $a_1, \dots, a_n$ 是已知实数或复数。直观判断标准：**变量只以一次幂出现，没有相乘、没有根号**。

是线性的（可整理成标准形式）：

$$
4x_1 - 5x_2 + 2 = x_1 \ \Rightarrow\  3x_1 - 5x_2 = -2,
$$

$$
x_2 = 2(\sqrt{6} - x_1) + x_3 \ \Rightarrow\  2x_1 + x_2 - x_3 = 2\sqrt{6}
$$

不是线性的：

$$
4x_1 - 5x_2 = x_1x_2 \quad (\text{出现 } x_1x_2 \text{ 乘积}),
$$

$$
x_2 = 2\sqrt{x_1} - 6 \quad (\text{出现根号})
$$

**英文原文**：
> A **linear equation** in the variables $x_1, \dots, x_n$ is an equation that can be written in the form $a_1x_1 + a_2x_2 + \cdots + a_nx_n = b$, where $b$ and the **coefficients** $a_1, \dots, a_n$ are real or complex numbers, usually known in advance.

### 什么是线性方程组？

**中文讲解**：**线性方程组**（也叫线性系统 linear system）就是**同一组变量**满足的一个或多个线性方程的集合。例如：

$$
2x_1 - x_2 + 1.5x_3 = 8
$$

$$
x_1 \quad\  - 4x_3 = -7
$$

- **解（solution）**：一组数 $(s_1, s_2, \dots, s_n)$，代入后使每个方程都成立。例如 $(5, 6.5, 3)$ 就是上面方程组的一个解。
- **解集（solution set）**：所有解的集合。
- **等价（equivalent）**：两个方程组解集完全相同，就称它们等价。

### 两个方程、两个变量的几何意义

求两个二元线性方程的解，就是找**两条直线的交点**：

$$
x_1 - 2x_2 = -1
$$

$$
-x_1 + 3x_2 = 3
$$

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

$$
x_1 - 2x_2 + \ x_3 = \ 0
$$

$$
2x_2 - 8x_3 = \ 8
$$

$$
5x_1 \quad\  - 5x_3 = 10
$$

系数矩阵与增广矩阵：

$$
\begin{bmatrix} 1 & -2 & 1 \\
0 & 2 & -8 \\
5 & 0 & -5 \end{bmatrix} \qquad \left[\begin{array}{rrrr} 1 & -2 & 1 & 0 \\
0 & 2 & -8 & 8 \\
5 & 0 & -5 & 10 \end{array}\right]
$$

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

$$
\begin{aligned} x_1 - 2x_2 + \ x_3 &= 0 \\
2x_2 - 8x_3 &= 8 \\
5x_1 \quad\  - 5x_3 &= 10 \end{aligned}
$$

**步骤（左边是方程，右边是对应增广矩阵，教材并排列出）**：

第1步：用第1方程消去其他方程的 $x_1$。第3行 += (−5)×第1行：

$$
\left[\begin{array}{rrrr} 1 & -2 & 1 & 0 \\
0 & 2 & -8 & 8 \\
0 & 10 & -10 & 10 \end{array}\right]
$$

第2步：第2行 × $\frac{1}{2}$，使 $x_2$ 系数变为 1：

$$
\left[\begin{array}{rrrr} 1 & -2 & 1 & 0 \\
0 & 1 & -4 & 4 \\
0 & 10 & -10 & 10 \end{array}\right]
$$

第3步：用第2行消去第3行的 $10x_2$。第3行 += (−10)×第2行，得 $30x_3 = -30$，再 × $\frac{1}{30}$：

$$
\left[\begin{array}{rrrr} 1 & -2 & 1 & 0 \\
0 & 1 & -4 & 4 \\
0 & 0 & 1 & -1 \end{array}\right] \quad \text{（三角形式 triangular form）}
$$

第4步（**回代 backward**，从下往上）：用 $x_3$ 消去上方方程中的 $x_3$ 项；再用 $x_2$ 消去第1方程中的 $-2x_2$。最终：

$$
\left[\begin{array}{rrrr} 1 & 0 & 0 & 1 \\
0 & 1 & 0 & 0 \\
0 & 0 & 1 & -1 \end{array}\right] \ \Rightarrow\  x_1 = 1,\quad x_2 = 0,\quad x_3 = -1
$$

**验证**：$1(1) - 2(0) + 1(-1) = 0$；$2(0) - 8(-1) = 8$；$5(1) - 5(-1) = 10$

**几何意义**：每个方程在三维空间中确定一个平面，三个平面交于唯一点 $(1, 0, -1)$：

![三平面交于一点](img/fig03_planes_one_point.png)

### 存在性与唯一性问题 Existence and Uniqueness Questions

**TWO FUNDAMENTAL QUESTIONS ABOUT A LINEAR SYSTEM 关于线性系统的两个基本问题**：
1. 系统相容吗？即解**存在**吗？（Is the system consistent?）
2. 如果解存在，它是唯一的吗？（Is the solution unique?）

**例2（Example 2）**：例1 的系统已化为三角形式。由 $x_3 = -1$ 可唯一确定 $x_2$，再唯一确定 $x_1$——**解存在且唯一，系统相容**。

**例3（Example 3）判断系统是否相容**：

$$
x_2 - \ 4x_3 = 8
$$

$$
2x_1 - 3x_2 + \ 2x_3 = 1
$$

$$
4x_1 - 8x_2 + 12x_3 = 1
$$

增广矩阵行变换：交换1、2行 → 第3行 += (−2)×第1行 → 第3行 += 2×第2行：

$$
\left[\begin{array}{rrrr} 2 & -3 & 2 & 1 \\
0 & 1 & -4 & 8 \\
0 & 0 & 0 & 15 \end{array}\right] \ \Rightarrow\  0 = 15 \ \text{矛盾！}
$$

最后一行 $0 = 15$ 是 $0x_1+0x_2+0x_3=15$ 的简写，永远不成立。**系统不相容（无解）**。

![三平面无公共点](img/fig04_planes_no_point.png)

> 记住：三角形式中若出现 $[0\ \cdots\ 0\ b]$（$b\neq0$）这样的行，系统就不相容。

**数值注记（Numerical Note）**：现实中计算机用**浮点算术（floating point arithmetic）**，数字以 ±.d1...dp×10^r 存储，会产生舍入误差（roundoff error），但通常不会造成问题。

### 练习示例（1.1 Exercises 精选）

**判断 $(3,4,-2)$ 是否为下列系统的解**（Practice Problem 3）：

$$
\begin{aligned} 5x_1 - \ x_2 + 2x_3 &= 7 \\
-2x_1 + 6x_2 + 9x_3 &= 0 \\
-7x_1 + 5x_2 - 3x_3 &= -7 \end{aligned}
$$

代入 $x_1=3, x_2=4, x_3=-2$：
- $5(3)-(4)+2(-2) = 15-4-4 = 7$ 成立
- $-2(3)+6(4)+9(-2) = -6+24-18 = 0$ 成立
- $-7(3)+5(4)-3(-2) = -21+20+6 = 5 \neq -7$ 不成立

**不是解**。几何上：该点在前两个平面的交线上，但不在第三个平面上：

![点(3,4,-2)只满足前两个方程](img/fig05_planes_partial.png)

**$h, k$ 取何值时系统相容？**（Practice Problem 4）：
$$
2x_1 - x_2 = h
$$

$$
-6x_1 + 3x_2 = k
$$
第2行 += 3×第1行，得 $0 = k + 3h$。**当且仅当 $k + 3h = 0$ 时系统相容**。

**练习题 1**：
$$
x_1 + 5x_2 = 7
$$

$$
-2x_1 - 7x_2 = -5
$$
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

$$
A = \begin{bmatrix} 0 & -3 & -6 & 4 & 9 \\
-1 & -2 & -1 & 3 & 1 \\
-2 & -3 & 0 & 3 & -1 \\
1 & 4 & 5 & -9 & -7 \end{bmatrix}
$$

最左非零列是主元列，把非零主元换到顶部（交换1、4行）→ 主元下方消零 → 对子矩阵重复……最终：

$$
\begin{bmatrix} 1 & 4 & 5 & -9 & -7 \\
0 & 2 & 4 & -6 & -6 \\
0 & 0 & 0 & -5 & 0 \\
0 & 0 & 0 & 0 & 0 \end{bmatrix} \quad\text{一般形式：}\  \begin{bmatrix} \blacksquare & * & * & * & * \\
0 & \blacksquare & * & * & * \\
0 & 0 & 0 & \blacksquare & * \\
0 & 0 & 0 & 0 & 0 \end{bmatrix}
$$

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

$$
\begin{bmatrix} 0 & 3 & -6 & 6 & 4 & -5 \\
3 & -7 & 8 & -5 & 8 & 9 \\
3 & -9 & 12 & -9 & 6 & 15 \end{bmatrix}
$$

1. 最左非零列是第1列 → 交换1、3行使主元 3 到顶部
2. 第2行 += (−1)×第1行
3. 子矩阵中第2列是新主元列，主元 2；第3行 += (−3/2)×第2行：

$$
\begin{bmatrix} 3 & -9 & 12 & -9 & 6 & 15 \\
0 & 2 & -4 & 4 & 2 & -6 \\
0 & 0 & 0 & 0 & 1 & 4 \end{bmatrix} \ \text{（阶梯形）}
$$

4. 后向阶段：第1行 += (−6)×第3行，第2行 += (−2)×第3行；第2行 ×1/2；第1行 += 9×第2行；第1行 ×1/3：

$$
\begin{bmatrix} 1 & 0 & -2 & 3 & 0 & -24 \\
0 & 1 & -2 & 2 & 0 & -7 \\
0 & 0 & 0 & 0 & 1 & 4 \end{bmatrix} \ \text{（简化阶梯形，唯一）}
$$

> **数值注记**：计算机通常选列中绝对值最大的元素作主元，这叫**部分主元法（partial pivoting）**，可减小舍入误差。

### 线性系统的解 Solutions of Linear Systems

**基本变量与自由变量（核心概念！）**：
- **基本变量（basic variables）**：对应于主元列的变量
- **自由变量（free variables）**：其余变量，可以任意取值

例如增广矩阵化简为：

$$
\left[\begin{array}{rrrr} 1 & 0 & -5 & 1 \\
0 & 1 & 1 & 4 \\
0 & 0 & 0 & 0 \end{array}\right] \ \Rightarrow\  \begin{aligned} x_1 \quad\  - 5x_3 &= 1 \\
x_2 + \ x_3 &= 4 \\
0 &= 0 \end{aligned}
$$

$x_1, x_2$ 是基本变量，$x_3$ 是自由变量。把基本变量用自由变量表示：

$$
\begin{cases} x_1 = 1 + 5x_3 \\
x_2 = 4 - x_3 \\
x_3 \text{ 自由（free）} \end{cases}
$$

$x_3$ 每取一个值就得到一个解：$x_3=0 \Rightarrow (1,4,0)$；$x_3=1 \Rightarrow (6,3,1)$。**不同的 $x_3$ 选择给出不同的解，每个解都由 $x_3$ 的取值决定**——这就是"无穷多解"的具体样子。这种写法叫解集的**参数描述（parametric description）**。

**例4（Example 4）**：增广矩阵（已化到阶梯形）：

$$
\left[\begin{array}{rrrrrr} 1 & 6 & 2 & -5 & -2 & -4 \\
0 & 0 & 2 & -8 & -1 & 3 \\
0 & 0 & 0 & 0 & 1 & 7 \end{array}\right]
$$

继续化成简化阶梯形：

$$
\sim \left[\begin{array}{rrrrrr} 1 & 6 & 0 & 3 & 0 & 0 \\
0 & 0 & 1 & -4 & 0 & 5 \\
0 & 0 & 0 & 0 & 1 & 7 \end{array}\right]
$$

主元列是 1、3、5，基本变量 $x_1, x_3, x_5$；自由变量 $x_2, x_4$。通解：

$$
\begin{cases} x_1 = -6x_2 - 3x_4 \\
x_2 \text{ 自由} \\
x_3 = 5 + 4x_4 \\
x_4 \text{ 自由} \\
x_5 = 7 \end{cases}
$$

> 写通解时必须**每个变量都写出**，并明确标注哪些是自由变量。$\{x_1 = -3+8x_3,\  x_2 = -1+x_3,\  x_3 = 1+x_2\}$ 是**错误**的写法，因为它错误地暗示 $x_2$ 和 $x_3$ 都是自由变量。

**回代（Back-Substitution）**：计算机程序通常只化到（非简化）阶梯形，再用回代求解：从最下方的方程解出最后一个基本变量，逐层向上代入。手算时教材建议**统一使用简化阶梯形**，矩阵格式能大幅减少出错。

**数值注记**：化到阶梯形约需 $2n^3/3$ 次浮点运算（flops），而继续化到简化阶梯形最多只需 $n^2$ 次——前向阶段远比后向阶段耗时。

### 存在性与唯一性定理 THEOREM 2

**Existence and Uniqueness Theorem 存在唯一性定理**：
> 线性系统**相容**，当且仅当增广矩阵的最右列**不是**主元列——即其阶梯形中没有 $[0\ \cdots\ 0\ b]$（$b\neq 0$）这样的行。
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

$$
\begin{bmatrix} 1 & 2 & 3 & 4 \\
4 & 5 & 6 & 7 \\
6 & 7 & 8 & 9 \end{bmatrix} \Rightarrow \begin{bmatrix} 1 & 0 & -1 & -2 \\
0 & 1 & 2 & 3 \\
0 & 0 & 0 & 0 \end{bmatrix}
$$

主元在第1、2列（原矩阵中对应位置圈出）。

**练习 19**：$x_1 + hx_2 = 2,\  4x_1 + 8x_2 = k$
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

$$
\mathbf{u} = \begin{bmatrix} 3 \\
-1 \end{bmatrix}, \qquad \mathbf{v} = \begin{bmatrix} .2 \\
.3 \end{bmatrix}, \qquad \mathbf{w} = \begin{bmatrix} w_1 \\
w_2 \end{bmatrix}
$$

- 所有"两个分量的向量"组成的集合记作 $\mathbb{R}^2$（r-two）
- **相等**：对应分量相等才算相等。$\begin{bmatrix} 4 \\\\ 7 \end{bmatrix} \neq \begin{bmatrix} 7 \\\\ 4 \end{bmatrix}$——向量是**有序**的！

### 向量的几何表示：点还是箭头？

**中文讲解**：向量有两种看图方式，都常用：

1. **看作点**：$\begin{bmatrix} -2 \\\\ 1 \end{bmatrix}$ 就是平面上坐标为 $(-2, 1)$ 的点
2. **看作箭头**：从原点 $(0,0)$ 指向点 $(-2,1)$ 的有向线段

![FIGURE 1 向量即点](img/fig13_1_vectors_as_points.png)

*FIGURE 1 Vectors as points. 图1 把向量看作平面上的点*

![FIGURE 2 带箭头的向量](img/fig13_2_vectors_with_arrows.png)

*FIGURE 2 Vectors with arrows. 图2 从原点出发的箭头表示向量*

> 英文原文：Individual numbers in a vector are called **entries** or **components**. Vectors in $\mathbb{R}^2$ can be visualized either as **points** or as **arrows** (directed line segments from the origin).

### 向量的加法与数乘（几何意义）

**中文讲解**：两个向量相加，用**平行四边形法则（parallelogram rule）**：把 $\mathbf{u}$ 和 $\mathbf{v}$ 的起点都放在原点，以它们为邻边作平行四边形，从原点到对角顶点的对角线就是 $\mathbf{u}+\mathbf{v}$。等价地，也可以把 $\mathbf{v}$ 的"尾巴"接到 $\mathbf{u}$ 的"箭头"上（三角形法则）。

![FIGURE 3 平行四边形法则](img/fig13_3_parallelogram_rule.png)

*FIGURE 3 The parallelogram rule. 图3 向量加法的平行四边形法则*

**例1（Example 1）**：

$$
\mathbf{u} = \begin{bmatrix} 1 \\
2 \end{bmatrix}, \quad \mathbf{v} = \begin{bmatrix} 2 \\
-2 \end{bmatrix} \ \Rightarrow\  \mathbf{u} + \mathbf{v} = \begin{bmatrix} 1+2 \\
2+(-2) \end{bmatrix} = \begin{bmatrix} 3 \\
0 \end{bmatrix}
$$

![FIGURE 4 例1的向量加法](img/fig13_4_example2.png)

**数乘（Scalar multiples）**：$c\mathbf{u}$ 就是把 $\mathbf{u}$ 的长度放大 $|c|$ 倍；$c>0$ 方向不变，$c<0$ 方向反转。

![FIGURE 5 向量的数乘](img/fig13_5_scalar_multiples.png)

*FIGURE 5 Scalar multiples. 图5 数乘的几何效果：拉伸与反向*

![FIGURE 6 三维空间中的数乘](img/fig13_6_scalar_multiples_R3.png)

*FIGURE 6 Scalar multiples in $\mathbb{R}^3$. 图6 三维空间中同样成立*

### n 维向量与代数性质

**中文讲解**：完全类似，所有 $n$ 个分量的向量组成 $\mathbb{R}^n$。两个向量**相等**当且仅当对应分量全相等；向量与**标量（scalar）**相乘就是每个分量都乘这个数。

**代数性质 THEOREM 1（Algebraic Properties of $\mathbb{R}^n$）**：对所有 $\mathbf{u},\mathbf{v},\mathbf{w} \in \mathbb{R}^n$ 和标量 $c,d$：

1. $\mathbf{u}+\mathbf{v} = \mathbf{v}+\mathbf{u}$（加法交换律）
2. $(\mathbf{u}+\mathbf{v})+\mathbf{w} = \mathbf{u}+(\mathbf{v}+\mathbf{w})$（加法结合律）
3. $\mathbf{u}+\mathbf{0} = \mathbf{0}+\mathbf{u} = \mathbf{u}$（零向量）
4. $\mathbf{u}+(-\mathbf{u}) = \mathbf{0}$（加法逆元）
5. $c(\mathbf{u}+\mathbf{v}) = c\mathbf{u}+c\mathbf{v}$（分配律）
6. $(c+d)\mathbf{u} = c\mathbf{u}+d\mathbf{u}$（分配律）
7. $c(d\mathbf{u}) = (cd)\mathbf{u}$（数乘结合律）
8. $1\mathbf{u} = \mathbf{u}$（单位元）

**减法**：$\mathbf{u}-\mathbf{v} = \mathbf{u}+(-\mathbf{v})$，几何上就是从 $\mathbf{v}$ 的终点指向 $\mathbf{u}$ 终点的箭头：

![FIGURE 7 向量减法](img/fig13_7_vector_subtraction.png)

*FIGURE 7 Vector subtraction. 图7 向量减法 $\mathbf{u}-\mathbf{v}$ 的几何意义*

### 线性组合 Linear Combinations（本节核心概念！）

**中文讲解**：给定向量 $\mathbf{v}_1, \dots, \mathbf{v}_p \in \mathbb{R}^n$ 和标量 $c_1, \dots, c_p$，向量

$$
\mathbf{y} = c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_p\mathbf{v}_p
$$

称为 $\mathbf{v}_1, \dots, \mathbf{v}_p$ 的**线性组合（linear combination）**，其中 $c_1, \dots, c_p$ 叫**权（weights）**。

> 英文原文：Given vectors $\mathbf{v}_1, \dots, \mathbf{v}_p$ in $\mathbb{R}^n$ and given scalars $c_1, \dots, c_p$, the vector $\mathbf{y} = c_1\mathbf{v}_1 + \cdots + c_p\mathbf{v}_p$ is called a **linear combination** of $\mathbf{v}_1, \dots, \mathbf{v}_p$ with **weights** $c_1, \dots, c_p$. Note that weights in a linear combination can be any real numbers, including zero.

**例2（Example 2）**：设 $\mathbf{v}_1 = \begin{bmatrix} -1 \\\\ 1 \end{bmatrix}$, $\mathbf{v}_2 = \begin{bmatrix} 2 \\\\ 1 \end{bmatrix}$，求 $\frac{1}{2}\mathbf{v}_1 + \mathbf{v}_2$、$0\mathbf{v}_1$、$1.5\mathbf{v}_1 - 2\mathbf{v}_2$。

$$
\frac{1}{2}\mathbf{v}_1 + \mathbf{v}_2 = \begin{bmatrix} -1/2+2 \\
1/2+1 \end{bmatrix} = \begin{bmatrix} 3/2 \\
3/2 \end{bmatrix}, \qquad 0\mathbf{v}_1 = \begin{bmatrix} 0 \\
0 \end{bmatrix}, \qquad 1.5\mathbf{v}_1 - 2\mathbf{v}_2 = \begin{bmatrix} -5.5 \\
-0.5 \end{bmatrix}
$$

**例3（Example 3）图释**：两个不平行的向量 $\mathbf{v}_1, \mathbf{v}_2$ 的所有整数线性组合，铺满了整个平面的网格——每个格点都是一个线性组合：

![FIGURE 8 线性组合网格](img/fig13_8_linear_combinations_grid.png)

*FIGURE 8 Linear combinations of $\mathbf{v}_1$ and $\mathbf{v}_2$. 图8 $\mathbf{v}_1$ 与 $\mathbf{v}_2$ 的线性组合铺满整个平面网格*

### 向量方程与线性方程组的联系（关键桥梁）

**中文讲解**：一个基本而重要的问题：

$$
\mathbf{b} \text{ 能否写成 } \mathbf{a}_1, \dots, \mathbf{a}_n \text{ 的线性组合？}
$$

即是否存在 $x_1, \dots, x_n$ 使

$$
x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n = \mathbf{b}
$$

这个**向量方程（vector equation）**有解，当且仅当对应的**线性方程组**有解——把每个 $\mathbf{a}_j$ 的分量按列排开，$x_j$ 就是方程组的未知量！

> 英文原文：A vector equation $x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n = \mathbf{b}$ has the same solution set as the linear system whose augmented matrix is $[\mathbf{a}_1\ \mathbf{a}_2\ \cdots\ \mathbf{a}_n\ \mathbf{b}]$. In particular, $\mathbf{b}$ can be generated by a linear combination of $\mathbf{a}_1, \dots, \mathbf{a}_n$ if and only if there exists a solution to the linear system corresponding to the matrix.

**例4（Example 4）**：设 $\mathbf{a}_1 = \begin{bmatrix} 1 \\\\ -2 \\\\ -5 \end{bmatrix}$, $\mathbf{a}_2 = \begin{bmatrix} 2 \\\\ 5 \\\\ 6 \end{bmatrix}$, $\mathbf{b} = \begin{bmatrix} 7 \\\\ 4 \\\\ -3 \end{bmatrix}$。判断 $\mathbf{b}$ 能否写成 $\mathbf{a}_1, \mathbf{a}_2$ 的线性组合，即解：

$$
x_1\begin{bmatrix} 1 \\
-2 \\
-5 \end{bmatrix} + x_2\begin{bmatrix} 2 \\
5 \\
6 \end{bmatrix} = \begin{bmatrix} 7 \\
4 \\
-3 \end{bmatrix}
$$

对应增广矩阵行化简：

$$
\left[\begin{array}{rrr} 1 & 2 & 7 \\
-2 & 5 & 4 \\
-5 & 6 & -3 \end{array}\right] \sim \left[\begin{array}{rrr} 1 & 2 & 7 \\
0 & 9 & 18 \\
0 & 16 & 32 \end{array}\right] \sim \left[\begin{array}{rrr} 1 & 2 & 7 \\
0 & 1 & 2 \\
0 & 0 & 0 \end{array}\right]
$$

得 $x_2 = 2$, $x_1 = 3$。**$\mathbf{b} = 3\mathbf{a}_1 + 2\mathbf{a}_2$**，向量方程有唯一解。

### 张成 Span（生成集）

**中文讲解**：$\mathbf{v}_1, \dots, \mathbf{v}_p$ 的**所有线性组合组成的集合**记作 $\operatorname{Span}\{\mathbf{v}_1, \dots, \mathbf{v}_p\}$，读作"$\mathbf{v}_1$ 到 $\mathbf{v}_p$ 张成的集合"。

> 英文原文：If $\mathbf{v}_1, \dots, \mathbf{v}_p$ are in $\mathbb{R}^n$, then the set of all linear combinations of $\mathbf{v}_1, \dots, \mathbf{v}_p$ is denoted by $\operatorname{Span}\{\mathbf{v}_1, \dots, \mathbf{v}_p\}$ and is called the **subset of $\mathbb{R}^n$ spanned** (or generated) by $\mathbf{v}_1, \dots, \mathbf{v}_p$. That is, $\operatorname{Span}\{\mathbf{v}_1, \dots, \mathbf{v}_p\}$ is the collection of all vectors that can be written in the form $c_1\mathbf{v}_1 + \cdots + c_p\mathbf{v}_p$ with $c_1, \dots, c_p$ scalars.

**几何直觉**：
- $\operatorname{Span}\{\mathbf{v}\}$：一条过原点的直线（$\mathbf{v}$ 的所有倍数）
- $\operatorname{Span}\{\mathbf{u}, \mathbf{v}\}$（$\mathbf{u},\mathbf{v}$ 不平行）：过原点的整张平面
- 问"$\mathbf{b}$ 是否在 $\operatorname{Span}\{\mathbf{v}_1, \dots, \mathbf{v}_p\}$ 中" = 问"$\mathbf{b}$ 能否由它们线性表示" = 问"对应的方程组是否有解"

**例5（Example 5）**：设 $\mathbf{v}_1 = \begin{bmatrix} 1 \\\\ -2 \\\\ 3 \end{bmatrix}$, $\mathbf{v}_2 = \begin{bmatrix} 5 \\\\ -13 \\\\ -3 \end{bmatrix}$, $\mathbf{b} = \begin{bmatrix} -3 \\\\ 8 \\\\ 1 \end{bmatrix}$。判断 $\mathbf{b}$ 是否在 $\operatorname{Span}\{\mathbf{v}_1, \mathbf{v}_2\}$ 中：

$$
\left[\begin{array}{rrr} 1 & 5 & -3 \\
-2 & -13 & 8 \\
3 & -3 & 1 \end{array}\right] \sim \left[\begin{array}{rrr} 1 & 5 & -3 \\
0 & -3 & 2 \\
0 & -18 & 10 \end{array}\right] \sim \left[\begin{array}{rrr} 1 & 5 & -3 \\
0 & -3 & 2 \\
0 & 0 & -2 \end{array}\right]
$$

第三行 $0 = -2$ 矛盾，**无解**——$\mathbf{b}$ 不在 $\operatorname{Span}\{\mathbf{v}_1, \mathbf{v}_2\}$ 中。

**例6（Example 6，几何判断）**：设 $\mathbf{u}, \mathbf{v}$ 是 $\mathbb{R}^3$ 中向量，判断 $\operatorname{Span}\{\mathbf{u}, \mathbf{v}\}$ 是什么样的集合。

若 $\mathbf{v}$ 是 $\mathbf{u}$ 的倍数，则 $\operatorname{Span}\{\mathbf{u},\mathbf{v}\}$ 是过原点的直线；否则是过原点的平面——即 $\mathbf{u}$ 与 $\mathbf{v}$ 确定的平面上的所有向量。

![FIGURE 9 例5的平行四边形](img/fig13_9_parallelogram_w.png)

*FIGURE 9 $\mathbf{w}$ 不是 $\mathbf{u}$ 与 $\mathbf{v}$ 的线性组合。 图9 不在张成集中的向量*

### 1.3 小结

| 概念 | 一句话理解 |
|---|---|
| 向量 vector | 有序数列，几何上即点或从原点出发的箭头 |
| 线性组合 linear combination | $c_1\mathbf{v}_1 + \cdots + c_p\mathbf{v}_p$，系数 $c_i$ 叫权 weights |
| 张成 Span | 所有线性组合的集合：一条直线、一个平面…… |
| 向量方程 = 方程组 | $\mathbf{b}$ 能否线性表示 ⇔ 对应方程组是否有解 |

> **本节精髓**：线性代数把两个看似无关的问题统一了——"一个向量能否由另一些向量线性生成"（几何/向量语言）与"一个线性方程组是否有解"（代数语言），本质上是**同一个问题**。

---

## 1.4 矩阵方程 Ax = b The Matrix Equation Ax = b

**中文讲解**：线性代数的一个 fundamental idea（根本思想）：把向量的线性组合看成**矩阵与向量的乘积**。这个定义让我们能用全新方式重述 1.3 节的概念。

### 定义：矩阵乘向量 DEFINITION

若 $A$ 是 $m \times n$ 矩阵，列为 $\mathbf{a}_1, \dots, \mathbf{a}_n$，$\mathbf{x}$ 在 $\mathbb{R}^n$ 中，则 $A$ 与 $\mathbf{x}$ 的**乘积（product）**记作 $A\mathbf{x}$，是 $A$ 的各列以 $\mathbf{x}$ 中对应元素为权的**线性组合**：

$$
A\mathbf{x} = \begin{bmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \end{bmatrix} \begin{bmatrix} x_1 \\
x_2 \\
\vdots \\
x_n \end{bmatrix} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n
$$

> 注意：$A\mathbf{x}$ 只有当 $A$ 的**列数**等于 $\mathbf{x}$ 的**分量个数**时才有定义（即 $\mathbf{x} \in \mathbb{R}^n$，$A$ 是 $m\times n$）。

**例1（Example 1）**：

$$
\text{a. } \begin{bmatrix} 1 & 2 & -1 \\
0 & -5 & 3 \end{bmatrix} \begin{bmatrix} 4 \\
3 \\
7 \end{bmatrix} = 4\begin{bmatrix} 1 \\
0 \end{bmatrix} + 3\begin{bmatrix} 2 \\
-5 \end{bmatrix} + 7\begin{bmatrix} -1 \\
3 \end{bmatrix} = \begin{bmatrix} 4 \\
0 \end{bmatrix} + \begin{bmatrix} 6 \\
-15 \end{bmatrix} + \begin{bmatrix} -7 \\
21 \end{bmatrix} = \begin{bmatrix} 3 \\
6 \end{bmatrix}
$$

**例2（Example 2）**：把线性组合 $3\mathbf{v}_1 - 5\mathbf{v}_2 + 7\mathbf{v}_3$ 写成矩阵乘向量：设 $A = [\mathbf{v}_1\ \mathbf{v}_2\ \mathbf{v}_3]$，$\mathbf{x} = (3, -5, 7)$，则 $A\mathbf{x} = 3\mathbf{v}_1 - 5\mathbf{v}_2 + 7\mathbf{v}_3$。

### THEOREM 3（矩阵方程、向量方程、方程组三合一）

**中文讲解**：下面三个问题有**完全相同的解集**——这是全书的"交通枢纽"：
1. 矩阵方程 $A\mathbf{x} = \mathbf{b}$
2. 向量方程 $x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n = \mathbf{b}$
3. 以 $[\mathbf{a}_1\ \mathbf{a}_2\ \cdots\ \mathbf{a}_n\ \mathbf{b}]$ 为增广矩阵的线性方程组

> **Existence of Solutions 解的存在性**：
> 方程 $A\mathbf{x} = \mathbf{b}$ 有解，当且仅当 $\mathbf{b}$ 是 $A$ 各列的线性组合。

**例3（Example 3）**：设 $A = \begin{bmatrix} 1 & 5 & -3 \\\\ -4 & 1 & 3 \\\\ -2 & -7 & 0 \end{bmatrix}$, $\mathbf{b} = \begin{bmatrix} -4 \\\\ 3 \\\\ h \end{bmatrix}$。$h$ 取何值时方程 $A\mathbf{x}=\mathbf{b}$ 有解？

增广矩阵行化简：

$$
\left[\begin{array}{rrrr} 1 & 5 & -3 & -4 \\
-4 & 1 & 3 & 3 \\
-2 & -7 & 0 & h \end{array}\right] \sim \left[\begin{array}{rrrr} 1 & 5 & -3 & -4 \\
0 & 1 & -2 & -1 \\
0 & 3 & -6 & h-8 \end{array}\right] \sim \left[\begin{array}{rrrr} 1 & 5 & -3 & -4 \\
0 & 1 & -2 & -1 \\
0 & 0 & 0 & h-5 \end{array}\right]
$$

当且仅当第4列没有主元时系统相容，即 **$h = 5$**。此时 $\mathbf{b} \in \operatorname{Span}\{\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3\}$。

> **Remember 切记**：系统中**有自由变量并不能保证系统相容**！自由变量只说明"如果"有解则有无穷多解。

![A的列张成一个过原点的平面](img/fig14_1_columns_span_plane.png)

*FIGURE 1 The columns of $A$ span a plane through $\mathbf{0}$. 图1 $A$ 的各列张成一个过原点的平面——只有当 $\mathbf{b}$ 在这个平面上时，$A\mathbf{x}=\mathbf{b}$ 才有解。*

### THEOREM 4（$A\mathbf{x}=\mathbf{b}$ 有解的等价条件）

设 $A$ 是 $m \times n$ 矩阵，则下列命题**逻辑等价**（同真同假）：
1. 对每个 $\mathbf{b} \in \mathbb{R}^m$，$A\mathbf{x} = \mathbf{b}$ 都有解
2. $\mathbb{R}^m$ 中每个 $\mathbf{b}$ 都是 $A$ 各列的线性组合
3. $A$ 的各列**张成** $\mathbb{R}^m$（$\operatorname{Span}\{\mathbf{a}_1,\dots,\mathbf{a}_n\} = \mathbb{R}^m$）
4. $A$ **每一行都有主元位置**

> **直观理解**：要"覆盖"整个 $\mathbb{R}^m$，$A$ 必须有足够多的"方向"——行数 $m$ 不能太多（否则某些 $\mathbf{b}$ 够不着）。特别地，若 $n < m$（列数少于行数），命题一定为假。

### 计算 Ax 的另一种方法：行向量法则 Row-Vector Rule

若 $A$ 的各**行**为 $\mathbf{r}_1, \dots, \mathbf{r}_m$，则：

$$
A\mathbf{x} = \begin{bmatrix} \mathbf{r}_1 \\
\vdots \\
\mathbf{r}_m \end{bmatrix} \mathbf{x} = \begin{bmatrix} \mathbf{r}_1 \cdot \mathbf{x} \\
\vdots \\
\mathbf{r}_m \cdot \mathbf{x} \end{bmatrix} = \begin{bmatrix} a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n \\
\vdots \\
a_{m1}x_1 + a_{m2}x_2 + \cdots + a_{mn}x_n \end{bmatrix}
$$

即 **$A\mathbf{x}$ 的每个分量 = $A$ 的对应行与 $\mathbf{x}$ 的点积（dot product）**。

**例4（Example 4）**：

$$
\begin{bmatrix} 1 & 2 & -1 \\
0 & -5 & 3 \end{bmatrix} \begin{bmatrix} 4 \\
3 \\
7 \end{bmatrix} = \begin{bmatrix} 1\cdot4 + 2\cdot3 + (-1)\cdot7 \\
0\cdot4 + (-5)\cdot3 + 3\cdot7 \end{bmatrix} = \begin{bmatrix} 3 \\
6 \end{bmatrix}
$$

两种算法（列的线性组合 / 行向量法则）结果一致，互为验证。

### 单位矩阵 Identity Matrix

主对角线（从左上到右下）全为 1、其余全为 0 的 $n\times n$ 矩阵叫**单位矩阵（identity matrix）**，记作 $I_n$。例如：

$$
I_3 = \begin{bmatrix} 1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1 \end{bmatrix}
$$

**关键性质**：对任何 $\mathbf{x} \in \mathbb{R}^3$，有 $I_3\mathbf{x} = \mathbf{x}$。单位矩阵是矩阵乘法中的"1"。

### THEOREM 5（矩阵乘法的代数性质）

设 $A$ 是 $m\times n$ 矩阵，$\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$，$c$ 是标量：
1. $A(\mathbf{u}+\mathbf{v}) = A\mathbf{u} + A\mathbf{v}$（分配律）
2. $A(c\mathbf{u}) = c(A\mathbf{u})$（数乘交换）

> **证明思路（小白版）**：把 $\mathbf{u}+\mathbf{v}$ 写成分量形式代入定义即可。例如第1条：$A(\mathbf{u}+\mathbf{v})$ 的第 $i$ 个分量是 $\sum_j a_{ij}(u_j+v_j) = \sum_j a_{ij}u_j + \sum_j a_{ij}v_j$，正是 $A\mathbf{u}+A\mathbf{v}$ 的第 $i$ 个分量。

### 1.4 小结

| 概念 | 一句话理解 |
|---|---|
| $A\mathbf{x}$ | $A$ 各列以 $\mathbf{x}$ 分量为权的线性组合；也等于"行·点积" |
| Theorem 3 | $A\mathbf{x}=\mathbf{b}$ = 向量方程 = 方程组，**三者同解** |
| Theorem 4 | 对所有 $\mathbf{b}$ 有解 ⇔ $A$ 的列张成 $\mathbb{R}^m$ ⇔ 每行有主元 |
| Theorem 5 | 矩阵乘法保持加法和数乘（线性性的雏形） |

> **本节精髓**：$A\mathbf{x}=\mathbf{b}$ 把"方程组求解"（1.1–1.2）、"线性组合/Span"（1.3）统一成一个简洁的矩阵语言。下一节 1.8–1.9 的**线性变换（linear transformation）**将揭示 $A\mathbf{x}$ 的几何本质：它是对空间的**旋转、拉伸、投影**。
## 1.5 线性方程组的解集 Solution Sets of Linear Systems

**中文讲解**：1.2 节我们学会了判断解的"个数"，这一节把解集**完整地、几何地**写出来。核心工具：把通解改写成**参数向量形式（parametric vector form）**。

### 齐次线性方程组 Homogeneous Linear Equations

**中文讲解**：形如 $A\mathbf{x} = \mathbf{0}$ 的方程组叫**齐次方程组**。它一定相容——因为 $\mathbf{x} = \mathbf{0}$ 永远是解（叫**平凡解 trivial solution**）。关键问题是：它有没有**非平凡解（nontrivial solution）**（即某个分量非零的解）？

> **基本事实**：齐次方程 $A\mathbf{x} = \mathbf{0}$ 有非平凡解，**当且仅当**方程至少有一个自由变量。（由 1.2 节定理 2 直接推出）

**例1（Example 1）**：判断下列齐次系统是否有非平凡解，并描述其解集：

$$
3x_1 + 5x_2 - 4x_3 = 0
$$

$$
-3x_1 - 2x_2 + 4x_3 = 0
$$

$$
6x_1 + \ x_2 - 8x_3 = 0
$$

增广矩阵（右端全 0，可省略）行化简到阶梯形再化到简化阶梯形：

$$
\begin{bmatrix} 3 & 5 & -4 & 0 \\
-3 & -2 & 4 & 0 \\
6 & 1 & -8 & 0 \end{bmatrix} \sim \begin{bmatrix} 3 & 5 & -4 & 0 \\
0 & 3 & 0 & 0 \\
0 & 0 & 0 & 0 \end{bmatrix} \sim \begin{bmatrix} 1 & 0 & -\frac{4}{3} & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 \end{bmatrix}
$$

$x_3$ 自由，$x_1 = \frac{4}{3}x_3$, $x_2 = 0$。通解写成向量形式并把 $x_3$ 提出来：

$$
\mathbf{x} = \begin{bmatrix} x_1 \\
x_2 \\
x_3 \end{bmatrix} = \begin{bmatrix} \frac{4}{3}x_3 \\
0 \\
x_3 \end{bmatrix} = x_3 \begin{bmatrix} \frac{4}{3} \\
0 \\
1 \end{bmatrix} = x_3\mathbf{v}, \quad \mathbf{v} = \begin{bmatrix} 4/3 \\
0 \\
1 \end{bmatrix}
$$

**每个解都是 $\mathbf{v}$ 的倍数，解集 = $\operatorname{Span}\{\mathbf{v}\}$**——几何上是 $\mathbb{R}^3$ 中过原点的一条直线（除 $\mathbf{x}=\mathbf{0}$ 外的点全是非平凡解）：

![FIGURE 1 齐次解集是一条过原点的直线](img/fig15_2_span_v_line.png)

*FIGURE 1 The solution set of $A\mathbf{x}=\mathbf{0}$ is a line through $\mathbf{0}$. 图1 $A\mathbf{x}=\mathbf{0}$ 的解集是过原点的直线*

![FIGURE 2 齐次解集是过原点的平面](img/fig15_1_homogeneous_plane.png)

*FIGURE 2 The solution set of $A\mathbf{x}=\mathbf{0}$ is a plane through $\mathbf{0}$（两个自由变量时）. 图2 若有两个自由变量，解集是过原点的平面（$\operatorname{Span}\{\mathbf{u},\mathbf{v}\}$）*

> **注意**：非平凡解 $\mathbf{x}$ 允许某些分量是 0，只要**不全为 0** 即可。

### 参数向量形式 Parametric Vector Form

**中文讲解**：上例中"把自由变量作为公因子提出来"的技巧可以推广。若齐次系统通解为

$$
\mathbf{x} = s\mathbf{u} + t\mathbf{v}
$$

（$s,t$ 自由），解集就是 $\operatorname{Span}\{\mathbf{u}, \mathbf{v}\}$——过原点的平面。这种把解写成自由变量的向量和的形式，叫解集的**参数向量形式**。

### 非齐次方程组的解集：平移！THEOREM 6

**中文讲解**：对非齐次方程 $A\mathbf{x} = \mathbf{b}$（$\mathbf{b} \neq \mathbf{0}$），解集不再是过原点的子空间，而是**把 $A\mathbf{x}=\mathbf{0}$ 的解集沿向量 $\mathbf{p}$ 平移**，其中 $\mathbf{p}$ 是 $A\mathbf{x}=\mathbf{b}$ 的任一特解。

> **THEOREM 6**：设 $A\mathbf{x}=\mathbf{b}$ 相容，$\mathbf{p}$ 是一个特解。则 $A\mathbf{x}=\mathbf{b}$ 的解集 = 满足 $\mathbf{w} = \mathbf{p} + \mathbf{v}_h$ 的所有向量 $\mathbf{w}$ 的集合，其中 $\mathbf{v}_h$ 是 $A\mathbf{x}=\mathbf{0}$ 的任意解。
> 换句话说：**非齐次通解 = 特解 + 齐次通解**（$\mathbf{x} = \mathbf{p} + \mathbf{v}_h$）

![FIGURE 3 平移](img/fig15_3_translation.png)

*FIGURE 3 Adding $\mathbf{p}$ to $\mathbf{v}$ translates $\mathbf{v}$ to $\mathbf{v}+\mathbf{p}$. 图3 把每个向量都加上 $\mathbf{p}$，等于把整条线平移*

**例2（Example 2）**：描述 $A\mathbf{x}=\mathbf{b}$ 的解集，其中 $A$ 为例1 的系数矩阵，$\mathbf{b} = \begin{bmatrix} 7 \\\\ -4 \\\\ -4 \end{bmatrix}$。

增广矩阵行化简：

$$
\left[\begin{array}{rrrr} 3 & 5 & -4 & 7 \\
-3 & -2 & 4 & -4 \\
6 & 1 & -8 & -4 \end{array}\right] \sim \left[\begin{array}{rrrr} 1 & 0 & -\frac{4}{3} & -1 \\
0 & 1 & 0 & 2 \\
0 & 0 & 0 & 0 \end{array}\right]
$$

$x_3$ 自由，通解：

$$
\mathbf{x} = \begin{bmatrix} -1 + \frac{4}{3}x_3 \\
2 \\
x_3 \end{bmatrix} = \begin{bmatrix} -1 \\
2 \\
0 \end{bmatrix} + x_3 \begin{bmatrix} 4/3 \\
0 \\
1 \end{bmatrix} = \mathbf{p} + x_3\mathbf{v}
$$

即 **$\mathbf{x} = \mathbf{p} + t\mathbf{v}$（$t \in \mathbb{R}$）**——解集是过 $\mathbf{p}$、平行于 $\mathbf{v}$ 的直线（不过原点！）：

![FIGURE 4 平移后的直线](img/fig15_4_translated_line.png)

*FIGURE 4 The translated line. 图4 非齐次解集 = 齐次解集平移后的直线*

**例3（Example 3）**：单个方程 $10x_1 - 3x_2 - 2x_3 = 0$ 的解集。把 $x_1$ 用自由变量表示：$x_1 = .3x_2 + .2x_3$，

$$
\mathbf{x} = \begin{bmatrix} .3x_2 + .2x_3 \\
x_2 \\
x_3 \end{bmatrix} = x_2\begin{bmatrix} .3 \\
1 \\
0 \end{bmatrix} + x_3\begin{bmatrix} .2 \\
0 \\
1 \end{bmatrix} = x_2\mathbf{u} + x_3\mathbf{v}
$$

解集是 $\operatorname{Span}\{\mathbf{u}, \mathbf{v}\}$（过原点的平面）。对非齐次版本 $10x_1 - 3x_2 - 2x_3 = b$（$b\neq 0$），解集是**平行平面** $\mathbf{p} + \operatorname{Span}\{\mathbf{u}, \mathbf{v}\}$。

**几何关系总结**：

![FIGURE 5 平行解集](img/fig15_5_parallel_lines.png)

*FIGURE 5 Parallel solution sets of $A\mathbf{x}=\mathbf{b}$ and $A\mathbf{x}=\mathbf{0}$. 图5 $A\mathbf{x}=\mathbf{b}$ 与 $A\mathbf{x}=\mathbf{0}$ 的解集互相平行*

![FIGURE 6 平行平面解集](img/fig15_6_parallel_solution_sets.png)

*FIGURE 6 Parallel solution sets of $A\mathbf{x}=\mathbf{b}$ and $A\mathbf{x}=\mathbf{0}$（两个自由变量时）. 图6 两自由变量时，解集是平行平面*

> **写作解集的步骤（WRITING A SOLUTION SET IN PARAMETRIC VECTOR FORM）**：
> 1. 把增广矩阵行化简到简化阶梯形
> 2. 把每个基本变量用自由变量表示
> 3. 把通解 $\mathbf{x}$ 写成自由变量是公因子的向量组合（参数向量形式）

### 1.5 小结

| 概念 | 一句话理解 |
|---|---|
| 齐次方程 $A\mathbf{x}=\mathbf{0}$ | 必有平凡解；有自由变量 ⇔ 有非平凡解 |
| 齐次解集 | 过原点的直线（1 个自由变量）/平面（2 个自由变量）= Span |
| THEOREM 6 | 非齐次通解 = 特解 $\mathbf{p}$ + 齐次通解 $\mathbf{v}_h$，几何上是**平移** |
| 参数向量形式 | 把自由变量提为系数，如 $\mathbf{x} = \mathbf{p} + s\mathbf{u} + t\mathbf{v}$ |

> **本节精髓**：齐次与非齐次的解集关系，就是"过原点"与"平移"的关系——这条几何直觉在 4.2 节会发展为严格的**向量空间（子空间）**理论。

---

## 第一章完 Chapter 1 Complete

至此第 1 章核心脉络完整：**解方程组（1.1–1.2）→ 线性组合与 Span（1.3）→ 矩阵方程 $A\mathbf{x}=\mathbf{b}$（1.4）→ 解集的结构：特解 + 齐次通解（1.5）**。接下来第 2 章将进入**矩阵代数**（逆矩阵、矩阵分解），第 1 章的一切都是它的地基。
---

## 2.2 逆矩阵 The Inverse of a Matrix

**中文讲解**：数的"倒数"概念推广到矩阵。回忆：5 的乘法逆元是 $5^{-1} = \frac{1}{5}$，满足 $5^{-1} \cdot 5 = 1$ 且 $5 \cdot 5^{-1} = 1$。矩阵推广需要**两个**方程都成立（因为矩阵乘法**不交换**），而且只有**方阵**才谈逆矩阵。

### 定义：可逆矩阵 DEFINITION

$n \times n$ 矩阵 $A$ 称为**可逆的（invertible）**，若存在 $n \times n$ 矩阵 $C$ 使

$$CA = I \quad \text{且} \quad AC = I$$

其中 $I = I_n$ 是 $n \times n$ 单位矩阵。这个 $C$ 由 $A$ **唯一决定**（证明：若 $B$ 也是逆，则 $B = BI = B(AC) = (BA)C = IC = C$），称为 $A$ 的逆，记作 $A^{-1}$：

$$A^{-1}A = I \quad \text{且} \quad AA^{-1} = I$$

- **奇异矩阵（singular）** = 不可逆矩阵；**非奇异矩阵（nonsingular）** = 可逆矩阵

**例1（Example 1）**：验证 $C = A^{-1}$

$$A = \begin{bmatrix} 2 & 5 \\
-3 & -7 \end{bmatrix}, \qquad C = \begin{bmatrix} -7 & -5 \\
3 & 2 \end{bmatrix}$$

$$AC = \begin{bmatrix} 1 & 0 \\
0 & 1 \end{bmatrix} = I \quad \text{且} \quad CA = \begin{bmatrix} 1 & 0 \\
0 & 1 \end{bmatrix} = I$$

两个方向都成立，所以 $C = A^{-1}$。

### THEOREM 4（2×2 矩阵的逆公式 + 判别法）

设

$$A = \begin{bmatrix} a & b \\
c & d \end{bmatrix}$$

若 $ad - bc \neq 0$，则 $A$ 可逆，且

$$A^{-1} = \frac{1}{ad-bc}\begin{bmatrix} d & -b \\
-c & a \end{bmatrix}$$

若 $ad - bc = 0$，则 $A$ 不可逆。

**行列式登场**：量 $ad - bc$ 称为 $A$ 的**行列式（determinant）**，记作 $\det A = ad - bc$。定理 4 说：**2×2 矩阵 $A$ 可逆 ⟺ $\det A \neq 0$**。

**例2（Example 2）**：求 $A = \begin{bmatrix} 3 & 4 \\
5 & 6 \end{bmatrix}$ 的逆。

$$\det A = 3(6) - 4(5) = -2 \neq 0 \Rightarrow A \text{ 可逆}$$

$$A^{-1} = \frac{1}{-2}\begin{bmatrix} 6 & -4 \\
-5 & 3 \end{bmatrix} = \begin{bmatrix} -3 & 2 \\
5/2 & -3/2 \end{bmatrix}$$

### THEOREM 5（可逆矩阵 ⇒ 方程组唯一解）

**中文讲解**：可逆矩阵在解方程时的作用，相当于"数的倒数解一元一次方程"。

> 若 $A$ 是可逆 $n \times n$ 矩阵，则对每个 $\mathbf{b} \in \mathbb{R}^n$，方程 $A\mathbf{x} = \mathbf{b}$ 有**唯一解** $\mathbf{x} = A^{-1}\mathbf{b}$。

**证明（小白版）**：存在性：代入 $\mathbf{x} = A^{-1}\mathbf{b}$，得 $A(A^{-1}\mathbf{b}) = (AA^{-1})\mathbf{b} = I\mathbf{b} = \mathbf{b}$ ✓。唯一性：若 $A\mathbf{u} = \mathbf{b}$，两边左乘 $A^{-1}$：$A^{-1}A\mathbf{u} = A^{-1}\mathbf{b} \Rightarrow I\mathbf{u} = A^{-1}\mathbf{b} \Rightarrow \mathbf{u} = A^{-1}\mathbf{b}$，即任何解都必等于 $A^{-1}\mathbf{b}$。

**例3（Example 3，工程应用）**：水平弹性梁两端支撑，在 1、2、3 点受力（FIGURE 1）。设 $\mathbf{f} \in \mathbb{R}^3$ 记录三点的力，$\mathbf{y} \in \mathbb{R}^3$ 记录三点的**挠度（deflection）**（梁的弯曲位移）。由物理学的胡克定律：

$$\mathbf{y} = D\mathbf{f}$$

$D$ 叫**柔度矩阵（flexibility matrix）**，其逆 $D^{-1}$ 叫**刚度矩阵（stiffness matrix）**。问：$D$ 和 $D^{-1}$ 的各列各有什么物理意义？

![FIGURE 1 弹性梁的挠度](img/fig22_1_elastic_beam.png)

*FIGURE 1 Deflection of an elastic beam. 图1 弹性梁在三点受力产生挠度*

**解答**：把 $I_3$ 按列拆开：$D = DI_3 = \begin{bmatrix} D\mathbf{e}_1 & D\mathbf{e}_2 & D\mathbf{e}_3 \end{bmatrix}$。

- $D$ 的第 1 列 $D\mathbf{e}_1$：在点 1 施加**单位力**（1 磅，向下）时，三个点的挠度各是多少。第 2、3 列同理。
- 再看 $D^{-1} = D^{-1}I_3 = \begin{bmatrix} D^{-1}\mathbf{e}_1 & D^{-1}\mathbf{e}_2 & D^{-1}\mathbf{e}_3 \end{bmatrix}$：$D^{-1}\mathbf{e}_1$ 是使点 1 产生**单位挠度**、其余两点挠度为零所需的三个力。
- 若柔度单位是"英寸/磅"，则刚度矩阵的单位是"磅/英寸"。

> ⚠️ **注意**：定理 5 的公式 $\mathbf{x} = A^{-1}\mathbf{b}$ 在**数值计算中很少用**——对 $[A\ \ \mathbf{b}]$ 做行化简几乎总是更快、更准（避免舍入误差累积）。唯一例外是 2×2 情形，用心算公式有时更快（例4）。

**例4（Example 4）**：用例2 的逆矩阵解方程组

$$
\begin{aligned}
3x_1 + 4x_2 &= 3 \\
5x_1 + 6x_2 &= 7
\end{aligned}
$$

即 $A\mathbf{x} = \mathbf{b}$，$\mathbf{b} = \begin{bmatrix} 3 \\
7 \end{bmatrix}$，所以

$$\mathbf{x} = A^{-1}\mathbf{b} = \begin{bmatrix} -3 & 2 \\
5/2 & -3/2 \end{bmatrix}\begin{bmatrix} 3 \\
7 \end{bmatrix} = \begin{bmatrix} 5 \\
-3 \end{bmatrix}$$

**验证**：$3(5) + 4(-3) = 3$ ✓；$5(5) + 6(-3) = 7$ ✓

### THEOREM 6（逆矩阵的三条运算律）

a. 若 $A$ 可逆，则 $A^{-1}$ 也可逆，且 $(A^{-1})^{-1} = A$

b. 若 $A, B$ 都是 $n \times n$ 可逆矩阵，则 $AB$ 也可逆，且**逆的顺序反过来**：

$$(AB)^{-1} = B^{-1}A^{-1}$$

c. 若 $A$ 可逆，则 $A^T$ 也可逆，且 $(A^T)^{-1} = (A^{-1})^T$

> **直觉记忆 (b)**：穿衣与脱衣——先穿袜子再穿鞋，脱的时候必须先鞋后袜（顺序反转）。这条"**穿脱衣法则**"在全书反复出现。

推广：任意 $n \times n$ 可逆矩阵的乘积仍可逆，逆 = 各逆矩阵按**相反顺序**的乘积。

### 初等矩阵 Elementary Matrices

**定义**：对**单位矩阵**做一次初等行变换得到的矩阵，叫**初等矩阵（elementary matrix）**。三种行变换对应三种初等矩阵。

**例5（Example 5）**：

$$E_1 = \begin{bmatrix} 1 & 0 & 0 \\
0 & 1 & 0 \\
-4 & 0 & 1 \end{bmatrix}, \qquad E_2 = \begin{bmatrix} 0 & 1 & 0 \\
1 & 0 & 0 \\
0 & 0 & 1 \end{bmatrix}, \qquad E_3 = \begin{bmatrix} 1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 5 \end{bmatrix}$$

$$A = \begin{bmatrix} a & b & c \\
d & e & f \\
g & h & i \end{bmatrix}$$

计算 $E_1A$、$E_2A$、$E_3A$：

$$E_1A = \begin{bmatrix} a & b & c \\
d & e & f \\
g-4a & h-4b & i-4c \end{bmatrix}, \qquad E_2A = \begin{bmatrix} d & e & f \\
a & b & c \\
g & h & i \end{bmatrix}, \qquad E_3A = \begin{bmatrix} a & b & c \\
d & e & f \\
5g & 5h & 5i \end{bmatrix}$$

**发现（重要！）**：左乘 $E_1$ = "第3行加 −4 倍第1行"；左乘 $E_2$ = "交换1、2行"；左乘 $E_3$ = "第3行乘 5"。**左乘初等矩阵的效果 = 对 A 做同样的行变换！**

### THEOREM 7（初等矩阵与行变换的桥梁）

> 若对 $m \times n$ 矩阵 $A$ 施行一次初等行变换，所得矩阵等于 $EA$，其中 $E$ 是对 $I_m$ 施行**同一行变换**得到的初等矩阵。

**推论**：每个初等矩阵 $E$ 都可逆，且 $E^{-1}$ 就是把 $E$ 变回 $I$ 的那个行变换对应的初等矩阵（行变换可逆 ⇒ 初等矩阵可逆）。

**例6（Example 6）**：

$$E_1 = \begin{bmatrix} 1 & 0 & 0 \\
0 & 1 & 0 \\
-4 & 0 & 1 \end{bmatrix} \quad \Rightarrow \quad E_1^{-1} = \begin{bmatrix} 1 & 0 & 0 \\
0 & 1 & 0 \\
4 & 0 & 1 \end{bmatrix}$$

（验证：$E_1^{-1}E_1 = I$ ✓）

### THEOREM 8（可逆 ⟺ 行等价于单位矩阵）

> $n \times n$ 矩阵 $A$ 可逆，**当且仅当** $A$ 行等价于 $I_n$。此时，任何把 $A$ 化简到 $I_n$ 的初等行变换序列，同时把 $I_n$ 变成 $A^{-1}$。

**这就是求逆矩阵的算法来源**：把 $A$ 和 $I$ 并排，对 $[A\ \ I]$ 行化简——左边变成 $I$ 的同时，右边自动变成 $A^{-1}$。

### 求 $A^{-1}$ 的算法 ALGORITHM

> 对增广矩阵 $[A\ \ I]$ 做行化简。若 $A$ 行等价于 $I$，则 $[A\ \ I]$ 行等价于 $[I\ \ A^{-1}]$；否则 $A$ 不可逆。

**例7（Example 7）**：求 $A = \begin{bmatrix} 0 & 1 & 2 \\
1 & 0 & 3 \\
4 & -3 & 8 \end{bmatrix}$ 的逆（若存在）。

$$
\left[\begin{array}{ccc|ccc} 0 & 1 & 2 & 1 & 0 & 0 \\
1 & 0 & 3 & 0 & 1 & 0 \\
4 & -3 & 8 & 0 & 0 & 1 \end{array}\right] \sim \left[\begin{array}{ccc|ccc} 1 & 0 & 3 & 0 & 1 & 0 \\
0 & 1 & 2 & 1 & 0 & 0 \\
4 & -3 & 8 & 0 & 0 & 1 \end{array}\right] \sim \left[\begin{array}{ccc|ccc} 1 & 0 & 3 & 0 & 1 & 0 \\
0 & 1 & 2 & 1 & 0 & 0 \\
0 & -3 & -4 & 0 & -4 & 1 \end{array}\right]
$$

（交换1、2行；第3行 += −4×第1行）

$$
\sim \left[\begin{array}{ccc|ccc} 1 & 0 & 3 & 0 & 1 & 0 \\
0 & 1 & 2 & 1 & 0 & 0 \\
0 & 0 & 2 & 3 & -4 & 1 \end{array}\right] \sim \left[\begin{array}{ccc|ccc} 1 & 0 & 3 & 0 & 1 & 0 \\
0 & 1 & 2 & 1 & 0 & 0 \\
0 & 0 & 1 & 3/2 & -2 & 1/2 \end{array}\right]
$$

（第3行 += 3×第2行；第3行 × 1/2）

$$
\sim \left[\begin{array}{ccc|ccc} 1 & 0 & 0 & -9/2 & 7 & -3/2 \\
0 & 1 & 0 & -2 & 4 & -1 \\
0 & 0 & 1 & 3/2 & -2 & 1/2 \end{array}\right] = \left[\begin{array}{ccc|ccc} I & A^{-1} \end{array}\right]
$$

（第1行 += −3×第3行；第2行 += −2×第3行）

$$A^{-1} = \begin{bmatrix} -9/2 & 7 & -3/2 \\
-2 & 4 & -1 \\
3/2 & -2 & 1/2 \end{bmatrix}$$

**验算**（好习惯）：

$$AA^{-1} = \begin{bmatrix} 0 & 1 & 2 \\
1 & 0 & 3 \\
4 & -3 & 8 \end{bmatrix}\begin{bmatrix} -9/2 & 7 & -3/2 \\
-2 & 4 & -1 \\
3/2 & -2 & 1/2 \end{bmatrix} = \begin{bmatrix} 1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1 \end{bmatrix} = I \checkmark$$

### 求逆的另一种视角 Another View of Matrix Inversion

把 $I_n$ 的列记为 $\mathbf{e}_1, \dots, \mathbf{e}_n$。则 $[A\ \ I] \sim [I\ \ A^{-1}]$ 可看作**同时解 $n$ 个方程组**：

$$A\mathbf{x} = \mathbf{e}_1, \qquad A\mathbf{x} = \mathbf{e}_2, \qquad \cdots, \qquad A\mathbf{x} = \mathbf{e}_n$$

这 $n$ 个方程组的解正好就是 $A^{-1}$ 的各列——因为 $AA^{-1} = I$ 按列看就是 $A \times (A^{-1}\text{ 的第 } j \text{ 列}) = \mathbf{e}_j$。

### 2.2 小结

| 概念 | 一句话理解 |
|---|---|
| $A^{-1}$ | 满足 $A^{-1}A = AA^{-1} = I$ 的唯一矩阵（方阵才有） |
| $\det A = ad-bc$ | 2×2 可逆 ⟺ 行列式非零 |
| Theorem 5 | $A$ 可逆 ⟹ $A\mathbf{x}=\mathbf{b}$ 唯一解 $\mathbf{x}=A^{-1}\mathbf{b}$ |
| Theorem 6(b) | $(AB)^{-1} = B^{-1}A^{-1}$（**穿脱衣法则**：顺序反转） |
| 初等矩阵 | 对 $I$ 做一次行变换；左乘它 = 对 $A$ 做同样变换 |
| Theorem 8 + 算法 | $A$ 可逆 ⟺ $A \sim I$；$[A\ \ I] \sim [I\ \ A^{-1}]$ 求逆 |

> **本节精髓**：逆矩阵是"倒数"思想的矩阵化。行变换与初等矩阵的对应（Theorem 7）把"解方程"与"矩阵乘法"焊接在一起，直接催生出行化简求逆算法——这是线性代数"**计算**"与"**代数结构**"两条主线的重要交汇点。2.3 节的**可逆矩阵定理（IMT）**将给出 20+ 个等价条件，把全章知识一网打尽。

---
---

## 2.3 可逆矩阵的特征 Characterizations of Invertible Matrices

**中文讲解**：本节把前 4 章的关键概念串成一条"项链"——**可逆矩阵定理（Invertible Matrix Theorem, IMT）**。它把 12 个看似无关的命题连成等价关系：验证其中**任何一个**，其余全部自动成立。这是全书被引用最多的定理。

### THE INVERTIBLE MATRIX THEOREM（可逆矩阵定理）

> 设 $A$ 是 $n \times n$ 矩阵，则下列命题**逻辑等价**（同真同假）：
>
> a. $A$ 是**可逆矩阵**
>
> b. $A$ **行等价于** $n \times n$ 单位矩阵 $I_n$
>
> c. $A$ 有 $n$ 个**主元位置**
>
> d. 方程 $A\mathbf{x} = \mathbf{0}$ 仅有**平凡解**
>
> e. $A$ 的各列**线性无关**
>
> f. 线性变换 $\mathbf{x} \mapsto A\mathbf{x}$ 是**一对一（one-to-one）**的
>
> g. 方程 $A\mathbf{x} = \mathbf{b}$ 对每个 $\mathbf{b} \in \mathbb{R}^n$ **至少有一个解**
>
> h. $A$ 的各列**张成** $\mathbb{R}^n$
>
> i. 线性变换 $\mathbf{x} \mapsto A\mathbf{x}$ 是**映上（onto）**的
>
> j. 存在 $n \times n$ 矩阵 $C$ 使 $CA = I$
>
> k. 存在 $n \times n$ 矩阵 $D$ 使 $AD = I$
>
> l. $A^T$ 是**可逆矩阵**

**证明思路（教材 Figure 1）**：关键链条是 (c) ⇒ (d) ⇒ (j) ⇒ (a) ⇒ (b) ⇒ (c)——绕一圈把所有命题串起来，再加上各种"直路"连接。其中 (a) ⇒ (l) 因为由 2.1 定理 6(c) 知 $(A^T)^{-1} = (A^{-1})^T$。

![FIGURE 1 IMT 证明的逻辑链条](img/fig23_1_imt_circle.png)

*FIGURE 1 The connections between the statements of the Invertible Matrix Theorem. 图1 IMT 各命题间的等价关系：左边圆圈是证明主干 (c)⇒(d)⇒(j)⇒(a)⇒(b)⇒(c)*

> **小白使用指南**：IMT 的威力在于"一证全证"。做题时只需验证**最容易验证**的一条（通常是对 $A$ 行化简看主元数），就能调用其余 11 条的全部结论。注意 **IMT 只对 $n \times n$ 方阵成立**！

**警告**：2.4 节定理 9 会把 (g)、(h)、(i) 推广到**任何矩阵**（不必是方阵）——那是"存在性定理"，但 (d)、(e)、(f) 的对应推广**不成立**。

### 例8（Example 8，IMT 实战）

**题目**：用 IMT 判断 $A$ 是否可逆：

$$A = \begin{bmatrix} 1 & 0 & -2 \\
3 & 1 & -2 \\
-5 & -1 & 9 \end{bmatrix}$$

**解**：对 $A$ 行化简（验证命题 (b)：$A$ 是否行等价于 $I$）：

$$A \sim \begin{bmatrix} 1 & 0 & -2 \\
0 & 1 & 4 \\
0 & 0 & 3 \end{bmatrix}$$

A 有 **3 个主元**（命题 (c) 成立）⇒ **A 可逆**，且 IMT 其余命题全部成立：方程 $A\mathbf{x}=\mathbf{b}$ 对每个 $\mathbf{b}$ 有唯一解、$A\mathbf{x}=\mathbf{0}$ 只有平凡解、列线性无关且张成 $\mathbb{R}^3$、变换既一对一又映上……

> **效率秘诀**：例8 只把 $A$ 化到**阶梯形**（不必到简化阶梯形）就能数主元、下结论——这是 IMT 实战的标准操作。

### 推论（THEOREM 9 预告）

若 $A, B$ 是方阵且 $AB = I$，则 $A, B$ 都可逆，且 $B = A^{-1}$、$A = B^{-1}$。

> **直觉**：IMT 的 (j)、(k) 说：对方阵而言，**单侧**逆就是双侧逆！这和"$AB=I \nRightarrow BA=I$"（2.1 节的警示）**不矛盾**——2.1 节说的反例都是**非方阵**。对方阵，单边成立 ⇒ 双边自动成立。

### 深入：IMT 的"扩展表"（定理 2-10 的完整版）

把后面章节的结果并入，IMT 实际可用的等价条件超过 **20 条**。核心几组对照（小白先混个眼熟，后面章节会逐一展开）：

| 视角 | 可逆（全部成立） | 不可逆（全部不成立） |
|---|---|---|
| 行化简 | $A \sim I_n$（b） | 主元数 < $n$ |
| 方程组 | $A\mathbf{x}=\mathbf{0}$ 只有平凡解（d） | 存在非平凡解 |
| 列向量 | 列线性无关（e）且张成 $\mathbb{R}^n$（h） | 列相关或张不满 |
| 线性变换 | 既一对一（f）又映上（i） | 至少破坏其一 |
| 行列式 | $\det A \neq 0$（第3章） | $\det A = 0$ |
| 特征值 | $0$ **不是**特征值（第5章） | $0$ 是特征值 |

### 例（补充练习，选自 2.3 Exercises 风格）

**判断题**："若 $A, B$ 是 $n \times n$ 矩阵且 $AB$ 可逆，则 $A$ 可逆。"——**True**。

**证明**：$AB$ 可逆 ⇒ 存在 $C$ 使 $C(AB) = I$ ⇒ $(CA)B = I$。由定理 9 推论，$B$ 可逆且 $B^{-1} = CA$；再由 2.1 定理 6，$A = C^{-1}B^{-1}$…… 更简洁：$AB$ 可逆 ⇒ $AB$ 有 $n$ 个主元（IMT (c)）⇒ $A$ 必有 $n$ 个主元（$A$ 的列数 = $n$，行化简 $A$ 时不可能比 $AB$ 的秩更多……用 $(AB)x=0$ 只有平凡解 ⇒ $Bx=0$ 时 $A(Bx)=0$ ⇒ $Bx=0$…… 故 $B$ 列无关 ⇒ $B$ 可逆 ⇒ $A = (AB)B^{-1}$ 可逆 ✓）。

### 2.3 小结

| 概念 | 一句话理解 |
|---|---|
| **IMT** | 12 条等价命题"一荣俱荣、一损俱损"，验证最易的一条即可调用全部 |
| 验证首选 | 行化简到阶梯形，数主元是否 = $n$（命题 (c)） |
| 方阵特权 | 单侧逆 $CA=I$ ⇒ 双侧逆，$C = A^{-1}$ |
| 使用前提 | **必须是 $n \times n$ 方阵**！非方阵请用 1.4 定理 4 与 2.4 定理 9 |

> **本节精髓**：IMT 是线性代数的"**总开关**"——它把"解方程"（d, g）、"向量关系"（e, h）、"变换性质"（f, i）、"矩阵运算"（a, j, k, l）、"行化简"（b, c）五大主题统一在"可逆"这一个概念之下。后续学习每章新概念（行列式、秩、特征值、对角化），请养成习惯问自己："这给 IMT 的清单添了哪一条等价条件？"

---
---

## 2.4 分块矩阵 Partitioned Matrices

**中文讲解**：把矩阵看作**列向量列表**的观点（1.4 节的核心）大获成功，自然推广：用水平、竖直虚线把矩阵切成**子矩阵（submatrix / block）**，按"块"参与运算。这就是**模块化思想**——像乐高积木一样拼拆矩阵。现代应用（飞机设计、电路板、大规模线性规划）中分块记号无处不在。

### 例1（Example 1）：分块记号

$$
A = \begin{bmatrix} 3 & 0 & -1 & 5 & 9 & -2 \\
-5 & 2 & 4 & 0 & -3 & 1 \\
-8 & -6 & 3 & 1 & 7 & -4 \end{bmatrix} = \begin{bmatrix} A_{11} & A_{12} & A_{13} \\
A_{21} & A_{22} & A_{23} \end{bmatrix}
$$

按 $2 \times 3$ 分块，各子矩阵为：

$$A_{11} = \begin{bmatrix} 3 & 0 & -1 \\
-5 & 2 & 4 \end{bmatrix}, \qquad A_{12} = \begin{bmatrix} 5 & 9 \\
0 & -3 \end{bmatrix}, \qquad A_{13} = \begin{bmatrix} -2 \\
1 \end{bmatrix}$$

$$A_{21} = \begin{bmatrix} -8 & -6 & 3 \end{bmatrix}, \qquad A_{22} = \begin{bmatrix} 1 & 7 \end{bmatrix}, \qquad A_{23} = \begin{bmatrix} -4 \end{bmatrix}$$

**例2（Example 2，应用背景）**：若微机电路板由三块 VLSI 芯片组成，电路板矩阵自然分块为

$$
A = \begin{bmatrix} A_{11} & A_{12} & A_{13} \\
A_{21} & A_{22} & A_{23} \\
A_{31} & A_{32} & A_{33} \end{bmatrix}
$$

每块对应一颗芯片的信息——**分块让矩阵的"物理结构"一目了然**。

### 分块矩阵的加法与数乘

**规则（与数的运算完全类似）**：两个**分块方式相同**的矩阵相加，对应块相加；数乘则每个块分别乘。

### 分块矩阵的乘法 THE KEY RULE

**中文讲解**：分块矩阵可以像普通矩阵一样用**行-列法则**相乘——**把每个块当作"数"**！但有前提：

> **相容条件（conformable）**：对乘积 $AB$，$A$ 的**列分法**必须与 $B$ 的**行分法**一致。

**例3（Example 3，完整演示）**：

$$
A = \begin{bmatrix} 2 & -3 & 1 & 0 & -4 \\
1 & 5 & -2 & 3 & -1 \\
0 & -4 & 2 & 7 & -1 \end{bmatrix} = \begin{bmatrix} A_{11} & A_{12} \\
A_{21} & A_{22} \end{bmatrix}, \qquad B = \begin{bmatrix} 6 & 4 \\
-2 & 1 \\
-3 & 7 \\
-1 & 3 \\
5 & 2 \end{bmatrix} = \begin{bmatrix} B_1 \\
B_2 \end{bmatrix}
$$

$A$ 的 5 列分成 $3+2$，$B$ 的 5 行也分成 $3+2$——**分法匹配**，可以分块相乘：

$$AB = \begin{bmatrix} A_{11}B_1 + A_{12}B_2 \\
A_{21}B_1 + A_{22}B_2 \end{bmatrix} = \begin{bmatrix} -5 & 4 \\
-6 & 2 \\
2 & 1 \end{bmatrix}$$

逐个验证（注意**子矩阵必须 A 在左、B 在右**，乘法不交换！）：

$$A_{11}B_1 = \begin{bmatrix} 2 & -3 & 1 \\
1 & 5 & -2 \end{bmatrix}\begin{bmatrix} 6 & 4 \\
-2 & 1 \\
-3 & 7 \end{bmatrix} = \begin{bmatrix} 15 & 12 \\
2 & -5 \end{bmatrix}$$

$$A_{12}B_2 = \begin{bmatrix} 0 & -4 \\
3 & -1 \end{bmatrix}\begin{bmatrix} -1 & 3 \\
5 & 2 \end{bmatrix} = \begin{bmatrix} -20 & -8 \\
-8 & 7 \end{bmatrix}$$

$$A_{11}B_1 + A_{12}B_2 = \begin{bmatrix} 15 & 12 \\
2 & -5 \end{bmatrix} + \begin{bmatrix} -20 & -8 \\
-8 & 7 \end{bmatrix} = \begin{bmatrix} -5 & 4 \\
-6 & 2 \end{bmatrix}$$

> **视角统一**：分块乘法的行-列法则是矩阵乘法**最一般**的看法。前面学过的都是特例：(1) $A\mathbf{x}$ 的列定义（$A$ 按列分块）；(2) $AB$ 的列定义（$B$ 按列分块）；(3) 行视角。

### THEOREM 10（AB 的列-行展开 Column-Row Expansion）

**中文讲解**：$A$ 按**列**分、$B$ 按**行**分，得到 $AB$ 的又一种表达——**外积之和**：

$$AB = \begin{bmatrix} \text{col}_1(A) & \text{col}_2(A) & \cdots & \text{col}_n(A) \end{bmatrix} \begin{bmatrix} \text{row}_1(B) \\
\text{row}_2(B) \\
\vdots \\
\text{row}_n(B) \end{bmatrix} = \text{col}_1(A)\,\text{row}_1(B) + \cdots + \text{col}_n(A)\,\text{row}_n(B)$$

**证明（小白版）**：$\text{col}_k(A)\,\text{row}_k(B)$ 的 $(i,j)$ 元是 $a_{ik}b_{kj}$。对 $k$ 求和：$a_{i1}b_{1j} + a_{i2}b_{2j} + \cdots + a_{in}b_{nj}$——这正是行-列法则给出的 $AB$ 的 $(i,j)$ 元 ✓

**数值注记**：计算机中 AB 的列-行展开把 $A$ 和 $B$ 切成小块，只把小矩阵送进高速内存相乘——这是处理**超大规模矩阵**的标准技术（见本节末尾注记：某线性规划问题矩阵有 837 行、51 块列、每块约 25 万列！）。

### 分块矩阵求逆（例5，重点公式）

**中文讲解**：形如

$$A = \begin{bmatrix} A_{11} & A_{12} \\
0 & A_{22} \end{bmatrix}$$

的**分块上三角矩阵**（$A_{11}$ 是 $p \times p$，$A_{22}$ 是 $q \times q$），求逆有简洁公式。设逆矩阵也按同法分块 $B = \begin{bmatrix} B_{11} & B_{12} \\
B_{21} & B_{22} \end{bmatrix}$，则

$$AB = \begin{bmatrix} A_{11} & A_{12} \\
0 & A_{22} \end{bmatrix}\begin{bmatrix} B_{11} & B_{12} \\
B_{21} & B_{22} \end{bmatrix} = \begin{bmatrix} I_p & 0 \\
0 & I_q \end{bmatrix}$$

分块相乘得到**四个方程**，逐一破解：

$$
\begin{aligned}
A_{11}B_{11} + A_{12}B_{21} &= I_p \\
A_{11}B_{12} + A_{12}B_{22} &= 0 \\
A_{22}B_{21} &= 0 \\
A_{22}B_{22} &= I_q
\end{aligned}
$$

**求解顺序**（从右下往左上，逐个击破）：

1. **(6)**：$A_{22}$ 是方阵，由 IMT 知 $A_{22}$ 可逆，$B_{22} = A_{22}^{-1}$
2. **(5)**：两边左乘 $A_{22}^{-1}$：$B_{21} = A_{22}^{-1} \cdot 0 = 0$
3. **(3)**：简化为 $A_{11}B_{11} = I_p$ ⇒ $A_{11}$ 可逆，$B_{11} = A_{11}^{-1}$
4. **(4)**：$A_{11}B_{12} = -A_{12}B_{22} = -A_{12}A_{22}^{-1}$ ⇒ $B_{12} = -A_{11}^{-1}A_{12}A_{22}^{-1}$

**结论公式（务必记住）**：

$$\begin{bmatrix} A_{11} & A_{12} \\
0 & A_{22} \end{bmatrix}^{-1} = \begin{bmatrix} A_{11}^{-1} & -A_{11}^{-1}A_{12}A_{22}^{-1} \\
0 & A_{22}^{-1} \end{bmatrix}$$

> **对照 2×2 三角矩阵**：$\begin{bmatrix} a & b \\
0 & d \end{bmatrix}^{-1} = \frac{1}{ad}\begin{bmatrix} d & -b \\
0 & a \end{bmatrix}$——形式完全一致，分块公式是它的"矩阵版"！

**分块对角矩阵（block diagonal）**：主对角块之外全是零块。它可逆 ⟺ **每个对角块都可逆**（练习 13、14）。

### 数值注记（为什么分块重要）

1. **内存管理**：矩阵太大装不进高速内存时，分块让计算机每次只取两三个子矩阵工作
2. **并行架构**：向量流水线计算机用分块算法效率更高
3. **工业软件**：LAPACK（高性能数值线性代数库）大量使用分块计算

### 2.4 小结

| 概念 | 一句话理解 |
|---|---|
| 分块矩阵 | 用横竖虚线切出子矩阵，按"块"运算 |
| 分块乘法 | 块当作数用行-列法则；**A 的列分法 = B 的行分法** |
| Theorem 10 | $AB$ = 各列-行外积之和（超大规模计算的基础） |
| 分块上三角求逆 | 右下角先求逆，逐块上推；$B_{12} = -A_{11}^{-1}A_{12}A_{22}^{-1}$ |
| 分块对角 | 可逆 ⟺ 每个对角块可逆 |

> **本节精髓**：分块是线性代数的"**结构化思维**"——好分法让复杂矩阵的结构显形（芯片对应块、稀疏结构对应零块）。Theorem 10 的外积展开更是连接"第 4 章秩的理论"与"第 7 章 SVD"的暗线，先在此留个印象。

---
---

## 2.5 矩阵分解 Matrix Factorizations

**中文讲解**：矩阵分解 = 把矩阵 $A$ 写成几个**更简单矩阵的乘积**。本节主角：**LU 分解** $A = LU$（$L$ 下三角、$U$ 阶梯形）。它是工程计算的核心技术——把 $A$ 的"消元工作"存进 $L$，之后解任何 $A\mathbf{x} = \mathbf{b}$ 只需**回代**，不必重复消元。

### 为什么需要 LU 分解？（动机）

设 $A$ 是 $n \times n$ 可逆矩阵。解单个 $A\mathbf{x} = \mathbf{b}$ 用行化简；但工程中常见**同一 $A$、25 个不同右端 $\mathbf{b}_1, \dots, \mathbf{b}_{25}$** 的场景（波音论文：机翼设计迭代）。做法：

$$[A\ \ \mathbf{b}_1\ \ \mathbf{b}_2\ \ \cdots\ \ \mathbf{b}_{25}] \sim [I\ \ \mathbf{x}_1\ \ \mathbf{x}_2\ \ \cdots\ \ \mathbf{x}_{25}]$$

虽然比逐个求逆快得多，但若有新的 $\mathbf{b}_{26}$ 进来，又得重新化简整张大矩阵。

**LU 分解的妙招**：把消元过程"存档"。解 $A\mathbf{x} = \mathbf{b}$ 分三步（设 $A = LU$）：

1. 解 $L\mathbf{y} = \mathbf{b}$（$L$ 下三角，**前代** forward substitution，快！）
2. 解 $U\mathbf{x} = \mathbf{y}$（$U$ 上三角，**回代** back substitution，快！）
3. 得 $A\mathbf{x} = LU\mathbf{x} = L\mathbf{y} = \mathbf{b}$ ✓

**运算量对比**：$n \times n$ 矩阵行化简约 $2n^3/3$ flops；而三角求解只需 $n^2$ flops。25 个右端项：行化简方案 ≈ $25 \times 2n^3/3$，LU 方案 ≈ $2n^3/3 + 25 \times 2n^2$——**省出一个数量级**！

### LU 分解的存在性 THEOREM（关键限制）

> **基本定理**：若 $n \times n$ 矩阵 $A$ 可以**只用行替换**（row replacement）**不做行交换**就化为阶梯形 $U$，则 $A = LU$ 成立，其中 $L$ 是**单位下三角矩阵**（对角线全 1、对角线以下存"消元倍数"）。
>
> **一句话：LU 分解存在 ⟺ 消元过程中不需要换行。**

**为什么行交换是障碍？** 行交换对应置换矩阵 $P$，需要分解推广为 $PA = LU$（2.5 节练习涉及）。

### 构造 L 的原理（理解性推导）

只用行替换：$E_p \cdots E_1 A = U$，每个 $E_k$ 都是**单位下三角**初等矩阵（"某行 += c×另一行"型）。

**关键事实**：单位下三角矩阵的乘积、逆仍是单位下三角（练习 19）。令

$$L = (E_p \cdots E_1)^{-1} = E_1^{-1} \cdots E_p^{-1}$$

则 $A = (E_p \cdots E_1)^{-1} U = LU$，且 $L$ 单位下三角。

**妙处**：把 $A$ 化到 $U$ 的**同一串行变换**，作用在 $L$ 上恰好把它变回 $I$（因为 $E_p \cdots E_1 L = I$）——所以 **$L$ 中存的正是每次行变换的倍数**（对角线以下）。

### LU 分解算法 ALGORITHM FOR AN LU FACTORIZATION

1. **把 $A$ 只用行替换化为阶梯形 $U$**（若必须换行则基本 LU 不存在）
2. **填 $L$**：在每个主元列，把该列中**高亮条目（确定行变换的数）除以主元**，商放进 $L$ 的对应位置；$L$ 对角线全填 1

### 例2（Example 2，完整实战）

求 $A$ 的 LU 分解：

$$A = \begin{bmatrix} 2 & 4 & -1 & 5 & -2 \\
-4 & -5 & 3 & -8 & 1 \\
2 & -5 & -4 & 1 & 8 \\
-6 & 0 & 7 & -3 & 1 \end{bmatrix}$$

**第1列**：$A$ 的第一列 $[2, -4, 2, -6]^T$ 除以主元 2，得 $L$ 的第一列 $[1, -2, 1, -3]^T$。

**行变换**（消第一列）：R2 += 2R1，R3 += −R1，R4 += 3R1：

$$A \sim \begin{bmatrix} 2 & 4 & -1 & 5 & -2 \\
0 & \mathbf{3} & 1 & 2 & -3 \\
0 & -9 & -3 & -4 & 10 \\
0 & 12 & 4 & 12 & -5 \end{bmatrix} = A_1$$

（高亮主元 3；高亮列下方条目 $-9, 12$ 将用于填 $L$ 第二列）

**第2列**：$(-9, 12)/3 = (-3, 4)$，填入 $L(3,2) = -3$，$L(4,2) = 4$。

**行变换**（消第二列）：R3 += 3R2，R4 += −4R2：

$$A_1 \sim \begin{bmatrix} 2 & 4 & -1 & 5 & -2 \\
0 & 3 & 1 & 2 & -3 \\
0 & 0 & 0 & \mathbf{2} & 1 \\
0 & 0 & 0 & 4 & 7 \end{bmatrix} = A_2$$

**第3列（主元列为第4列）**：主元 2 下方条目 $4$，$4/2 = 2$，填入 $L(4,3) = 2$。

**行变换**：R4 += −2R3：

$$A_2 \sim \begin{bmatrix} 2 & 4 & -1 & 5 & -2 \\
0 & 3 & 1 & 2 & -3 \\
0 & 0 & 0 & 2 & 1 \\
0 & 0 & 0 & 0 & 5 \end{bmatrix} = U$$

**组装结果**：

$$L = \begin{bmatrix} 1 & 0 & 0 & 0 \\
-2 & 1 & 0 & 0 \\
1 & -3 & 1 & 0 \\
-3 & 4 & 2 & 1 \end{bmatrix}, \qquad U = \begin{bmatrix} 2 & 4 & -1 & 5 & -2 \\
0 & 3 & 1 & 2 & -3 \\
0 & 0 & 0 & 2 & 1 \\
0 & 0 & 0 & 0 & 5 \end{bmatrix}$$

**验算**：$LU = A$ ✓（读者可按行-列法则抽查第 2 行：$-2 \times [2\ 4\ -1\ 5\ -2] + [0\ 3\ 1\ 2\ -3] = [-4\ -5\ 3\ -8\ 1]$ ✓）

> **小白技巧**：填 $L$ 的口诀——"**主元除高亮，对角填 1**"。只需盯着行化简过程中每个主元列的高亮条目，除以主元即得 $L$ 对应列。

### LU 分解解方程组实战（练习题风格）

解 $A\mathbf{x} = \mathbf{b}$，其中 $A = LU$ 如上，$\mathbf{b} = [2, -4, 5, -6]^T$：

**第1步 前代** $L\mathbf{y} = \mathbf{b}$（从上到下逐行解）：

$$
\begin{aligned}
y_1 &= 2 \\
-2y_1 + y_2 &= -4 & \Rightarrow\ y_2 &= 0 \\
y_1 - 3y_2 + y_3 &= 5 & \Rightarrow\ y_3 &= 3 \\
-3y_1 + 4y_2 + 2y_3 + y_4 &= -6 & \Rightarrow\ y_4 &= 4
\end{aligned}
$$

**第2步 回代** $U\mathbf{x} = \mathbf{y}$（从下到上）：

$$
\begin{aligned}
5x_5 &= 4 & \Rightarrow\ x_5 &= \tfrac{4}{5} \\
2x_4 + x_5 &= 3 & \Rightarrow\ x_4 &= \tfrac{11}{10} \\
&\ \vdots
\end{aligned}
$$

（逐行向上即可解出全部 $x_i$——每个方程只含一个未知量！）

### 2.5 小结

| 概念 | 一句话理解 |
|---|---|
| $A = LU$ | 消元过程存档：$L$ 存倍数、$U$ 存阶梯形 |
| 存在条件 | 消元**不需换行**（否则要 $PA = LU$） |
| 解方程三步 | 前代 $L\mathbf{y} = \mathbf{b}$ → 回代 $U\mathbf{x} = \mathbf{y}$ |
| 效率 | $n^2$ vs $2n^3/3$，多右端问题时完胜行化简 |
| 填 $L$ 口诀 | 主元除高亮，对角填 1 |

> **本节精髓**：LU 分解把"一次昂贵的消元"与"多次廉价的求解"分离——这正是数值线性代数"**分解优先于求逆**"原则的第一个实例。第 7 章的 QR 分解、SVD 都是同一思想的升级版。

---
---

## 2.6 列昂惕夫投入产出模型 The Leontief Input-Output Model

**中文讲解**：回到第 1 章开头的列昂惕夫（1973 年诺贝尔经济学奖）！他的投入产出模型是第一个真正大规模应用线性代数的经济学模型，至今仍在世界各国使用。本节用 $(I-C)\mathbf{x} = \mathbf{d}$ 回答一个实际经济问题：**各部门各生产多少，才能恰好满足全社会的需求？**

### 模型设定 Setup

把国民经济分成 $n$ 个生产部门（制造、农业、服务业……），另有一个**开放部门（open sector）**只消费不生产（消费者、政府、出口）。

- **生产向量** $\mathbf{x} \in \mathbb{R}^n$：各部门一年的产出
- **最终需求向量** $\mathbf{d} \in \mathbb{R}^n$：开放部门对各部门产品的需求
- **单位消费向量** $\mathbf{c}_j$：第 $j$ 部门每生产 **1 单位**产品，需要从各部门"订购"多少投入

**平衡方程（核心）**：

$$\underbrace{\mathbf{x}}_{\text{生产量}} = \underbrace{C\mathbf{x}}_{\text{中间需求}} + \underbrace{\mathbf{d}}_{\text{最终需求}} \tag{4}$$

- 第 $j$ 部门计划生产 $x_j$ 单位，就要消耗 $x_j\mathbf{c}_j$（自己的投入需求）
- 所有部门的中间需求之和：$x_1\mathbf{c}_1 + x_2\mathbf{c}_2 + x_3\mathbf{c}_3 = C\mathbf{x}$，其中

$$C = \begin{bmatrix} \mathbf{c}_1 & \mathbf{c}_2 & \mathbf{c}_3 \end{bmatrix}$$

叫**消费矩阵（consumption matrix）**——注意 $\mathbf{c}_j$ 是 $C$ 的第 $j$ **列**（列 j 记录"生产 1 单位 j 产品要消耗什么"，按"从谁那里购买"分行）。

**例1（Example 1）**：制造业生产 100 单位要消耗什么？

$$100\mathbf{c}_1 = 100\begin{bmatrix} .50 \\
.20 \\
.10 \end{bmatrix} = \begin{bmatrix} 50 \\
20 \\
10 \end{bmatrix}$$

即从制造业内部订购 50、从农业订购 20、从服务业订购 10。

### 教材数值例（Table 1 + Example 2 完整求解）

三个部门的消费矩阵（**列 = 生产 1 单位该部门产品的投入来源**，表 1）：

| 投入来源 ↓ / 产品 → | 制造业 | 农业 | 服务业 |
|---|---|---|---|
| 制造业 | .50 | .40 | .20 |
| 农业 | .20 | .30 | .10 |
| 服务业 | .10 | .10 | .30 |

$$C = \begin{bmatrix} .50 & .40 & .20 \\
.20 & .30 & .10 \\
.10 & .10 & .30 \end{bmatrix}$$

**例2（Example 2）**：最终需求 $\mathbf{d} = \begin{bmatrix} 50 \\
30 \\
20 \end{bmatrix}$，求生产水平 $\mathbf{x}$。

方程 $(I - C)\mathbf{x} = \mathbf{d}$，系数矩阵：

$$I - C = \begin{bmatrix} 1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1 \end{bmatrix} - \begin{bmatrix} .50 & .40 & .20 \\
.20 & .30 & .10 \\
.10 & .10 & .30 \end{bmatrix} = \begin{bmatrix} .5 & -.4 & -.2 \\
-.2 & .7 & -.1 \\
-.1 & -.1 & .7 \end{bmatrix}$$

增广矩阵行化简（×10 消小数）：

$$\left[\begin{array}{rrrr} 5 & -4 & -2 & 500 \\
-2 & 7 & -1 & 300 \\
-1 & -1 & 7 & 200 \end{array}\right] \sim \cdots \sim \left[\begin{array}{rrrr} 1 & 0 & 0 & 226 \\
0 & 1 & 0 & 119 \\
0 & 0 & 1 & 78 \end{array}\right]$$

**答案**：制造业约 **226 单位**、农业约 **119 单位**、服务业约 **78 单位**（末列取整）。

> **观察**：最终需求只要 50+30+20 = 100 单位，生产却要做 226+119+78 = 423 单位！**中间需求（部门互相购买）是最终需求的 3 倍多**——这正是模型的价值：直观算不清的链式反应，线性代数一步搞定。

### THEOREM 11（存在性 + 经济可行性）

> 设 $C$ 是消费矩阵、$\mathbf{d}$ 是最终需求。若 $C$ 和 $\mathbf{d}$ 的元素都非负，且 **$C$ 的每列元素之和都小于 1**（每部门生产 1 单位的投入成本小于 1 单位——经济学常识），则 $(I-C)^{-1}$ 存在，且

$$\mathbf{x} = (I-C)^{-1}\mathbf{d}$$

元素全非负，是 $\mathbf{x} = C\mathbf{x} + \mathbf{d}$ 的唯一解。

### $(I-C)^{-1}$ 的级数公式（美的发现）

**中文讲解**：直接求 $(I-C)^{-1}$ 很贵。换个思路：年初先按 $\mathbf{x} = \mathbf{d}$ 生产，各部门为生产 $\mathbf{d}$ 又要订购 $C\mathbf{d}$，为满足 $C\mathbf{d}$ 又要生产 $C^2\mathbf{d}$……需求一轮轮传导：

| 需要满足的 | 为达此需求又产生的 |
|---|---|
| 最终需求 $\mathbf{d}$ | $C\mathbf{d}$ |
| 第 1 轮中间需求 $C\mathbf{d}$ | $C(C\mathbf{d}) = C^2\mathbf{d}$ |
| 第 2 轮 $C^2\mathbf{d}$ | $C^3\mathbf{d}$ |
| 第 3 轮 $C^3\mathbf{d}$ | $C^4\mathbf{d}$ |
| ⋮ | ⋮ |

全部需求之和：

$$\mathbf{x} = \mathbf{d} + C\mathbf{d} + C^2\mathbf{d} + C^3\mathbf{d} + \cdots = (I + C + C^2 + C^3 + \cdots)\mathbf{d} \tag{6}$$

**代数恒等式**（可验证）：

$$(I - C)(I + C + C^2 + \cdots + C^m) = I - C^{m+1} \tag{7}$$

当 $C$ 的列和都 < 1 时，$C^m \to 0$（类比：$|t|<1$ 时 $t^m \to 0$），所以 $I - C^{m+1} \to I$，于是

$$(I-C)^{-1} \approx I + C + C^2 + C^3 + \cdots + C^m \tag{8}$$

$m$ 越大越精确。**实际模型中 $C^m$ 收敛极快**（部门链条几轮就衰减没了），所以 (8) 是实用的计算方法！(6) 还自动保证：$C,\mathbf{d}$ 非负 ⇒ $\mathbf{x}$ 非负。

### $(I-C)^{-1}$ 元素的经济意义

$(I-C)^{-1}$ 的第 $j$ **列** = 当第 $j$ 部门的最终需求**增加 1 单位**时，各部门需要**追加**的产量。它是经济预测的敏感分析工具（见练习 8）。

### 例3（Example 3，建模练习）

设两部门经济体：产品（goods）与服务（services）。

- 生产 1 单位产品需要：0.2 单位产品投入 + 0.5 单位服务投入
- 生产 1 单位服务需要：0.4 单位产品投入 + 0.3 单位服务投入
- 最终需求：产品 20 单位，服务 30 单位

**建立模型**：

$$C = \begin{bmatrix} .2 & .4 \\
.5 & .3 \end{bmatrix}, \qquad \mathbf{d} = \begin{bmatrix} 20 \\
30 \end{bmatrix}$$

$$(I - C)\mathbf{x} = \mathbf{d} \quad \text{即} \quad \begin{bmatrix} .8 & -.4 \\
-.5 & .7 \end{bmatrix}\begin{bmatrix} x_1 \\
x_2 \end{bmatrix} = \begin{bmatrix} 20 \\
30 \end{bmatrix}$$

（求解作为练习，答案：$x_1 = 50$，$x_2 = 50$。）

### 数值注记：条件数

任何应用问题（不限于经济）中，$A\mathbf{x} = \mathbf{b}$ 的解对 $A$ 或 $\mathbf{b}$ 的小扰动的敏感程度由**条件数（condition number）**刻画。病态（ill-conditioned）方程组：输入数据的测量误差（或舍入误差）会被放大到解中。奇异值分解（第 7 章）给出条件数的严格定义。

### 2.6 小结

| 概念 | 一句话理解 |
|---|---|
| 消费矩阵 $C$ | 列 $j$ = 生产 1 单位产品 $j$ 的全部投入 |
| 模型方程 | $\mathbf{x} = C\mathbf{x} + \mathbf{d}$，即 $(I-C)\mathbf{x} = \mathbf{d}$ |
| Theorem 11 | 列和 < 1 ⇒ 唯一非负解（经济可行） |
| 级数公式 | $(I-C)^{-1} \approx I + C + C^2 + \cdots$（需求传导的数学化身） |
| $(I-C)^{-1}$ 的列 | 第 $j$ 部门需求 +1 ⇒ 各部门增产多少 |

> **本节精髓**：列昂惕夫模型展示了线性代数"**建模**"的一面：把经济系统抽象成矩阵方程，平衡条件写成 $(I-C)\mathbf{x} = \mathbf{d}$，求解即得生产计划。级数公式 (8) 更是把"经济直觉"（需求层层传导）翻译成了"矩阵级数"——这个 $(I - \text{算子})^{-1}$ 的展开思想在马尔可夫链（4.9）、差分方程（4.8）中还会反复出现。

---
---

## 2.7 计算机图形学的应用 Applications to Computer Graphics

**中文讲解**：现代计算机图形学（CAD、游戏、动画）的底层数学，就是**矩阵乘法**。物体 = 顶点列表，变换 = 矩阵，渲染 = 算矩阵乘积。本节看线性代数如何让图形"动"起来。

### 数据矩阵 Data Matrix

屏幕上的图形由顶点（vertices）决定。把 $n$ 个顶点的坐标存进**数据矩阵** $D$（2×n，第一行 x 坐标、第二行 y 坐标）。顶点被变换后，新图形的数据矩阵为 $AD$（矩阵乘法的**列定义**在此大显身手：$AD$ 的每列 = 对应顶点变换后的像）。

**例1（Example 1）**：大写字母 N（FIGURE 1）由 8 个顶点确定，其数据矩阵为：

$$D = \begin{bmatrix} 0 & .5 & .5 & 6 & 6 & 5.5 & 5.5 & 0 \\
0 & 0 & 6.42 & 0 & 8 & 8 & 1.58 & 8 \end{bmatrix}$$

> **几何观察**：左竖条宽 0.5（x 从 0 到 .5），右竖条宽 0.5（x 从 5.5 到 6）；对角线从左竖条上端 $(.5, 6.42)$ 连到右竖条下端 $(5.5, 1.58)$——注意 $6.42 + 1.58 = 8$，对角线端点关于中线 $y = 4$ 对称，这是精心设计的标准 N。

![FIGURE 1 原来的 N](img/fig27_1_letter_N.png)

*FIGURE 1 Regular N. 图1 原来的大写字母 N（8 个顶点）*

用剪切矩阵 $A = \begin{bmatrix} 1 & .25 \\
0 & 1 \end{bmatrix}$ 左乘 $D$，得剪切后的 N（FIGURE 2）：

$$AD = \begin{bmatrix} 1 & .25 \\
0 & 1 \end{bmatrix}\begin{bmatrix} 0 & .5 & .5 & 6 & 6 & 5.5 & 5.5 & 0 \\
0 & 0 & 6.42 & 0 & 8 & 8 & 1.58 & 8 \end{bmatrix} = \begin{bmatrix} 0 & .5 & 2.105 & 6 & 8 & 7.5 & 5.895 & 2 \\
0 & 0 & 6.42 & 0 & 8 & 8 & 1.58 & 8 \end{bmatrix}$$

**逐个验证**（剪切 = x 坐标加 0.25 倍 y 坐标）：$.5 + .25(6.42) = 2.105$ ✓；$6 + .25(8) = 8$ ✓；$5.5 + .25(8) = 7.5$ ✓；$5.5 + .25(1.58) = 5.895$ ✓；$0 + .25(8) = 2$ ✓

![FIGURE 2 剪切后的 N](img/fig27_2_slanted_N.png)

*FIGURE 2 Slanted N. 图2 剪切后的斜体 N*

**例2（Example 2，复合变换）**：FIGURE 2 的 N 有点太宽，补偿办法：再做一个**缩放变换** $S = \begin{bmatrix} .75 & 0 \\
0 & 1 \end{bmatrix}$（x 坐标乘 .75）。复合变换矩阵 = **乘积**（顺序：先剪切、后缩放 = $S$ 在左）：

$$SA = \begin{bmatrix} .75 & 0 \\
0 & 1 \end{bmatrix}\begin{bmatrix} 1 & .25 \\
0 & 1 \end{bmatrix} = \begin{bmatrix} .75 & .1875 \\
0 & 1 \end{bmatrix}$$

一个矩阵完成"先剪切、再缩放"两步（FIGURE 3）：

![FIGURE 3 复合变换后的 N](img/fig27_3_composite_N.png)

*FIGURE 3 Composite transformation of N. 图3 复合变换（剪切 + 缩放）后的 N*

> **硬件优势**：数据矩阵的乘法天然并行（每列独立），这正是 GPU 架构的数学模型——现代显卡每秒完成数十亿次矩阵乘法。

### 齐次坐标 Homogeneous Coordinates（核心技巧）

**中文讲解**：平移 $(x, y) \mapsto (x+h, y+k)$ **不是线性变换**（它把原点挪走了），无法用 2×2 矩阵表示。标准技巧：把平面点"升级"成三维**齐次坐标**：

$$(x, y) \;\longleftrightarrow\; (x, y, 1)$$

平移就能用 3×3 矩阵实现了：

$$\begin{bmatrix} 1 & 0 & h \\
0 & 1 & k \\
0 & 0 & 1 \end{bmatrix}\begin{bmatrix} x \\
y \\
1 \end{bmatrix} = \begin{bmatrix} x+h \\
y+k \\
1 \end{bmatrix}$$

**旋转 + 缩放 + 平移任意复合** = 对应 3×3 矩阵连乘，一次乘法完成整个动画帧。CAD 与游戏引擎的全部几何运算，都在这个 3×3 世界里完成。

### 三维旋转 3D Rotations

三维空间中，图形绕坐标轴旋转的标准矩阵（飞机姿态 yaw/pitch/roll 的数学本体）：

$$R_x(\alpha) = \begin{bmatrix} 1 & 0 & 0 \\
0 & \cos\alpha & -\sin\alpha \\
0 & \sin\alpha & \cos\alpha \end{bmatrix}$$

$$R_y(\beta) = \begin{bmatrix} \cos\beta & 0 & \sin\beta \\
0 & 1 & 0 \\
-\sin\beta & 0 & \cos\beta \end{bmatrix}$$

$$R_z(\gamma) = \begin{bmatrix} \cos\gamma & -\sin\gamma & 0 \\
\sin\gamma & \cos\gamma & 0 \\
0 & 0 & 1 \end{bmatrix}$$

任意朝向 = 三个基本旋转的乘积（欧拉角分解）。

### 透视投影 Perspective Projection

**中文讲解**：三维场景显示到二维屏幕 = 投影。近大远小的透视效果，用**齐次坐标 + 除以深度**实现。标准透视投影矩阵（投影中心在 z 轴上、屏幕为 xy 平面）作用后，齐次坐标的第 4 个分量变为 $1 - z/d$，最后"归一化"（各分量同除以它）：

$$\begin{bmatrix} 1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 \\
0 & 0 & -1/d & 1 \end{bmatrix}\begin{bmatrix} x \\
y \\
z \\
1 \end{bmatrix} = \begin{bmatrix} x \\
y \\
0 \\
1 - z/d \end{bmatrix} \;\longrightarrow\; \left(\frac{x}{1-z/d},\ \frac{y}{1-z/d},\ 0\right)$$

**直觉**：$z$ 越大（越远），除数 $1 - z/d$ 越大，投影点越靠近原点——**近大远小**自动实现。计算机图形学中"摄像机在内、投影在外"的奇观（例：1997 年北卡罗来纳大学虚拟现实项目，观察者移动、投影实时跟随），都是这条公式的工程实现。

### 2.7 小结

| 概念 | 一句话理解 |
|---|---|
| 数据矩阵 D | 图形顶点坐标存成 2×n 矩阵 |
| 图形变换 | 左乘矩阵；复合 = 矩阵连乘（顺序从右到左） |
| 齐次坐标 | (x,y) → (x,y,1)，把平移"伪装"成线性变换 |
| 3D 旋转 | 绕 x/y/z 轴三个基本矩阵，欧拉角的积木 |
| 透视投影 | 齐次坐标 + 除以深度 = 近大远小 |

> **本节精髓**：图形学把矩阵乘法的"**变换**"本质演绎到极致——矩阵是运动，乘法是复合，逆矩阵是"撤销"。下一节的子空间理论，则把矩阵的"**结构**"本质（列空间、零空间）推向中心舞台。

---
