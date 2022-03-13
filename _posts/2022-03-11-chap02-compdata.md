---
layout: post
excerpt_separator: <!--more-->
title: 《成分数据分析及其应用》笔记：第二章 成分数据简介
permalink: compositional-data-analysis/chap02
date: 2022-03-11
---

## 2.1   成分数据的定义及运算

### 2.1.1   成分数据空间

***定义2.1.1（$D$维成分数据空间）***   称满足
$$
S^D = \left\{\mathbf{x} = [x_1, x_2, \cdots, x_D]: x_i > 0,\ i = 1, 2, \cdots, D; \sum_{i=1}^D x_i = c\right\}
$$
的空间为$D$维成分数据空间，这里$c$是任意常数。

$S^D$中的元素是$D$维行向量，但由于成分和为定值，所以是$D-1$维向量空间。

成分数据所在的空间$S^D$称为单形空间，单形空间的几何结构在研究成分数据时有着很重要的作用。

### 2.1.2   $S^D$中的运算

***定义2.1.2（加法运算）***   对于任意的$\mathbf{x, y} \in S^D$，定义$\mathbf{x, y}$的加法运算$\oplus$为：
$$
\mathbf{x} \oplus \mathbf{y} = C [x_1 y_1, x_2 y_2, \cdots, x_D y_D]
$$
$C(\cdot)$表示对各个分量先除以各分量之和再乘以$c$，即：
$$
\mathbf{x} \oplus \mathbf{y} = \left[\frac{c x_1 y_1}{\sum_{i=1}^D x_i y_i}, \frac{c x_2 y_2}{\sum_{i=1}^D x_i y_i}, \cdots, \frac{c x_D y_D}{\sum_{i=1}^D x_i y_i}\right]
$$
由于加法运算的定义，可以得到单位元：$\mathbf{e} = \left[\dfrac{1}{D}, \dfrac{1}{D}, \cdots, \dfrac{1}{D}\right]$，且可见，对于任意的成分数据$\mathbf{x} \in S^D$，有$\mathbf{x} = \mathbf{x} \oplus \mathbf{e}$成立。

***定义2.1.3（数乘运算）***   对于任意的$\mathbf{x} \in S^D, \alpha \in \mathbb{R}$，定义数乘运算$\otimes$为：
$$
\begin{aligned} \alpha \otimes \mathbf{x} &= C[x_1^\alpha, x_2^\alpha, \cdots, x_D^\alpha] \\ &= \Bigg[\frac{c x_1^\alpha}{\sum_{i=1}^D x_i^\alpha}, \frac{c x_2^\alpha}{\sum_{i=1}^D x_i^\alpha}, \cdots, \frac{c x_D^\alpha}{\sum_{i=1}^D x_i^\alpha}\Bigg] \end{aligned}
$$
***定义2.1.4（减法运算）***   对于任意的$\mathbf{x, y} \in S^D$，定义$\mathbf{x, y}$的减法运算$\ominus$为：
$$
\mathbf{x} \ominus \mathbf{y} = \mathbf{x} \oplus ((-1) \otimes \mathbf{y}) = \mathbf{x} \oplus \mathbf{y}^{-1}
$$

### 2.1.3   $S^D$中运算的性质

***性质2.1.1***   $(S^D, \oplus)$是可交换群，且满足：

1. 交换律：$\mathbf{x} \oplus \mathbf{y} = \mathbf{y} \oplus \mathbf{x}$
2. 结合律：$(\mathbf{x} \oplus \mathbf{y}) \oplus \mathbf{z} = \mathbf{x} \oplus (\mathbf{y} \oplus \mathbf{z})$
3. 存在零元：$\mathbf{e} = C[1, 1, \cdots, 1]$
4. 对于$S^D$中的每一个元素$\mathbf{x}$都存在逆元：$\mathbf{x}^{-1} = C[x_1^{-1}, x_2^{-1}, \cdots, x_D^{-1}]$。且满足：$\mathbf{x} \oplus \mathbf{x}^{-1} = \mathbf{x}^{-1} \oplus \mathbf{x} = \mathbf{e}$

***性质2.1.2***   乘法运算的性质：

1. 结合律：$\alpha \otimes (\beta \otimes \mathbf{x}) = (\alpha \cdot \beta) \otimes \mathbf{x}$
2. 分配律1：$\alpha \otimes (\mathbf{x} \oplus \mathbf{y}) = (\alpha \otimes \mathbf{x}) \oplus (\alpha \otimes \mathbf{y})$
3. 分配律2：$(\alpha + \beta) \otimes \mathbf{x} = (\alpha \otimes \mathbf{x}) \oplus (\beta \otimes \mathbf{x})$
4. 存在单位元：$1 \otimes \mathbf{x} = \mathbf{x}$

### 2.1.4   Aitchison距离及其性质

***定义2.1.5（Aitchison距离）***   对于任意的$\mathbf{x, y} \in S^D$，定义$\mathbf{x, y}$之间的Aitchison距离为：
$$
d_a(\mathbf{x, y}) = \left[\frac{1}{D} \sum_{i<j} \left(\ln \frac{x_i}{x_j} - \ln \frac{y_i}{y_j}\right)^2 \right]^\frac{1}{2} = \left[\sum_{i=1}^D \left(\ln \frac{x_i}{g(\mathbf{x})} - \ln \frac{y_i}{g(\mathbf{y})}\right)^2 \right]^\frac{1}{2}
$$
其中：
$$
g(\mathbf{x}) = \sqrt[D]{\prod_{j=1}^D x_j},\ g(\mathbf{y}) = \sqrt[D]{\prod_{j=1}^D y_j}
$$
Aitchison距离反映的是<u>成分比例之间的差距</u>，而不是成分值的差距。

***定理2.1.1***   Aitchison距离具有加法不变性，即：
$$
d_a(\mathbf{x, y}) = d_a(\mathbf{z} \oplus \mathbf{x},\ \mathbf{z} \oplus \mathbf{y}),\quad \mathbf{x, y, z} \in S^D
$$
因此，有
$$
d_a(\mathbf{x, y}) = d_a(\mathbf{x} \oplus \mathbf{y}^{-1}, \mathbf{e})
$$
对于$S^D$中的任意$\mathbf{x, y}$，其Aitchison距离等于$\mathbf{x}$与$\mathbf{y}^{-1}$的和与单位元的Aitchison距离。

***定理2.1.2***   对任意的$\mathbf{x, y} \in S^D,\ \alpha \in \mathbb{R}$，有$d_a(\alpha \otimes \mathbf{x}, \alpha \otimes \mathbf{y}) = \mid\alpha\mid \cdot\ d_a(\mathbf{x,y})$。

***定理2.1.3***   在向量成比例时，Aitchison距离不变。

***定义2.1.6（成分向量的幂运算）***   定义成分数据的幂运算为：任意的$\mathbf{x} \in S^D,\ b>0$，有
$$
b^\mathbf{x} = C[b^{x_1}, b^{x_2}, \cdots, b^{x_D}]
$$
***定义2.1.7（成分向量的内积与模长）***   对于成分向量$\mathbf{x, y} \in S^D$，定义它们之间的内积与模长分别为：
$$
\langle\mathbf{x, y}\rangle_a = \frac{1}{D} \sum_{i<j} \ln \left(\frac{x_i}{x_j}\right) \ln \left(\frac{y_i}{y_j}\right) = \sum_{i=1}^D \ln \left(\frac{x_i}{g(\mathbf{x})}\right) \ln \left(\frac{y_i}{g(\mathbf{y})}\right) \\ \|\mathbf{x}\|^2_a = \langle\mathbf{x, x} \rangle_a,\ d_a(\mathbf{x, y}) = \|\mathbf{x} \ominus \mathbf{y}\|_a
$$
***定义2.1.8（正交）***   若成分向量$\mathbf{x, y} \in S^D$，且$\langle\mathbf{x,y}\rangle_a = 0$，则称$\mathbf{x, y}$正交。

## 2.2   成分数据的变换及其性质

### 2.2.1   非对称对数比（alr）变换

***定义2.2.1（非对称对数比变换）***   设$\mathbf{x} = [x_1, x_2, \cdots, x_D]$是成分向量，令
$$
y_i = \ln \frac{x_i}{x_D},\quad i=1,2,\cdots,D-1
$$
称此变换为非对称对数比变换。

如果变换后的$\mathbf{y} = (y_1, y_2, \cdots, y_{D-1})$服从多元正态分布，则称$\mathbf{x} = [x_1, x_2, \cdots, x_D]$服从加法逻辑正态分布。

非对称对数比变换的特点：

- 优势
  1. 可以克服“定和限制”
  2. 可以部分消除成分间的完全相关性，以便运用最小二乘法
  3. 变换后的数据在$(-\infty, +\infty)$内取值，便于进行模型选择
- 缺陷
  1. 经对数比变换得到的新变量完全不能和原始变量相对应，使得模型的解释性不强

在非对称对数比变换中，如果用$y_i$来表示$x_i$，有如下表达式：
$$
\begin{cases}x_i = \dfrac{e^{y_i}}{1 + \sum\limits_{j=1}^{D-1} e^{y_j}},\quad i=1,2,\cdots,D-1 \\ x_D = \dfrac{1}{1 + \sum\limits_{j=1}^{D-1} e^{y_j}}\end{cases}
$$
假设$\mathbf{y} = (y_1, y_2, \cdots, y_{D-1})\ \sim\ N(\mathbf{\mu, \Sigma})$，则$\mathbf{y}$的密度函数为
$$
\left(\frac{1}{\sqrt{2\pi}}\right)^{D-1} |\mathbf{\Sigma}|^{-1/2} \exp\left\{-\frac{1}{2}(\mathbf{y-\mu})^\top\mathbf{\Sigma}^{-1}(\mathbf{y-\mu})\right\}
$$
$\mathbf{x}$的密度函数为
$$
\left(\frac{1}{\sqrt{2\pi}}\right)^{D-1} |\mathbf{\Sigma}|^{-1/2} \prod_{j=1}^D x_j^{-1} \exp\left\{-\frac{1}{2}Q\right\}
$$
其中
$$
Q = (F \ln \mathbf{x - \mu})^\top |\mathbf\Sigma|^{-1/2} (F \ln \mathbf{x - \mu}),\quad F = (-\mathbf{1}, \mathbf{I}_{D-1})
$$
***定义2.2.2***   假定成分向量$\mathbf{x} = [x_1, x_2, \cdots, x_D]$的函数：
$$
\mathbf{y} = (y_1, y_2, \cdots, y_D),\ y_i = \ln \left(\frac{x_i}{1 - \sum_{j=1}^i x_j}\right),\quad i=1,2,\cdots,D
$$
服从多元正态分布$N(\mathbf{\mu,\Sigma})$，则称成分$\mathbf{x}$服从乘法逻辑正态分布。

由上式可知，如果用$y_i$来表示$x_i$，有如下表达式<span class="sidenote-number"></span> 
<span class="sidenote">可用数学归纳法与定和限制$\sum_{i=1}^D x_i = 1$加以证明。原书此处$x_D$的表达式似乎有误。</span>：
$$
\begin{cases}x_i = \dfrac{e^{y_i}}{\prod\limits_{j=1}^i (1 + e^{y_j})},\quad i=1,2,\cdots,D-1 \\ x_D = \dfrac{1}{\prod\limits_{j=1}^{D-1} (1 + e^{y_j})}\end{cases}
$$

### 2.2.2   对称对数比（clr）变换

***定义2.2.3（对称对数比变换）***   设成分向量$\mathbf{x} = [x_1, x_2, \cdots, x_D]$，做变换：
$$
y_i = \ln \frac{x_i}{\sqrt[D]{\prod\limits_{j=1}^D x_j}},\quad i=1,2,\cdots,D
$$
这种变换称为对称对数比变换。

如果用$y_i$来表示$x_i$，有如下表达式：
$$
x_i = \frac{e^{y_i}}{\sum\limits_{j=1}^D e^{y_j}},\quad i=1,2,\cdots,D
$$
对称对数比变换的特点：

- 优势：
  1. 相对于alr变换而言，有效解决了变换后分量不能与原有分量一一对应的不足，增强了解释性
- 缺陷：
  1. 变换后的各个分量之和为零
  2. 与之相对应的协方差阵是奇异的

### 2.2.3   等距对数比（ilr）变换

***定义2.2.4（等距对数比变换）***   设$\mathbf{x} = [x_1, x_2, \cdots, x_D]$是成分向量，做如下变换：
$$
z_i = \sqrt{\frac{D-i}{D-i+1}} \ln \frac{\sqrt[D-i]{\prod\limits_{l=i+1}^D x_l}}{x_i},\quad i=1,2,\cdots,D-1
$$
则可将$D$维向量$\mathbf{x} = [x_1, x_2, \cdots, x_D]$等距对数比变换成$D-1$维向量$\mathbf{z} = (z_1, z_2, \cdots, z_{D-1})$。

等距变换使得成分数据从单形空间映射到欧氏空间中，并且可以保证单形空间中的两组成分数据之间的Aitchison距离与通过等距变换后的两组数据的欧氏距离相等。

***定理2.2.1***   令$\mathbf{u}_i \in \mathbb{R}^D,\ i=1,2,\cdots,D-1$，且
$$
\mathbf{u}_i = \sqrt{\frac{i}{i+1}} \Big(\underbrace{\frac{1}{i}, \cdots, \frac{1}{i}}_{i个}, -1, 0, \cdots, 0\Big)
$$
则$\mathbf{u}_i$在欧氏空间中是正交的，且组成了$D-1$维子空间。

***定理2.2.2***   令$\mathbf{e}_i,\ i=1,2,\cdots,D-1$是$S^D$空间中的向量，且<span class="sidenote-number"></span> 
<span class="sidenote">原书$\mathbf{e}_i$的表达式似乎有误。</span>
$$
\mathbf{e}_i = C(\exp(\mathbf{u}_i)) = C\Big(\exp\Big(\underbrace{\sqrt{\frac{1}{i(i+1)}}, \cdots, \sqrt{\frac{1}{i(i+1)}}}_{i个}, -\sqrt{\frac{i}{i+1}}, 0, \cdots, 0\Big)\Big)
$$
则成分向量$\mathbf{e}_i,\ i=1,2,\cdots,D-1$在$S^D$空间以Aitchison内积尺度正交，组成了$S^D$中的正交基。

**alr变换，clr变换和ilr变换的性质及相互之间的关系如下所述：**

***性质2.2.1***   ilr变换把$S^D$空间的向量映射到$\mathbb{R}^{D-1}$空间中，对于任意的$\mathbf{x} \in S^D$都有：
$$
\mathbf{y} = \mathrm{ilr}(\mathbf{x}) = (\langle\mathbf{x, e}_1\rangle_a, \langle\mathbf{x, e}_2\rangle_a, \cdots, \langle\mathbf{x, e}_{D-1}\rangle_a)
$$
其中：
$$
\mathbf{e}_i = C(\exp(\mathbf{u}_i)) = C\Big(\exp\Big(\underbrace{\sqrt{\frac{1}{i(i+1)}}, \cdots, \sqrt{\frac{1}{i(i+1)}}}_{i个}, -\sqrt{\frac{i}{i+1}}, 0, \cdots, 0\Big)\Big)，\quad i=1,2,\cdots,D-1
$$
***性质2.2.2***   对于任意的$\mathbf{x}_1, \mathbf{x}_2 \in S^D$及任意的$\alpha \in \mathbb{R}$，且$\mathbf{y}_1 = \mathrm{ilr}(\mathbf{x}_1), \mathbf{y}_2 = \mathrm{ilr}(\mathbf{x}_2)$，有如下结论成立：
$$
\mathrm{ilr}(\mathbf{x}_1 \oplus \mathbf{x}_2) = \mathbf{y}_1 + \mathbf{y}_2 \\ \mathrm{ilr}(\alpha \otimes \mathbf{x}_i) = \alpha \mathbf{y}_i,\quad i=1,2
$$
***性质2.2.3***   对于任意的$\alpha, \beta \in \mathbb{R}$及$\mathbf{x}_1, \mathbf{x}_2 \in S^D, \mathbf{e}_i \in S^D$和$\mathbf{u}_i \in \mathbb{R}^{D-1}$如定理2.2.1和定理2.2.2所述，clr变换有如下结论成立：
$$
\mathrm{clr}(\mathbf{e}_i) = \mathbf{u}_i \\ \mathrm{clr}((\alpha \otimes \mathbf{x}_1) \oplus (\beta \otimes \mathbf{x}_2)) = \alpha \mathrm{clr}(\mathbf{x}_1) + \beta \mathrm{clr}(\mathbf{x}_2)
$$
***定理2.2.3***   对任意的$\mathbf{x} \in S^D$，考虑正交基$\mathbf{e}_i = C(\exp(\mathbf{u}_i)),\ i=1,2,\cdots,D-1$，其中$\mathbf{e}_i \in S^D, \mathbf{u}_i \in \mathbb{R}^{D-1}$如定理2.2.1和定理2.2.2所述。则对于alr变换、clr变换、ilr变换有如下关系：
$$
\mathrm{clr}(\mathbf{x}) = \mathrm{ilr}(\mathbf{x}) \mathbf{U},\ \mathrm{alr}(\mathbf{x}) = \mathrm{ilr}(\mathbf{x}) \mathbf{UF},\ \mathrm{ilr}(\mathbf{x}) = \mathrm{clr}(\mathbf{x}) \mathbf{U}^\top = \mathrm{alr}(\mathbf{x}) \mathbf{AU}^\top
$$
其中：

$\mathbf{U}$是$(D-1) \times D$维矩阵，$\mathbf{u}_i$是它的行向量

$\mathbf{F} = [\mathbf{I}_{D-1},$ $\mathbf{1}_{D-1}^\top]^\top$

$\mathbf{A} = \dfrac{1}{D} \begin{bmatrix}D-1 & -1 & -1 & \cdots & -1 & -1 \\ -1 & D-1 & -1 & \cdots & -1 & -1 \\ -1 & -1 & D-1 & \cdots & -1 & -1 \\ \vdots & \vdots & \vdots & & \vdots & \vdots \\ -1 & -1 & -1 & \cdots & D-1 & -1\end{bmatrix}_{(D-1) \times D}$

### 2.2.4   球坐标变换

球坐标变换可以解决对数比变换中零成分变换的困境，<span class="sidenote-number"></span><span class="sidenote">王惠文,刘强. (2002). 成分数据预测模型及其在中国产业结构趋势分析中的应用. *中外管理导报*(05), 27-29.</span>提供了一种变换方法，但它并没有具体的理论指导。

***定义2.2.5*** 对于成分向量$\mathbf{x} = [x_1, x_2, \cdots, x_D]$，满足定和限制$\sum\limits_{i=1}^D x_i = 1$，首先对成分向量的各分量开根号得：$y_i = \sqrt{x_i}$，此时有$\sum\limits_{i=1}^D y_i^2 = 1$，那么向量$\mathbf{y} = (y_1, y_2, \cdots, y_D)$可以看成超球面上的点，球坐标变换把$D$维向量$\mathbf{y} = (y_1, y_2, \cdots, y_D) \in \mathbb{R}^D$映射到超球面$(r, \theta_2, \cdots, \theta_D) \in \Theta^D$，这里$r^2 = \|\mathbf{y}\|^2 = 1$，且有
$$
\begin{cases}y_1 = \sin\theta_2 \sin\theta_3 \sin\theta_4 \cdots \sin\theta_D \\ y_2 = \cos\theta_2 \sin\theta_3 \sin\theta_4 \cdots \sin\theta_D \\ y_3 = \cos\theta_3 \sin\theta_4 \cdots \sin\theta_D \\ \cdots \\ y_{D-2} = \cos\theta_{D-2} \sin\theta_{D-1} \sin\theta_D \\ y_{D-1} = \cos\theta_{D-1} \sin\theta_D \\ y_D = \cos\theta_D\end{cases}
$$




做反变换可得：
$$
\begin{cases}\theta_D = \arccos y_D \\ \theta_{D-1} = \arccos\left(\dfrac{y_{D-1}}{\sin\theta_D}\right) \\ \theta_{D-2} = \arccos\left(\dfrac{y_{D-2}}{\sin\theta_D \sin\theta_{D-1}}\right) \\ \cdots \\ \theta_2 = \arccos\left(\dfrac{y_2}{\sin\theta_D \sin\theta_{D-1} \cdots \sin\theta_3}\right)\end{cases}
$$

## 2.3   成分数据分析的困难

### 2.3.1   高维数据困难

对于高维困难，最常用的一种解决方式是投影。关于高维成分数据的投影，有以下三个特点：

1. 部分分析方法：投影方法将注意力集中于子成分的选择，而不是整个成分。
2. 图像分析不一定能解决高维成分数据的问题，因为不能保证单形空间上的图形与一般实数空间上的图形有相同的解释。一个可供选择的单形空间中的投影是：在所有成分中选择一个方差变化最大的成分，其余$D-1$个成分分别与该成分进行画图。
3. 借助多元统计分析技术可以实现成分数据作为一个整体的分析，而不必求助于部分分析。

### 2.3.2   缺乏可解释的方差结构

***定义2.3.1（成分数据的协方差结构）***   成分向量$\mathbf{x} = [x_1, x_2, \cdots, x_D]$的原始协方差结构是指所有的
$$
k_{ij} = \mathrm{Cov}(x_i, x_j),\quad i,j=1,2,\cdots,D
$$
组成的集，$D \times D$原始协方差矩阵为：
$$
K= [k_{ij}:\ i,j=1,2,\cdots,D]
$$
原始相关系数为：
$$
\rho_{ij} = \mathrm{Corr}(x_i, x_j) = \frac{k_{ij}}{\sqrt{k_{ii} k_{jj}}},\quad i=1,\cdots,D;\ j=i+1,\cdots,D
$$
相关系数矩阵由$\dfrac{D(D-1)}{2}$个相关系数决定。

缺乏可解释的协方差结构可以主要表述为以下几个方面：

1. 负偏困难
2. 子成分困难
3. 基困难
4. 零相关困难







