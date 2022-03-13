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

对数比变换的特点：

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


