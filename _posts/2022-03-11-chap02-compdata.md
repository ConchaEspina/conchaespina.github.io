---
layout: post
excerpt_separator: <!--more-->
title: 《成分数据分析及其应用》笔记：第二章 成分数据简介
permalink: compositional-data-analysis/chap02
date: 2022-03-11
---

## 2.1   成分数据的定义及运算

***定义2.1.1（$D$维成分数据空间）***   称满足
$$
S^D = \left\{\mathbf{x} = [x_1, x_2, \cdots, x_D]: x_i > 0,\ i = 1, 2, \cdots, D; \sum_{i=1}^D x_i = c\right\}
$$
的空间为**$D$维成分数据空间**，这里$c$是任意常数。

$S^D$中的元素是$D$维行向量，但由于成分和为定值，所以是$D-1$维向量空间。

成分数据所在的空间$S^D$称为单形空间，单形空间的几何结构在研究成分数据时有着很重要的作用。

### 单形空间上的运算

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
\alpha \otimes \mathbf{x} = C[x_1^\alpha, x_2^\alpha, \cdots, x_D^\alpha]
$$




