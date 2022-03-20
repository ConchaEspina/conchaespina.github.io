---
layout: post
excerpt_separator: <!--more-->
title: 《成分数据分析及其应用》笔记：第三章 成分数据缺失数据的处理
permalink: compositional-data-analysis/chap03
date: 2022-03-20
categories: review
tags: 
  - compositional data analysis
---

成分数据的零值：<span class="sidenote-number"></span><span class="sidenote">Martín-Fernández, J. A., Hron, K., Templ, M., Filzmoser, P., & Palarea-Albaladejo, J. (2012). Model-based replacement of rounded zeros in compositional data: classical and robust approaches. *Computational Statistics & Data Analysis*, *56*(9), 2688-2704.</span>

- **真实零值（essential zeros）：**成分数据不含某一成分，在这一成分的观测值为0，是不可剔除的。
- **近似零值（rounded zeros）：**成分数据含有某一成分，但是观测值非常小，低于仪器特定的观测限制，常指大小非常接近0的数据。

零值存在时，对成分数据做对数比变换会出现无穷大或无意义。

## 3.1   非参数方法

设$D$维成分数据$\mathbf{x} \in S^D$有$Z$个零值，替换后的新成分数据用$\mathbf{r} = [r_1, r_2, \cdots, r_D]$表示。

### 3.1.1   加法替换方法（Additive Replacement Strategy)

$$
r _j = \begin{cases}\dfrac{\delta(Z+1)(D-Z)}{D^2},\quad x_j = 0 \\ x_j - \dfrac{\delta(Z+1)Z}{D^2},\quad x_j>0\end{cases}
$$

其中，$j=1, 2, \cdots, D$，$\delta$是小于给定阈值的数。

加法替换方法具有以下性质：

1. 该方法

### 3.1.2   简单替换方法（Simple Replacement Strategy）





### 3.1.3   乘法替换方法（Multiplication Replacement Strategy）













