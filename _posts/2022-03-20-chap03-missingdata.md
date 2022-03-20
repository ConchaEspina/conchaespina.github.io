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

