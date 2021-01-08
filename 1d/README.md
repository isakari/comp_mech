$$
\newcommand{\bs}[1]{\boldsymbol{#1}}
\newcommand{\dfrac}[2]{\displaystyle\frac{\text{d}{#1}}{\text{d}{#2}}}
\newcommand{\ddfrac}[2]{\displaystyle\frac{\text{d}^2{#1}}{\text{d}{#2}^2}}
$$

# 1. イントロダクション: 1次元問題に対する有限要素法

本章では, 有限要素法 (finite element method; FEM) の概念を説明することを目的とし, $x\in[0,1]$ で定義された二階常微分方程式

$$
\frac{\text{d}^2 u(x)}{\text{d} x^2} = -f(x)
$${#eq:governing_eq}

を対象とした有限要素法について述べる. ここに, $u(x)$ は未知のスカラー関数, $f$ は既知のスカラー関数である. 

## 境界値問題の定式化

はじめに, 弦のたわみが式 ([@eq:governing_eq]) で記述される物理現象 (の一つ) であることを確認しよう. 今, [@fig:setting] のように, 両端が固定された長さ 1 の弦に, 下向き分布荷重 $q(x)$ が作用した時, 点 $x$ におけるたわみ $v(x)$ を求める問題を考えよう. なお, 荷重は下向きに作用するので, 縦軸は下向きを正にとる. 

![分布荷重を受けてたわむ弦](./figs/setting.jpeg){#fig:setting}

微小領域 $[x,x+\text{d}x]$ における弦に作用する力の釣り合いを考える ([@fig:equiv]). 

![微小領域における力のつりあい](./figs/equiv.jpeg){#fig:equiv}

弦に働く張力を $T$ と書き, 変形が微小であることを仮定すれば $\sin\theta\simeq\tan\theta=\displaystyle\frac{\text{d}u(x)}{\text{d}x}$ だから, 左端 $x$ に作用する張力の鉛直方向成分は $-T\sin\theta\simeq-T\dfrac{u(x)}{x}$ と書ける. 同様にして, 右端 $x+\text{d}x$ に作用する力の鉛直方向成分は

$$
T\dfrac{u(x+\text{d}x)}{x}
\simeq T\left(\dfrac{u(x)}{x}+\ddfrac{u(x)}{x}\text{d}x\right)
$$

である. さらに, 与えられた分布荷重が微小領域において一定値 $q(x)$ を取るとすれば, 力の釣り合いは

$$
-T\dfrac{u(x)}{x}+T\left(\dfrac{u(x)}{x}+\ddfrac{u(x)}{x}\text{d}x\right)+q(x)\text{d}x=0
$$

と書ける. これを整理して, $f(x)=q(x)/T$ とおけば式 ([@eq:governing_eq]) を得る. 

今, 弦は両端で固定されているので, 境界条件は 

$$
u(0)=u(1)=0
$${#eq:bc}

である. 微分方程式 ([@eq:governing_eq]) と境界条件 ([@eq:bc]) の組みを**境界値問題 (boudnary value problem; BVP)** と呼ぶ. 有限要素法は境界値問題の数値解法のひとつである. なお, 式 ([@eq:bc]) のように, 未知関数そのものの境界値を与える境界条件を Dirichlet 境界条件という. 一方で, 未知関数の境界値での微分値を与える境界条件は Neumann 境界条件, Dirichlet 条件と Neumann 条件の線形結合を与える境界条件は Robin 境界条件と呼ばれる. 

> **問**: $f(x)=1$ のとき, 上の Dirichlet 境界値問題の解を解析的に (=手計算で積分し) 求めよ. 特に, $x=1/2$ における $u$ の値はいくらか?

> **解**: $u(x)=-\displaystyle\frac{1}{2}x^2+\displaystyle\frac{1}{2}x$, \ $u(1/2)=1/8$

## 重み付き残差法 (weighted residual method; WRM) と Galerkin 法

先の「問」においては, 解析的に微分方程式を解くことができた. しかし, 一般にはこれは難しい. そこで $u$ をなんらかの意味で近似する関数 $\tilde{u}$ を求めることが必要となる. 以下, $\tilde{u}$ を求める方法の一つである重み付き残差法について解説する. 

まず,  $\tilde{u}$ が $u$ が満たすべき境界条件 ([@eq:bc]) を満足することを要請しよう. すると, $g_i(0)=g_i(1)=0$ を満たす関数 $g_i(x)$ (例えば$　g_i(x)=x^{i+1}(1-x)$ など) を用いて構成される

$$
\sum_{i=0}^{n-1} a_i g_i(x)
$${#eq:aigi}

は $\tilde{u}$ の候補となる. ここに, $a_i\in\mathbb{R}$ は未知の係数である. 未知係数 $a_i$ を求めるため, ここでは**重み付き残差法**を採用する. その他の方法としては, 微分方程式の定義された領域内部の有限個の点 $x_i \ (i=0,\cdots, n-1)$ において微分方程式を満たすように $a_i$ を決める**選点法**が代表的である.

適当な重み関数 $v_i(x)$ を $i=0, \cdots, n-1$ の $n$ 種類準備し, 微分方程式 ([@eq:governing_eq]) の両辺にこれを乗じ, 微分方程式の定義された領域 $[0, 1]$ において積分すると以下を得る. 

$$
\int_0^1 v_i \left( \ddfrac{\tilde{u}(x)}{x}+f(x) \right) \text{d}x=0
$${#eq:wre}

式 ([@eq:wre]) を**重み付き残差式**と呼ぶ. これを満たすように未知係数 $a_i$ を定める方法を**重み付き残差法**である.

> **例**: $f(x)=1$の場合を考える. $n=1$, $g_i(x)=\sin i \pi x$ とし, 

> - **v_i=1**

> - **v_i=g_i**

## 有限要素法

---

[../](../index.html)
