$$
\newcommand{\bs}[1]{\boldsymbol{#1}}
$$

# Einstein の総和規約

本章では, 力学や電磁気学及びその数値解法を記述する際に便利な **Einstein の総和規約**について説明する. なお, 本稿では $\mathbb{R}^d$ のベクトルを太字のアルファベットで $\bs{a}$ などと書く. また, ベクトル $\bs{a}$ の第 $i$ 成分を細字のアルファベットを用いて $a_i \ (i=1,\cdots,d)$ などと書く. すなわち, 

$$
\bs{a}=\begin{pmatrix}a_1\\a_2\\a_3\end{pmatrix}
$$
である. ここに, $d$ は対象とする空間の次元である. 本章では $d=3$ の場合について考える. 


## ベクトルの内積

はじめに, ベクトルの内積を例にとり, Einstein の総和規約について説明する. ベクトル $\bs{a}, \bs{b} \in \mathbb{R}^3$ に対して, 

$$
\bs{a}\cdot\bs{b} := \sum_{i=1}^3 a_i b_i
$${#eq:naiseki_orig}

をベクトル $\bs{a}, \bs{b} \in \mathbb{R}^3$ の**内積 (inner product)** というのであった. $\bs{a}\cdot\bs{b}$ を $(\bs{a},\bs{b})$ と書くこともある. また, ベクトル $\bs{a}$ の転置 $\bs{a}^t$ を用いて $\bs{a}^t \bs{b}$ と書いても良い. さて, 次元 $d$ が前後の文脈などから明らかなときには式([@eq:naiseki_orig])において, $\sum_{i=1}^3$ を省略しても差し支えないように思える. そこで, 以降, 

$$
\bs{a}\cdot\bs{b} = a_i b_i
$${#eq:naiseki}

と書くことにする. つまり, **同じ項において添え字が二度現れる場合にはその添字について和をとると約束**する。式([@eq:naiseki])の $i$ のように, 一つの項において二度出現する添字を**ダミーインデックス (dummy index) **と呼ぶ. 一方, 一つの項において一度だけ出現する添字を**フリーインデックス (free index) **と呼ぶ.

> **例**: ダミーインデックスとフリーインデックス

> ベクトル $\bs{b}$ を単位ベクトル $\bs{a}$ が定める直線に正射影したベクトル $\bs{c}=\displaystyle\frac{\bs{a}\cdot\bs{b}}{|\bs{a}|^2} \bs{a}$ の第 $i$ 成分は

> $$
> c_i=\frac{a_jb_j}{a_ka_k}a_i \ \ \left(=\frac{a_1b_1+a_2b_2+a_3b_3}{a_1a_1+a_2a_2+a_3a_3}a_i\right)
> $${#eq:projection}

> と書ける. ここで, $j$ と $k$ はダミーインデックス, $i$ はフリーインデックスである. なお, 式([@eq:projection]) において, ダミーインデックスとして $j$ や $k$ を採用したことは本質的でなく, 例えば

> $$
> c_i=\frac{a_\ell b_\ell}{a_ma_m}a_i
> $$

> などと書いても構わない. ただし, (一つの項において) 添字が 3 回以上現れる場合についてはここでは定義しないので, 

> $$
> c_i=\frac{a_j b_j}{a_ja_j}a_i
> $$

> と書くことは許されない. 

このような記法は Einstein が初めて用いたことから, **Einstein の総和規約**と呼ばれる.

## Kronecker のデルタ $\delta_{ij}$ と交代記号 $e_{ijk}$

ここで、2 つの記号を導入する。1 つめは Kronecker のデルタであり、

$$
\delta_{ij}=\begin{cases} 1 & \text{if} \ \ i=j \\ 0 & \text{otherwise}\end{cases}
$${#eq:kronecker}

と定義される. すなわち, Kronecker のデルタは, その二つの添字が等しいときに 1 , そうでないときに 0 である二階のテンソルである.

> **問**: $\delta_{ii}$ の値は?

> **解**: 前述の Einstein の総和規約を適用すれば, $\delta_{ii} \ \left(=\delta_{11}+\delta_{22}+\delta_{33}\right)=3$ である. $\delta_{ii}=1$ としないように気をつけよう. 

ここで, Kronecker のデルタ $\delta_{ij}$ を成分に持つ行列 $\mathsf{E}$ を考える. 行列 $\mathsf{E}$ を成分表示すれば, 

$$
\mathsf{E}=
\begin{pmatrix}
1 & 0 & 0\\
0 & 1 & 0\\
0 & 0 & 1
\end{pmatrix}
$$

であるから, これは $\mathbb{R}^{3\times 3}$ の単位行列に他ならない. 単位行列 $\mathsf{E}$ と適当なベクトル $\bs{a}\in\mathbb{R}^3$ の掛け算が $\mathsf{E}\bs{a}=\bs{a}$ となることから, 

$$
\delta_{ij}a_j = a_i
$${#eq:deltaijaj}

となることが分かるであろう. すなわち, **$\delta_{ij}$ と $j$ を添字に持つベクトル (一般にはテンソル) の掛け算は, ベクトル (またはテンソル) の添字 $j$ を $i$ に置換したもの**となる. 同様に, $\delta_{ij}$ と $i$ を添字に持つベクトル (一般にはテンソル) の掛け算は, ベクトル (またはテンソル) の添字 $i$ を $j$ に置換したものとなる. 式 ([@eq:deltaijaj]) は今後しばしば使うことになる. 

2 つめは交代記号 (alternating symbol) と呼ばれ、

$$
e_{ijk}=
\begin{cases}
 1 & \text{if} \ (i,j,k)=(1,2,3) \ \text{or} \ (2,3,1) \ \text{or} \ (3,1,2)\\
-1 & \text{if} \ (i,j,k)=(2,1,3) \ \text{or} \ (3,2,1) \ \text{or} \ (1,3,2)\\
 0 & \text{otherwise}
\end{cases}
$${#eq:eijk}

と定義される. 交代記号は, 順列記号 (permutation symbol), Levi-Civita の記号などと呼ばれることもある. 

Kronecker のデルタと交代記号の間には以下の恒等式が成り立つ. 

$$
 e_{ijk}e_{ipq}=\delta_{jp}\delta_{kq}-\delta_{jq}\delta_{kp}
$${#eq:eedddd}

式 ([@eq:eedddd]) の証明は, 例えば $e_{ijk}$ が $\mathbb{R}^3$ の 3 つの標準基底のスカラー三重積 (=標準基底をならべた行列の行列式) であることと行列式に関するいくつかの公式を用いることでできるが, ここでは省略する. あるいは $j, k, p, q$ に 1 から 3 までの ($3^4=81$ 通りの) 数字の組を代入して確かめることもできる. 読者のみなさんには, いくつかの $j, k, p, q$ の組み合わせに対して式 ([@eq:eedddd]) が成り立つことを実際に確かめてみてほしい. **恒等式 ([@eq:eedddd]) は記憶する価値がある**. なぜなら, 後に見るように, この式さえ記憶していれば, ベクトル解析の種々の公式を導出できるからである.

## ベクトルの外積

次に, ベクトルの外積について考える. $\bs{a}, \bs{b} \in \mathbb{R}^3$ に対して, 以下で定義されるベクトル
$$
\bs{a}\times\bs{b} = 
\begin{pmatrix}
a_2b_3-a_3b_2\\
a_3b_1-a_1b_3\\
a_1b_2-a_2b_1
\end{pmatrix}
$${#eq:outerproduct}
をベクトル $\bs{a}$, $\bs{b}$ の**外積 (outer product, cross product)** という. 

さて, 先に導入した交代記号 $e_{ijk}$ を使うと, 式 ([@eq:outerproduct]) に示した外積 $\bs{a}\times\bs{b}$ の第 $i$ 成分は以下のように簡単に表現できる. 

$$
(\bs{a}\times\bs{b})_i = e_{ijk} a_j b_k
$${#eq:eijkajbk}

ここで, $i$ がフリーインデックス, $j$, $k$ がダミーインデックスであることに注意すれば式 ([@eq:eijkajbk]) が外積 (の成分) を表すことは明らかであろう. もし分からなければ, 例えば, 

$$
\left(\bs{a}\times\bs{b}\right)_1 
= e_{1jk} a_j b_k
　
\left( =\sum_{j=1}^3\sum_{k=1}^3e_{1jk}a_jb_k
=e_{111}a_1 b_1+e_{112}a_1 b_2+e_{113}a_1 b_3
+e_{121}a_2 b_1+e_{122}a_2 b_2+e_{123}a_2 b_3
+e_{131}a_3 b_1+e_{132}a_3 b_2+e_{133}a_3 b_3
\right)
$$

を計算してみると良いだろう. 

## 微分演算子 $\nabla$

微分演算子 $\nabla$ を, 偏微分作用素 $\partial/\partial x_i$ を成分にもつベクトルと解釈することにしよう. ただし, $\nabla$ は演算子 (あるいは作用素) なので, 作用する順序に注意する必要がある. 例えば, $\nabla\cdot\bs{a}$ と $\bs{a}\cdot\nabla$ の表すものは異なる (前者はスカラー, 後者は微分作用素). ここで, 例えばスカラー関数 $f$ の $x_i \ (i=1,2,3)$ による偏微分を, 
$$
f,_i:=\frac{\partial f}{\partial x_i}
$$
と書くことにしよう. この記法を用いれば, 勾配, 発散, 回転は以下のように簡潔に書くことができる. 

- スカラー値関数 $f$ の勾配 (の第 $i$ 成分):  $(\text{grad} f)_i:=(\nabla f)_i=f,_i$
- ベクトル値関数 $\bs{a}$ の発散: $\text{div} \ \bs{a}:=\nabla\cdot\bs{a}=a_{i,i}$
- ベクトル値関数$\bs{a}$の回転 (の第 $i$ 成分)： $(\text{rot} \ \bs{a})_i=(\text{curl} \ \bs{a})_i:=(\nabla\times\bs{a})_i=e_{ijk}a_{k,j}$

ここで, 回転が $e_{ijk}a_{j,k}$ ではないことに注意しよう。$k$, $j$ という順になることが分かりにくければ, $e_{ijk}\frac{\partial}{\partial x_j} a_k$ と書いてみると良いだろう. 

## ベクトル解析に現れる種々の公式の導出

ここで紹介した記法と恒等式 ([@eq:eedddd]) を用いると, ベクトル解析に現れる種々の公式を簡単に導出することができる. 

> **例**: $\bs{a}\times\left(\bs{b}\times\bs{c}\right)=(\bs{a}\cdot\bs{c})\bs{b}-(\bs{a}\cdot\bs{b})\bs{c}$

> 証明:
> $$
> \begin{align}
> \text{左辺の第} \ i \ \text{成分}
> &=e_{ijk}a_j(\bs{b}\times \bs{c})_k\\
> &=e_{ijk}a_je_{kpq}b_pc_q\\
> &=e_{kij}e_{kpq}a_jb_pc_q\\
> &=(\delta_{ip}\delta_{jq}-\delta_{iq}\delta_{jp})a_jb_pc_q\\
> &=a_jc_jb_i-a_jb_jc_i
> =\text{右辺の第} \ i \ \text{成分}
> \end{align}
> $$

> **例**: $(\bs{a}\times\bs{b})\cdot(\bs{c}\times\bs{d})=(\bs{a}\cdot\bs{c})(\bs{b}\cdot\bs{d})-(\bs{b}\cdot\bs{c})(\bs{a}\cdot\bs{d})$

> 証明: 
> $$
> \begin{align}
> \text{左辺}
> &=(e_{ijk}a_jb_k) (e_{ipq}c_pd_q)\\
> &=e_{ijk}e_{ipq} a_jb_kc_pd_q\\
> &=(\delta_{jp}\delta_{kq}-\delta_{jq}\delta_{kp})a_jb_kc_pd_q\\
> &=a_jb_kc_jd_k-a_jb_kc_kd_j
> =\text{右辺}
> \end{align}
> $$

> **例**:  $\nabla\times(\nabla\times\bs{a})=\nabla(\nabla\cdot \bs{a})-\nabla^2 \bs{a}$　(ここに, $\nabla^2=\partial^2/\partial x_i\partial x_i$ は Laplacian である. 注意: $i$ はダミーインデックス. )

> 証明: 
> $$
> \begin{align}
> \text{左辺の第} \ i \ \text{成分}
> &=e_{ijk}(\nabla\times\bs{a})_{k,j}\\
> &=e_{ijk}(e_{kpq}a_{q,p}),_{j}\\
> &=e_{kij}e_{kpq}a_{q,pj}\\
> &=(\delta_{ip}\delta_{jq}-\delta_{iq}\delta_{jp})a_{q,pj}\\
> &=a_{j,ij}-a_{i,jj}\\
> &=(a_{j,j})_{,i}-a_{i,jj}=\text{右辺の第} \ i \ \text{成分}
> \end{align}
> $$

## 積分公式

ベクトル解析におけるもっとも重要な定理として, Gauss の (発散) 定理は以下のようである. 

$$
\int\int\int_V \text{div} \bs{F} \text{d}V = \int\int_S \bs{F}\cdot\bs{n} \text{d}S
$${#eq:gauss_theorem}

ここに, $\bs{F}: \mathbb{R}^3\rightarrow\mathbb{R}^3$ はベクトル場, $V$ は $\mathbb{R}^3$ の有界領域, $S:=\partial V$ はその境界である. 単位法線ベクトル $\bs{n}$ は $V$ の外向きに定義した. ガウスの発散定理 ([@eq:gauss_theorem]) を添字を用いて表記すれば, 

$$
\int\int\int_V F_{i,i} \text{d}V = \int\int_S F_in_i \text{d}S
$${#eq:gauss_theorem_index}

である. Gauss の(発散) 定理と書いたのは, 式 ([@eq:gauss_theorem_index]) において, 発散を取る前の

$$
\int\int\int_V F_1 \text{d}V = \int\int_S F_1n_1 \text{d}S
$$

$$
\int\int\int_V F_2 \text{d}V = \int\int_S F_2n_2 \text{d}S
$$

$$
\int\int\int_V F_3 \text{d}V = \int\int_S F_3n_3 \text{d}S
$$

もそれぞれ成立することを強調するためである. なお, 一次元の Gauss の定理 (微分積分学の基本定理) は

$$
\int_a^b \frac{\text{d}f}{\text{d}x} \text{d}x = f(b)-f(a)
$$

と書ける. 右辺を $x=a \ (b)$ における「法線」$-1 \ (+1)$ を乗じた $f$ の「点積分」と見做せば, これが Gauss の定理に他ならないことが分かるであろう. 


## 演習問題

ここで導入した添字記法と Einstein の総和規約を用いることで, 力学や電磁気学の記述が簡単になり, 計算の見通しが格段に良くなるのでぜひ習得して欲しい. 以下に練習用の問題をいくつか挙げておく.

1. $\bs{a}$, $\bs{b}$, $\bs{c}\in\mathbb{R}^3\rightarrow\mathbb{R}^3$ をベクトル場とする. $(\bs{a}\times\bs{b})\times\bs{c}=(\bs{a}\cdot\bs{c})\bs{b}-(\bs{b}\cdot\bs{c})\bs{a}$ を示せ.

2. $\phi: \mathbb{R}^3\rightarrow\mathbb{R}$ をスカラー場, $\bs{u}$, $\bs{v}: \mathbb{R}^3\rightarrow\mathbb{R}^3$ をベクトル場とする. 以下の恒等式を証明せよ. 
   1. $\nabla\cdot(\phi \bs{v})=\nabla\phi\cdot\bs{v}+\phi \ \text{div} \ \bs{v}$
   2. $\nabla\times(\phi \bs{v})=\nabla\phi\times\bs{v}+\phi \ (\nabla\times\bs{v})$
   3. $\nabla\cdot(\bs{u}\times\bs{v})=\bs{v}\cdot(\nabla\times\bs{u})-\bs{u}\cdot(\nabla\times\bs{v})$
   4. $\nabla\times\nabla\phi = \bs{0}$
   5. $\nabla\cdot\nabla\times\bs{v} = 0$

3. $\bs{r}=(x_1, x_2, x_3)^t$, $r=|\bs{r}|$ (ただし $\bs{r}\neq\bs{0}$) とする. また, $\bs{a}$ を ($x_i$ に依らない) 定ベクトルとする. 以下を計算せよ. 
   1. $\nabla r$ \ \ (ヒント: $r^2=x_jx_j$ であり, また $x_{i,j}=\delta_{ij}$ である.)
   2. $\nabla\cdot\bs{r}$
   3. $\nabla\times\bs{r}$
   4. $\nabla \ (\bs{a}\cdot\bs{r})$
   5. $(\bs{a}\times\nabla)\cdot\bs{r}$
   6. $\nabla\times\left(\bs{a}\times\displaystyle\frac{\bs{r}}{r}\right)$
