$$
\newcommand{\bs}[1]{\boldsymbol{#1}}
\newcommand{\dfrac}[2]{\displaystyle\frac{\text{d}{#1}}{\text{d}{#2}}}
\newcommand{\ddfrac}[2]{\displaystyle\frac{\text{d}^2{#1}}{\text{d}{#2}^2}}
\newcommand{\dintegral}[3]{\int_{#1}^{#2} {#3} \mathrm{d}x}
$$

# 3. 連続体力学: 連続体の変形

連続体とは、運動する物体のモデルである. これまでに学んだ

- 質点: 物体を質量のみが定義された空間的な広がりを持たない点によりモデル化したもの. 物体の並進運動を記述することができる. 

- 剛体: 物体を質点の集合体としてモデル化したもの. 剛体内部の任意の二点間の距離が不変. 物体の並進運動に加え, 回転を記述することできる. 

の発展的なモデルとみなすことができ, 物体の運動に伴う **変形** を記述することができる. 本章の目標は, 変形を表す物理量として **ひずみ (strain)** を定義することである. なお, 本節ではテンソルを明確な定義なしに用いるが, その定義は次章において述べる. 本章においてはテンソルを行列と混同して構わない. 

## 変位

適当な座標系をもつ空間において, 境界を $S$ とする有限の大きさの物体 $V$ の内部に埋め込まれた点 $\bs{X}$ を考える. この物体が変形し, $V$, $S$, $\bs{X}$ がそれぞれ $v$, $s$, $\bs{x}$ に変化したとする ([@fig:displacement]). この時, 

$$
\bs{u}:=\bs{x}-\bs{X}
$${#eq:displacement}

を **変位 (displacement)** という. 変位 $\bs{u}$ は連続体内部の各点の移動量を表すが, 変位を眺めても物体 $V$ がどのように変形したかはよく分からない. 

以降の議論においては, 変形前の点 $\bs{X}$ を変形後の点 $\bs{x}$ にうつす写像 $\bs{x}=\bs{x}(\bs{X})$ が, 一対一であって (したがって逆写像 $\bs{x}=\bs{X}(\bs{x})$ が定義できる), かつ十分に滑らかである (必要なだけ微分できる) と仮定する.

![変形する連続体とその内部に埋め込まれた点](./figs/displacement.jpeg){#fig:displacement}

## 変形勾配テンソル

先の例においては変形する連続体内部の一点に注目したが、ここでは [@fig:strain] に示す連続体内部に埋め込まれた短いベクトル $\mathrm{d}\bs{X}:=\bs{Y}-\bs{X}$ を考え, 変形に伴いこれが $\mathrm{d}\bs{x}:=\bs{y}-\bs{x}$ となったとする. 

![変形する連続体とその内部に埋め込まれたベクトル](./figs/strain.jpeg){#fig:strain}

先に定義した写像 $\bs{x}$ を用いて $\mathrm{d}\bs{x}$ の一次近似を評価すると

$$
\begin{align}
\mathrm{d}x_i
&=x_i(\bs{Y})-x_i(\bs{X})\\
&=x_i(\bs{X}+\mathrm{d}\bs{X})-x_i(\bs{X})\\
&\simeq x_i(\bs{X})+\frac{\partial x_i(\bs{X})}{\partial X_j}\mathrm{d}X_j-x_i(\bs{X})\\
&=\frac{\partial x_i(\bs{X})}{\partial X_j}\mathrm{d}X_j
\end{align}
$${#eq:F}

と書ける (式 ([@eq:F]) において Einstein の総和規約を用いたことに注意せよ. すなわち, $j$ はダミーインデックスであり, $j$ についての和の記号を省略している). ここで, ([@eq:F]) 右辺は, $\partial x_i/\partial X_j$ を $ij$ 成分にもつ行列 $\mathsf{F}$ とベクトル $\mathrm{
d}\bs{X}$ の積 (の第 $i$ 成分) に他ならない. ここで, $\mathsf{F}$ を **変形勾配テンソル** と呼ぶ. 

> **問 1**: 変形勾配テンソルの「行列式」 (=$\text{det }\mathsf{F}$) が何を意味するか述べよ. また, $\text{det } \mathsf{F}>0$ すなわち $F$ が正定であることを示せ. ヒント: ベクトル $\mathrm{d}\bs{X}$, $\mathrm{d}\bs{Y}$, $\mathrm{d}\bs{Z}$ からなる平行六面体の体積が, $\mathsf{F}$ により特徴づけられる変形に伴ってどのように変化するかを考えよ.

> **問 2**: 変形勾配テンソルの成分を Kronecker のデルタと変位を用いて表わせ. 

## Green のひずみテンソルと Almansi のひずみテンソル

$$
\mathrm{d}\bs{x}=\mathsf{F}\mathrm{d}\bs{X}
$$

に現れる変形勾配 $\mathsf{F}$ は, ベクトルの変形させる作用のみならずこれを回転させるも含んでいる. そこで, 以降, $\mathsf{F}$ を回転を表す成分と変形を表す成分に分解することを目指す.

はじめに, $\mathsf{F}$ が正定値行列であるとき, $\mathsf{F}^T\mathsf{F}$ は正定値対称行列である. 

## 微小ひずみ

---

[../](../index.html)
