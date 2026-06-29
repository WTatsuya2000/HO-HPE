# src — 不変量バケツ方式の高速分類【Wilson / Polyakov ともに検証済】

旧 `src/` の分類（`is_circshift_equal` の総当たり＋予選決勝マージ）は L≥10 で O(n²) マージが破綻する。
ここでは **各ループに真の不変量を計算してバケツ分けし、同じバケツ内だけで完全同値判定**する方式で
高速化・破綻回避する。**Wilson は W⁰(n)、Polyakov は L₁⁰(Nt,n) を論文と一致**確認済み。

## ★ 推奨実装: 不変量バケツ方式（`invariant_*.jl`）

- `invariant_core.jl`: 真の不変量キー `invariant_key`（軸使用回数 / 隣接リンク straight·backtrack·turn /
  距離 k=2,3 のリンク相関 / 頂点訪問分布 / 無向リンク通過分布）と `classify_invariant`。
  - 分類は **バケツ内を `is_circshift_equal` の連結成分(union-find)で reduce（並列）→ 残った代表に大域連結成分マージ**。
  - **連結成分（貪欲でなく）にするのが要点**: `is_circshift_equal` は非推移的（特に Polyakov）なので、
    貪欲だと順序依存で過剰分割する。連結成分なら順序非依存・物理的に正しい分割（推移閉包）になる。
- `invariant_wilson.jl`: `classify_loop_inv(L, type_N)`（Wilson, 24軸置換）。
- `invariant_polyakov.jl`: `classify_loop_inv_poly(L, Nt)`（Polyakov, 6軸置換, **winding 対応生成で `-d` 込み**）。

各不変量は同値ループで必ず一致するので**クラスがバケツをまたいで分断されない**。バケツ内だけ
完全判定するので、旧 `classify_loop`（global 貪欲）と**厳密に同じ分類**を返しつつ、O(n²) を
O(n·L) ＋ Σ(小バケツ²) に落とす。バケツ間は非同値なのでマージ不要 → スレッド並列が素直に効く。

使い方:
```julia
include("src/invariant_wilson.jl")
classify_loop_inv(8, 4)   # Wilson L=8 type4 → 7本
```

## Polyakov（`invariant_polyakov.jl`）— 検証済み

**winding 対応生成**（巻き数 m=1: a,b,c は net 0、d は net=Nt、`-d` 許容）＋**真の群の正準形**で分類。
`net(d)=Nt` 固定なので巻き数 m≥2 は混入しない。
winding+1 集合の実効対称群は **{巡回シフト × a,b,c の6置換 × x,y,z の8反射}**（48·L 個, 閉じた真の群）で、
**軌道の辞書順最小 = 完全不変量** → これでグループ化するだけ（O(n·48L)・並列・総当たり不要）。

`L₁⁰(Nt, n) = 6·ΣM·D/LLL` が論文 Table 2 と一致:

| (Nt, n) | クラス | 計算 | 論文 |
|---|---|---|---|
| (4,6)  | 2     | 1,728      | 1,728 ✅ |
| (4,8)  | 24    | 45,792     | 45,792 ✅ |
| (4,10) | 500   | 645,120    | 645,120 ✅ |
| (6,10) | 64    | 716,544    | 716,544 ✅ |
| (6,12) | 2,093 | 35,891,712 | 35,891,712 ✅ |

使い方:
```julia
include("src/invariant_polyakov.jl")
classify_loop_inv_poly(8, 4)   # Nt=4, L=8 → 24クラス
```

## 生成ヘルパ（スクリプトが使用）

- `fast_wilson.jl`: 正準形+Booth の Wilson 高速分類 `classify_loop_fastw` と `wilson_mult`。
  `src/loop_core_wilson.jl` を include。→ `scripts/record_wilson.jl` が使用。
- `gen_canon_polyakov.jl`: winding 対応の正準生成 `gen_canon_poly`, `tosym_p`。
  `invariant_polyakov.jl` を include。→ `scripts/record_polyakov.jl` が使用。
