# HQQCD-HPE — Heavy-quark determinant HPE measurement (Nt = 4)

格子ゲージ理論における重い Wilson クォークの行列式 `ln det M(κ)` を
**hopping parameter expansion (HPE)** で展開し、各次数の κⁿ 係数を配位ごとに測定する
self-contained パッケージです。

```
ln det M(κ) = N_site · Σ_n D_n κⁿ
```
`D_n` は長さ n の閉ループ和で、**Wilson 型**（時間方向に巻かないループ）と
**Polyakov 型**（時間方向に m 回巻くループ、巻き数 m 別）に分かれます。
本パッケージは **LO〜N3LO（長さ 4 / 6 / 8 / 10）** を測定します
（LO=plaquette / Polyakov loop は別枠で出力）。

| 次数 | 長さ | Wilson | Polyakov |
|---|---|---|---|
| NLO  | 6  | W6  | P6 (w1) |
| NNLO | 8  | W8  | P8 (w1, w2) |
| N3LO | 10 | W10 | P10 (w1, w2) |

---

## 測定方法（オーダー別 mixed）

各次数で最適／正準な手法を使い分けます（既定ビルド）。**どの手法でも数値は完全一致**します。

| 次数 | 長さ | 手法 | 実装 |
|---|---|---|---|
| NLO  | 6  | **直接積**（ループを直接トレース） | `util/rectangle.hpp` `chair.hpp` `crown.hpp` `polyakov_loop_bend.hpp` |
| NNLO | 8  | **staple**（SingleStaple を combine） | `NNLO_human/nnlo_w.cpp` `NNLO_W` / `NNLO_human/nnlo_p4.cpp` `NNLO_P4L1,L2` |
| N3LO | 10 | **trie**（meet-in-the-middle, 前後半共有） | `src/hpe_meas.cpp` |

- NLO の直接積→単一値への合成係数は trie に対し**機械精度で校正**済:
  - `W6 = 1152·rec_s + 576·rec_t1 + 576·rec_t2 + 1152·cha_s + 2304·cha_t1 + 1152·cha_t2 + 384·cro_s + 1152·cro_t`
  - `P6 = 1152·Polyakov_bend1 + 576·Polyakov_bend2`
- NNLO の staple（W8 / P8w1 / P8w2）は trie と直接一致（フィット不要）。
- **`-DLOW_ORDER_TRIE` で全セクターを trie に切り替え可**（数値同一・クロスチェック用、下記 `make trie`）。

---

## ディレクトリ構成

```
HPE/
├── src/              measure_n3lo.cpp (本番ドライバ) , hpe_meas.cpp (trie エンジン)
├── NNLO_human/       人手 staple 実装
│                       nnlo.cpp = single_staple / local_elem / nnlo_w / nnlo_p4 を束ねる
│                       nnlo_safe.cpp = W8traj/W10traj/MakeTrajList/cyclic_symmetry
│                       calc_hpe.cpp
├── trajectory_list/  ループ方向列リスト nnlo_traj_list.cpp , n3lo_traj_list.cpp
│   └── tables/        係数テーブル *.dat と生成器 *_table.cpp（W6/8/10, P6/8/10, P12, poly_68）
├── util/             rectangle/chair/crown/polyakov_loop(_bend)/staple_pol/gauge_io_mult/
│                       gauge_update_routines_pol/stopwatch
├── loop_classification/  長さ L ごとの Wilson/Polyakov 閉ループ列挙・分類・γトレース (Julia)
│                       src/ src_v2/ scripts/ Project.toml,Manifest.toml
├── TML/              格子ライブラリ (mdp/v2606, CART_HALO=2)
├── jobs/             Slurm 投入スクリプト (sub_openmp.sh, sub_mpi.sh)
├── input
├── makefile
└── conf/ data/ log/  実行時生成
```

> **実行時のテーブル参照**: `a.out` は係数テーブルを `trajectory_list/tables/` から **CWD 相対**で読みます。
> 必ずパッケージのルートから起動してください（`jobs/*.sh` は `cd $SLURM_SUBMIT_DIR` 済み）。

## 必要環境
- Intel oneAPI（Intel MPI + OpenMP）— `mpiicpx`
- ビルド／実行の前に: `bash -l` の上で `module load intel/2025.3.1`
  （Slurm ジョブは `#!/bin/bash -l` でログインシェルにし、ジョブ内でも `module load` する）

## ビルド
実行ファイル（`*.out`）は**同梱しない**。下記でその場でビルドする（`make clean` で消える build 生成物）。

```sh
make            # a.out          : mixed, 単一ランク + OpenMP (div=1,1,1,1)   ← 推奨
make mpi        # a_mpi.out      : mixed, 8ランク MPI        (div=1,2,2,2)
make trie       # a_trie.out     : 全 trie (-DLOW_ORDER_TRIE), 単一ランク
make mpi_trie   # a_mpi_trie.out : 全 trie, 8ランク
make clean      # *.out を削除
```

各バイナリの違いは**並列化の方式だけ**で、測る中身（mixed 測定）も出力数値も同一:

| バイナリ | make | 並列 | 用途 |
|---|---|---|---|
| `a.out`          | `make`      | 単一ランク + OpenMP スレッド          | 小さめ格子・手軽に（推奨） |
| `a_mpi.out`      | `make mpi`  | **8ランク MPI**（空間を 2×2×2 分割）+ OpenMP | 大きい格子をノード分散（本番 48³ はこちら） |
| `a_trie.out`     | `make trie` | 単一ランク（全セクター trie）          | クロスチェック／性能比較 |
| `a_mpi_trie.out` | `make mpi_trie` | 8ランク MPI（全 trie）             | 同上・MPI |

- `-DLOW_ORDER_TRIE`（`trie`/`mpi_trie`）= 全セクターを trie 経路で測る。mixed と**数値完全一致**（検証用）。
- `-DDIV1=2 -DDIV2=2 -DDIV3=2`（`mpi`/`mpi_trie`）= 空間3方向を2分割 → 計8ランク。時間方向は分割しない。
- 測定は OpenMP スレッド並列で安全（数値は 13 桁再現）。
- 配位更新（heat bath / over-relaxation）は単一プロセス想定。MPI 版は領域分割で並列。

## input の形式（8 行・順番固定。各行の先頭トークンが値、残りはコメント）
```
4 48 48 48   lattice size (Nt Nx Ny Nz)
5.8899       beta
0.0014       lambda   (LO Polyakov 結合; 内部で lambda/=beta)
conf         configuration directory
data         data storing directory
5000         number of measurements   (=0 で測定のみ・保存なし)
10           sweeps between measurements
10           HPE max loop length (10=N3LO, 偶数)
```
- 8 行目 = 取り込む最大ループ長。`6`/`8`/`10` を指定可。`12`(N4LO) は `CART_HALO=3` で TML 再ビルドが必要。

## 実行（計算ノードへ Slurm 投入。ログインノードで直接走らせない）
```sh
sbatch jobs/sub_openmp.sh    # 単一ランク OpenMP (a.out, input)
sbatch jobs/sub_mpi.sh       # 8ランク MPI      (a_mpi.out, input)  ※先に make mpi
```
- cold start から走らせると、初回（更新前）の測定行が canonical 値（下記）になる。
  既存配位だけ測りたい場合は input 6 行目（測定回数）を `0`（測定のみ・更新/保存なし）にする。
- MPI 版の注意: 分割する空間方向の局所サブ格子 ≥ 2·halo（N3LO は reach=2 → 局所 ≥ 4）。
  起動時にガードが検査し、足りなければ FATAL で停止。

## 出力（`data/hpe_<param>_<counter>`、1 配位 1 行）
```
# plaq_s plaq_t Re_pol Im_pol  W6 W8 W10  Re_L1n6 Im_L1n6  Re_L1n8 Im_L1n8  Re_L2n8 Im_L2n8  Re_L1n10 Im_L1n10  Re_L2n10 Im_L2n10
```
- Wilson（W6/8/10）は実 1 列、Polyakov（`L<m>n<n>` = 巻き数 m・長さ n）は Re/Im 2 列。
- 列の順序・本数は手法（mixed / 全 trie）に依らず**同一**。

## 検証（再現性）
- **cold (U=1) 値 = 論文値**（本パッケージのビルドで確認済）:
  `W6=8448, W8=245952, W10=7372800, L1n6=1728, L1n8=45792, L2n8=−384, L1n10=645120, L2n10=−64512`。
- **mixed == 全 trie**: 同一の熱化配位で全列が機械精度一致（cold・熱化配位とも確認）。
- **単一 ↔ MPI(halo=2)**: 非自明配位で 12 桁一致。
- **OpenMP ↔ MPI**: 全列ビット一致（4×16³, 8 ランク）。

## ループ分類（`loop_classification/`）

格子測定とは独立に、HPE に現れる**閉ループ（Wilson 型 / Polyakov 型）を長さ L ごとに列挙・分類し、
その γ 行列トレースを計算**する Julia コード群。Wakabayashi 論文 [arXiv:2112.06340](https://arxiv.org/abs/2112.06340)
の Table 1（W⁰(n)）/ Table 2（L_m⁰(Nₜ,n)）を再現し、`trajectory_list/tables/` の係数テーブルの出所となる。

| ディレクトリ | 役割 |
|---|---|
| `src/`       | 計算コア（分類・対称性・γトレース; `classify_*`, `gamma_trace_*`, `loop_core_*`, `symmetry_*`） |
| `src_v2/`    | 高速版（正準形+Booth / 賢い生成）← 現行の主力 |
| `scripts/`   | 実行ドライバ（分類→値計算・CSV 出力・論文値突き合わせ。例 `writeout_tables.jl` で係数表生成） |

- 環境: `Project.toml`/`Manifest.toml`（`Combinatorics`, `JLD2`）。`julia --project=. scripts/<name>.jl` で実行（`--project=.` 必須）。
- 重い計算（n≥12 等）の本番リモートは `kitazaway:~/hpe_calc`（`src/ src_v2/` は本同梱版と一致）。
- 使い方は同梱 `loop_classification/README.md`。Wilson と Polyakov は別プロセスで実行（`xyzt_sym` の置換版衝突回避）。

## 実装メモ（正しさの肝）
- **trie + 開始位置最適化 (reach_min)**: 各ループを「空間 reach 最小」の開始位置に正準化し、
  N3LO の必要 halo を 2 に（naive な L/2=5 ではない）。trace の巡回不変性＋全サイト和で値は不変。
  → `TML/mdp/v2606/mdp_lattice.h` の `CART_HALO = 2`。
- **MPI: 測定前に `U.update()`**: 配位更新は depth-1 の袖（staple 用）しか同期しないが、trie は
  コーナー等の深い袖を読む。`src/measure_n3lo.cpp` の MC ループで測定直前に全袖を同期する
  （単一ランクでは no-op、MPI で必須）。
- **MeasureHPE は全ランクで HPE 評価**: collective（内部 `mpi.add`）を `if(ME==0)` の中に置くと
  MPI デッドロックするため、値は全ランクで計算し書き出しのみ rank 0。

## 既知の制限
- **N4LO（長さ 12）は別ビルド**: reach=3 なので `CART_HALO=3` に TML を再ビルドし、`trajectory_list/tables/`
  に長さ 12 を追加する必要がある（既定パッケージは N3LO まで）。
- 配位の保存（`GaugeSave_SU3_mult`）は MPI 領域分割でも動作確認済み。測定のみで回す場合は
  input 6 行目を `0` にすると保存をスキップできる。

---

物理・規約の詳細は同梱ノート、ならびに arXiv:2112.06340 を参照。
