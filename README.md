# HO-HPE — Heavy-quark QCD hopping parameter expansion

Code for the **hopping parameter expansion (HPE)** of the heavy Wilson-quark determinant
`ln det M(κ)` in lattice gauge theory:

```
ln det M(κ) = N_site · Σ_n C_n κⁿ
```

`C_n` is a sum over closed trajectories of length `n`, split into **Wilson type** `W(n)`
(trajectories that do not wind in time) and **Polyakov type** `L_l(Nt,n)` (winding `l`
times in time). Trajectories are grouped into **shapes** — classes that coincide under
cyclic permutation, reversal, and the cubic rotations/reflections of the hypercubic group —
and each shape `j` contributes with its multiplicity `M_j`, symmetry factor `S_j` and
Dirac coefficient `D_j`:

```
W(n)       = -2 N_c Σ_j (M_j D_j / S_j) Re W_j(n)
L_l(Nt,n)  = -(-1)^l N_c Σ_j (M_j D_j / S_j) L_j(Nt,n,l)
```

| Order | Length | Wilson | Polyakov (Nt = 4) |
|---|---|---|---|
| LO   | 4  | W4  | L1(4,4) |
| NLO  | 6  | W6  | L1(4,6) |
| NNLO | 8  | W8  | L1, L2 (4,8) |
| N3LO | 10 | W10 | L1, L2 (4,10) |
| N4LO | 12 | W12 | L1, L2, L3 (4,12) |

Everything here covers the expansion up to **N4LO**. Tables are also provided for
`Nt = 6, 8, 10`; note that N4LO means length `Nt + 8` for the Polyakov type, so the
N4LO Polyakov term at `Nt = 6` is `L(6,14)`.

---

## Contents

```
loop_classification/   trajectory enumeration, shape classification & γ-trace  (Julia, runnable)
  src/ scripts/ README.md

shape_tables/          the complete shape lists: representative trajectory, M_j, S_j, D_j
  wilson/w<n>.csv                n = 4 ... 12
  polyakov/l_nt<Nt>_n<n>.csv     Nt = 4 (n ≤ 12), Nt = 6 (n ≤ 14)

src/
  hpe_pool.cpp         trie engine: shared-pool meet-in-the-middle, lengths 6 ... 12
  hpe_cost.cpp         cost counter — builds and runs with no lattice library
  lattice_stub.hpp     minimal declarations that let the cost counter compile

trajectory_list/
  nnlo_traj_list.cpp n3lo_traj_list.cpp   trajectory direction lists (lengths 8, 10)
  tables/              coefficient tables w<n>.dat, p<n>_nt<Nt>_w<l>.dat  (lengths 6 ... 12),
                       plus *_table.cpp, compile-time variants of the short tables that the
                       current engine does not use

NNLO_human/            human-developed staple implementation and its tuning history v1 ... v6
  readme.md            which version corresponds to which arXiv revision
  NNLO_list.pdf        the W(8) shape list with the human-readable shape names

misc/
  AGENTS.md            instruction file used for the initial AI-assisted development of the
                       N3LO calculation (in Japanese)
  nnnlo_v3.cpp         first output from Codex/ChatGPT 5.5 — the ancestor of src/hpe_pool.cpp
```

## Measurement method

| Order | Length | Method | Implementation |
|---|---|---|---|
| NNLO | 8 | staple: `SingleStaple` combinations chosen shape by shape | `NNLO_human/nnlo_w.cpp` (`NNLO_W`), `NNLO_human/nnlo_p4.cpp` (`NNLO_P4L1,L2`) |
| NLO … N4LO | 6 … 12 | trie (meet-in-the-middle, shared pool) | `src/hpe_pool.cpp` |

The trie algorithm cuts every trajectory at a chosen starting point into a forward and a
backward subtrajectory of equal length, reverses and conjugates the backward one so that
both are ordinary link products, and stores them in a single trie. Matrix products shared
by several subtrajectories are then evaluated once. The trajectory is closed by
`Tr(P·S†)` over the two leaves.

The shape-based approach remains the better one at N2LO, where the number of trajectories
is still small enough to optimise each shape by hand; the trie takes over beyond that,
where the trajectory count grows explosively.

Reference values on a unit (`U = 1`) configuration, reproduced both by the C++ engine and
by `loop_classification/`:

```
W(4)=288   W(6)=8448   W(8)=245952   W(10)=7372800   W(12)=225232896
L1(4,6)=1728
L1(4,8)=45792        L2(4,8)=-384
L1(4,10)=645120      L2(4,10)=-64512
L1(4,12)=-26224128   L2(4,12)=-6842880   L3(4,12)=4096
L1(6,6)=128   L1(6,8)=11520   L1(6,10)=716544   L1(6,12)=35891712   L2(6,12)=-4096
```

## Cost of the trie algorithm

Every trie node other than the root is one matrix multiplication (MM) per lattice site, so
the cost is directly countable from the tables alone — no gauge field required.
`src/hpe_cost.cpp` does this and needs nothing but a C++17 compiler:

```bash
c++ -std=gnu++17 -O2 -fpermissive -I src -o hpe_cost src/hpe_cost.cpp
./hpe_cost 4 12          # Nt=4, up to length 12, one trie per term
```

Costs are quoted in units of one `SingleStaple` evaluation = 48 MM/site. Measured with the
default anchor rule, one pool per term:

| term | MM/site | staples | traces | reach |
|---|---|---|---|---|
| W(6)  | 93    | 1.9   | 76        | 1 |
| W(8)  | 694   | 14.5  | 1716      | 2 |
| W(10) | 4815  | 100.3 | 38040     | 3 |
| W(12) | 34217 | 712.9 | 1031830   | 3 |
| L1(4,10) | 4414  | 92.0  | 9093   | 3 |
| L1(4,12) | 38834 | 809.0 | 361497 | 3 |

`traces` counts the `Re Tr[XY]` operations closing the trajectories; for `SU(N_c)` each is
worth `1/(2N_c)` MM. `reach` is the spatial halo depth the sector needs under MPI domain
decomposition — the choice of starting point bounds it well below the naive `L/2`.

Letting all terms share a single trie reduces the total for the twelve sectors at `Nt = 4`,
`Lmax = 12` from 115396 MM/site (one pool per term) to **51596** MM/site:

```bash
./hpe_cost 4 12 1 2      # anchor = split_min, all sectors share one pool
```

`src/hpe_cost.cpp` also implements the alternative anchor rules, so their cost can be
compared directly (`./hpe_cost 4 12 <anchor> <pool>`; see the header of the file).

> **Note — the C++ part is implementation sources, not a standalone build.**
> The production driver, `makefile`, the lattice library, batch job scripts and input files
> are maintained in the full measurement environment and are **not** included here. The
> `.cpp` files are meant to be `#include`d by that driver. The one exception is
> `src/hpe_cost.cpp`, which builds on its own via `src/lattice_stub.hpp` — a set of
> placeholder declarations that computes nothing and exists only so that the cost counter
> and the symmetry orbits can be exercised without the lattice library.

## Quick start — classification

```bash
cd loop_classification
julia -e 'import Pkg; Pkg.add(["Combinatorics","SpecialFunctions","CSV","DataFrames"])'   # once

julia -t 8 scripts/record_wilson.jl      # W(4) ... W(12) against the reference values
julia    scripts/table1.jl               # shape and trajectory counts, from shape_tables/
```

See [`loop_classification/README.md`](loop_classification/README.md) and
[`shape_tables/README.md`](shape_tables/README.md) for the full usage.

---

For physics and conventions, see arXiv:2112.06340 (Wakabayashi, Ejiri, Kanaya, Kitazawa).
