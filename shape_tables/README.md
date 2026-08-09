# shape_tables — complete shape classification

The trajectories contributing to each order of the hopping parameter expansion are
grouped into **shapes**: classes of trajectories that become identical under cyclic
permutation (CyP), reversal, and the cubic rotations/reflections (CuRR) of the
hypercubic group — H(4) for the Wilson type, H(3) for the Polyakov type, where the
time axis is fixed and reversal is not a symmetry because it flips the winding number.

Each shape `j` carries three numbers:

| | meaning |
|---|---|
| `M` | multiplicity: distinct trajectories per lattice site belonging to shape `j` |
| `S` | symmetry factor: cyclic permutations leaving the trajectory invariant |
| `D` | Dirac coefficient, `tr_D[(1-γ^{μ₁})(1-γ^{μ₂})···]` |

from which the expansion coefficients follow:

```
W(n)        = -2 N_c Σ_j (M_j D_j / S_j) Re W_j(n)
L_l(Nt,n)   = -(-1)^l N_c Σ_j (M_j D_j / S_j) L_j(Nt,n,l)
```

## Files

```
wilson/w<n>.csv            n = 4, 6, 8, 10, 12       columns: j,trajectory,M,S,D
polyakov/l_nt<Nt>_n<n>.csv Nt = 4 (n ≤ 12), 6 (n ≤ 14)  columns: j,l,trajectory,M,S,D
```

`l` is the temporal winding number. `trajectory` is the representative, written as
signed axis indices (`1 2 3 4` = `x y z t`, negative = traversed backwards) — the same
convention as `trajectory_list/tables/*.dat`.

Shapes with `D = 0` **are listed here** (they appear in the paper's appendix too) but
contribute nothing, and the shape counts below exclude them. The `.dat` coefficient
tables used by the measurement engine omit them entirely, which is why a `.dat` file
has fewer entries than the corresponding `.csv`.

Human-readable names for the W(8) shapes (`U=U`, `PP`, `ULU`, `crown + U`, `4d-crown`, …)
are not reproduced here; see `NNLO_human/NNLO_list.pdf` and Table IV of the paper.

## Shape and trajectory counts

`N_shape` counts shapes with `D ≠ 0`; `N_traj = Σ_j M_j / S_j` is the number of
different trajectories per lattice site. Polyakov rows sum over positive `l`.

| | N_shape | N_traj | | | N_shape | N_traj |
|---|---|---|---|---|---|---|
| W(4)  | 1    | 6      | | L(6,6)  | 1     | 1/6    |
| W(6)  | 3    | 76     | | L(6,8)  | 3     | 15     |
| W(8)  | 24   | 1.71e3 | | L(6,10) | 56    | 870    |
| W(10) | 189  | 3.80e4 | | L(6,12) | 1477  | 4.19e4 |
| W(12) | 3701 | 1.03e6 | | L(6,14) | 56395 | 2.06e6 |

Regenerate this table from the files with `julia loop_classification/scripts/table1.jl`.

## Regenerating

Run from `loop_classification/`. Wilson and Polyakov must run in **separate
processes** — their `xyzt_sym` (24- vs 6-permutation) definitions collide.

```bash
julia -t 8 scripts/writeout_shapes.jl wilson   0 12   # W(4) ... W(12)
julia -t 8 scripts/writeout_shapes.jl polyakov 4 12   # L(4,4) ... L(4,12)
julia -t 8 scripts/writeout_shapes.jl polyakov 6 14   # L(6,6) ... L(6,14)
```

Timings on a laptop (8 threads): the whole Wilson side including W(12) takes about
25 s, `L(4,12)` about 95 s, and `L(6,14)` — the largest case, 85832 shapes listed —
about 12 minutes.

The classification is checked against the known values of the expansion coefficients
on a unit gauge configuration, which the scripts print as they go:

```
W(4)=288  W(6)=8448  W(8)=245952  W(10)=7372800  W(12)=225232896
L1(4,6)=1728   L1(4,8)=45792   L2(4,8)=-384
L1(4,10)=645120   L2(4,10)=-64512
L1(4,12)=-26224128   L2(4,12)=-6842880   L3(4,12)=4096
L1(6,6)=128   L1(6,8)=11520   L1(6,10)=716544
L1(6,12)=35891712   L2(6,12)=-4096
```
