# HQQCD-HPE — Heavy-quark determinant HPE measurement (Nt = 4)

A self-contained package that expands the heavy Wilson-quark determinant `ln det M(κ)` of
lattice gauge theory via the **hopping parameter expansion (HPE)** and measures the κⁿ
coefficient of each order, configuration by configuration.

```
ln det M(κ) = N_site · Σ_n D_n κⁿ
```
`D_n` is a sum over closed loops of length n, split into **Wilson type** (loops that do not
wind in the time direction) and **Polyakov type** (loops winding m times in time, separated by
winding number m). This package measures **LO–N3LO (lengths 4 / 6 / 8 / 10)**
(LO = plaquette / Polyakov loop is output separately).

| Order | Length | Wilson | Polyakov |
|---|---|---|---|
| NLO  | 6  | W6  | P6 (w1) |
| NNLO | 8  | W8  | P8 (w1, w2) |
| N3LO | 10 | W10 | P10 (w1, w2) |

---

## Measurement method (per-order, mixed)

Each order uses the most efficient / canonical method (default build). **All methods give
numerically identical results.**

| Order | Length | Method | Implementation |
|---|---|---|---|
| NLO  | 6  | **direct product** (trace the loops directly) | `util/rectangle.hpp` `chair.hpp` `crown.hpp` `polyakov_loop_bend.hpp` |
| NNLO | 8  | **staple** (combine SingleStaple) | `NNLO_human/nnlo_w.cpp` `NNLO_W` / `NNLO_human/nnlo_p4.cpp` `NNLO_P4L1,L2` |
| N3LO | 10 | **trie** (meet-in-the-middle, shared halves) | `src/hpe_meas.cpp` |

- The NLO direct-product → single-value combination coefficients are **calibrated to machine
  precision** against the trie:
  - `W6 = 1152·rec_s + 576·rec_t1 + 576·rec_t2 + 1152·cha_s + 2304·cha_t1 + 1152·cha_t2 + 384·cro_s + 1152·cro_t`
  - `P6 = 1152·Polyakov_bend1 + 576·Polyakov_bend2`
- The NNLO staples (W8 / P8w1 / P8w2) match the trie directly (no fit needed).
- **`-DLOW_ORDER_TRIE` switches every sector to the trie path** (identical numbers; for
  cross-checking, see `make trie` below).

---

## Directory layout

```
HPE/
├── src/              measure_n3lo.cpp (production driver) , hpe_meas.cpp (trie engine)
├── NNLO_human/       human-written staple implementation
│                       nnlo.cpp = bundles single_staple / local_elem / nnlo_w / nnlo_p4
│                       nnlo_safe.cpp = W8traj/W10traj/MakeTrajList/cyclic_symmetry
│                       calc_hpe.cpp
├── trajectory_list/  loop direction lists nnlo_traj_list.cpp , n3lo_traj_list.cpp
│   └── tables/        coefficient tables *.dat and generators *_table.cpp (W6/8/10, P6/8/10, P12, poly_68)
├── util/             rectangle/chair/crown/polyakov_loop(_bend)/staple_pol/gauge_io_mult/
│                       gauge_update_routines_pol/stopwatch
├── loop_classification/  per-length Wilson/Polyakov closed-loop enumeration, classification & γ-trace (Julia)
│                       src/ scripts/ (see its own README)
├── TML/              lattice library (mdp/v2606, CART_HALO=2)
├── jobs/             Slurm submission scripts (sub_openmp.sh, sub_mpi.sh)
├── input
├── makefile
└── conf/ data/ log/  generated at runtime
```

> **Runtime table lookup**: `a.out` reads the coefficient tables from `trajectory_list/tables/`
> **relative to the CWD**, so always launch from the package root (`jobs/*.sh` already `cd
> $SLURM_SUBMIT_DIR`).

## Requirements
- Intel oneAPI (Intel MPI + OpenMP) — `mpiicpx`
- Before building/running: `module load intel/2025.3.1` under `bash -l`
  (Slurm jobs use `#!/bin/bash -l` for a login shell and `module load` inside the job).

## Build
Executables (`*.out`) are **not shipped**; build them in place (build products, removed by `make clean`).

```sh
make            # a.out          : mixed, single rank + OpenMP (div=1,1,1,1)   <- recommended
make mpi        # a_mpi.out      : mixed, 8-rank MPI          (div=1,2,2,2)
make trie       # a_trie.out     : all-trie (-DLOW_ORDER_TRIE), single rank
make mpi_trie   # a_mpi_trie.out : all-trie, 8-rank
make clean      # remove *.out
```

The binaries differ **only in parallelization**; the measured content (mixed measurement) and
the output numbers are identical:

| Binary | make | Parallelism | Use |
|---|---|---|---|
| `a.out`          | `make`      | single rank + OpenMP threads          | small lattices, quick (recommended) |
| `a_mpi.out`      | `make mpi`  | **8-rank MPI** (space split 2×2×2) + OpenMP | large lattices across nodes (use this for the 48³ production) |
| `a_trie.out`     | `make trie` | single rank (all sectors via trie)    | cross-check / performance comparison |
| `a_mpi_trie.out` | `make mpi_trie` | 8-rank MPI (all-trie)             | same, MPI |

- `-DLOW_ORDER_TRIE` (`trie`/`mpi_trie`) = measure every sector via the trie path; **numerically
  identical** to mixed (for verification).
- `-DDIV1=2 -DDIV2=2 -DDIV3=2` (`mpi`/`mpi_trie`) = split the 3 spatial directions in 2 → 8 ranks
  total; the time direction is not split.
- Measurement is safe under OpenMP threading (13-digit reproducible).
- Gauge update (heat bath / over-relaxation) assumes a single process; the MPI build parallelizes
  via domain decomposition.

## Input format (8 lines, fixed order; the first token of each line is the value, the rest is a comment)
```
4 48 48 48   lattice size (Nt Nx Ny Nz)
5.8899       beta
0.0014       lambda   (LO Polyakov coupling; internally lambda/=beta)
conf         configuration directory
data         data storing directory
5000         number of measurements   (=0 -> measure only, no save)
10           sweeps between measurements
10           HPE max loop length (10=N3LO, even)
```
- Line 8 = maximum loop length to include. `6`/`8`/`10` are allowed. `12` (N4LO) requires
  rebuilding TML with `CART_HALO=3`.

## Run (submit to compute nodes via Slurm; do not run directly on the login node)
```sh
sbatch jobs/sub_openmp.sh    # single rank OpenMP (a.out, input)
sbatch jobs/sub_mpi.sh       # 8-rank MPI        (a_mpi.out, input)  (run make mpi first)
```
- Starting from a cold start, the first measurement row (before any update) holds the canonical
  values (below). To measure existing configurations only, set input line 6 (number of
  measurements) to `0` (measure only; no update/save).
- MPI note: the local sublattice along each split spatial direction must be ≥ 2·halo (N3LO has
  reach=2 → local ≥ 4). A startup guard checks this and aborts FATAL if not satisfied.

## Output (`data/hpe_<param>_<counter>`, one line per configuration)
```
# plaq_s plaq_t Re_pol Im_pol  W6 W8 W10  Re_L1n6 Im_L1n6  Re_L1n8 Im_L1n8  Re_L2n8 Im_L2n8  Re_L1n10 Im_L1n10  Re_L2n10 Im_L2n10
```
- Wilson (W6/8/10) is one real column; Polyakov (`L<m>n<n>` = winding m, length n) is Re/Im (two columns).
- The column order and count are **identical** regardless of method (mixed / all-trie).

## Verification (reproducibility)
- **cold (U=1) values = paper values** (confirmed with this package's build):
  `W6=8448, W8=245952, W10=7372800, L1n6=1728, L1n8=45792, L2n8=−384, L1n10=645120, L2n10=−64512`.
- **mixed == all-trie**: all columns match to machine precision on the same thermalized
  configurations (verified for both cold and thermalized).
- **single ↔ MPI (halo=2)**: 12-digit agreement on a non-trivial configuration.
- **OpenMP ↔ MPI**: bit-identical for all columns (4×16³, 8 ranks).

## Loop classification (`loop_classification/`)

Independent of the lattice measurement, this is Julia code that **enumerates and classifies the
closed loops (Wilson type / Polyakov type) of the HPE per length L and computes their γ-matrix
traces**. It reproduces Table 1 (W⁰(n)) / Table 2 (L_m⁰(Nₜ,n)) of Wakabayashi et al.
[arXiv:2112.06340](https://arxiv.org/abs/2112.06340), and is the origin of the coefficient tables
in `trajectory_list/tables/`.

| Directory | Role |
|---|---|
| `src/`     | computation core (classification, symmetry, γ-trace, and the fast invariant engine) |
| `scripts/` | drivers (classify → compute values, CSV output, paper-value cross-check) |

- Verified: `W⁰(4)=288, W⁰(6)=8448, W⁰(8)=245952, W⁰(10)=7372800` and `L₁⁰(4,8)=45792`.
- See `loop_classification/README.md` for usage. Run Wilson and Polyakov in separate processes
  (their `xyzt_sym` permutation variants collide otherwise).

## Implementation notes (the crux of correctness)
- **trie + start-position optimization (reach_min)**: canonicalize each loop to the start with
  minimal spatial reach, bringing N3LO's required halo down to 2 (not the naive L/2=5). The value
  is unchanged thanks to cyclic invariance of the trace plus the full-site sum.
  → `CART_HALO = 2` in `TML/mdp/v2606/mdp_lattice.h`.
- **MPI: `U.update()` before measuring**: the gauge update only syncs the depth-1 halo (for
  staples), but the trie reads deeper halos (corners etc.). The MC loop in `src/measure_n3lo.cpp`
  syncs all halos just before measuring (no-op for single rank; required for MPI).
- **MeasureHPE evaluates the HPE on all ranks**: putting the collective (internal `mpi.add`)
  inside `if(ME==0)` would deadlock MPI, so all ranks compute the values and only rank 0 writes.

## Known limitations
- **N4LO (length 12) is a separate build**: reach=3, so TML must be rebuilt with `CART_HALO=3`
  and length 12 added to `trajectory_list/tables/` (the default package goes up to N3LO).
- Configuration saving (`GaugeSave_SU3_mult`) is verified to work under MPI domain decomposition.
  For measure-only runs, set input line 6 to `0` to skip saving.

---

For physics/convention details, see the bundled notes and arXiv:2112.06340.
