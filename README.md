# HO-HPE — Heavy-quark QCD hopping parameter expansion

Code for the **hopping parameter expansion (HPE)** of the heavy Wilson-quark determinant
`ln det M(κ)` in lattice gauge theory:

```
ln det M(κ) = N_site · Σ_n D_n κⁿ
```
`D_n` is a sum over closed loops of length n, split into **Wilson type** (loops that do not wind
in time) and **Polyakov type** (loops winding m times in time, by winding number m).

This repository contains two parts:

1. **`loop_classification/`** — a self-contained, runnable Julia tool that enumerates and
   classifies the closed loops per length L and computes their Dirac (γ) traces, reproducing
   Table 1 (W⁰(n)) / Table 2 (L_m⁰(Nₜ,n)) of Wakabayashi et al.
   [arXiv:2112.06340](https://arxiv.org/abs/2112.06340). See its own
   [README](loop_classification/README.md).
2. **NNLO measurement sources** (C++) — the operator implementations, trie engine, and
   coefficient tables used to measure W6/W8/W10 and the Polyakov sectors on gauge configurations.
   These are **implementation sources only** (see the note below).

| Order | Length | Wilson | Polyakov |
|---|---|---|---|
| NLO  | 6  | W6  | P6 (w1) |
| NNLO | 8  | W8  | P8 (w1, w2) |
| N3LO | 10 | W10 | P10 (w1, w2) |

---

## Contents

```
loop_classification/   per-length Wilson/Polyakov loop enumeration, classification & γ-trace (Julia, runnable)
  src/ scripts/ README.md

NNLO_human/            Human-developed history (C++)

misc/
  AGENTS.md            Instruction file used for the initial AI-assisted development of the N^3LO calculation (Japanese language)
  nnnlo_v3.cpp         First output from Codex/ChatGPT 5.5

src/
  hpe_meas.cpp         trie engine (HPE_W6/8/10, HPE_P6/8/10; meet-in-the-middle)

trajectory_list/
  nnlo_traj_list.cpp n3lo_traj_list.cpp   loop direction lists
  tables/            coefficient tables *.dat and generators *_table.cpp (W6/8/10, P6/8/10, poly_68)
```

### Measurement method (per-order, mixed)

Each order uses the most efficient / canonical method; all methods give numerically identical results.

| Order | Length | Method | Implementation |
|---|---|---|---|
| NLO  | 6  | direct product | `NNLO_human/rectangle.hpp` `chair.hpp` `crown.hpp` `polyakov_loop_bend.hpp` |
| NNLO | 8  | staple | `NNLO_human/nnlo_w.cpp` (`NNLO_W`) / `NNLO_human/nnlo_p4.cpp` (`NNLO_P4L1,L2`) |
| N3LO | 10 | trie (meet-in-the-middle) | `src/hpe_meas.cpp` |

Canonical cold (U=1) values reproduced by these sources / by `loop_classification/`:
`W6=8448, W8=245952, W10=7372800, L1n6=1728, L1n8=45792, L2n8=−384, L1n10=645120, L2n10=−64512`
(and `W⁰(4)=288`).

> **Note — the C++ part is implementation sources, not a standalone build.**
> The production driver (`measure_n3lo.cpp`), `makefile`, the TML lattice library, Slurm job
> scripts, and the input files are maintained in the full measurement environment and are **not**
> included here. The `.cpp/.hpp` files above are meant to be `#include`d by that driver; they do
> not compile on their own in this repository.

## loop_classification quick start

```bash
cd loop_classification
import_pkg='import Pkg; Pkg.add(["Combinatorics","SpecialFunctions","CSV","DataFrames"])'
julia -e "$import_pkg"            # one-time, into your global environment
julia scripts/record_wilson.jl   # -> W⁰(4)=288, W⁰(6)=8448, W⁰(8)=245952, W⁰(10)=7372800
```
See [`loop_classification/README.md`](loop_classification/README.md) for the full usage.

---

For physics/convention details, see arXiv:2112.06340.
