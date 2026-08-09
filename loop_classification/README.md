# loop_classification — trajectory enumeration, shape classification & γ-trace

Julia code that **enumerates the closed trajectories (Wilson-type / Polyakov-type)
appearing in the hopping parameter expansion (HPE) of lattice QCD, length by length,
groups them into shapes, and computes their Dirac (γ-matrix) traces**.
It reproduces Table 1 (W⁰(n)) and Table 2 (L_m⁰(Nₜ,n)) of Wakabayashi et al.,
[arXiv:2112.06340](https://arxiv.org/abs/2112.06340), and extends them to N4LO.

Two independent classification engines are kept, and agreement between them is the
check on the classification itself:

* **constructive generation** (`gen_canon_*`) — the symmetry is folded into the
  generation, so only near-canonical sequences are produced and the orbit is closed by
  Booth's minimal-rotation canonical form. This is what makes N4LO tractable.
* **invariant buckets** (`invariant_*`) — a genuine symmetry invariant is computed,
  candidates are bucketed by it and merged with union-find inside each bucket. Slower,
  used to cross-check the above at NNLO/N3LO.

## Dependencies (install into your global Julia environment)
```julia
import Pkg; Pkg.add(["Combinatorics", "SpecialFunctions", "CSV", "DataFrames"])
```

## Run

```bash
julia -t 8 scripts/record_wilson.jl        # W⁰(n), n = 4 ... 12       -> records/wilson_*.csv
julia -t 8 scripts/record_polyakov.jl      # L_m⁰(Nt,n)                -> records/polyakov_*.csv
julia -t 8 scripts/writeout_shapes.jl wilson       0 12   # -> ../shape_tables/wilson/
julia -t 8 scripts/writeout_shapes.jl polyakov     6 14   # -> ../shape_tables/polyakov/
julia      scripts/table1.jl               # shape / trajectory counts, read back from shape_tables/
julia -t 8 scripts/writeout_nnlo_nnnlo.jl  # independent cross-check at NNLO / N3LO
```

Verified output of `record_wilson.jl`:
`W⁰(4)=288, W⁰(6)=8448, W⁰(8)=245952, W⁰(10)=7372800, W⁰(12)=225232896`.
`record_polyakov.jl` walks an in-script `grid` of `(n, Nt)` pairs and checks winding 1
against the published values; edit the grid to narrow the range.

**Run Wilson and Polyakov in separate processes** — their `xyzt_sym` (24- vs
6-permutation versions) collide if loaded into the same session.

Timings on a laptop with 8 threads: the entire Wilson side through W(12) takes about
25 s. On the Polyakov side `L(4,12)` takes about 95 s and `L(6,14)` — the N4LO term for
`Nt = 6`, 85832 shapes listed — about 12 minutes.

## Layout
```
src/   (shared primitives)
  gamma_trace_{wilson,polyakov}.jl   γ-trace, character2gamma  (includes loop_core)
  loop_core_{wilson,polyakov}.jl     trajectory generation & symmetry (includes symmetry)
  symmetry_{wilson,polyakov}.jl      axis permutations, reflections, xyzt_sym
  classify_polyakov.jl               is_circshift_equal & generation helpers

  fast_wilson.jl                     canonical form + Booth; classify_loop_fastw, wilson_mult
  gen_canon_wilson.jl                classify_loop_gen       constructive generation (Wilson)
  gen_canon_polyakov.jl              gen_canon_poly          constructive generation (Polyakov)

  invariant_core.jl                  invariant key + bucket/union-find classification
  invariant_wilson.jl                classify_loop_inv       (24 axis permutations)
  invariant_polyakov.jl              classify_loop_inv_poly  (winding-aware canonical form)
scripts/
  record_wilson.jl  record_polyakov.jl     coefficient values, per length
  writeout_shapes.jl                       the complete shape tables -> ../shape_tables/
  table1.jl                                shape / trajectory counts
  writeout_nnlo_nnnlo.jl                   independent cross-check (invariant engine)
```

## Symmetry groups

| | axis permutations | reflections | reversal | translation |
|---|---|---|---|---|
| Wilson | all 4 axes, 24 | all 4 axes, 16 signs | yes | cyclic shift |
| Polyakov | spatial only, 6 | spatial only, 8 signs | no — it flips the winding number | cyclic shift |

## Conventions

Directions `:a :b :c :d` = x y z t; reversed links are `Symbol("-a")` etc. (integers 1..4
and their negatives). A trajectory is a `Vector{Symbol}`.

## Notes
- Generated `records/` output is not committed; the shape tables in `../shape_tables/` are.
- Heavy cases are best given a dedicated machine, but everything through N4LO runs on a
  laptop.
