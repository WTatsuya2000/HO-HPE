# loop_classification — closed-loop enumeration, classification & γ-trace

Julia code that **enumerates and classifies the closed loops (Wilson-type / Polyakov-type)
appearing in the hopping parameter expansion (HPE) of lattice QCD, length L by length L,
and computes their Dirac (γ-matrix) traces**.
It reproduces Table 1 (W⁰(n)) and Table 2 (L_m⁰(Nₜ,n)) of Wakabayashi et al.,
[arXiv:2112.06340](https://arxiv.org/abs/2112.06340).

This is a minimal, working subset: only the current, verified method is kept
(the old brute-force `classify`/`canonical` variants, notebooks, regression `verify_*`/`data`,
and `Project.toml`/`Manifest.toml` have been removed).

## Dependencies (install into your global Julia environment)
```julia
import Pkg; Pkg.add(["Combinatorics", "SpecialFunctions", "CSV", "DataFrames"])
```

## Run (verified)
```bash
julia scripts/record_wilson.jl        # Wilson   W⁰(n)        -> records/wilson_*.csv
julia scripts/record_polyakov.jl      # Polyakov L_m⁰(Nt,n)   -> records/polyakov_*.csv   (heavy; edit the `grid` to narrow the range)
julia scripts/writeout_nnlo_nnnlo.jl  # NNLO/N3LO representative loops & coefficient rows
```
Verified output of `record_wilson.jl`: `W⁰(4)=288, W⁰(6)=8448, W⁰(8)=245952, W⁰(10)=7372800`
(matches the paper). Polyakov (invariant method): `L₁⁰(4,6)=1728, L₁⁰(4,8)=45792` (matches the paper).

## Layout
```
src/   (shared primitives)
  gamma_trace_{wilson,polyakov}.jl   γ-trace, character2gamma, etc. (includes loop_core)
  loop_core_{wilson,polyakov}.jl     loop generation & symmetry    (includes symmetry)
  symmetry_{wilson,polyakov}.jl      axis permutations, reflections, xyzt_sym
  classify_polyakov.jl               is_circshift_equal & generation helpers (used by invariant_polyakov)
       (fast classification engine)
  invariant_core.jl                  invariant key + classify_invariant (bucket + union-find)
  invariant_wilson.jl                classify_loop_inv      (Wilson, 24 axis permutations)
  invariant_polyakov.jl              classify_loop_inv_poly (Polyakov, winding-aware canonical form)
  fast_wilson.jl                     classify_loop_fastw    (canonical form + Booth; used by record_wilson)
  gen_canon_polyakov.jl              gen_canon_poly         (used by record_polyakov -> invariant_polyakov)
scripts/
  record_wilson.jl  record_polyakov.jl  writeout_nnlo_nnnlo.jl
```

## Conventions
Directions `:a :b :c :d` = x y z t; reversed links are `Symbol("-a")` etc. (integers 1..4 / negatives).
A loop is a `Vector{Symbol}`. **Run Wilson and Polyakov in separate processes** — their `xyzt_sym`
(24- vs 6-permutation versions) collide if loaded in the same session.

## Notes
- Heavy cases (n ≥ 12, etc.) are best run on a dedicated machine; the code is self-contained.
- `record_polyakov.jl` has an in-script `grid` (target (n, Nt)) you can edit to limit the range.
