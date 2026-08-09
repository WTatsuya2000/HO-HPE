# writeout_shapes.jl — write the complete shape classification tables to shape_tables/.
#
# For every shape (symmetry class of closed trajectories) it records the representative
# trajectory together with the three quantities the paper defines:
#
#   M_j  multiplicity : number of distinct trajectories per lattice site belonging to shape j
#   S_j  symmetry factor : number of cyclic permutations that leave the trajectory invariant
#   D_j  Dirac coefficient : tr_D[(1-gamma^{mu_1})(1-gamma^{mu_2}) ... ]
#
# so that   W(n)        = -2 N_c Sum_j (M_j D_j / S_j) Re W_j(n)
#           L_l(Nt,n)   = -(-1)^l N_c Sum_j (M_j D_j / S_j) L_j(Nt,n,l)
#
# Shapes with D_j = 0 are kept in the tables (they are listed in the paper's appendix too)
# but do not contribute; the shape counts of Table I exclude them.
#
# Wilson and Polyakov MUST run in separate processes: their xyzt_sym (24- vs 6-permutation
# versions) collide when loaded into the same session.
#
#   julia -t 8 scripts/writeout_shapes.jl wilson            # -> shape_tables/wilson/w<n>.csv
#   julia -t 8 scripts/writeout_shapes.jl polyakov 6        # -> shape_tables/polyakov/l_nt6_n<n>.csv
#
# Optional third argument caps the loop length (default 12; the Polyakov N4LO row for
# Nt = 6 needs 14 and takes a few minutes).

cd(joinpath(@__DIR__, ".."))
using Base.Threads

const OUTDIR = joinpath("..", "shape_tables")

# Trajectory as a space-separated list of signed axis indices, 1..4 = x,y,z,t,
# negative = traversed backwards. Same convention as trajectory_list/tables/*.dat.
const CODE = Dict{Symbol,Int}(:a=>1, :b=>2, :c=>3, :d=>4,
    Symbol("-a")=>-1, Symbol("-b")=>-2, Symbol("-c")=>-3, Symbol("-d")=>-4)
trajstr(w) = join((string(CODE[s]) for s in w), " ")

# number of cyclic shifts that leave the trajectory invariant
symfactor(w) = (s = 1; for j in 1:length(w)-1; w == circshift(w, j) && (s += 1); end; s)

mode = length(ARGS) >= 1 ? ARGS[1] : "wilson"
Lmax = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 12

if mode == "wilson"
    include(joinpath(@__DIR__, "..", "src", "gamma_trace_wilson.jl"))
    include(joinpath(@__DIR__, "..", "src", "gen_canon_wilson.jl"))
    dir = joinpath(OUTDIR, "wilson"); mkpath(dir)
    classify_loop_gen(4, 1)   # warm up the JIT so the reported timings are meaningful
    for n in 4:2:Lmax
        t = @elapsed begin
            io = open(joinpath(dir, "w$n.csv"), "w")
            println(io, "j,trajectory,M,S,D")
            j = 0; nshape = 0; ntraj = 0//1; omega = 0//1
            for ty in 1:length(count_combinations(n)), w in classify_loop_gen(n, ty)
                D = Int(real(character2gamma(w)))
                M = wilson_mult(w); S = symfactor(w)
                j += 1
                println(io, "$j,\"$(trajstr(w))\",$M,$S,$D")
                if D != 0
                    nshape += 1; ntraj += M//S; omega += M*D//S
                end
            end
            close(io)
            global RES = (j, nshape, ntraj, -6 * omega)
        end
        println("W($n): $(RES[1]) shapes listed, N_shape=$(RES[2]) (D!=0), " *
                "N_traj=$(RES[3]), W0=$(RES[4])   [$(round(t, digits=1)) s]")
    end

elseif mode == "polyakov"
    Nt = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6
    include(joinpath(@__DIR__, "..", "src", "gamma_trace_polyakov.jl"))
    include(joinpath(@__DIR__, "..", "src", "gen_canon_polyakov.jl"))
    dir = joinpath(OUTDIR, "polyakov"); mkpath(dir)

    # multiplicity under the Polyakov symmetry group (6 spatial axis permutations x
    # 8 spatial reflections; the time axis is fixed and reversal is not an element
    # because it flips the winding number)
    function poly_mult(w)
        P = Vector{Vector{Symbol}}()
        for f in transforms, i in 1:6
            g = f(xyzt_sym(w, i)); (g in P) || push!(P, g)
        end
        length(exclude_same(P))
    end

    gen_canon_poly(Nt, Nt, 1)   # warm up
    for n in Nt:2:Lmax
        t = @elapsed begin
            io = open(joinpath(dir, "l_nt$(Nt)_n$(n).csv"), "w")
            println(io, "j,l,trajectory,M,S,D")
            j = 0; nshape = 0; ntraj = 0//1; vals = Dict{Int,Rational{BigInt}}()
            for m in 1:div(n, Nt)
                (n - m*Nt >= 0 && (n - m*Nt) % 2 == 0) || continue
                reps, _ = gen_canon_poly(n, Nt, m)
                acc = 0//1
                for v in reps
                    w = tosym_p(v)
                    D = Int(real(character2gamma(w)))
                    M = poly_mult(w); S = symfactor(w)
                    j += 1
                    println(io, "$j,$m,\"$(trajstr(w))\",$M,$S,$D")
                    if D != 0
                        nshape += 1; ntraj += M//S; acc += M*D//S
                    end
                end
                vals[m] = (-1)^(m-1) * 6 * acc
            end
            close(io)
            global RES = (j, nshape, ntraj, vals)
        end
        vs = join(("L$m=$(RES[4][m])" for m in sort(collect(keys(RES[4])))), " ")
        println("L($Nt,$n): $(RES[1]) shapes listed, N_shape=$(RES[2]) (D!=0), " *
                "N_traj=$(RES[3])   $vs   [$(round(t, digits=1)) s]")
    end

else
    error("unknown mode: $mode  (use 'wilson' or 'polyakov [Nt] [Lmax]')")
end
