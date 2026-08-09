# table1.jl — reproduce Table I of the paper from the committed shape tables.
#
#   N_shape  number of distinct shapes with D_j != 0
#   N_traj   number of different trajectories per lattice site,
#            N_traj = sum_j M_j / S_j     (summed over positive winding for L)
#
# Reads shape_tables/, which is produced by scripts/writeout_shapes.jl, so this
# runs in a fraction of a second and needs no classification.
#
#   julia scripts/table1.jl
#
# Pass a Polyakov Nt as the first argument to tabulate a different one
# (the paper shows Nt = 6).

cd(joinpath(@__DIR__, ".."))
const SHAPES = joinpath("..", "shape_tables")

Nt = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 6

# read "j,[l,]trajectory,M,S,D" and return (N_shape, N_traj) over D != 0
function tally(path)
    isfile(path) || return nothing
    nshape = 0; ntraj = 0//1
    open(path) do io
        header = split(strip(readline(io)), ',')
        iM = findfirst(==("M"), header); iS = findfirst(==("S"), header)
        iD = findfirst(==("D"), header)
        for line in eachline(io)
            isempty(strip(line)) && continue
            # the trajectory field is quoted and contains no comma
            f = split(line, ',')
            D = parse(Int, f[iD]); D == 0 && continue
            nshape += 1
            ntraj += parse(Int, f[iM]) // parse(Int, f[iS])
        end
    end
    (nshape, ntraj)
end

# three significant digits, as in the paper; small non-integers stay exact
function sci(v::Real)
    abs(v) < 1000 && return isinteger(v) ? string(Integer(v)) : ""
    e = floor(Int, log10(abs(v)))
    string(round(v / 10.0^e, digits=2), "e", e)
end
function sci(x::Rational)
    denominator(x) == 1 && return sci(numerator(x))
    abs(float(x)) < 1000 && return string(numerator(x), "/", denominator(x))
    sci(float(x))
end

println("TABLE I. Number of distinct shapes N_shape and total number of")
println("different trajectories per lattice site N_traj = sum_j M_j/S_j.")
println("Shapes with D_j = 0 are excluded. Polyakov rows sum over positive l.")
println()
println(rpad("", 12), lpad("N_shape", 10), lpad("N_traj", 14))

for n in 4:2:12
    r = tally(joinpath(SHAPES, "wilson", "w$n.csv"))
    r === nothing && continue
    println(rpad("W($n)", 12), lpad(r[1], 10), lpad(sci(r[2]), 14))
end
println()
for n in Nt:2:(Nt+8)
    r = tally(joinpath(SHAPES, "polyakov", "l_nt$(Nt)_n$(n).csv"))
    r === nothing && continue
    println(rpad("L($Nt,$n)", 12), lpad(r[1], 10), lpad(sci(r[2]), 14))
end
