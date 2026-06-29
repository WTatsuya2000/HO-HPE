# writeout_nnlo_nnnlo.jl
# NNLO (length L=8) と NNNLO (length L=10) のループを書き下す。
# 完全版の正準形分類 src を用いて各代表ループの多重度 M・γトレース D・寄与を出力する。
# (data/ の旧分類は L=10 type4 が欠落しているため使わない)
#
# 使い方（プロジェクトルートから・スレッド指定推奨・別プロセスで）:
#   julia -t 8 --project=. scripts/writeout_nnlo_nnnlo.jl wilson
#   julia -t 8 --project=. scripts/writeout_nnlo_nnnlo.jl polyakov 4
#
# Wilson:   W^0(n)     = -6 * Σ (M/LLL) * D        (24軸置換, src classify_loop_inv)
# Polyakov: L1^0(Nt,n) = +6 * Σ (M/LLL) * D        (6軸置換, src classify_loop_inv_poly, 巻き数 m=1)

cd(joinpath(@__DIR__, ".."))
using Base.Threads

# 記号を読みやすい文字列に: a->x, b->y, c->z, d->t（- 接頭辞は逆向き）
const NAME = Dict("a"=>"x","b"=>"y","c"=>"z","d"=>"t")
function pretty(arr)
    parts = String[]
    for s in String.(arr)
        neg = startswith(s, "-")
        base = neg ? s[2:end] : s
        nm = get(NAME, base, base)
        push!(parts, neg ? "-$nm" : " $nm")
    end
    return "[" * join(parts, ",") * "]"
end

# 1ループの (M, D, LLL) を計算。nsym = 軸置換の数（Wilson 24 / Polyakov 6）
function loop_MDL(arr2, nsym)
    Path1_pre = Vector{Any}()
    for f in transforms, j = 1:nsym
        g = f(xyzt_sym(arr2, j))
        any(x -> x == g, Path1_pre) || push!(Path1_pre, g)
    end
    Morbit = length(exclude_same(Path1_pre))
    LLL = 1
    for j = 1:length(arr2)-1
        arr2 == circshift(arr2, j) ? (LLL += 1) : nothing
    end
    D = Int(real(character2gamma(arr2)))
    return Morbit, D, LLL
end

mode = length(ARGS) >= 1 ? ARGS[1] : "wilson"
println("# nthreads = ", nthreads())

if mode == "wilson"
    include(joinpath(@__DIR__, "..", "src", "invariant_wilson.jl"))
    include(joinpath(@__DIR__, "..", "src", "gamma_trace_wilson.jl"))
    for L in (8, 10)
        order = L == 8 ? "NNLO" : "NNNLO"
        ntype = length(count_combinations(L))
        println("\n########## Wilson  L=$L  ($order)  ##########")
        println("# idx | type | loop | M | D | LLL | contrib(-6*(M/LLL)*D)")
        total = 0; idx = 0
        for ty = 1:ntype
            reps = classify_loop_inv(L, ty)
            for arr2 in reps
                Morbit, D, LLL = loop_MDL(arr2, 24)
                M = Int(Morbit ÷ LLL)
                contrib = -6 * M * D
                total += contrib; idx += 1
                println("$idx | t$ty | $(pretty(arr2)) | $M | $D | $LLL | $contrib")
            end
        end
        println("# representatives: $idx")
        println("# W^0($L) = $total")
    end

elseif mode == "polyakov"
    Nt = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
    include(joinpath(@__DIR__, "..", "src", "invariant_polyakov.jl"))
    include(joinpath(@__DIR__, "..", "src", "gamma_trace_polyakov.jl"))
    for L in (8, 10)
        order = L == 8 ? "NNLO" : "NNNLO"
        println("\n########## Polyakov  Nt=$Nt  L=$L  ($order, winding m=1)  ##########")
        println("# idx | loop | M | D | LLL | contrib(+6*(M/LLL)*D)")
        reps = classify_loop_inv_poly(L, Nt)
        total = 0; idx = 0
        for arr2 in reps
            Morbit, D, LLL = loop_MDL(arr2, 6)
            contrib = 6 * Morbit * D ÷ LLL
            total += contrib; idx += 1
            println("$idx | $(pretty(arr2)) | $Morbit | $D | $LLL | $contrib")
        end
        println("# representatives: $idx")
        println("# L1^0(Nt=$Nt, $L) = $total")
    end
else
    error("unknown mode: $mode  (use 'wilson' or 'polyakov [Nt]')")
end
