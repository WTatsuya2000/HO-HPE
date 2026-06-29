# Polyakov ループ分類（不変量バケツ方式・winding 対応生成, -d 込み）
#
# 生成方針:
#   巻き数 m=1 の Polyakov ループ＝「a,b,c は net 0、d は net = Nt」。
#   k = (L-Nt)/2 を「空間ペア数 na+nb+nc」と「余分な (+d,-d) ペア数 q」に分配:
#     na+nb+nc+q = k,   #(+d)=Nt+q, #(-d)=q  →  net(d)=Nt, 長さ=2k+Nt=L。
#   q>=1 のとき -d を含む（時間方向の巻き戻り）ループが生成される。
#   net(d)=Nt 固定なので巻き数は常に 1（m>=2 の混入なし）。
#
# 分類は不変量バケツ方式（invariant_core.jl）。同値判定は Polyakov 版 is_circshift_equal
# （a,b,c の 6 置換, d 固定）を使用。
include(joinpath(@__DIR__, "classify_polyakov.jl"))  # is_circshift_equal(6置換), 生成補助
include(joinpath(@__DIR__, "invariant_core.jl"))

# total を parts 個の非負整数和に分ける全パターン（弱組成）
function weak_compositions(total::Int, parts::Int)
    parts == 1 && return [[total]]
    res = Vector{Vector{Int}}()
    for i in 0:total, tail in weak_compositions(total - i, parts - 1)
        push!(res, vcat(i, tail))
    end
    return res
end

"""
    gen_polyakov_loops(L, Nt) -> Vector{Vector{Symbol}}

巻き数 1 の Polyakov ループ候補をすべて生成（-d 込み）。
"""
function gen_polyakov_loops(L::Int, Nt::Int)
    @assert L >= Nt && iseven(L - Nt) "L>=Nt かつ L-Nt が偶数であること"
    k = div(L - Nt, 2)
    cands = Vector{Vector{Symbol}}()
    for comp in weak_compositions(k, 4)          # (na, nb, nc, q)
        na, nb, nc, q = comp
        bag = Symbol[]
        for (cnt, sym) in ((na,:a),(nb,:b),(nc,:c))
            append!(bag, fill(sym, cnt))
            append!(bag, fill(Symbol("-"*String(sym)), cnt))
        end
        append!(bag, fill(:d, Nt + q))           # +d を Nt+q 個
        append!(bag, fill(Symbol("-d"), q))      # -d を q 個
        for p in remove_adjacent_pairs_circular.(unique_permutations(bag))
            length(p) == L && push!(cands, p)
        end
    end
    return cands
end

# Polyakov winding+1 の実効対称群 = {巡回シフト × a,b,c の6置換 × x,y,z の8反射}。
# Dagger は d 符号を反転し winding を −1 にするので winding+1 集合内では効かない。
# この群は閉じている（真の群）ので、軌道の辞書順最小 = 完全不変量。
# → グループ化だけで分類完了（O(n·48L)、is_circshift_equal の総当たり不要）。
const _REFL8 = Function[identity, xReflection, yReflection, zReflection,
                        xyReflection, xzReflection, yzReflection, xyzReflection]

"""
    canonical_poly(loop) -> Vector{Int}

winding+1 Polyakov 群（6置換×8反射×巡回シフト）の軌道の辞書順最小（整数コード列）。完全不変量。
"""
function canonical_poly(loop::Vector{Symbol})
    n = length(loop); best = nothing
    for i in 1:6, f in _REFL8
        base = to_int(f(xyzt_sym(loop, i)))
        for j in 0:n-1
            c = circshift(base, j)
            if best === nothing || c < best; best = c; end
        end
    end
    return best::Vector{Int}
end

"""
    classify_loop_inv_poly(L, Nt) -> Vector{Vector{Symbol}}

長さ L・時間サイズ Nt・巻き数 1 の独立 Polyakov ループ代表（-d 込み）。
winding 対応生成 → 正準形 `canonical_poly` でグループ化（並列）。
"""
function classify_loop_inv_poly(L::Int, Nt::Int)
    cands = gen_polyakov_loops(L, Nt)
    isempty(cands) && return Vector{Vector{Symbol}}()
    nthr  = nthreads()
    chunk = cld(length(cands), nthr)
    locals = Vector{Dict{Vector{Int},Vector{Symbol}}}(undef, nthr)
    @threads for t in 1:nthr
        d = Dict{Vector{Int},Vector{Symbol}}()
        for idx in ((t-1)*chunk+1):min(t*chunk, length(cands))
            k = canonical_poly(cands[idx])
            haskey(d, k) || (d[k] = cands[idx])
        end
        locals[t] = d
    end
    merged = Dict{Vector{Int},Vector{Symbol}}()
    for t in 1:nthr, (k,v) in locals[t]; haskey(merged,k) || (merged[k]=v); end
    return collect(values(merged))
end
