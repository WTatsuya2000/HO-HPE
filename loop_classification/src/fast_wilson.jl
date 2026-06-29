# Wilson 分類の高速版: 閉群の正準形（完全不変量）でグループ化（O(n)・並列）
#
# Wilson の対称群 = 符号付き置換 H(4)=(Z2)^4⋊S4（位数 2^4·4!=384） × 二面体 D_L（シフト+反転, 2L）
#                  = 768·L 個。これは閉じた真の群なので「軌道の辞書順最小 = 完全不変量」。
# is_circshift_equal の総当たり（O(バケツ×クラス)）を、正準キーの Dict 集約（O(n)）に置き換える。
include(joinpath(@__DIR__, "loop_core_wilson.jl"))
using Base.Threads, Combinatorics

const WCODE = Dict{Symbol,Int}(:a=>1,:b=>2,:c=>3,:d=>4,
    Symbol("-a")=>-1,Symbol("-b")=>-2,Symbol("-c")=>-3,Symbol("-d")=>-4)
toint(w::Vector{Symbol}) = Int[WCODE[s] for s in w]

# 384 個の符号付き置換（各 = 長さ4: +μ の行き先（符号付き整数））
function build_wops()
    ops = Vector{Vector{Int}}()
    for p in permutations(1:4), sg in Iterators.product((1,-1),(1,-1),(1,-1),(1,-1))
        push!(ops, Int[sg[μ]*p[μ] for μ in 1:4])
    end
    return ops          # length 384
end
const WOPS = build_wops()

# Booth 法: 系列 a の辞書順最小回転の開始位置（0-based）を O(n) で返す
function booth(a::Vector{Int})
    n = length(a); s = vcat(a, a); f = fill(-1, 2n); k = 0
    @inbounds for j in 1:(2n-1)
        i = f[j-k]
        while i != -1 && s[j+1] != s[k+i+2]
            s[j+1] < s[k+i+2] && (k = j-i-1); i = f[i+1]
        end
        if i == -1 && s[j+1] != s[k+1]
            s[j+1] < s[k+1] && (k = j); f[j-k+1] = -1
        else
            f[j-k+1] = i+1
        end
    end
    return k
end
minrot(a::Vector{Int}) = (k = booth(a); vcat(@view(a[k+1:end]), @view(a[1:k])))

# 閉群（384 軸操作 × {恒等,反転} × 巡回シフト）の軌道最小 = 完全不変量
# 巡回シフトの最小化は Booth 法 O(n)（旧: 全 n 回転を都度 circshift で確保していた O(n²)）。
function canon_w(wi::Vector{Int})
    n = length(wi); best = nothing
    rel = Vector{Int}(undef, n)
    @inbounds for tab in WOPS
        for t in 1:n
            s = wi[t]; rel[t] = s > 0 ? tab[s] : -tab[-s]
        end
        for seq in (rel, reverse(rel))
            mr = minrot(seq)
            if best === nothing || mr < best; best = mr; end
        end
    end
    return best::Vector{Int}
end

# 高速 多重度 M:
#   軌道 P = {f(xyzt_sym(w,i)) : f∈transforms, i∈1:24} を Set(整数符号化)で構築（O(|P|)）→
#   {巡回シフト, Dagger∘巡回シフト} の正準形（Booth 最小回転）で distinct 数を数える（O(|P|·L)）。
#   旧来の O(|P|²)（線形 `g in P` ＋ exclude_same の二重ループ）を置き換え。同じ M を返す。
#   ※ transforms / xyzt_sym（gamma_trace_wilson.jl の 24軸置換版）が読み込まれている前提。
_dagi(g::Vector{Int}) = Int[-x for x in Iterators.reverse(g)]
function wilson_mult(w::Vector{Symbol})
    Pset = Set{Vector{Int}}()
    for f in transforms, i in 1:24
        push!(Pset, toint(f(xyzt_sym(w, i))))
    end
    cset = Set{Vector{Int}}()
    for g in Pset
        push!(cset, min(minrot(g), minrot(_dagi(g))))
    end
    return length(cset)
end

"""
    classify_loop_fastw(L, type_N) -> Vector{Vector{Symbol}}

長さ L・指定 type の独立 Wilson ループ代表（閉群正準形でグループ化）。`classify_loop_inv` と同値・高速。
"""
function classify_loop_fastw(L::Int, type_N::Int)
    combos = count_combinations(L)
    (1 <= type_N <= length(combos)) || error("type_N=$type_N 範囲外")
    target = combos[type_N]
    bag1 = nothing
    for arr in add_minus.(all_combinations_multi_types(L))
        if sort([count(x->x==s, arr) for s in (:a,:b,:c,:d)]) == target
            bag1 = arr; break
        end
    end
    bag1 === nothing && return Vector{Vector{Symbol}}()
    perms2 = [p for p in remove_adjacent_pairs_circular.(unique_permutations(bag1)) if length(p)==L]
    isempty(perms2) && return Vector{Vector{Symbol}}()
    ints = map(toint, perms2)
    nthr = nthreads(); chunk = cld(length(ints), nthr)
    locals = Vector{Dict{Vector{Int},Vector{Symbol}}}(undef, nthr)
    @threads for t in 1:nthr
        d = Dict{Vector{Int},Vector{Symbol}}()
        for idx in ((t-1)*chunk+1):min(t*chunk, length(ints))
            k = canon_w(ints[idx])
            haskey(d, k) || (d[k] = perms2[idx])
        end
        locals[t] = d
    end
    merged = Dict{Vector{Int},Vector{Symbol}}()
    for t in 1:nthr, (k,v) in locals[t]; haskey(merged,k) || (merged[k]=v); end
    return collect(values(merged))
end
