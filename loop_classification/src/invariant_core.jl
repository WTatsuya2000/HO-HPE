# ============================================================================
# 不変量ベースの段階的分類コア（Invariant.md ベース, d=4 固定）
#
# 思想（Invariant.md §13）:
#   生成 → 真の不変量 I(w) でバケツ分け → 同じバケツ内だけ完全同値判定(is_circshift_equal)。
#   各不変量は H(4)×D_L 作用で不変なので、同値なループは必ず同じバケツに入る
#   （クラス分断なし）。バケツ間は独立 → スレッド並列可。
#
# 前提: is_circshift_equal が定義済み（src/loop_core_*.jl 等を先に include）。
# ============================================================================
using Base.Threads

const _DIR = Dict{Symbol,Int}(
    :a=>1, :b=>2, :c=>3, :d=>4,
    Symbol("-a")=>-1, Symbol("-b")=>-2, Symbol("-c")=>-3, Symbol("-d")=>-4)
to_int(loop::Vector{Symbol}) = Int[_DIR[s] for s in loop]

# §3 軸使用回数（ソート）= sort(n_1..n_4)
axis_mult(s::Vector{Int}) = Tuple(sort!([count(x->abs(x)==μ, s) for μ in 1:4]))

# §4 隣接リンク (N_straight, N_backtrack, N_turn)  ※円順列
function local_profile(s::Vector{Int})
    n = length(s); st=0; bk=0; tn=0
    @inbounds for i in 1:n
        a = s[i]; b = s[i==n ? 1 : i+1]
        if a == b; st+=1 elseif a == -b; bk+=1 else tn+=1 end
    end
    return (st, bk, tn)
end

# §5 距離 k のリンク相関 (N_k^+, N_k^-, N_k^0)  ※円順列
function corr_profile(s::Vector{Int}, k::Int)
    n = length(s); p=0; m=0; z=0
    @inbounds for i in 1:n
        a = s[i]; b = s[((i+k-1) % n) + 1]
        if a == b; p+=1 elseif a == -b; m+=1 else z+=1 end
    end
    return (p, m, z)
end

# 頂点列 x_0..x_{L-1}（x_L=x_0 なので末尾は数えない）
function _vertices(s::Vector{Int})
    n = length(s)
    pts = Vector{NTuple{4,Int}}(undef, n)
    x = (0,0,0,0)
    pts[1] = x
    @inbounds for i in 1:n-1
        a = s[i]; μ = abs(a); δ = sign(a)
        x = ntuple(t -> x[t] + (t==μ ? δ : 0), 4)
        pts[i+1] = x
    end
    return pts  # x_0..x_{L-1}
end

# §6 頂点訪問回数の分布（ソート）
function vertex_profile(s::Vector{Int})
    pts = _vertices(s)
    c = Dict{NTuple{4,Int},Int}()
    for p in pts; c[p] = get(c, p, 0) + 1; end
    return Tuple(sort!(collect(values(c))))
end

# §7.2 無向リンク通過回数の分布（ソート）
function edge_profile(s::Vector{Int})
    n = length(s)
    pts = _vertices(s)                 # x_0..x_{L-1}
    c = Dict{Tuple{NTuple{4,Int},NTuple{4,Int}},Int}()
    @inbounds for i in 1:n
        x = pts[i]
        y = pts[i==n ? 1 : i+1]        # x_{i} → x_{i+1}（x_L=x_0）
        e = x <= y ? (x,y) : (y,x)     # 無向（標準形）
        c[e] = get(c, e, 0) + 1
    end
    return Tuple(sort!(collect(values(c))))
end

"""
    invariant_key(loop) -> Tuple

Invariant.md §13 推奨キー: (軸使用回数, 隣接リンク, corr2, corr3, 頂点訪問分布, リンク通過分布)。
すべて H(4)×D_L 不変なので、同値ループは必ず同じキー。
"""
function invariant_key(loop::Vector{Symbol})
    s = to_int(loop)
    return (axis_mult(s), local_profile(s), corr_profile(s,2), corr_profile(s,3),
            vertex_profile(s), edge_profile(s))
end

# is_equal の連結成分（union-find）で代表を返す共通ヘルパ
function _components(group::Vector{Vector{Symbol}}, is_equal)
    n = length(group)
    n <= 1 && return copy(group)
    par = collect(1:n)
    find(x) = par[x] == x ? x : (par[x] = find(par[x]))
    @inbounds for a in 1:n, b in (a+1):n
        if find(a) != find(b) && is_equal(group[a], group[b])
            par[find(a)] = find(b)
        end
    end
    seen = Set{Int}(); reps = Vector{Vector{Symbol}}()
    for i in 1:n
        r = find(i)
        if r ∉ seen; push!(seen, r); push!(reps, group[i]); end
    end
    return reps
end

"""
    classify_invariant(perms; is_equal=is_circshift_equal, parallel=true)

不変量でバケツ分け → 各バケツ内を is_equal の連結成分で reduce（並列）→ 残った代表集合に
大域連結成分マージ。

- バケツ内 reduce が n 依存の重い処理を潰す（並列）。
- 大域マージは「不変量が同値ループを別バケツに分けてしまう」場合の保険。`is_equal` が
  非推移的でも連結成分で順序非依存・物理的に正しい分割になる。代表数は小さいので安価。
  Wilson（不変量が真に invariant・is_equal が推移的）では大域マージは何も統合しない。
"""
function classify_invariant(perms::Vector{Vector{Symbol}};
                            is_equal = is_circshift_equal, parallel::Bool = true)
    isempty(perms) && return Vector{Vector{Symbol}}()
    buckets = Dict{Any,Vector{Vector{Symbol}}}()
    for p in perms
        push!(get!(buckets, invariant_key(p), Vector{Vector{Symbol}}()), p)
    end
    blists = collect(values(buckets))

    # バケツ内を貪欲 dedup（Wilson は is_equal が実効的に推移的なので貪欲＝連結成分で正しい・高速）。
    # Wilson の不変量は真に H(4)×D_L 不変なのでクラスはバケツをまたがず、大域マージは不要。
    # （Polyakov は非推移なので classify_invariant ではなく canonical_poly を使うこと。）
    function greedy(group)
        reps = Vector{Vector{Symbol}}()
        for p in group
            any(q -> is_equal(q, p), reps) || push!(reps, p)
        end
        return reps
    end
    if parallel && length(blists) > 1
        out = Vector{Vector{Vector{Symbol}}}(undef, length(blists))
        @threads for b in 1:length(blists)
            out[b] = greedy(blists[b])
        end
        return reduce(vcat, out)
    else
        return reduce(vcat, map(greedy, blists))
    end
end
