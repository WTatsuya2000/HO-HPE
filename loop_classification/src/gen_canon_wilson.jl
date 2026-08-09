# Wilson ループの「賢い生成」（レベル1: 標準形を直接構成）
#
# 現方式の無駄: 全順列(L=14 最密で 13.6億)を生成してから 768L 倍の重複を dedup。
# ここでは対称性を生成に織り込み、各クラスの代表(に近いもの)だけを直接作る:
#   構成的バックトラッキングで以下を満たす列のみ生成
#     ・軸ラベルは初出順（新しい軸は必ず次の番号 1→2→3→4）
#     ・各軸の初出は + のみ（最初の向きを正に統一）
#     ・非後退（隣接で逆向き禁止, 円としても）
#     ・各軸の +/- 本数が type に一致（閉路）
#   → H(4)=軸置換24×符号16 の自由度が消え、先頭=+1 で巡回シフトも大幅固定。
#   残る自由度（回転の一部・向き反転 Dagger・同数軸の入替）は最後に canon_w(Booth) で dedup。
#
# 生成数 ≈ #クラス × (小さな残余因子) になり、十数億 → 数百万オーダーに激減する。
include(joinpath(@__DIR__, "fast_wilson.jl"))   # canon_w, minrot, WOPS, toint, WCODE
using Base.Threads

# type_counts: (Na,Nb,Nc,Nd)（0 含む）。代表ループ(Vector{Int})の配列と、生成した葉の数を返す。
function gen_canon_wilson(type_counts::Vector{Int})
    cc = sort(filter(>(0), type_counts))      # 非ゼロ counts（昇順）
    L  = 2 * sum(cc)
    reps = Dict{Vector{Int},Vector{Int}}()
    seq = Int[]; rem_plus = Int[]; rem_minus = Int[]
    leaves = Ref(0)

    function rec(pool::Vector{Int}, last::Int)
        n = length(seq)
        if n == L
            seq[end] == -seq[1] && return          # 円として後退ならNG
            leaves[] += 1
            k = canon_w(seq)
            haskey(reps, k) || (reps[k] = copy(seq))
            return
        end
        nactive = length(rem_plus)
        # (A) 新しい軸を導入（初出 = +(nactive+1)）。pool の distinct 値のみ試す
        seen = Set{Int}()
        for c in pool
            (c in seen) && continue; push!(seen, c)
            lab = nactive + 1
            push!(rem_plus, c-1); push!(rem_minus, c); push!(seq, lab)
            np = copy(pool); deleteat!(np, findfirst(==(c), np))
            rec(np, lab)
            pop!(seq); pop!(rem_plus); pop!(rem_minus)
        end
        # (B) 既出の軸 μ を ± で配置
        @inbounds for μ in 1:nactive, s in (1, -1)
            d = s*μ
            r = s > 0 ? rem_plus[μ] : rem_minus[μ]
            (r > 0 && d != -last) || continue
            s > 0 ? (rem_plus[μ] -= 1) : (rem_minus[μ] -= 1)
            push!(seq, d); rec(pool, d); pop!(seq)
            s > 0 ? (rem_plus[μ] += 1) : (rem_minus[μ] += 1)
        end
    end

    rec(cc, 0)                                     # last=0（番兵: -0=0 は実方向と一致しない）
    return collect(values(reps)), leaves[]
end

# 既存インタフェース互換: int 代表を Symbol に戻す
const _DECW = Dict{Int,Symbol}(1=>:a,2=>:b,3=>:c,4=>:d,-1=>Symbol("-a"),-2=>Symbol("-b"),-3=>Symbol("-c"),-4=>Symbol("-d"))
tosym(v::Vector{Int}) = Symbol[_DECW[x] for x in v]

"長さ L・type の独立 Wilson ループ代表（賢い生成版）。classify_loop_fastw と同値。"
function classify_loop_gen(L::Int, type_N::Int)
    combos = count_combinations(L)
    (1 <= type_N <= length(combos)) || error("type_N=$type_N 範囲外")
    reps, _ = gen_canon_wilson(combos[type_N])
    return [tosym(v) for v in reps]
end
