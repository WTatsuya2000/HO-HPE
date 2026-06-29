# Polyakov ループの「賢い生成 + Booth」（レベル1）
#
# Polyakov(巻き数 m, winding+1 系) の実効対称群は閉じた
#   {巡回シフト × a,b,c の6置換 × x,y,z の8反射}（d は固定・反射なし・Dagger なし）。
# 生成: 空間軸(a,b,c)を Wilson と同様に「初出順・初出+」で構成的に作り、
#       時間軸 d は特別扱い（±d を随時配置, net(d)=m·Nt）。
#       → 軸置換6×反射8=48 の自由度が消える。残る巡回シフト等は最後に Booth canon で dedup。
include(joinpath(@__DIR__, "invariant_polyakov.jl"))   # symmetry_polyakov(xyzt_sym 6, 反射), to_int, weak_compositions
using Base.Threads

# --- Booth 法（最小回転） ---
function _booth(a::Vector{Int})
    n=length(a); s=vcat(a,a); f=fill(-1,2n); k=0
    @inbounds for j in 1:(2n-1)
        i=f[j-k]
        while i!=-1 && s[j+1]!=s[k+i+2]; s[j+1]<s[k+i+2]&&(k=j-i-1); i=f[i+1]; end
        if i==-1 && s[j+1]!=s[k+1]; s[j+1]<s[k+1]&&(k=j); f[j-k+1]=-1 else f[j-k+1]=i+1 end
    end
    k
end
_minrot(a::Vector{Int}) = (k=_booth(a); vcat(@view(a[k+1:end]), @view(a[1:k])))

# Polyakov 軸操作 48 個（a,b,c の6置換 × x,y,z の8符号, d は固定）
function _build_pops()
    ops=Vector{Vector{Int}}()
    for p in permutations(1:3), sg in Iterators.product((1,-1),(1,-1),(1,-1))
        push!(ops, Int[sg[1]*p[1], sg[2]*p[2], sg[3]*p[3], 4])  # d(4)→4 固定
    end
    ops
end
const POPS = _build_pops()   # 48

# 閉群（48 軸操作 × 巡回シフト）の軌道最小 = 完全不変量（Polyakov, winding 固定）
function canon_poly_b(wi::Vector{Int})
    n=length(wi); best=nothing; rel=Vector{Int}(undef,n)
    @inbounds for tab in POPS
        for t in 1:n; s=wi[t]; rel[t]= s>0 ? tab[s] : -tab[-s]; end
        mr=_minrot(rel); (best===nothing||mr<best)&&(best=mr)
    end
    best::Vector{Int}
end

# 巻き数 m の Polyakov ループを構成的に生成し、canon_poly_b で dedup。代表(Vector{Int})と葉数を返す。
function gen_canon_poly(L::Int, Nt::Int, m::Int=1)
    k = div(L - m*Nt, 2)                  # = na+nb+nc+q
    (L >= m*Nt && (L-m*Nt) % 2 == 0) || return (Vector{Vector{Int}}(), 0)
    reps = Dict{Vector{Int},Vector{Int}}(); leaves = Ref(0)
    seq = Int[]; rem_plus=Int[]; rem_minus=Int[]
    dplus = Ref(0); dminus = Ref(0)       # 残り #(+d), #(-d)

    function rec(pool::Vector{Int}, last::Int)
        if length(seq) == L
            seq[end] == -seq[1] && return
            leaves[] += 1
            key = canon_poly_b(seq)
            haskey(reps,key) || (reps[key]=copy(seq))
            return
        end
        nactive = length(rem_plus)
        # (A) 新しい空間軸 a,b,c を導入（初出 = +(nactive+1), label≤3）
        if nactive < 3
            seen=Set{Int}()
            for c in pool
                (c in seen)&&continue; push!(seen,c)
                lab=nactive+1
                push!(rem_plus,c-1); push!(rem_minus,c); push!(seq,lab)
                np=copy(pool); deleteat!(np,findfirst(==(c),np)); rec(np,lab)
                pop!(seq); pop!(rem_plus); pop!(rem_minus)
            end
        end
        # (B) 既出 空間軸 μ を ±
        @inbounds for μ in 1:nactive, sgn in (1,-1)
            d=sgn*μ; r= sgn>0 ? rem_plus[μ] : rem_minus[μ]
            (r>0 && d != -last) || continue
            sgn>0 ? (rem_plus[μ]-=1) : (rem_minus[μ]-=1)
            push!(seq,d); rec(pool,d); pop!(seq)
            sgn>0 ? (rem_plus[μ]+=1) : (rem_minus[μ]+=1)
        end
        # (C) +d
        if dplus[]>0 && last != -4
            dplus[]-=1; push!(seq,4); rec(pool,4); pop!(seq); dplus[]+=1
        end
        # (D) -d
        if dminus[]>0 && last != 4
            dminus[]-=1; push!(seq,-4); rec(pool,-4); pop!(seq); dminus[]+=1
        end
    end

    # (na,nb,nc,q) を列挙: 空間 s=na+nb+nc, q=k-s, #(+d)=m·Nt+q, #(-d)=q
    for q in 0:k
        s = k - q
        for comp in weak_compositions(s, 3)     # (na,nb,nc) 順序付き（重複は dedup で吸収）
            spool = sort(filter(>(0), comp))     # 非ゼロ空間 counts
            empty!(seq); empty!(rem_plus); empty!(rem_minus)
            dplus[] = m*Nt + q; dminus[] = q
            rec(spool, 0)
        end
    end
    return collect(values(reps)), leaves[]
end

const _DECP = Dict{Int,Symbol}(1=>:a,2=>:b,3=>:c,4=>:d,-1=>Symbol("-a"),-2=>Symbol("-b"),-3=>Symbol("-c"),-4=>Symbol("-d"))
tosym_p(v::Vector{Int}) = Symbol[_DECP[x] for x in v]

"巻き数 m=1 の独立 Polyakov ループ代表（賢い生成版）。classify_loop_inv_poly と同値。"
function classify_loop_gen_poly(L::Int, Nt::Int)
    reps, _ = gen_canon_poly(L, Nt, 1)
    return [tosym_p(v) for v in reps]
end
