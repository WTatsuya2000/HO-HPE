# Polyakov ループを巻き数 m ごとに分類し、D・M・寄与を records/*.csv に追記
# 実行: julia -t 8 --project=. scripts/record_polyakov.jl   （record_wilson.jl の後）
# 分類は src/gen_canon_polyakov.jl の「賢い生成 + Booth」（gen_canon_poly）を使用（旧 classify_w より高速）。
cd(joinpath(@__DIR__, ".."))
using Base.Threads
include(joinpath(@__DIR__, "..", "src", "gamma_trace_polyakov.jl"))
include(joinpath(@__DIR__, "..", "src", "gen_canon_polyakov.jl"))  # gen_canon_poly, tosym_p, canon_poly_b

const XYZT = Dict{Symbol,String}(:a=>"x",:b=>"y",:c=>"z",:d=>"t",
    Symbol("-a")=>"-x",Symbol("-b")=>"-y",Symbol("-c")=>"-z",Symbol("-d")=>"-t")
loopstr(w) = join((XYZT[s] for s in w), " ")  # 方向は xyzt 表記

# 巻き数 m の独立 Polyakov ループ代表（賢い生成版）
function classify_w(L, Nt, m)
    reps, _ = gen_canon_poly(L, Nt, m)
    return [tosym_p(v) for v in reps]
end
function dml(w)
    P=Vector{Vector{Symbol}}(); for f in transforms, i in 1:6; g=f(xyzt_sym(w,i)); (g in P)||push!(P,g); end
    M=length(exclude_same(P)); D=real(character2gamma(w))
    LLL=1; for i in 1:length(w)-1; w==circshift(w,i)&&(LLL+=1); end
    (D,M,LLL)
end

# 論文 Table2（巻き1のみ。m>=2 は (-1)^(m-1) L1^0(mNt,n) と照合可能）
paper1 = Dict((4,4)=>48,(4,6)=>1728,(4,8)=>45792,(4,10)=>645120,(4,12)=>-26224128,(4,14)=>-3201067008,
              (6,6)=>128,(6,8)=>11520,(6,10)=>716544,(6,12)=>35891712,
              (8,8)=>384,(8,10)=>64512,(8,12)=>6842880,
              (10,10)=>1228.8,(10,12)=>331776,(12,12)=>4096)

mkpath("records")
# Nt ごとに別ファイル。各ファイル内で通し番号（id）を独立採番。
LIO = Dict{Int,IOStream}(); SIO = Dict{Int,IOStream}()
GID = Dict{Int,Int}(); SID = Dict{Int,Int}()
function files(Nt)
    if !haskey(LIO, Nt)
        l = open("records/polyakov_Nt$(Nt)_loops.csv", "w")
        s = open("records/polyakov_Nt$(Nt)_summary.csv", "w")
        println(l, "id,kind,Nt,n,m,index,loop,D,M,LLL,contrib")
        println(s, "id,kind,Nt,n,m,classes,sum_contrib,coeff,value,paper")
        LIO[Nt]=l; SIO[Nt]=s; GID[Nt]=0; SID[Nt]=0
    end
    (LIO[Nt], SIO[Nt])
end
classify_w(4,4,1)  # warmup

grid = [(4,4),(6,4),(8,4),(10,4),(12,4),(14,4),(6,6),(8,6),(10,6),(12,6),(8,8),(10,8),(12,8),(10,10),(12,10),(12,12)]
for (n,Nt) in grid
    for m in 1:div(n,Nt)
        (n - m*Nt >= 0 && (n - m*Nt) % 2 == 0) || continue
        cls = classify_w(n, Nt, m)
        isempty(cls) && continue
        lio, sio = files(Nt)
        sub = 0.0; idx = 0
        for w in cls
            D, M, LLL = dml(w); c = M*D/LLL; idx += 1; GID[Nt] += 1
            println(lio, "$(GID[Nt]),P,$Nt,$n,$m,$idx,\"$(loopstr(w))\",$(Int(D)),$M,$LLL,$c")
            sub += c
        end
        flush(lio)
        coeff = (-1)^(m-1) * 6
        val = coeff * sub
        p = (m == 1) ? get(paper1, (Nt, n), "") : ""
        SID[Nt] += 1
        println(sio, "$(SID[Nt]),P,$Nt,$n,$m,$idx,$sub,$coeff,$(val),$p"); flush(sio)
        println(stderr, "P Nt=$Nt n=$n m=$m: $idx classes, L_$m^0=$(val)")
    end
end
for io in values(LIO); close(io); end
for io in values(SIO); close(io); end
