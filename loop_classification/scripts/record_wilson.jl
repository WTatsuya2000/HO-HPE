# Wilson ループを分類し、D・M・寄与を records/*.csv に記録（m=0）
# 実行: julia -t 8 --project=. scripts/record_wilson.jl
cd(joinpath(@__DIR__, ".."))
using Base.Threads
include(joinpath(@__DIR__, "..", "src", "gamma_trace_wilson.jl"))
include(joinpath(@__DIR__, "..", "src", "fast_wilson.jl"))  # classify_loop_fastw（正準形+Booth）, wilson_mult

const XYZT = Dict{Symbol,String}(:a=>"x",:b=>"y",:c=>"z",:d=>"t",
    Symbol("-a")=>"-x",Symbol("-b")=>"-y",Symbol("-c")=>"-z",Symbol("-d")=>"-t")
loopstr(w) = join((XYZT[s] for s in w), " ")  # 方向は xyzt 表記

# 各代表ループの (D, M, LLL)  ※ M は wilson_mult（Set+Booth, 旧 exclude_same と同値・高速）
function dml(w)
    M = wilson_mult(w); D = real(character2gamma(w))
    LLL = 1; for i in 1:length(w)-1; w == circshift(w,i) && (LLL+=1); end
    (D, M, LLL)
end

paper = Dict(4=>288, 6=>8448, 8=>245952, 10=>7372800)
mkpath("records")
lio = open("records/wilson_loops.csv", "w")
sio = open("records/wilson_summary.csv", "w")
println(lio, "id,kind,Nt,n,m,index,loop,D,M,LLL,contrib")
println(sio, "id,kind,Nt,n,m,classes,sum_contrib,coeff,value,paper")

gid = 0; sid = 0   # 通し番号（loops / summary）
classify_loop_fastw(4,1)  # warmup
for n in (4, 6, 8, 10)
    sub = 0.0; idx = 0
    for t in 1:length(count_combinations(n)), w in classify_loop_fastw(n, t)
        D, M, LLL = dml(w); c = M*D/LLL; idx += 1; global gid += 1
        println(lio, "$gid,W,0,$n,0,$idx,\"$(loopstr(w))\",$(Int(D)),$M,$LLL,$c")
        sub += c
    end
    flush(lio)
    val = -6*sub; p = get(paper, n, ""); global sid += 1
    println(sio, "$sid,W,0,$n,0,$idx,$sub,-6,$(val),$p"); flush(sio)
    println(stderr, "W n=$n: $idx classes, W^0=$(val)")
end
close(lio); close(sio)
