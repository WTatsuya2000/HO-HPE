using Combinatorics  # combinations を使うために必要です

"""
    all_combinations_multi_types(l::Int)

長さ l は偶数であることを要求します。
プール ["a","b","c"] から 2 種類以上を選び，
合計 k = l÷2 個を各種類に１個以上割り当てた
ベクトルのリストを返します。
"""
function all_combinations_multi_types(l::Int)
    @assert iseven(l) "Error: l must be even (got l=$l)."
    k = div(l,2)

    pool = [:a,:b,:c,:d]
    result = Vector{Vector{Symbol}}()

    # プールから 2 種類以上を選ぶ
    for m in 2:length(pool)
        for types in combinations(pool, m)
            # 選んだ m 種類に対して，合計 k かつ各 1 以上となる分割（組成）を列挙
            for parts in compositions(k, m)
                # parts[i] 個だけ types[i] を並べたベクトルを連結
                push!(result, vcat([ fill(types[i], parts[i]) for i in 1:m ]...))
            end
        end
    end

    return result
end

"""
    compositions(total::Int, m::Int)

total を m 個の正の整数の和として表すすべての組成 (ordered partitions) を返す。
例えば compositions(4,2) == [[1,3],[2,2],[3,1]] となる。
"""
function compositions(total::Int, m::Int)
    if m == 1
        return [[total]]
    end
    result = Vector{Vector{Int}}()
    # i を最初の要素としたとき，残り total-i を m-1 個で分割
    for i in 1:(total - (m-1))
        for tail in compositions(total - i, m-1)
            push!(result, vcat(i, tail))
        end
    end
    return result
end


"""
    add_minus(strs::Vector{String}) -> Vector{String}

文字列のベクトル `strs` を受け取り、
各要素 s について s と "-"*s の両方を含む新しいベクトルを返す。
"""
function add_minus(strs::Vector{Symbol})
    result = Vector{Symbol}(undef, 2*length(strs))
    for (i, s) in enumerate(strs)
        result[2*i-1] = s
        result[2*i]   = Symbol("-" * String(s))
    end
    return result
end
"""
    unique_permutations(arr::Vector{T}) where T

"""
function unique_permutations(arr::Vector{T}) where T
    # 入力をソートしておく（重複スキップのため）
    sorted = sort(arr)
    n = length(arr)
    used = falses(n)
    current = Vector{T}()
    result = Vector{Vector{T}}()

    function backtrack()
        if length(current) == n
            push!(result, copy(current))
            return
        end
        for i in 1:n
            # すでに使った要素はスキップ
            if used[i]
                continue
            end
            # 同じ値の要素の並び替えで重複を避ける
            if i > 1 && sorted[i] == sorted[i-1] && !used[i-1]
                continue
            end
            # 選択
            push!(current, sorted[i])
            used[i] = true
            backtrack()
            # 後片付け
            used[i] = false
            pop!(current)
        end
    end

    backtrack()
    return result
end

function remove_adjacent_pairs_circular(arr::Vector{Symbol})
    # キャンセル可能なペアがなくなるまで繰り返す
    while true
        n = length(arr)
        # 要素数が1以下なら終了
        if n < 2
            break
        end
        found_pair = false
        # 1からnまでの各インデックスについて、次の要素（円順列なので、nの次は1）と比較
        for i in 1:n
            j = (i % n) + 1  # 次のインデックス（iがnならjは1）
            if (arr[i] == :a && arr[j] == Symbol("-a")) || (arr[i] == Symbol("-a") && arr[j] == :a) ||
               (arr[i] == :b && arr[j] == Symbol("-b")) || (arr[i] == Symbol("-b") && arr[j] == :b) ||
               (arr[i] == :c && arr[j] == Symbol("-c")) || (arr[i] == Symbol("-c") && arr[j] == :c) ||
               (arr[i] == :d && arr[j] == Symbol("-d")) || (arr[i] == Symbol("-d") && arr[j] == :d) 
                # 2つの要素を削除する。インデックスがずれないように大きい方から削除
                if i < j
                    deleteat!(arr, j)
                    deleteat!(arr, i)
                else
                    # i > j の場合、i は末尾、j は先頭に相当
                    deleteat!(arr, i)
                    deleteat!(arr, j)
                end
                found_pair = true
                break  # ペアを削除したら再度先頭からチェック
            end
        end
        # どのペアも削除されなかったらループ終了
        if !found_pair
            break
        end
    end
    return arr
end

# 対称性演算子についての関数
include("symmetry_polyakov.jl")

function is_circshift_equal(arr1::Vector{Symbol}, arr2::Vector{Symbol})
    counts_abc_arr1 = [
        count(x -> x == :a, arr1),
        count(x -> x == :b, arr1),
        count(x -> x == :c, arr1),
        count(x -> x == :d, arr1),
    ]
    counts_abc_arr2 = [
        count(x -> x == :a, arr2),
        count(x -> x == :b, arr2),
        count(x -> x == :c, arr2),
        count(x -> x == :d, arr2)
    ]
    
    if sort(counts_abc_arr1) != sort(counts_abc_arr2)
        return false
    end

    # 長さが異なればすぐに false
    if length(arr1) != length(arr2)
    end
    # もしどちらも空ならば同じとみなす
    if isempty(arr1) && isempty(arr2)
        return true
    end
    n = length(arr1)
    transforms = [
        identity,
        xReflection,
        yReflection,
        zReflection,
        # tReflection,
        xyReflection,
        xzReflection,
        yzReflection,
        # xtReflection,
        # ytReflection,
        # ztReflection,
        xyzReflection,
        # yztReflection,
        # xztReflection,
        # xytReflection,
        # xyztReflection,
        
        x -> Dagger(identity(x)),
        x -> Dagger(xReflection(x)),
        x -> Dagger(yReflection(x)),
        x -> Dagger(zReflection(x)),
        # x -> Dagger(tReflection(x)),
        x -> Dagger(xyReflection(x)),
        x -> Dagger(xzReflection(x)),
        x -> Dagger(yzReflection(x)),
        # x -> Dagger(xtReflection(x)),
        # x -> Dagger(ytReflection(x)),
        # x -> Dagger(ztReflection(x)),
        x -> Dagger(xyzReflection(x)),
        # x -> Dagger(yztReflection(x)),
        # x -> Dagger(xztReflection(x)),
        # x -> Dagger(xytReflection(x)),
        # x -> Dagger(xyztReflection(x)),
    ]

    for f in transforms
        for i = 1:6
            deformed = f(xyzt_sym(arr2,i))
            
            counts_abc_deformed = [
                count(x -> x == :a, deformed),
                count(x -> x == :b, deformed),
                count(x -> x == :c, deformed),
                count(x -> x == :d, deformed)
            ]
            if counts_abc_arr1 == counts_abc_deformed
                for j in 0:n-1
                    shifted = circshift(arr1, j)
                    if shifted == deformed
                        return true
                    end
                end
            end
        end
    end
    return false
end

"""
    count_combinations(l::Int) -> Vector{NTuple{4,Int}}
    長さ l（偶数）を指定し，
    k = l÷2 となる (Na,Nb,Nc,Nd) の全組み合わせを返す。
    ただし 0 が３つ以上含まれる組は除外（0 の個数 ≤ 2）。
"""
function count_combinations(l::Int)
    @assert iseven(l) "Error: l must be even (got l=$l)."
    k = div(l, 2)
    result = Vector{Vector{Int}}()

    for Na in 0:k, Nb in 0:k, Nc in 0:k
        Nd = k - (Na + Nb + Nc)
        if Nd < 0
            continue
        end
        # 0 の個数を数える
        zero_count = (Na==0) + (Nb==0) + (Nc==0) + (Nd==0)
        if zero_count <= 2 && Na <= Nb && Nb <= Nc && Nc <= Nd
            push!(result, [Na, Nb, Nc, Nd])
        end
    end

    return result
end

"""
    classify_loop(L::Int, type_N::Int)
    L の長さの配列を受け取り、指定されたタイプ type_N に基づいて分類します。
    並列化のために、タイプごとに配列を処理します。
"""
function classify_loop(L,type_N::Int)
    # L = 8
    CL = add_minus.(all_combinations_multi_types(L))
    NCL = length(CL)
    CL_set = []
    for arr1 in CL
        counts_abc_arr1 = sort([
            count(x -> x == :a, arr1),
            count(x -> x == :b, arr1),
            count(x -> x == :c, arr1),
            count(x -> x == :d, arr1),
        ])
        push!(CL_set, [arr1,counts_abc_arr1])
    end
    # CL_set

    Type_loop = Dict()
    for i in 1:length(count_combinations(L))
        Type_loop[Symbol(i)] = []
    end

    for i in 1:length(count_combinations(L)), j in 1:length(CL_set)
        if CL_set[j][2] == count_combinations(L)[i]
            push!(Type_loop[Symbol(i)], CL_set[j][1])
        end
    end
    perms_result = []
    is_new = true
    c = 0
    dir_check("./log/L$L")
    log = open("./log/L$L/log$(L)_$type_N.txt","w")          
        println(log, "start", )
    close(log)
    
    Num_all = 1

    @time for p1 in Type_loop[Symbol(type_N::Int)]
        input = p1
        # すべての並び方についての処理を行う
        perms = unique_permutations(input)
        # 隣接するペアを削除
        perms_new = remove_adjacent_pairs_circular.(perms)
        perms_new2 = []
        for p in perms_new
            if length(p) == L # 長さが n でないものは除外
                push!(perms_new2,p)
            end
        end
        
        c==0 ? Num_all = length(perms_new2) * length(Type_loop[Symbol(type_N::Int)]) : nothing
        # ここでは同じものから作られるもののみを判定している。
        # 結果格納用
        perms_new3 = Vector{Vector{Symbol}}()
        for p in perms_new2
            log = open("log$(L)_$type_N.txt","a")            
                c==0 ? println(log,Num_all)  : nothing
                c+=1
                # c % 100 == 0  ? println(log, "step:$c/$Num_all", ) : nothing
                c % div(Num_all,1000)  == 0  ? println(log, "step:$(c)/$(Num_all),$(c/Num_all*100)%,", ) : nothing
            close(log)
            # p が既に何かと等しいかどうかフラグ
            is_new = true
            for q in perms_new3
                flag = is_circshift_equal(q, p)
                # display(flag)
                if flag
                    is_new = false
                    break
                end
            end
            # どの既存要素とも等しくなければ追加
            if is_new
                push!(perms_new3, p)
            end
        end
        # 独立なものだけをリスト化
        for p2 in perms_new3
            push!(perms_result,p2)
        end
    end
    # 重複排除後の結果格納用を空で初期化
    perms_result2 = Vector{Vector{Symbol}}()

    for p in perms_result
        # p が既に登録済みかどうかを判定するフラグ
        is_new = true
        for q in perms_result2
            if is_circshift_equal(q, p)
                is_new = false
                break
            end
        end
        # どの登録済み要素とも同一でなければ追加
        if is_new
            push!(perms_result2, p)
        end
    end
    return perms_result2
end

using Base.Threads

# ─────────────────────────────────────────────────────────────────────────────
# グローバルなログ同期用ロック
const log_lock = ReentrantLock()

# ログ書き込み関数
function write_log(log_path::String, msg::String)
    lock(log_lock) do
        open(log_path, "a") do io
            println(io, msg)
        end
    end
end

# ── 並列ユニーク化＋ログ出力 ──────────────────────────────────────────────────
function parallel_unique(perms::Vector{Vector{Symbol}}, L::Int, log_path::String)
    n_threads = nthreads()
    @show n_threads 
    n = length(perms)   # 完全分類のため全件処理（旧: min(10000, length(perms)) — L>=10 で取りこぼす危険があった）
    chunk_size = ceil(Int, n/n_threads)
    local_results = Vector{Vector{Vector{Symbol}}}(undef, n_threads)

    write_log(log_path, ">> [parallel_unique] start: total=$n, threads=$n_threads")

    @threads for t in 1:n_threads
        start_idx = (t-1)*chunk_size + 1
        end_idx   = min(t*chunk_size, n)
        uniq = Vector{Vector{Symbol}}()
        for i in start_idx:end_idx
            p = perms[i]
            # ローカル uniq のみと比較
            is_new = true
            for q in uniq
                if is_circshift_equal(q, p)
                    is_new = false
                    break
                end
            end
            is_new && push!(uniq, p)
        end
        local_results[t] = uniq
        write_log(log_path, ">> Thread $t processed indices [$start_idx:$end_idx], uniq=$(length(uniq))")
    end

    # マージフェーズ
    merged = Vector{Vector{Symbol}}()
    for uniq in local_results, p in uniq
        if all(!is_circshift_equal(q, p) for q in merged)
            push!(merged, p)
        end
    end
    write_log(log_path, ">> [parallel_unique] merge complete, merged=$(length(merged))")

    return merged
end

# ── classify_loop 本体（ログ統合版）────────────────────────────────────────
function classify_loop(L::Int, type_N::Int)
    # ログファイル名
    log_path = "./log/log_$(L)_$(type_N).txt"
    # ログファイル初期化
    open(log_path, "w") do io
        println(io, "=== classify_loop start: L=$(L), type_N=$(type_N) ===")
    end

    # 組み合わせ→タイプ分類
    CL = add_minus.(all_combinations_multi_types(L))
    CL_set = [ [arr, sort([ count(x->x==s,arr) for s in (:a,:b,:c,:d) ])] for arr in CL ]
    combos = count_combinations(L)

    Type_loop = Dict( Symbol(i)=>Vector{Vector{Symbol}}() for i in 1:length(combos) )
    for (arr, counts) in CL_set
        idx = findfirst(x-> x==counts, combos)
        push!(Type_loop[Symbol(idx)], arr)
    end
    write_log(log_path, ">> Type_loop prepared, types=$(keys(Type_loop))")

    perms_result = Vector{Vector{Symbol}}()
    total_groups = length(Type_loop[Symbol(type_N)])
    write_log(log_path, ">> Entering main loop: total groups = $total_groups")

    @time for (grp_i, p1) in enumerate(Type_loop[Symbol(type_N)])
        # 進捗ログ
        write_log(log_path, ">> Group $(grp_i)/$(total_groups): base=$(p1)")

        # (1) 全順列→ペア削除→長さ絞り
        perms      = unique_permutations(p1)
        perms_new  = remove_adjacent_pairs_circular.(perms)
        perms_new2 = [ p for p in perms_new if length(p)==L ]

        # (2) 並列ユニーク化＋ログ
        perms_new3 = parallel_unique(perms_new2, L, log_path)

        # (3) 結果集約
        append!(perms_result, perms_new3)
        write_log(log_path, ">> Group $(grp_i) result count=$(length(perms_new3)), cumulative=$(length(perms_result))")
    end

    # 最終マージ（円環シフト同値排除）
    final = Vector{Vector{Symbol}}()
    for p in perms_result
        if all(!is_circshift_equal(q,p) for q in final)
            push!(final, p)
        end
    end
    write_log(log_path, "=== classify_loop end: final count=$(length(final)) ===")

    return final
end

# チェック・作成するディレクトリのパス
dir = "path/to/your/dir"
function dir_check(dir::String)
    # ディレクトリが存在しなければ作成
    if !isdir(dir)
        mkpath(dir)   # 指定したパスの途中のディレクトリもまとめて作成
        println("Created directory: $dir")
    else
        println("Directory already exists: $dir")
    end
end