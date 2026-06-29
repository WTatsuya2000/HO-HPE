# Wilson ループ分類（不変量バケツ方式・Invariant.md ベース）
include(joinpath(@__DIR__, "loop_core_wilson.jl"))
include(joinpath(@__DIR__, "invariant_core.jl"))

"""
    classify_loop_inv(L, type_N) -> Vector{Vector{Symbol}}

長さ L・指定 type の独立 Wilson ループ代表（不変量バケツ＋バケツ内 is_circshift_equal）。
"""
function classify_loop_inv(L::Int, type_N::Int)
    combos = count_combinations(L)
    (1 <= type_N <= length(combos)) || error("type_N=$type_N は範囲外 (1..$(length(combos)))")
    target = combos[type_N]
    bag1 = nothing
    for arr in add_minus.(all_combinations_multi_types(L))
        if sort([count(x->x==s, arr) for s in (:a,:b,:c,:d)]) == target
            bag1 = arr; break
        end
    end
    bag1 === nothing && return Vector{Vector{Symbol}}()
    perms  = unique_permutations(bag1)
    perms2 = [p for p in remove_adjacent_pairs_circular.(perms) if length(p) == L]
    return classify_invariant(perms2)
end
