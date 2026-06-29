# 2種　4C2 = 6
#=
function x2y(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :b,Symbol("-a") => Symbol("-b"),
        :b  => :a,Symbol("-b") => Symbol("-a"),
        :c  => :c,Symbol("-c") => Symbol("-c"),
        :d  => :d,Symbol("-d") => Symbol("-d")
        )
    return [mapping[x] for x in arr]
end

function x2z(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :c,Symbol("-a") => Symbol("-c"),
        :b  => :b,Symbol("-b") => Symbol("-b"),
        :c  => :a,Symbol("-c") => Symbol("-a"),
        :d  => :d,Symbol("-d") => Symbol("-d")
        )
    return [mapping[x] for x in arr]
end

function x2t(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :d,Symbol("-a") => Symbol("-d"),
        :b  => :b,Symbol("-b") => Symbol("-b"),
        :c  => :c,Symbol("-c") => Symbol("-c"),
        :d  => :a,Symbol("-d") => Symbol("-a")
        )
    return [mapping[x] for x in arr]
end

function y2z(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :a,Symbol("-a") => Symbol("-a"),
        :b  => :c,Symbol("-b") => Symbol("-c"),
        :c  => :b,Symbol("-c") => Symbol("-b"),
        :d  => :d,Symbol("-d") => Symbol("-d")
        )
    return [mapping[x] for x in arr]
end
function y2t(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :a,Symbol("-a") => Symbol("-a"),
        :b  => :d,Symbol("-b") => Symbol("-d"),
        :c  => :c,Symbol("-c") => Symbol("-c"),
        :d  => :b,Symbol("-d") => Symbol("-b")
        )
    return [mapping[x] for x in arr]
end

function z2t(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :a,Symbol("-a") => Symbol("-a"),
        :b  => :b,Symbol("-b") => Symbol("-b"),
        :c  => :d,Symbol("-c") => Symbol("-d"),
        :d  => :c,Symbol("-d") => Symbol("-c")
        )
    return [mapping[x] for x in arr]
end



# 3種　4C3 * 2 = 4 * 2  = 8
# tぬけ
function x2y2z(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :b,Symbol("-a") => Symbol("-b"),
        :b  => :c,Symbol("-b") => Symbol("-c"),
        :c  => :a,Symbol("-c") => Symbol("-a"),
        :d  => :d,Symbol("-d") => Symbol("-d")
        )
    return [mapping[x] for x in arr]
end
function x2z2y(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :c,Symbol("-a") => Symbol("-c"),
        :b  => :a,Symbol("-b") => Symbol("-a"),
        :c  => :b,Symbol("-c") => Symbol("-b")
        )
    return [mapping[x] for x in arr]
end
# z抜け
function x2y2t(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :b,Symbol("-a") => Symbol("-b"),
        :b  => :d,Symbol("-b") => Symbol("-d"),
        :d  => :a,Symbol("-d") => Symbol("-a"),
        :c  => :c,Symbol("-c") => Symbol("-c"),
        )
    return [mapping[x] for x in arr]
end
function x2t2y(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :d,Symbol("-a") => Symbol("-d"),
        :b  => :a,Symbol("-b") => Symbol("-a"),
        :d  => :b,Symbol("-d") => Symbol("-b"),
        :c  => :c,Symbol("-c") => Symbol("-c"),
        )
    return [mapping[x] for x in arr]
end

# y抜け
function x2z2t(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :c,Symbol("-a") => Symbol("-c"),
        :c  => :d,Symbol("-c") => Symbol("-d"),
        :d  => :a,Symbol("-d") => Symbol("-a"),
        :b  => :b,Symbol("-b") => Symbol("-b"),
        )
    return [mapping[x] for x in arr]
end
function x2t2z(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :d,Symbol("-a") => Symbol("-d"),
        :c  => :a,Symbol("-c") => Symbol("-a"),
        :d  => :c,Symbol("-d") => Symbol("-c"),
        :b  => :b,Symbol("-b") => Symbol("-b"),
        )
    return [mapping[x] for x in arr]
end
# x抜け
function y2z2t(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :a,Symbol("-a") => Symbol("-a"),
        :b  => :c,Symbol("-b") => Symbol("-c"),
        :c  => :d,Symbol("-c") => Symbol("-d"),
        :d  => :b,Symbol("-d") => Symbol("-b"),
        )
    return [mapping[x] for x in arr]
end
function y2t2z(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :a,Symbol("-a") => Symbol("-a"),
        :b  => :d,Symbol("-b") => Symbol("-d"),
        :c  => :b,Symbol("-c") => Symbol("-b"),
        :d  => :c,Symbol("-d") => Symbol("-c"),
        )
    return [mapping[x] for x in arr]
end

# 4種　4C4 * 2 = 1 * 2 = 2
function x2y2z2t(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :b,Symbol("-a") => Symbol("-b"),
        :b  => :c,Symbol("-b") => Symbol("-c"),
        :c  => :d,Symbol("-c") => Symbol("-d"),
        :d  => :a,Symbol("-d") => Symbol("-a"),
        )
    return [mapping[x] for x in arr]
end
function x2t2z2y(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :d,Symbol("-a") => Symbol("-d"),
        :b  => :a,Symbol("-b") => Symbol("-a"),
        :c  => :b,Symbol("-c") => Symbol("-b"),
        :d  => :c,Symbol("-d") => Symbol("-c"),
        )
    return [mapping[x] for x in arr]
end
function x2t2z2y(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :c,Symbol("-a") => Symbol("-c"),
        :b  => :d,Symbol("-b") => Symbol("-d"),
        :c  => :a,Symbol("-c") => Symbol("-a"),
        :d  => :b,Symbol("-d") => Symbol("-b"),
        )
    return [mapping[x] for x in arr]
end

 ----------------------- =#
 using Combinatorics
using Combinatorics

# 1. a,b,c の順列に基づいて mapping を作る関数（d は固定）
function generate_mapping_abc(p::Vector{Symbol})
    mapping = Dict(
        :a  => p[1],  Symbol("-a") => Symbol("-" * String(p[1])),
        :b  => p[2],  Symbol("-b") => Symbol("-" * String(p[2])),
        :c  => p[3],  Symbol("-c") => Symbol("-" * String(p[3])),
        :d  => :d,    Symbol("-d") => Symbol("-d")
    )
    return mapping
end

# 2. a,b,c の全順列に対して写像を生成
function generate_all_mappings_abc()
    basis_abc = [:a, :b, :c]
    perms = collect(permutations(basis_abc))  # 6通り
    return [generate_mapping_abc(p) for p in perms]
end

# 3. 一般の変換関数
function x2y(arr::Vector{Symbol}, mapping::Dict{Symbol, Symbol})
    return [mapping[x] for x in arr]
end

# 4. すべての写像を生成
all_mappings_abc = generate_all_mappings_abc()
@show length(all_mappings_abc)  # 6
# 5. 特定の順列番号に対応する変換
function xyzt_sym(arr::Vector{Symbol}, i::Int)
    mapping = all_mappings_abc[i]
    return [mapping[x] for x in arr]
end

# --- 使用例 ---
# arr = [:a, Symbol("-b"), :c, :d, Symbol("-a"), :b]

# for (i, mapping) in enumerate(all_mappings_abc)
#     println("Permutation $(i): ", mapping)
#     println("Transformed array: ", x2y(arr, mapping))
#     println("------")
# end


function Dagger(arr::Vector{Symbol})
    mapping = Dict(
        :a => Symbol("-a"), 
        Symbol("-a") => :a, 
        :b => Symbol("-b"), 
        Symbol("-b") => :b,
        :c => Symbol("-c"), 
        Symbol("-c") => :c,
        :d => Symbol("-d"), 
        Symbol("-d") => :d
        )
    return [mapping[x] for x in reverse(arr)]
end

function xReflection(arr::Vector{Symbol})
    mapping = Dict(
        :a  => Symbol("-a"),
        Symbol("-a") => :a,
        :b  => :b,
        Symbol("-b") => Symbol("-b"),
        :c  => :c,
        Symbol("-c") => Symbol("-c"),
        :d => :d, 
        Symbol("-d") => Symbol("-d")
        )
    return [mapping[x] for x in arr]
end
function yReflection(arr::Vector{Symbol})
    mapping = Dict(
        :a  => :a,
        Symbol("-a") => Symbol("-a"),
        :b  => Symbol("-b"),
        Symbol("-b") => :b,
        :c => :c,
        Symbol("-c") => Symbol("-c"),
        :d => :d, 
        Symbol("-d") => Symbol("-d")
    )
    return [mapping[x] for x in arr]
end
function zReflection(arr::Vector{Symbol})
    mapping = Dict(
        :a => :a,
        Symbol("-a") => Symbol("-a"),
        :b => :b,
        Symbol("-b") => Symbol("-b"),
        :c  => Symbol("-c"),
        Symbol("-c") => :c,
        :d => :d, 
        Symbol("-d") => Symbol("-d")
    )
    return [mapping[x] for x in arr]
end
function tReflection(arr::Vector{Symbol})
    mapping = Dict(
        :a => :a,
        Symbol("-a") => Symbol("-a"),
        :b => :b,
        Symbol("-b") => Symbol("-b"),
        :c  => :c,
        Symbol("-c") => Symbol("-c"),
        :d  => Symbol("-d"),
        Symbol("-d") => :d,
    )
    return [mapping[x] for x in arr]
end
function xyReflection(arr::Vector{Symbol})
    arr = xReflection(arr)
    arr = yReflection(arr)
    return arr
end
function xzReflection(arr::Vector{Symbol})
    arr = xReflection(arr)
    arr = zReflection(arr)
    return arr
end
function yzReflection(arr::Vector{Symbol})
    arr = yReflection(arr)
    arr = zReflection(arr)
    return arr
end
function xtReflection(arr::Vector{Symbol})
    arr = xReflection(arr)
    arr = tReflection(arr)
    return arr
end
function ytReflection(arr::Vector{Symbol})
    arr = yReflection(arr)
    arr = tReflection(arr)
    return arr
end
function ztReflection(arr::Vector{Symbol})
    arr = tReflection(arr)
    arr = zReflection(arr)
    return arr
end
function xyzReflection(arr::Vector{Symbol})
    arr = xReflection(arr)
    arr = yReflection(arr)
    arr = zReflection(arr)
    return arr
end
function yztReflection(arr::Vector{Symbol})
    arr = yReflection(arr)
    arr = zReflection(arr)
    arr = tReflection(arr)
    return arr
end
function xztReflection(arr::Vector{Symbol})
    arr = xReflection(arr)
    arr = zReflection(arr)
    arr = tReflection(arr)
    return arr
end
function xytReflection(arr::Vector{Symbol})
    arr = xReflection(arr)
    arr = yReflection(arr)
    arr = tReflection(arr)
    return arr
end
function xyztReflection(arr::Vector{Symbol})
    arr = xReflection(arr)
    arr = yReflection(arr)
    arr = zReflection(arr)
    arr = tReflection(arr)
    return arr
end