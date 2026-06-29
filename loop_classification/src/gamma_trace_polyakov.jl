using DelimitedFiles

include("loop_core_polyakov.jl")

transforms = [
        identity,
        xReflection,
        yReflection,
        zReflection,
        xyReflection,
        xzReflection,
        yzReflection,
        xyzReflection,
    ]
c=1

function exclude_same(data)
    l = length(data)
    L = length(data[1])
    result = []
    push!(result , data[1])
    for i in 2:l
        judge = true 
        for test in result
            for j= 1:L
                if data[i] == circshift(test, j)
                    judge = false
                end
                if data[i] == Dagger(circshift(test, j))
                    judge = false
                end
            end
        end
        judge ? push!(result, data[i]) : nothing
    end
    return result
end



using LinearAlgebra

# 2×2 のパウリ行列
const σ1 = [0 1;
            1 0]
const σ2 = [0 -im;
            im  0]
const σ3 = [1  0;
            0 -1]

# 2×2 の単位行列と零行列（複素型）
const I2    = Matrix{ComplexF64}(I, 2, 2)
const I4    = Matrix{ComplexF64}(I, 4, 4)
const zero2 = zeros(ComplexF64, 2, 2)

# Euclid 空間での γ_μ （μ = 1,2,3,4）
# anticommutator {γ_μ, γ_ν} = 2 δ_{μν} を満たし，
# すべて Hermitian: γ_μ† = γ_μ

const γ1 = -im*[ zero2        σ1;
            -σ1        zero2 ]

const γ2 = -im*[ zero2        σ2;
             -σ2        zero2 ]

const γ3 = -im*[ zero2        σ3;
            -σ3        zero2 ]

# Euclid 時間方向（Minkowski では γ^0 に相当）
const γ4 = [ I2        zero2;
             zero2    -I2]

# 　chiral演算子 γ5
const γ5 = γ1 * γ2 * γ3 * γ4
# 確認：反交換子が 2δμν になるかチェック
for μ in 1:4, ν in 1:4
    anticomm = γ1 == 0 ? nothing : nothing  # placeholder
end

# println("γ1 =\n", γ1)
# println("γ2 =\n", γ2)
# println("γ3 =\n", γ3)
# println("γ4 =\n", γ4)
# println("γ5 =\n", γ5)

# 【チェック用関数】
function check_gamma(γs)
    # 抗交換子 {γμ, γν} == 2δ_{μν}I をチェック
    println("=== Anticommutation ===")
    for μ in 1:4, ν in 1:4
        A = γs[μ]*γs[ν] + γs[ν]*γs[μ]
        target = 2*(μ==ν ? I4 : zeros(ComplexF64,4,4))
        if norm(A - target) > 1e-12
            # println("ng at (",μ,",",ν,")  norm=", norm(A-target))
        else
            # println("$A,Noooo! ng at (",μ,",",ν,")  norm=", norm(A-target))
            # println(norm(A))
        end
    end

    # エルミート性 γμ† == γμ をチェック
    println("\n=== Hermiticity ===")
    for μ in 1:4
        diff = γs[μ] - γs[μ]'
        println("γ",μ, " Hermitian?  ", norm(diff) < 1e-12)
    end
end

# 4×4 単位行列
const I4 = Matrix{ComplexF64}(I, 4, 4)

# 実行
check_gamma((γ1, γ2, γ3, γ4))

function character2gamma(arr::Vector{Symbol})
    mapping = Dict(
        :a => (I4 -γ1),
        :b => (I4 -γ2),
        :c => (I4 -γ3),
        :d => (I4 -γ4), 
        Symbol("-a") => (I4 +γ1),
        Symbol("-b") => (I4 +γ2),
        Symbol("-c") => (I4 +γ3),
        Symbol("-d") => (I4 +γ4) 
    )
    matrix = 1.0
    for x in arr
        matrix *= mapping[x]
    end
    # return (matrix)
    return tr(matrix)
end