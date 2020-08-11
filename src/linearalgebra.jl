for F in (:+, :-, :*, :/, :\)
  @eval begin
    $F(x::Vector{Robust32}, y::Vector{Robust32}) = reinterpret(Robust32)($F(reinterpret(Float64, x), reinterpret(Float64, y)))
    $F(x::Vector{Robust32}, y::Vector{Float64})  = reinterpret(Robust32)($F(reinterpret(Float64, x), y))
    $F(x::Vector{Float64}, y::Vector{Robust32}) = reinterpret(Robust32)($F(x, reinterpret(Float64, y)))
    $F(x::Vector{Robust32}, y::Vector{Float32})  = reinterpret(Robust32)($F(reinterpret(Float64, x), y))
    $F(x::Vector{Float32}, y::Vector{Robust32}) = reinterpret(Robust32)($F(x, reinterpret(Float64, y)))

    $F(x::Matrix{Robust32}, y::Matrix{Robust32}) = reinterpret(Robust32)($F(reinterpret(Float64, x), reinterpret(Float64, y)))
    $F(x::Matrix{Robust32}, y::Matrix{Float64})  = reinterpret(Robust32)($F(reinterpret(Float64, x), y))
    $F(x::Matrix{Float64}, y::Matrix{Robust32}) = reinterpret(Robust32)($F(x, reinterpret(Float64, y)))
    $F(x::Matrix{Robust32}, y::Matrix{Float32})  = reinterpret(Robust32)($F(reinterpret(Float64, x), y))
    $F(x::Matrix{Float32}, y::Matrix{Robust32}) = reinterpret(Robust32)($F(x, reinterpret(Float64, y)))
  end
end

for F in (:tr, :det)
    @eval LinearAlgebra.$F(x::Matrix{Robust32}) = Rob32($F(reinterpret(Float64, x))
end

for F in (:isdiag, :ishermitian, :isposdef, :isposdef!, :issuccess, :issymmetric, :istril, :istriu)
  @eval LinearAlgebra.$F(x::Matrix{Robust32}) = $F(reinterpret(Float64, x))
end

for F in (:inv, :sqrt, :exp, :log, 
          :sin, :cos, :tan, :csc, :sec, :cot, :asin, :acos, :atan, :acsc, :asec, :acot,
          :sinh, :cosh, :tanh, :csch, :sech, :coth, :asinh, :acosh, :atanh, :acsch, :asech, :acoth)
    @eval LinearAlgebra.$F(x::Matrix{Robust32}) = reinterpret(Robust32, $F(reinterpret(Float64, x)))
end

LinearAlgebra.dot(x::Array{N,Robust32}) where {N} = Rob32(dot(reinterpret(Float64, x)))

LinearAlgebra.adjoint(x::Matrix{Robust32}) = Adjoint(x)
LinearAlgebra.transpose(x::Matrix{Robust32}) = Transpose(x)

# diag diagm diagind
LinearAlgebra.diag(x::Matrix{Robust32}) = rewrap(diag(rewrap(Float64, x)))
                                
LinearAlgebra.svdvals(A::Matrix{Robust32}; kw...) = rewrap(svdvals(rewrap(A); kw...))
LinearAlgebra.eigvals(A::Matrix{Robust32}; kw...) = rewrap(eigvals(rewrap(A); kw...))
LinearAlgebra.svdvals!(A::Matrix{Robust32}) = rewrap(svdvals!(rewrap(A)))
LinearAlgebra.eigvals!(A::Matrix{Robust32}) = rewrap(eigvals!(rewrap(A)))

function LinearAlgebra.svd(x::Matrix{Robust32})
    u, s, v = svd(rewrap(x))
    U = rewrap(u)
    S = rewrap(x)
    V = adjoint(rewrap(adjoint(v)))
    return U, S, V
end
                                    
function LinearAlgebra.svd!(x::Matrix{Robust32}; kw...)
    u, s, v = svd!(rewrap(x); kw...)
    U = rewrap(u)
    S = rewrap(x)
    V = adjoint(rewrap(adjoint(v)))
    return U, S, V
end

function LinearAlgebra.eigen(x::Matrix{Robust32}; kw...)
    v, m = eigen(rewrap(x); kw...)
    V = rewrap(v)
    M = rewrap(m)
    return Eigen{ComplexR32,ComplexR32,Matrix{ComplexR32},Vector{ComplexR32}}(V,M)
end

function LinearAlgebra.eigen!(x::Matrix{Robust32}; kw...)
    v, m = eigen!(rewrap(x); kw...)
    V = rewrap(v)
    M = rewrap(m)
    return Eigen{ComplexR32,ComplexR32,Matrix{ComplexR32},Vector{ComplexR32}}(V,M)
end
