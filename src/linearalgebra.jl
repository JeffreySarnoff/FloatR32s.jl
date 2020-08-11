#=
[  "/", "\\", 
   "adjoint", "adjoint!", "axpby!", "axpy!", "bunchkaufman", "bunchkaufman!", "cholesky", "cholesky!", "cond", "condskeel", 
   "copy_transpose!", "copyto!", "cross", 
   # "det", "diag", 
   "diagind", "diagm", 
   # "dot", 
   # "eigen", "eigen!", 
   "eigmax", "eigmin", 
   # "eigvals", "eigvals!", 
    "eigvecs", 
   "factorize", "givens", "hessenberg", "hessenberg!", 
   # "isdiag", "ishermitian", "isposdef", "isposdef!", "issuccess", "issymmetric", "istril", "istriu", 
   "kron", "ldiv!", "ldlt", "ldlt!", "lmul!", "logabsdet", "logdet", 
   "lowrankdowndate", "lowrankdowndate!", "lowrankupdate", "lowrankupdate!", 
   "lq", "lq!", "lu", "lu!", "lyap", 
   "mul!", "norm", "normalize", "normalize!", "nullspace", "opnorm", 
   "ordschur", "ordschur!", "pinv", "qr", "qr!", 
   "rank", "rdiv!", "reflect!", "rmul!", "rotate!", "schur", "schur!", 
   # "svd", "svd!", "svdvals", "svdvals!", 
   "sylvester", 
   # "tr", "transpose", "transpose!", 
   "tril", "tril!", "triu", "triu!", "×", "⋅" ]
=#

import LinearAlgebra: isdiag, ishermitian, isposdef, isposdef!, issuccess, issymmetric, istril, istriu,
     tr, det, adjoint, adjoint!, transpose, transpose!, diag, diagm, diagind, 
     svdvals, svdvals!, svd, svd!, eigvals, eigvals!, eigen, eigen!

for F in (:+, :-, :*, :/, :\)
  @eval begin
    $F(x::Vector{FloatR32}, y::Vector{FloatR32}) = reinterpret(FloatR32)($F(reinterpret(Float64, x), reinterpret(Float64, y)))
    $F(x::Vector{FloatR32}, y::Vector{Float64})  = reinterpret(FloatR32)($F(reinterpret(Float64, x), y))
    $F(x::Vector{Float64}, y::Vector{FloatR32}) = reinterpret(FloatR32)($F(x, reinterpret(Float64, y)))
    $F(x::Vector{FloatR32}, y::Vector{Float32})  = reinterpret(FloatR32)($F(reinterpret(Float64, x), y))
    $F(x::Vector{Float32}, y::Vector{FloatR32}) = reinterpret(FloatR32)($F(x, reinterpret(Float64, y)))

    $F(x::Matrix{FloatR32}, y::Matrix{FloatR32}) = reinterpret(FloatR32)($F(reinterpret(Float64, x), reinterpret(Float64, y)))
    $F(x::Matrix{FloatR32}, y::Matrix{Float64})  = reinterpret(FloatR32)($F(reinterpret(Float64, x), y))
    $F(x::Matrix{Float64}, y::Matrix{FloatR32}) = reinterpret(FloatR32)($F(x, reinterpret(Float64, y)))
    $F(x::Matrix{FloatR32}, y::Matrix{Float32})  = reinterpret(FloatR32)($F(reinterpret(Float64, x), y))
    $F(x::Matrix{Float32}, y::Matrix{FloatR32}) = reinterpret(FloatR32)($F(x, reinterpret(Float64, y)))
  end
end

for F in (:tr, :det)
    @eval LinearAlgebra.$F(x::Matrix{FloatR32}) = Rob32($F(reinterpret(Float64, x)))
end

for F in (:isdiag, :ishermitian, :isposdef, :isposdef!, :issuccess, :issymmetric, :istril, :istriu)
  @eval LinearAlgebra.$F(x::Matrix{FloatR32}) = $F(reinterpret(Float64, x))
end

LinearAlgebra.abs2(x::Array{FloatR32,N}) where {N} = Rob32(dot(x, x))
LinearAlgebra.abs(x::Array{FloatR32,N}) where {N} = Rob32(sqrt(dot(x, x)))
LinearAlgebra.dot(x::Array{FloatR32,N}, y::Array{FloatR32,N}) where {N} = Rob32(dot(rewrap(x), rewrap(y)))

for F in (:inv, :sqrt, :exp, :log, 
          :sin, :cos, :tan, :csc, :sec, :cot, :asin, :acos, :atan, :acsc, :asec, :acot,
          :sinh, :cosh, :tanh, :csch, :sech, :coth, :asinh, :acosh, :atanh, :acsch, :asech, :acoth)
    @eval LinearAlgebra.$F(x::Matrix{FloatR32}) = reinterpret(FloatR32, LinearAlgebra.$F(reinterpret(Float64, x)))
end

LinearAlgebra.adjoint(x::Matrix{FloatR32}) = Adjoint(x)
LinearAlgebra.transpose(x::Matrix{FloatR32}) = Transpose(x)

# diag diagm diagind
LinearAlgebra.diag(x::Matrix{FloatR32}) = rewrap(diag(rewrap(Float64, x)))
                                
LinearAlgebra.svdvals(A::Matrix{FloatR32}; kw...) = rewrap(svdvals(rewrap(A); kw...))
LinearAlgebra.eigvals(A::Matrix{FloatR32}; kw...) = rewrap(eigvals(rewrap(A); kw...))
LinearAlgebra.svdvals!(A::Matrix{FloatR32}) = rewrap(svdvals!(rewrap(A)))
LinearAlgebra.eigvals!(A::Matrix{FloatR32}) = rewrap(eigvals!(rewrap(A)))

function LinearAlgebra.svd(x::Matrix{FloatR32})
    u, s, v = svd(rewrap(x))
    U = rewrap(u)
    S = rewrap(s)
    V = adjoint(rewrap(adjoint(v)))
    return SVD{FloatR32,FloatR32,Matrix{FloatR32}}(U, S, V)
end
                                    
function LinearAlgebra.svd!(x::Matrix{FloatR32}; kw...)
    u, s, v = svd!(rewrap(x); kw...)
    U = rewrap(u)
    S = rewrap(s)
    V = adjoint(rewrap(adjoint(v)))
    return SVD{FloatR32,FloatR32,Matrix{FloatR32}}(U, S, V)
end

function LinearAlgebra.eigen(x::Matrix{FloatR32}; kw...)
    v, m = eigen(rewrap(x); kw...)
    V = rewrap(v)
    M = rewrap(m)
    return Eigen{ComplexR32,ComplexR32,Matrix{ComplexR32},Vector{ComplexR32}}(V,M)
end

function LinearAlgebra.eigen!(x::Matrix{FloatR32}; kw...)
    v, m = eigen!(rewrap(x); kw...)
    V = rewrap(v)
    M = rewrap(m)
    return Eigen{ComplexR32,ComplexR32,Matrix{ComplexR32},Vector{ComplexR32}}(V,M)
end
