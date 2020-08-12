
#=
   functions to enfold:
   organized by applicable signature.
=#

const scalar_functions_of_one_arg = 
    (abs2, acos, acosd, acosh, acot, acotd, acoth, 
     acsc, acscd, acsch, asec, asecd, asech, 
     asin, asind, asinh, atan, atand, atanh,
     cos, cosc, cosd, cosh, cospi, cot, cotd, coth,
     csc, cscd, csch, deg2rad, evalpoly, exp, exp10, exp2, expm1,
     log, log10, log1p, log2, mod2pi, modf, rad2deg, rem2pi,
     sec, secd, sech, sin, sinc, sind,
     sinh, sinpi, tan, tand, tan)

# sincos, sincospi, sincosd, return 2 

const scalar_functions_of_two_args   = ()
const scalar_functions_of_three_args = ()

const vector_functions_of_one_arg  = ()
const vector_functions_of_two_args = ()

const vector_adjoint_functions_of_two_args = ()
const adjoint_vector_functions_of_two_args = ()

const matrix_vector_functions_of_two_args = ()
const vector_matrix_functions_of_two_args = ()

const matrix_functions_of_one_arg  = ()
const matrix_functions_of_two_args = ()

const matrix_adjoint_functions_of_two_args = ()
const adjoint_matrix_functions_of_two_args = ()

#=     
     provide(x) uses reinterpret
     
     provide(x::Vector{Robust32}), provide(x::Matrix{Robust32})
        presents an Array of Robust32s as an Array of Float64s

     provide(x::Vector{Float64}), provide(x::Matrix{Float64})
        presents an Array of Float64s as an Array of Robust32s

     provide(x::Vector{ComplexR32}), provide(x::Matrix{ComplexR32})
        presents an Array of ComplexR32s as an Array of ComplexF64s

     provide(x::Vector{ComplexF64}), provide(x::Matrix{ComplexF64})
        presents an Array of ComplexF64s as an Array of ComplexR32s
=#

for (R,F) in ((:Robust32, :Float64), (:ComplexR32, :ComplexF64))
  for A in (:Vector, :Matrix)
    @eval begin
        provide(x::$A{$R}) = reinterpret($F, x)
        provide(x::$A{$F}) = reinterpret($R, x)
    end
  end
end

#=
    enfold(fn, ..) simplifies implementing functions correctly
=#

for F in scalar_functions_of_one_arg
  @eval enfold($F, x::Robust32) = provide(fn(provide(x)))
end

for F in scalar_functions_of_two_args
  @eval enfold($F, x::Robust32, y::Robust32) = provide(fn(provide(x), provide(y)))
end

for F in scalar_functions_of_three_args
  @eval enfold($F, x::Robust32, y::Robust32, z::Robust32) = provide(fn(provide(x), provide(y), provide(z)))
end

for F in vector_functions_of_one_arg
  @eval enfold($F, x::Vector{Robust32}) = provide(fn(provide(x)))
end

for F in vector_functions_of_two_args
  @eval enfold($F, x::Vector{Robust32}, y::Vector{Robust32}) = provide(fn(provide(x), provide(y)))
end

for F in vector_adjoint_functions_of_two_args
  @eval enfold($F, x::Vector{Robust32}, y::Adjoint{Robust32, Vector{Robust32}}) = provide(fn(provide(x), provide(y)))
end

for F in adjoint_vector_functions_of_two_args
  @eval enfold($F, x::Adjoint{Robust32, Vector{Robust32}}, y::Vector{Robust32}) = provide(fn(provide(x), provide(y)))
end

for F in matrix_functions_of_one_arg
  @eval enfold($F, x::Matrix{Robust32}) = provide(fn(provide(x)))
end

for F in matrix_functions_of_two_args
  @eval enfold($F, x::Matrix{Robust32}, y::Matrix{Robust32}) = provide(fn(provide(x), provide(y)))
end

for F in vector_functions_of_one_arg
  @eval enfold($F, x::Vector{Robust32}) = provide(fn(provide(x)))
end

for F in matrix_vector_functions_of_two_args
  @eval enfold($F, x::Matrix{Robust32}, y::Vector{Robust32}) = provide(fn(provide(x), provide(y)))
end

for F in vector_matrix_functions_of_two_args
  @eval enfold($F, x::Vector{Robust32}, y::Matrix{Robust32}) = provide(fn(provide(x), provide(y)))
end

for F in matrix_adjoint_functions_of_two_args
  @eval enfold($F, x::Matrix{Robust32}, y::Adjoint{Robust32, Vector{Robust32}}) = provide(fn(provide(x), provide(y)))
end

for F in adjoint_matrix_functions_of_two_args
  @eval enfold($F, x::Adjoint{Robust32, Vector{Robust32}}, y::Matrix{Robust32}) = provide(fn(provide(x), provide(y)))
end
