module Robust32s

export Robust32

import Base: signbit, significand, exponent, sign, eps, inv, sqrt, cbrt, +, -, *, \, /, ^, hypot, clamp, clamp!,
             min, max, minmax, frexp, ldexp, abs, copysign, flipsign
             
import Base.Math: acos, acosd, acosh, acot, acotd, acoth, acsc, acscd, acsch, asec, asecd, asech, 
                  asin, asind, asinh, atan, atand, atanh, cos, cosc, cosd, cosh, cospi, cot, cotd, coth,
                  csc, cscd, csch, deg2rad, evalpoly, exp, exp10, exp2, expm1,
                  log, log10, log1p, log2, mod2pi, modf, rad2deg, rem2pi, sec, secd, sech,
                  sin, sinc, sincos, sincosd, sincospi, sind, sinh, sinpi, tan, tand, tanh

type Robust32
    val::Float64
end

value(x::Robust32) = x.val

Robust32(x::Robust32) = x

#=
  to maintain the package intent correctly
     explicit construction of a Float64 requires the target become a Float32
=#     
Base.Float64(x::Robust32) = Float64(Float32(x.val))
Base.Float32(x::Robust32) = Float32(x.val)
Base.Float16(x::Robust32) = Float16(x.val)

for T in (:Float32, :Float16, :Signed, :Unsigned)
   @eval Robust32(x::$T) = Roubust32(Float64(X))
end   

Base.show(io::IO, x::Roubust32) = print(io, Float32(x))

Base.promote_rule(::Type{Robust32}, ::Type{Base.IEEEFloat}) = Robust32
Base.promote_rule(::Type{Robust32}, ::Type{T}) where {T<:Signed} = Robust32
Base.promote_rule(::Type{Robust32}, ::Type{T}) where {T<:Unsigned} = Robust32

#=
  to maintain the package intent correctly
     implicit construction of a Float64 does not require the target become a Float32
=# 
Base.convert(::Type{Robust32}, x::Float64) = Robust32(x)
Base.convert(::Type{Robust32}, x::T) where {T<:Union{Float32, Float16}} = Robust32(x)
Base.convert(::Type{Robust32}, x::T) where {T<:Signed} = Robust32(x)
Base.convert(::Type{Robust32}, x::T) where {T<:Unsigned} = Robust32(x)

#=
  to maintain the package intent correctly
    some primitive operations must be taken with respect to Float32
=#
eps(x::Robust32) = eps(Float32(x))
significand(x::Robust32) = significand(Float32(x)))
exponent(x::Robust32) = exponent(Float32(x))

signbit(x::Robust32) = signbit(value(x))

for F in (:-, :abs, :sign, :inv, :sqrt, :cbrt)
  @eval $F(x::Robust32) = Robust32($F(value(x)))
end

for F in (:+, :-, :*, :/, :\, :hypot, :copysign, :flipsign)
  @eval begin
    $F(x::Robust32, y::Robust32) = Robust32($F(value(x), value(y)))
    $F(x::Robust32, y::Real) = $F(promote(x,y)...)
    $F(x::Real, y::Robust32) = $F(promote(x,y)...)
  end  
end    

for F in (:hypot, :clamp)
  @eval begin
    $F(x::Robust32, y::Robust32, z::Robust32) = Robust32($F(value(x), value(y), value(z)))
    $F(x::Robust32, y::Real, z::Real) = $F(promote(x,y,z)...)
    $F(x::Real, y::Robust32, z::Real) = $F(promote(x,y,z)...)
    $F(x::Real, y::Real, z::Robust32) = $F(promote(x,y,z)...)
  end  
end

for F in (:acos, :acosd, :acosh, :acot, :acotd, :acoth, :acsc, :acscd, :acsch, :asec,
          :asecd, :asech, :asin, :asind, :asinh, :atan, :atand, :atanh, :cos, :cosc,
          :cosd, :cosh, :cospi, :cot, :cotd, :coth, :csc, :cscd, :csch, :deg2rad,
          :exp, :exp10, :exp2, :expm1, :log, :log10, :log1p, :log2, :mod2pi, :modf,
          :rad2deg, :rem2pi, :sec, :secd, :sech, :sin, :sinc, :sincos, :sincosd,
          :sincospi, :sind, :sinh, :sinpi, :tan, :tand, :tanh)
    @eval $F(x::Robust32) = Robust32($F(value(x)))
end

function evalpoly(x::Robust32, p::NTuple{N, Robust32}) where {N}
    Robust32(evalpoly(value(x), map(value, p)))
end
function evalpoly(x::T, p::NTuple{N, Robust32}) where {T,N}
    Robust32(evalpoly(Float64(x), map(value, p)))
end
function evalpoly(x::Robust32, p::NTuple{N, T}) where {T,N}
    Robust32(evalpoly(value(x), p)))
end

end  # Robust32s
