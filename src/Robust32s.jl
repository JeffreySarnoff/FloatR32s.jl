module Robust32s

export Robust32, InfR32, NaNR32

import Base: ==, !=, <, <=, >, >=, isless, isequal, +, -, *, \, /, ^,
             signbit, significand, exponent, sign, eps, inv, sqrt, cbrt, hypot, clamp, clamp!,
             min, max, minmax, frexp, ldexp, abs, copysign, flipsign, zero, one, iszero, isone,
             isfinite, issubnormal, isinf, isnan
             
import Base.Math: abs2, acos, acosd, acosh, acot, acotd, acoth, acsc, acscd, acsch, asec, asecd, asech, 
                  asin, asind, asinh, atan, atand, atanh, cos, cosc, cosd, cosh, cospi, cot, cotd, coth,
                  csc, cscd, csch, deg2rad, evalpoly, exp, exp10, exp2, expm1,
                  log, log10, log1p, log2, mod2pi, modf, rad2deg, rem2pi, sec, secd, sech,
                  sin, sinc, sincos, sincosd, sind, sinh, sinpi, tan, tand, tanh
                  # sincospi,

# internal use only
struct As64 end 

struct Robust32 <: AbstractFloat
    val::Float64
  
    Robust32(x::Float64) = new(Float64(Float32(x)))
    Robust32(::Type{As64}, x::Float64) = new(x)
end

Rob32(x::Float64) = Robust32(As64, x)

Robust32(x::Robust32) = x

value64(x::Robust32) = x.val
value32(x::Robust32) = Float32(x.val)

Robust32(x::Float32) = Rob32(Float64(x))

const InfR32 = Rob32(Inf)
const NaNR32 = Rob32(NaN)

#=
  to maintain the package intent correctly
     explicit construction of T requires the target become a Float32
=#
Base.BigFloat(x::Robust32) = BigFloat(value32(x))
Base.Float64(x::Robust32) = Float64(value32(x))
Base.Float32(x::Robust32) = value32(x)
Base.Float16(x::Robust32) = Float16(value32(x))

for T in (:Float32, :Float16, :Signed, :Unsigned)
   @eval Robust32(x::$T) = Robust32(Float64(X))
end   

Base.show(io::IO, x::Robust32) = print(io, Float32(x))

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
significand(x::Robust32) = significand(Float32(x))
exponent(x::Robust32) = exponent(Float32(x))
iszero(x::Robust32) = iszero(Float32(x))
isone(x::Robust32) = isone(Float32(x))
isfinite(x::Robust32) = isfinite(Float32(x))
issubnormal(x::Robust32) = issubnormal(Float32(x))
isinf(x::Robust32) = isinf(Float32(x))
isnan(x::Robust32) = isnan(Float32(x))

signbit(x::Robust32) = signbit(value(x))

zero(::Type{Robust32}) = Robust32(0.0)
one(::Type{Robust32}) = Robust32(1.0)
zero(x::Robust32) = zero(Robust32)
one(x::Robust32) = one(Robust32)

frexp(x::Robust32) = frexp(Float32(x))

for F in (:sign, :exponent, :significand)
  @eval $F(x::Robust32) = $F(Float32(x))
end

for F in (:-, :abs, :sign, :inv, :sqrt, :cbrt)
  @eval $F(x::Robust32) = Rob32($F(value64(x)))
end

for F in (:(==), :(!=), :(<), :(<=), :(>), :(>=), :isless, :isequal)
  @eval begin
    $F(x::Robust32, y::Robust32) = $F(value32(x), value32(y))
    $F(x::Robust32, y::Real) = $F(promote(x,y)...)
    $F(x::Real, y::Robust32) = $F(promote(x,y)...)
  end  
end

for F in (:+, :-, :*, :/, :\, :hypot, :copysign, :flipsign)
  @eval begin
    $F(x::Robust32, y::Robust32) = Rob32($F(value64(x), value64(y)))
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

for F in (:abs2, :acos, :acosd, :acosh, :acot, :acotd, :acoth, :acsc, :acscd, :acsch, :asec,
          :asecd, :asech, :asin, :asind, :asinh, :atan, :atand, :atanh, :cos, :cosc,
          :cosd, :cosh, :cospi, :cot, :cotd, :coth, :csc, :cscd, :csch, :deg2rad,
          :exp, :exp10, :exp2, :expm1, :log, :log10, :log1p, :log2, :mod2pi,
          :rad2deg, :rem2pi, :sec, :secd, :sech, :sin, :sinc, :sind, :sinh,
          :sinpi, :tan, :tand, :tanh)
    @eval $F(x::Robust32) = Rob32($F(value64(x)))
end

for F in (:modf, :sincos, :sincosd) # , :sincospi)
  @eval function $F(x::Robust32)
            s, c = $F(value64(x))
            return Rob32(s), Rob32(c)
         end
end

function evalpoly(x::Robust32, p::NTuple{N, Robust32}) where {N}
    Rob32(evalpoly(value64(x), map(value64, p)))
end
function evalpoly(x::T, p::NTuple{N, Robust32}) where {T,N}
    Rob32(evalpoly(Float64(x), map(value64, p)))
end
function evalpoly(x::Robust32, p::NTuple{N, T}) where {T,N}
    Rob32(evalpoly(value64(x), p))
end

end  # Robust32s
