#=
  to maintain the package intent correctly
     implicit construction of a Float64 does not require the target become a Float32

  to maintain the package intent correctly
    some primitive operations must be taken with respect to Float32

  to maintain the package intent correctly
     explicit construction of T requires the target become a Float32
=#

"""
    Robust32s

A module for robust Float32 computation.

Exports: `Robust32`, `ComplexR32`
"""
module Robust32s

export Robust32, ComplexR32

import Base: convert, promote_rule, show, string,
             Float64, Float32,
             ==, !=, <, <=, >, >=, isless, isequal, +, -, *, \, /, ^,
             signbit, precision, significand, exponent, sign, eps, inv, sqrt, cbrt, hypot, clamp, clamp!,
             min, max, minmax, frexp, ldexp, abs, copysign, flipsign, zero, one, iszero, isone,
             isfinite, issubnormal, isinf, isnan, float, floatmin, floatmax, maxintfloat, typemax, typemin,
             evalpoly
             
import Base.Math: abs2, acos, acosd, acosh, acot, acotd, acoth, acsc, acscd, acsch, asec, asecd, asech, 
                  asin, asind, asinh, atan, atand, atanh, cos, cosc, cosd, cosh, cospi, cot, cotd, coth,
                  csc, cscd, csch, deg2rad, evalpoly, exp, exp10, exp2, expm1,
                  log, log10, log1p, log2, mod2pi, modf, rad2deg, rem2pi, sec, secd, sech,
                  sin, sinc, sincos, sincosd, sind, sinh, sinpi, tan, tand, tan #  sincospi

using LinearAlgebra
# using Gaius

struct As64 end # internal use only

struct Robust32 <: AbstractFloat
    val::Float64
  
    Robust32(x::Float64) = new(Float64(Float32(x)))
    Robust32(::Type{As64}, x::Float64) = new(x)
end

Robust32(x::Robust32) = x # idempotency

value64(x::Robust32) = x.val
value32(x::Robust32) = Float32(x.val)

const ComplexR32 = Complex{Robust32}

value64(x::ComplexR32) = (value64(x.re), value64(x.im))
value32(x::ComplexR32) = (value32(x.re), value32(x.im))

# internal use only
Rob32(x::Float64) = Robust32(As64, x)

Float64(x::Robust32) = Float64(value32(x))
convert(::Type{Float64}, x::Robust32) = Float64(x)
promote_rule(::Type{Robust32}, ::Type{Float64}) = Robust32
# Robust32(x::Float64) defined as an internal constructor

Float32(x::Robust32) = value32(x)
convert(::Type{Float32}, x::Robust32) = Float32(x)
promote_rule(::Type{Robust32}, ::Type{Float32}) = Robust32
Robust32(x::Float32) = Rob32(Float64(x))

show(io::IO, x::Robust32) = show(io, value32(x))
string(x::Robust32) = string(value32(x))

for T in (:BigFloat, :Float16)
  @eval begin
    Base.$T(x::Robust32) = $T(value32(x))
    Robust32(x::$T) = Rob32(Float64(Float32(x)))
    Base.convert(::Type{$T}, x::Robust32) = $T(x)
    Base.promote_rule(::Type{Robust32}, ::Type{$T}) = Robust32
  end
end

for T in (:BigInt, :Int128, :Int64, :Int32, :Int16, :Int8,
                   :UInt128, :UInt64, :UInt32, :UInt16, :UInt8)
  @eval begin
    Base.$T(x::Robust32) = $T(value32(x))
    Robust32(x::$T) = Rob32(Float64(Float32(x)))
    convert(::Type{$T}, x::Robust32) = $T(x)
    promote_rule(::Type{Robust32}, ::Type{$T}) = Robust32
  end
end

convert(::Type{Rational{T}}, x::Robust32) where {T} = convert(Rational{T}, value32(x))
convert(::Type{Rational}, x::Robust32) = convert(Rational{Int64}, x)
promote_rule(::Type{Robust32}, ::Type{Rational}) = Robust32

Robust32(x::T) where {T<:Irrational} = Robust32(Float64(x))
convert(::Type{Robust32}, x::T) where {T<:Irrational} = Rob32(Float64(x))
promote_rule(::Type{Robust32}, ::Type{Irrational})  = Robust32

const Robust32_0 = Rob32(0.0)
const Robust32_1 = Rob32(1.0)
const Robust32_2 = Rob32(2.0)

Robust32(x::Bool) = x ? Robust32_1 : Robust32_0

Base.hash(x::Robust32, h::UInt64) = Base.hash(value32(x), h)

Base.decompose(x::Robust32) = Base.decompose(value32(x))
Base.rtoldefault(x::Robust32) = Base.rtoldefault(Float32(x))

precision(::Type{Robust32}) = precision(Float32)

for F in (:floatmin, :floatmax, :maxintfloat, :typemax, :typemin)
  @eval $F(::Type{Robust32}) = Robust32($F(Float32))
end

eps(x::Robust32) = eps(value32(x))
significand(x::Robust32) = significand(value32(x))
exponent(x::Robust32) = exponent(value32(x))
sign(x::Robust32) = exponent(value32(x))
iszero(x::Robust32) = iszero(value32(x))
isone(x::Robust32) = isone(value32(x))
isfinite(x::Robust32) = isfinite(value32(x))
issubnormal(x::Robust32) = issubnormal(value32(x))
isinf(x::Robust32) = isinf(value64(x))
isnan(x::Robust32) = isnan(value64(x))

signbit(x::Robust32) = signbit(value32(x))

zero(::Type{Robust32}) = Robust32_0
one(::Type{Robust32}) = Robust32_1
two(::Type{Robust32}) = Robust32_2
zero(x::Robust32) = zero(Robust32)
one(x::Robust32) = one(Robust32)
two(x::Robust32) = two(Robust32)

frexp(x::Robust32) = frexp(value32(x))
ldexp(fr::Float32, ex::Int) = Rob32(ldexp(Float64(fr), ex))

for F in (:-, :abs, :inv, :sqrt, :cbrt)
  @eval Base.$F(x::Robust32) = Rob32($F(value64(x)))
end

for F in (:(==), :(!=), :(<), :(<=), :(>), :(>=), :isless, :isequal)
  @eval begin
    Base.$F(x::Robust32, y::Robust32) = $F(value32(x), value32(y))
    Base.$F(x::Robust32, y::Real) = $F(promote(x,y)...)
    Base.$F(x::Real, y::Robust32) = $F(promote(x,y)...)
  end  
end

for F in (:+, :-, :*, :/, :\, :hypot, :copysign, :flipsign)
  @eval begin
    Base.$F(x::Robust32, y::Robust32) = Rob32($F(value64(x), value64(y)))
    Base.$F(x::Robust32, y::Real) = $F(promote(x,y)...)
    Base.$F(x::Real, y::Robust32) = $F(promote(x,y)...)
  end  
end

for F in (:hypot, :clamp)
  @eval begin
    Base.$F(x::Robust32, y::Robust32, z::Robust32) = Robust32($F(value(x), value(y), value(z)))
    Base.$F(x::Robust32, y::Real, z::Real) = $F(promote(x,y,z)...)
    Base.$F(x::Real, y::Robust32, z::Real) = $F(promote(x,y,z)...)
    Base.$F(x::Real, y::Real, z::Robust32) = $F(promote(x,y,z)...)
  end  
end

for F in (:abs2, :acos, :acosd, :acosh, :acot, :acotd, :acoth, :acsc, :acscd, :acsch, :asec,
          :asecd, :asech, :asin, :asind, :asinh, :atan, :atand, :atanh, :cos, :cosc,
          :cosd, :cosh, :cospi, :cot, :cotd, :coth, :csc, :cscd, :csch, :deg2rad,
          :exp, :exp10, :exp2, :expm1, :log, :log10, :log1p, :log2, :mod2pi,
          :rad2deg, :rem2pi, :sec, :secd, :sech, :sin, :sinc, :sind, :sinh,
          :sinpi, :tan, :tand, :tanh)
    @eval Base.Math.$F(x::Robust32) = Rob32($F(value64(x)))
end

Base.Math.atan(x::Robust32, y::Robust32) = Rob32(atan(value64(x), value64(y)))

for F in (:modf, :sincos, :sincosd) # , :sincospi)
  @eval function $F(x::Robust32)
            s, c = Base.Math.$F(value64(x))
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

#=     
     provide(x) uses reinterpret
     
     provide(x::Vector{Robust32}), provide(x::Matrix{Robust32})
        presents an Array of Robust32s as an Array of Float64s
     provide(x::Vector{Float64}), provide(x::Matrix{Float64})
        presents an Array of Float64s as an Array of Robust32s
     rewrap(x) uses unsafe_wrap
     
     rewrap(x::Vector{Robust32}, rewrap(x::Matrix{Robust32})
        obtains an Array of Robust32s as an Array of Float64s
     rewrap(x::Vector{Float64}, rewrap(x::Matrix{Float64})
        obtains an Array of Float64s  as an Array of Robust32s
=#

for (R,F) in ((:Robust32, :Float64), (:ComplexR32, :ComplexF64))
  for A in (:Vector, :Matrix)
    @eval begin
        provide(x::$A{$R}) = reinterpret($F, x)
        provide(x::$A{$F}) = reinterpret($R, x)

        rewrap(x::$A{$R}) = unsafe_wrap($A{$F}, cvtptr($F, x), size(x))             
        rewrap(x::$A{$F}) = unsafe_wrap($A{$R}, cvtptr($R, x), size(x))
    end
  end
end

@inline cvtptr(::Type{T}, m::Array{S,N}) where {N,T,S} =
    convert(Ptr{T}, pointer(m,1))

include("linearalgebra.jl")
  
end  # Robust32s
