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

Exports: `Robust32`, `ComplexR32`, `robust32`, `complexr32`, `float32`, `float64`
""" Robust32s

module Robust32s

export Robust32, ComplexR32,
       robust32, complexr32,
       float64, float32

import Base: convert, promote_rule, show, string,
             Float64, Float32, Float16, ComplexF64, ComplexF32, ComplexF16,
             float, prevfloat, nextfloat,
             ==, !=, <, <=, >, >=, isless, isequal, +, -, *, \, /, ^, fma, muladd,
             signbit, precision, significand, exponent, sign, eps, inv, sqrt, cbrt, hypot, clamp, clamp!,
             min, max, minmax, frexp, ldexp, abs, copysign, flipsign, zero, one, iszero, isone,
             isfinite, issubnormal, isinf, isnan, floatmin, floatmax, maxintfloat, typemax, typemin,
             mod, rem, div, fld, cld, divrem, fldmod,
             floor, ceil, trunc, round,
             rand, randn, evalpoly,
             sum, prod, cumsum, cumprod, accumulate, sum!, prod!, cumsum!, cumprod!, accumulate!

import Base.Math: abs2, acos, acosd, acosh, acot, acotd, acoth, acsc, acscd, acsch, asec, asecd, asech, 
                  asin, asind, asinh, atan, atand, atanh, cos, cosc, cosd, cosh, cospi, cot, cotd, coth,
                  csc, cscd, csch, deg2rad, evalpoly, exp, exp10, exp2, expm1,
                  log, log10, log1p, log2, mod2pi, modf, rad2deg, rem2pi, sec, secd, sech,
                  sin, sinc, sincos, sincosd, sind, sinh, sinpi, tan, tand, tan #  sincospi

using Random
using LinearAlgebra
# using Gaius

primitive type Robust32 <: AbstractFloat 64 end
const ComplexR32 = Complex{Robust32}

# idempotency
Robust32(x::Robust32) = x
ComplexR32(x::ComplexR32) = x

Float64(x::Robust32) = reinterpret(Float64, x)
Float32(x::Robust32) = Float32(reinterpret(Float64,x))
Robust32(x::Float64) = reinterpret(Robust32, x)
Robust32(x::Float32) = reinterpret(Robust32, Float64(x))
# process `divrem`, `fldmod` automatically 
Robust32(x::Tuple{Float64, Float64}) = (Robust32(x[1]), Robust32(x[2]))

ComplexF64(x::ComplexR32) = ComplexF64(reinterpret(Float64, x.re), reinterpret(Float64, x.im))
ComplexF32(x::ComplexR32) = ComplexF32(reinterpret(Float64, x.re), reinterpret(Float64, x.im))
ComplexR32(x::ComplexF64) = ComplexR32(reinterpret(Robust32, x.re), reinterpret(Robust32, x.im))
ComplexR32(x::ComplexF32) = 
    ComplexR32(reinterpret(Robust32, Float64(x.re)), reinterpret(Robust32, Float64(x.im)))

#=
   `robust32` is used to map input,  Float64 values into Robust32s
   `robust32` is used to map input,  Float32 values into Robust32s
   `float32`  is used to map output, Robust32 values into Float32s
   `float64`  is used to map output, Robust32 values into Float64s
=#
robust32(x::Float64)  = reinterpret(Robust32, Float64(Float32(x)))
# `numeric=false` gives digit memetic representation of Float32 in Float64
robust32(x::Float32; numeric::Bool=true) =
    numeric ? reinterpret(Robust32, Float64(x)) : Robust32(Meta.parse(string(x))
robust32(x::Robust32) = x
float64(x::Robust32)  = Float64(Float32(reinterpret(Float64, x)))
float64(x::Float64)   = x
float64(x::Float32)   = Float64(x)
float32(x::Robust32)  = Float32(reinterpret(Float64, x))
float32(x::Float32)   = x
float32(x::Float64)   = Float32(x)

#=
   `complexr32` is used to map input,  ComplexF64 values into ComplexR32s
   `complexr32` is used to map input,  ComplexF32 values into ComplexR32s
   `complexf32` is used to map output, ComplexR32 values into ComplexF32s
   `complexf64` is used to map output, ComplexR32 values into ComplexF64s
=#
complexR32(x::ComplexF64) =
    ComplexR32(convert(Robust32, x.re), convert(Robust32, x.im))
complexR32(x::ComplexF32) =
    ComplexR32(convert(Robust32, x.re), convert(Robust32, x.im))
complexF64(x::ComplexR32) =
    ComplexF64(convert(Float64, x.re), convert(Float64, x.im))
complexF32(x::ComplexR32) =
    ComplexF32(convert(Float32, x.re), convert(Float32, x.im))

# internal (64bit) conversion
convert(::Type{Robust32}, x::Float64) = reinterpret(Robust32, x)
convert(::Type{Robust32}, x::Float32) = reinterpret(Robust32, Float64(x))
convert(::Type{Float64}, x::Robust32) = reinterpret(Float64, x)
convert(::Type{Float32}, x::Robust32) = Float32(reinterpret(Float64, x))
# internal promotion
promote_rule(::Type{Robust32}, ::Type{Float64}) = Robust32
promote_rule(::Type{Robust32}, ::Type{Float32}) = Robust32

# internal (64bit) conversion
convert(::Type{ComplexR32}, x::ComplexF64) =
    ComplexR32(reinterpret(Robust32, x.re), reinterpret(Robust32, x.im))
convert(::Type{ComplexR32}, x::ComplexF32) =
    ComplexR32(reinterpret(Robust32, Float64(x.re)), reinterpret(Robust32, Float64(x.im)))
convert(::Type{ComplexF64}, x::ComplexR32) =
    ComplexF64(reinterpret(Float64, x.re), reinterpret(Float64, x.im))
convert(::Type{ComplexF32}, x::ComplexR32) =
    ComplexF32(Float32(reinterpret(Float64, x.re)), Float32(reinterpret(Float64, x.im)))
# internal promotion
promote_rule(::Type{ComplexR32}, ::Type{ComplexF64}) = ComplexR32
promote_rule(::Type{ComplexR32}, ::Type{ComplexF32}) = ComplexR32

string(x::Robust32) = string(float32(x))
string(x::ComplexR32) = string(complexf32(x))

show(io::IO, x::Robust32) = print(io, float32(x))
show(io::IO, x::ComplexR32) = print(io, complexf32(x))

for T in (:BigFloat, :Float16)
  @eval begin
    Base.$T(x::Robust32) = $T(value64(x))
    Robust32(x::$T) = Rob32(Float64(x))
    convert(::Type{$T}, x::Robust32) = $T(x)
    promote_rule(::Type{Robust32}, ::Type{$T}) = Robust32
    convert(::Type{Robust32}, x::$T) = Robust32(x)
  end
end

for T in (:BigInt, :Int128, :Int64, :Int32, :Int16, :Int8,
                   :UInt128, :UInt64, :UInt32, :UInt16, :UInt8)
  @eval begin
    Base.$T(x::Robust32) = $T(value32(x))
    Robust32(x::$T) = Rob32(Float64(Float32(x)))
    convert(::Type{$T}, x::Robust32) = $T(x)
    promote_rule(::Type{Robust32}, ::Type{$T}) = Robust32
    convert(::Type{Robust32}, x::$T) = Robust32(x)
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
Bool(x::Robust32) = iszero(x) ? false : true
Base.hash(x::Robust32, h::UInt64) = Base.hash(value32(x), h)

Base.decompose(x::Robust32) = Base.decompose(value32(x))
Base.rtoldefault(x::Robust32) = Base.rtoldefault(Float32(x))

precision(::Type{Robust32}) = precision(Float32)

for F in (:floatmin, :floatmax, :maxintfloat, :typemax, :typemin)
  @eval $F(::Type{Robust32}) = Robust32($F(Float32))
end

eps(::Type{Robust32}) = Robust32(eps(Float32))
eps(x::Robust32) = Robust32(eps(value32(x)))

significand(x::Robust32) = significand(value32(x))
exponent(x::Robust32) = exponent(value32(x))
sign(x::Robust32) = sign(value32(x))

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

nextfloat(x::Robust32, n::Int=1) = Rob32(nextfloat(value64(x),n))
prevfloat(x::Robust32, n::Int=1) = Rob32(prevfloat(value64(x),n))

#=
   >>> Important Implementation Note <<<

   Base.frexp(x::Robust32) = frexp64(x)
   Base.ldexp(fr::Robust32, xp::Int) = ldexp64(fr, xp)

   -----

   frexp32(::Robust32)
   ldexp32(::Robust32, ::Int)

   This pair is to be used with input data (`ldexp32`)
   and output values (`frexp32`) should the need arise. 

   For example, if one prefer to generate 
   initial input data values by combining 
   (significand, exponent) pairings,
   use `ldexp32.(significands, exponents)`.

   Similarly, if one prefer to store results
   in the form of (significand, exponent) pairings
   use `frexp32.(results)`.
  
   ----

   frexp64(::Robust32)
   ldexp64(::Robust32, ::Int)

   This pair works with 64-bit values directly.
   That design decision makes the following hold:
      `ldexp64(frexp64(x::Robust32)) == x`,
   
   _otherwise it would not hold much of the time_
   _consider_ `x = sqrt(Robust32(2.0))`.
=#

frexp(x::Robust32) = frexp64(x)
ldexp(fr::Robust32, xp::Int) = ldexp64(fr, xp)

function frexp64(x::Robust32)
  fr64, xp = frexp(value64(x))
  fr = Rob32(fr64)
  return (fr, xp)
end

function ldexp64(fr::Robust32, ex::Int)
  fr64 = value64(fr)
  return Rob32(ldexp(fr64, ex))
end
ldexp64(x::Tuple{Robust32, Int}) = ldexp64(x[1], x[2])

function frexp32(x::Robust32)
  fr64, xp = frexp(Float64(x))
  fr = Rob32(fr64)
  return (fr, xp)
end

function ldexp32(fr::Robust32, ex::Int)
  fr64 = Float64(fr)
  return Rob32(ldexp(fr64, ex))
end
ldexp32(x::Tuple{Robust32, Int}) = ldexp32(x[1], x[2])

for F in (:-, :abs, :inv, :sqrt, :cbrt)
  @eval $F(x::Robust32) = Rob32($F(value64(x)))
end

for F in (:(==), :(!=), :(<), :(<=), :(>), :(>=), :isless, :isequal)
  @eval begin
    $F(x::Robust32, y::Robust32) = $F(value64(x), value64(y))
    $F(x::Robust32, y::Float64) = $F(value64(x), y)
    $F(x::Robust32, y::Float32) = $F(value64(x), y)
    $F(x::Float64, y::Robust32) = $F(x, value64(y))
    $F(x::Float32, y::Robust32) = $F(x, value64(y))
    $F(x::AbstractFloat, y::Robust32) = $F(promote(x,y)...)
    $F(x::Robust32, y::AbstractFloat) = $F(promote(x,y)...)
    $F(x::Rational, y::Robust32) = $F(x, value64(y))
    $F(x::Robust32, y::Rational) = $F(value64(x), y)
    $F(x::Bool, y::Robust32) = $F(x, value64(y))
    $F(x::Robust32, y::Bool) = $F(value64(x), y)
  end  
end

^(x::Robust32, y::Robust32) = Rob32(value64(x)^value64(y))
^(x::Robust32, y::Integer) = Rob32(value64(x)^y)
^(x::Real, y::Robust32) = Rob32(Float64(x^value64(y)))

for F in (:+, :-, :*, :/, :\, :hypot, :copysign, :flipsign,
          :mod, :rem, :div, :fld, :cld)
  @eval begin
    $F(x::Robust32, y::Robust32) = Rob32($F(value64(x), value64(y)))
    $F(x::Robust32, y::Float32) = Rob32($F(value64(x), y))
    $F(x::Float32, y::Robust32) = Rob32($F(x, value64(y)))
    $F(x::Robust32, y::Float16) = Rob32($F(value64(x), Float64(y)))
    $F(x::Float16, y::Robust32) = Rob32($F(Float64(x), value64(y)))
    $F(x::Bool, y::Robust32) = Rob32($F(x, value64(y)))
    $F(x::Robust32, y::Bool) = Rob32($F(value64(x), y))
    $F(x::AbstractFloat, y::Robust32) = $F(promote(x,y)...)
    $F(x::Robust32, y::AbstractFloat) = $F(promote(x,y)...)
    $F(x::Rational, y::Robust32) = $F(promote(x,y)...)
    $F(x::Robust32, y::Rational) = $F(promote(x,y)...)
    $F(x::Robust32, y::T) where {T<:Real} = $F(x, Rob32(y))
    $F(x::T, y::Robust32) where {T<:Real} = $F(Rob32(x), y)
  end  
end

# divrem(x, y) = divrem(x, y, RoundToZero)
# fldmod(x,y) = divrem(x, y, RoundDown)
divrem(x::Robust32, y::Robust32) = Rob32(divrem(value64(x), value64(y), RoundToZero))
divrem(x::Robust32, y::Float64) = Rob32(divrem(value64(x), y, RoundToZero))
divrem(x::Robust32, y::Float32) = Rob32(divrem(value64(x), Float64(y), RoundToZero))
divrem(x::Float64, y::Robust32) = Rob32(divrem(x, value64(y), RoundToZero))
divrem(x::Float32, y::Robust32) = Rob32(divrem(Float64(x), value64(y), RoundToZero))

fldmod(x::Robust32, y::Robust32) = Rob32(fldmod(value64(x), value64(y), RoundDown))
fldmod(x::Robust32, y::Float64) = Rob32(fldmod(value64(x), y, RoundDown))
fldmod(x::Robust32, y::Real) = Rob32(fldmod(value64(x), Float64(y), RoundDown))
fldmod(x::Float64, y::Robust32) = Rob32(fldmod(x, value64(y), RoundDown))
fldmod(x::Real, y::Robust32) = Rob32(fldmod(Float64(x), value64(y), RoundDown))
 

for F in (:hypot, :clamp)
  @eval begin
    $F(x::Robust32, y::Robust32, z::Robust32) = Robust32($F(value64(x), value64(y), value64(z)))
    $F(x::Robust32, y::Robust32, z::T) where {T<:Real} = $F(promote(x,y,z)...)
    $F(x::T, y::Robust32, z::Robust32) where {T<:Real} = $F(promote(x,y,z)...)
    $F(x::Robust32, y::T, z::Robust32) where {T<:Real} = $F(promote(x,y,z)...)
    $F(x::Robust32, y::T, z::T) where {T<:Real} = $F(promote(x,y,z)...)
    $F(x::T, y::Robust32, z::T) where {T<:Real} = $F(promote(x,y,z)...)
    $F(x::T, y::T, z::Robust32) where {T<:Real} = $F(promote(x,y,z)...)
  end  
end

for F in (:floor, :ceil, :trunc)
  @eval Base.$F(x::Robust32) = Rob32($F(value64(x)))
    for T in (:BigInt, :Int128, :Int64, :Int32, :Int16, :Int8) 
      @eval Base.$F(::Type{$T}, x::Robust32) = $T(Float32($F(value64(x))))
    end  
end

Base.round(x::Robust32; digits=0, sigdigits=0, base=10) = 
    iszero(sigdigits) ? Rob32(round(value64(x), digits=digits, base=base)) : Rob32(round(value64(x), sigdigits=sigdigits, base=base))

for T in (:BigInt, :Int128, :Int64, :Int32, :Int16, :Int8) 
  @eval Base.round(::Type{$T}, x::Robust32) = round($T, Float64(x))
end  

include("minmax.jl")

for F in (:abs2, :acos, :acosd, :acosh, :acot, :acotd, :acoth, :acsc, :acscd, :acsch, :asec,
          :asecd, :asech, :asin, :asind, :asinh, :atan, :atand, :atanh, :cos, :cosc,
          :cosd, :cosh, :cospi, :cot, :cotd, :coth, :csc, :cscd, :csch, :deg2rad,
          :exp, :exp10, :exp2, :expm1, :log, :log10, :log1p, :log2, :mod2pi,
          :rad2deg, :rem2pi, :sec, :secd, :sech, :sin, :sinc, :sind, :sinh,
          :sinpi, :tan, :tand, :tanh)
    @eval $F(x::Robust32) = Rob32($F(value64(x)))
end

atan(x::Robust32, y::Robust32) = Rob32(atan(value64(x), value64(y)))

for F in (:modf, :sincos, :sincosd) # , :sincospi)
  @eval function $F(x::Robust32)
            s, c = $F(value64(x))
            return Rob32(s), Rob32(c)
         end
end

for F in (:fma, :muladd)
  @eval begin
    $F(x::Robust32, y::Robust32, z::Robust32) = Robust32($F(value64(x), value64(y), value64(z)))
    $F(x::Robust32, y::Robust32, z::T) where {T<:Real} = $F(promote(x,y,z)...)
    $F(x::T, y::Robust32, z::Robust32) where {T<:Real} = $F(promote(x,y,z)...)
    $F(x::Robust32, y::T, z::Robust32) where {T<:Real} = $F(promote(x,y,z)...)
    $F(x::Robust32, y::T, z::T) where {T<:Real} = $F(promote(x,y,z)...)
    $F(x::T, y::Robust32, z::T) where {T<:Real} = $F(promote(x,y,z)...)
    $F(x::T, y::T, z::Robust32) where {T<:Real} = $F(promote(x,y,z)...)
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

include("randnums.jl")
include("linearalgebra.jl")
include("specialfuncs.jl")

end  # Robust32s
