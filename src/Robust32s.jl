module Robust32s

export Robust32

import LinearAlgebra

struct As64 end # internal use only

struct Robust32 <: AbstractFloat
    val::Float64
  
    Robust32(x::Float64) = new(Float64(Float32(x)))
    Robust32(::Type{As64}, x::Float64) = new(x)
end

value64(x::Robust32) = x.val
value32(x::Robust32) = Float32(x.val)

# internal use only
Rob32(x::Float64) = Robust32(As64, x)
Flo64(x::Robust32) = x.val

Robust32(x::Robust32) = x # idempotency

Base.show(io::IO, x::Robust32) = print(io, value32(x))
Base.string(x::Robust32) = string(value32(x))

#=
  to maintain the package intent correctly
     explicit construction of T requires the target become a Float32
=#
float6432(x) = Float64(Float32(x))

for T in (:BigFloat, :Float64, :Float32, :Float16)
  @eval begin
    Base.$T(x::Robust32) = $T(value32(x))
    Robust32(x::$T) = Rob32(float6432(x))
    Base.convert(::Type{$T}, x::Robust32) = $T(x)
    Base.promote_rule(::Type{Robust32}, ::Type{$T}) = Robust32
  end
end

for T in (:BigInt, :Int128, :Int64, :Int32, :Int16, :Int8,
                   :UInt128, :UInt64, :UInt32, :UInt16, :UInt8)
  @eval begin
    Base.$T(x::Robust32) = $T(value32(x))
    Robust32(x::$T) = Rob32(float6432(x))
    Base.convert(::Type{$T}, x::Robust32) = $T(x)
    Base.promote_rule(::Type{Robust32}, ::Type{$T}) = Robust32
  end
end

const Robust32_0 = Rob32(0.0)
const Robust32_1 = Rob32(1.0)

Robust32(x::Bool) = x ? Robust32_1 : Robust32_0

#=
  to maintain the package intent correctly
     implicit construction of a Float64 does not require the target become a Float32
=# 
#=
Base.convert(::Type{Robust32}, x::Float64) = Robust32(x)
Base.convert(::Type{Robust32}, x::T) where {T<:Union{Float32, Float16}} = Robust32(x)
Base.convert(::Type{Robust32}, x::T) where {T<:Signed} = Robust32(x)
Base.convert(::Type{Robust32}, x::T) where {T<:Unsigned} = Robust32(x)

Base.convert(::Type{Float64}, x::Robust32) = Float64(x)
Base.convert(::Type{Float32}, x::Robust32) = Float32(x)
=#

#=
  to maintain the package intent correctly
    some primitive operations must be taken with respect to Float32
=#
Base.eps(x::Robust32) = eps(value32(x))
Base.significand(x::Robust32) = significand(value32(x))
Base.exponent(x::Robust32) = exponent(value32(x))
Base.sign(x::Robust32) = exponent(value32(x))
Base.iszero(x::Robust32) = iszero(value32(x))
Base.isone(x::Robust32) = isone(value32(x))
Base.isfinite(x::Robust32) = isfinite(value32(x))
Base.issubnormal(x::Robust32) = issubnormal(value32(x))
Base.isinf(x::Robust32) = isinf(value64(x))
Base.isnan(x::Robust32) = isnan(value64(x))

Base.signbit(x::Robust32) = signbit(value32(x))

Base.zero(::Type{Robust32}) = Robust32_0
Base.one(::Type{Robust32}) = Robust32_1
Base.zero(x::Robust32) = zero(Robust32)
Base.one(x::Robust32) = one(Robust32)

Base.frexp(x::Robust32) = frexp(value32(x)) # ??????????? and ldexp

for F in (:-, :abs, :sign, :inv, :sqrt, :cbrt)
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

for F in (:modf, :sincos, :sincosd) # , :sincospi)
  @eval function $F(x::Robust32)
            s, c = Base.Math.$F(value64(x))
            return Rob32(s), Rob32(c)
         end
end

# ?????? @evalpoly

function Base.evalpoly(x::Robust32, p::NTuple{N, Robust32}) where {N}
    Rob32(evalpoly(value64(x), map(value64, p)))
end
function Base.evalpoly(x::T, p::NTuple{N, Robust32}) where {T,N}
    Rob32(evalpoly(Float64(x), map(value64, p)))
end
function Base.evalpoly(x::Robust32, p::NTuple{N, T}) where {T,N}
    Rob32(evalpoly(value64(x), p))
end

function Base.Vector{Float64}(x::Vector{Robust32})
    n = length(x)
    res = Vector{Float64}(undef, n)
    for idx in 1:n
       @inbounds res[idx] = x[idx].val
    end
    return res
end

Base.Vector{Robust32}(x::Matrix{Float64}) = Rob32.(x)

function Base.Matrix{Float64}(x::Matrix{Robust32})
    rows, cols = size(x)
    res = Matrix{Float64}(undef, rows, cols)
    for r in 1:rows
        for c in 1:cols
           @inbounds res[r,c] = x[r,c].val
        end
    end
    return res
end

Base.Matrix{Robust32}(x::Matrix{Float64}) = Rob32.(x)

function Base.Array{N,Float64}(x::Array{N,Robust32}) where {N}
    dims = size(x)
    res = Array{Float64}(undef, dims)
    for idx in eachindex(x)
       @inbounds res[idx] = x[idx].val
    end
    return res
end

Base.Array{N,Robust32}(x::Array{N,Float64}) = Rob32.(x)

for F in (:tr, :det)
    @eval LinearAlgebra.$F(x::Matrix{Robust32}) = Rob32($F(Matrix{Float64}(x))
end

for F in (:inv, :sqrt, :exp, :log, 
          :sin, :tan, :csc, :cot, :asin, :atan, :acsc, :acot)
    @eval LinearAlgebra.$F(x::Matrix{Robust32}) = Rob32.($F(Matrix{Float64}(x))
end

LinearAlgebra.dot(x::Array{N,Robust32}) where {N} = Rob32(dot(Array{N,Float64}(x))

end  # Robust32s

#=
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
=#
