#=
  to maintain the package intent correctly
     implicit construction of a Float64 does not require the target become a Float32

  to maintain the package intent correctly
    some primitive operations must be taken with respect to Float32

  to maintain the package intent correctly
     explicit construction of T requires the target become a Float32
=#

"""
    FloatR32s

A module for robust Float32 computation.

Exports: `FloatR32`, `ComplexR32`
"""
module FloatR32s

export FloatR32, ComplexR32

using LinearAlgebra

struct As64 end # internal use only

struct FloatR32 <: AbstractFloat
    val::Float64
  
    FloatR32(x::Float64) = new(Float64(Float32(x)))
    FloatR32(::Type{As64}, x::Float64) = new(x)
end

value64(x::FloatR32) = x.val
value32(x::FloatR32) = Float32(x.val)

# internal use only
Rob32(x::Float64) = FloatR32(As64, x)

FloatR32(x::FloatR32) = x # idempotency

Base.show(io::IO, x::FloatR32) = print(io, value32(x))
Base.string(x::FloatR32) = string(value32(x))

const ComplexR32 = Complex{FloatR32}

value64(x::ComplexR32) = (value64(x.re), value64(x.im))
value32(x::ComplexR32) = (value32(x.re), value32(x.im))

Base.Float64(x::FloatR32) = Float64(value32(x))
Base.convert(::Type{Float64}, x::FloatR32) = Float64(x)
Base.promote_rule(::Type{FloatR32}, ::Type{Float64}) = FloatR32

for T in (:BigFloat, :Float32, :Float16)
  @eval begin
    Base.$T(x::FloatR32) = $T(value32(x))
    FloatR32(x::$T) = Rob32(Float64(Float32(x)))
    Base.convert(::Type{$T}, x::FloatR32) = $T(x)
    Base.promote_rule(::Type{FloatR32}, ::Type{$T}) = FloatR32
  end
end

for T in (:BigInt, :Int128, :Int64, :Int32, :Int16, :Int8,
                   :UInt128, :UInt64, :UInt32, :UInt16, :UInt8)
  @eval begin
    Base.$T(x::FloatR32) = $T(value32(x))
    FloatR32(x::$T) = Rob32(Float64(Float32(x)))
    Base.convert(::Type{$T}, x::FloatR32) = $T(x)
    Base.promote_rule(::Type{FloatR32}, ::Type{$T}) = FloatR32
  end
end

const FloatR32_0 = Rob32(0.0)
const FloatR32_1 = Rob32(1.0)

FloatR32(x::Bool) = x ? FloatR32_1 : FloatR32_0

Base.eps(x::FloatR32) = eps(value32(x))
Base.significand(x::FloatR32) = significand(value32(x))
Base.exponent(x::FloatR32) = exponent(value32(x))
Base.sign(x::FloatR32) = exponent(value32(x))
Base.iszero(x::FloatR32) = iszero(value32(x))
Base.isone(x::FloatR32) = isone(value32(x))
Base.isfinite(x::FloatR32) = isfinite(value32(x))
Base.issubnormal(x::FloatR32) = issubnormal(value32(x))
Base.isinf(x::FloatR32) = isinf(value64(x))
Base.isnan(x::FloatR32) = isnan(value64(x))

Base.signbit(x::FloatR32) = signbit(value32(x))

Base.zero(::Type{FloatR32}) = FloatR32_0
Base.one(::Type{FloatR32}) = FloatR32_1
Base.zero(x::FloatR32) = zero(FloatR32)
Base.one(x::FloatR32) = one(FloatR32)

Base.frexp(x::FloatR32) = frexp(value32(x)) # ??????????? and ldexp

for F in (:-, :abs, :inv, :sqrt, :cbrt)
  @eval Base.$F(x::FloatR32) = Rob32($F(value64(x)))
end

for F in (:(==), :(!=), :(<), :(<=), :(>), :(>=), :isless, :isequal)
  @eval begin
    Base.$F(x::FloatR32, y::FloatR32) = $F(value32(x), value32(y))
    Base.$F(x::FloatR32, y::Real) = $F(promote(x,y)...)
    Base.$F(x::Real, y::FloatR32) = $F(promote(x,y)...)
  end  
end

for F in (:+, :-, :*, :/, :\, :hypot, :copysign, :flipsign)
  @eval begin
    Base.$F(x::FloatR32, y::FloatR32) = Rob32($F(value64(x), value64(y)))
    Base.$F(x::FloatR32, y::Real) = $F(promote(x,y)...)
    Base.$F(x::Real, y::FloatR32) = $F(promote(x,y)...)
  end  
end

for F in (:hypot, :clamp)
  @eval begin
    Base.$F(x::FloatR32, y::FloatR32, z::FloatR32) = FloatR32($F(value(x), value(y), value(z)))
    Base.$F(x::FloatR32, y::Real, z::Real) = $F(promote(x,y,z)...)
    Base.$F(x::Real, y::FloatR32, z::Real) = $F(promote(x,y,z)...)
    Base.$F(x::Real, y::Real, z::FloatR32) = $F(promote(x,y,z)...)
  end  
end

for F in (:abs2, :acos, :acosd, :acosh, :acot, :acotd, :acoth, :acsc, :acscd, :acsch, :asec,
          :asecd, :asech, :asin, :asind, :asinh, :atan, :atand, :atanh, :cos, :cosc,
          :cosd, :cosh, :cospi, :cot, :cotd, :coth, :csc, :cscd, :csch, :deg2rad,
          :exp, :exp10, :exp2, :expm1, :log, :log10, :log1p, :log2, :mod2pi,
          :rad2deg, :rem2pi, :sec, :secd, :sech, :sin, :sinc, :sind, :sinh,
          :sinpi, :tan, :tand, :tanh)
    @eval Base.Math.$F(x::FloatR32) = Rob32($F(value64(x)))
end

for F in (:modf, :sincos, :sincosd) # , :sincospi)
  @eval function $F(x::FloatR32)
            s, c = Base.Math.$F(value64(x))
            return Rob32(s), Rob32(c)
         end
end

# ?????? @evalpoly

function Base.evalpoly(x::FloatR32, p::NTuple{N, FloatR32}) where {N}
    Rob32(evalpoly(value64(x), map(value64, p)))
end
function Base.evalpoly(x::T, p::NTuple{N, FloatR32}) where {T,N}
    Rob32(evalpoly(Float64(x), map(value64, p)))
end
function Base.evalpoly(x::FloatR32, p::NTuple{N, T}) where {T,N}
    Rob32(evalpoly(value64(x), p))
end

rewrap(m::Vector{Float64}) =
    unsafe_wrap(Array{FloatR32,1}, Ptr{FloatR32}(pointer(m,1)), length(m))
rewrap(m::Vector{FloatR32}) =
    unsafe_wrap(Array{Float64,1}, Ptr{Float64}(pointer(m,1)), length(m))
rewrap(m::Matrix{Float64}) =
    unsafe_wrap(Array{FloatR32,2}, Ptr{FloatR32}(pointer(m,1)), size(m))
rewrap(m::Matrix{FloatR32}) =
    unsafe_wrap(Array{Float64,2}, Ptr{Float64}(pointer(m,1)), size(m))

rewrap(m::Vector{ComplexF64}) =
    unsafe_wrap(Array{ComplexR32,1}, Ptr{ComplexR32}(pointer(m,1)), length(m))
rewrap(m::Vector{ComplexR32}) =
    unsafe_wrap(Array{ComplexF64,1}, Ptr{ComplexF64}(pointer(m,1)), length(m))
rewrap(m::Matrix{ComplexF64}) =
    unsafe_wrap(Array{ComplexR32,2}, Ptr{ComplexR32}(pointer(m,1)), size(m))
rewrap(m::Matrix{ComplexR32}) =
    unsafe_wrap(Array{ComplexF64,2}, Ptr{ComplexF64}(pointer(m,1)), size(m))

include("linearalgebra.jl")

end  # FloatR32s
