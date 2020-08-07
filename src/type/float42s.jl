"""
    Float42 <: AbstractFloat
    
A well supported floating point type
  -  75% of the performance relative to Float32
  - 125% of the accuracy is relative to Float32
""" Float42

struct Float42 <: AbstractFloat
    val::Float64
end

Float42(x::Float42) = x

value(x::Float42) = x.val

Float42(x::Float32) = Float42(Float64(x))
Base.convert(::Type{Float42}, x::Float32) = Float42(Float64(x))
Float32(x::Float42) = Float32(value(x))
Base.convert(::Type{Float32}, x::Float42) = Float32(value(x))

Float42(x::Float64) = Float42(x)
Base.convert(::Type{Float42}, x::Float64) = Float42(Float64(Float32(x)))
Float64(x::Float42) = Float64(Float32(value(x)))
Base.convert(::Type{Float64}, x::Float42) = value(x)

Base.BigFloat(x::Float42) = BigFloat(Float32(x))
Base.Float64(x::Float42)  = Float64(Float32(x))
Base.Float32(x::Float42)  = Float32(value(x))
Base.Float16(x::Float42)  = Float16(value(x))

Base.BigInt(x::Float42) = BigInt(Float32(x))
Base.Int128(x::Float42) = Int128(Float32(x))
Base.Int64(x::Float42)  = Int64(Float32(x))
Base.Int32(x::Float42)  = Int32(Float32(x))
Base.Int16(x::Float42)  = Int16(Float16(x))
Base.Int8(x::Float42)   = Int8(Float32(x))


Base.show(io::IO, x::Float42) = print(io, Float32(x))
Base.parse(::Type{Float42}, x::String) = Float42(parse(Float32, x))
