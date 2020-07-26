"""
    Float42 <: AbstractFloat
    
A well supported floating point type
  -  75% of the performance relative to Float32
  - 125% of the accuracy is relative to Float32
""" Float42

struct Float42 <: AbstractFloat
    val::Float64
end

value(x::Float42) = x.val

Base.Float64(x::Float42) = value(x)
Base.Float32(x::Float42) = Float32(value(x))
Base.Float16(x::Float42) = Float32(value(x))
