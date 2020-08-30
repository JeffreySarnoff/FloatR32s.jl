#=
   My faster rewrites of `min`, `max`, `minmax`, `maxmin`.
   
   These are processing isomorphs for all normal uses.
   Julia's NaN, NaN32 are handled properly.
   Negatively signed NaNs, which cannot appear unless
   intentionally generated, are not IEEE754 compliant.

   For use with `Float32`, `Float64`.
=#

const FastFloat = Union{Float32, Float64}

function fastermin(a::T, b::T) where {T<:FastFloat}
    a_b = a - b
    (signbit(a_b) || isnan(a)) ? a : b
end

function fastermax(a::T, b::T) where {T<:FastFloat}
    b_a = b - a
    (signbit(b_a) || isnan(a)) ? a : b
end

fasterminmax(a::T, b::T) where {T<:FastFloat} =
    (min(a, b), max(a, b))

fastermaxmin(a::T, b::T) where {T<:FastFloat} =
    (max(a, b), min(a, b))

#=
   their utilzation
=#

min(x::Robust32, y::Robust32) = Rob32(fastermin(value64(x), value64(y)))
max(x::Robust32, y::Robust32) = Rob32(fastermax(value64(x), value64(y)))
minmax(x::Robust32, y::Robust32) = Rob32(fasterminmax(value64(x), value64(y)))
maxmin(x::Robust32, y::Robust32) = Rob32(fastmaxmin(value64(x), value64(y)))

min(x::Robust32, y::Float64) = Rob32(fastermin(value64(x), y))
max(x::Robust32, y::Float64) = Rob32(fastermax(value64(x), y))
minmax(x::Robust32, y::Float64) = Rob32(fasterminmax(value64(x), y))
maxmin(x::Robust32, y::Float64) = Rob32(fasterminmax(value64(x), y))

min(x::Float64, y::Robust32) = Rob32(fastermin(x, value64(y)))
max(x::Float64, y::Robust32) = Rob32(fastermax(x, value64(y)))
minmax(x::Float64, y::Robust32) = Rob32(fasterminmax(x, value64(y)))
maxmin(x::Float64, y::Robust32) = Rob32(fasterminmax(x, value64(y)))

min(x::Robust32, y::Real) = Rob32(fastermin(value64(x), Float64(y)))
max(x::Robust32, y::Real) = Rob32(fastermax(value64(x), Float64(y)))
minmax(x::Robust32, y::Real) = Rob32(fasterminmax(value64(x), Float64(y)))
maxmin(x::Robust32, y::Real) = Rob32(fasterminmax(value64(x), Float64(y)))

min(x::Real, y::Robust32) = Rob32(fastermin(Float64(x), value64(y)))
max(x::Real, y::Robust32) = Rob32(fastermax(Float64(x), value64(y)))
minmax(x::Real, y::Robust32) = Rob32(fasterminmax(Float64(x), value64(y)))
maxmin(x::Real, y::Robust32) = Rob32(fasterminmax(Float64(x), value64(y)))
