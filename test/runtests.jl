using Robust32s, LinearAlgebra, Test

streq(a, b) = string(a) === string(b)

include("primitive.jl")
include("arithmetic.jl")
include("vector.jl")
include("matrix.jl")
include("elementary.jl")
