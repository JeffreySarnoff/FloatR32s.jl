#=
   My faster rewrites of `min`, `max`, `minmax`, `maxmin`.
   
   These are processing isomorphs. NaNs are handled properly.

   For use with `Float32`, `Float64`, and other IEEEFloats.
=#
   
jsmin(a, b) = signbit(a-b) ? a : b + (a-a)
jsmax(a, b) =  signbit(b-a) ? a : b + (a-a)
jsminmax(a, b) = (jsmin(a,b), jsmax(a,b))
jsmaminx(a, b) = (jsmax(a,b), jsmin(a,b))

#=
   their utilzation
=#

min(x::Robust32, y::Robust32) = Rob32(jsmin(value64(x), value64(y)))
max(x::Robust32, y::Robust32) = Rob32(jsmax(value64(x), value64(y)))
minmax(x::Robust32, y::Robust32) = Rob32(jsminmax(value64(x), value64(y)))
maxmin(x::Robust32, y::Robust32) = Rob32(jsminmax(value64(x), value64(y)))

min(x::Robust32, y::Float64) = Rob32(jsmin(value64(x), y))
max(x::Robust32, y::Float64) = Rob32(jsmax(value64(x), y))
minmax(x::Robust32, y::Float64) = Rob32(jsminmax(value64(x), y))
maxmin(x::Robust32, y::Float64) = Rob32(jsminmax(value64(x), y))

min(x::Float64, y::Robust32) = Rob32(jsmin(x, value64(y)))
max(x::Float64, y::Robust32) = Rob32(jsmax(x, value64(y)))
minmax(x::Float64, y::Robust32) = Rob32(jsminmax(x, value64(y)))
maxmin(x::Float64, y::Robust32) = Rob32(jsminmax(x, value64(y)))

min(x::Robust32, y::Real) = Rob32(jsmin(value64(x), Float64(y)))
max(x::Robust32, y::Real) = Rob32(jsmax(value64(x), Float64(y)))
minmax(x::Robust32, y::Real) = Rob32(jsminmax(value64(x), Float64(y)))
maxmin(x::Robust32, y::Real) = Rob32(jsminmax(value64(x), Float64(y)))

min(x::Real, y::Robust32) = Rob32(jsmin(Float64(x), value64(y)))
max(x::Real, y::Robust32) = Rob32(jsmax(Float64(x), value64(y)))
minmax(x::Real, y::Robust32) = Rob32(jsminmax(Float64(x), value64(y)))
maxmin(x::Real, y::Robust32) = Rob32(jsminmax(Float64(x), value64(y)))
