#=
   High performance work-alike rewrites of `min`, `max`, `minmax`, `maxmin`.
   
   For use with `Float32`, `Float64`, and other IEEEFloats.
=#
   
jsmin(a, b) = signbit(a-b) ? a : b + (a-a)
jsmax(a, b) =  signbit(b-a) ? a : b + (a-a)
jsminmax(a, b) = (jsmin(a,b), jsmax(a,b))
jsmaminx(a, b) = (jsmax(a,b), jsmin(a,b))

