# Robust32s

## A more robust Float32 that preserves `float` performance.

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://jeffreysarnoff.github.io/Robust32s.jl/dev)

----


> An old Rule-of-Thumb offers the simplest way to evade numerical embarrassment.
 Perform your compututations carrying somewhat more than twice the precision
 of your data and somewhat more than twice the precision of you seek in your results.
 This rule has long served statistics, optimization, root-finding, geometry,
 and differential equations. Rare exceptions exist, of course.
 > (adapted from W. Kahan)
 
 -----
 
 
 This package implements that Rule-of-Thumb in a highly performant manner.
 To offer the desired performance, this package works with `Float32` data
 and provides `Float32` results.  This is handled automatically; so as long
 as you are committed to working with `Float32` data and obtaining `Float32`
 results, all is well.
 
 ### exports
 
 - `Robust32`   a robust 32bit floating point type
 - `ComplexR32` a robust 32bit complex floating point type (named like `ComplexF32`)
 
 ### installation
 
 ```julia
 julia> using Pkg
 julia> Pkg.add("Robust32s")
 ```
 
## Basic Examples
 
 ```julia
using FloatR32s

julia> a, b = sqrt.(Float32.((2.0, 0.5)))
(1.4142135f0, 0.70710677f0)

julia> c = a * b    # product of Float32s
0.99999994f0

julia> a, b = sqrt.(Robust32.((2.0, 0.5)))
(1.4142135f0, 0.70710677f0)

julia> c = a * b    # product of FloatR32s
1.0f0
```


----

----

".. the simplest way to evade numerical embarrassment is to perform computation carrying extravagantly
more precision throughout than you think necessary, and pray that it is enough. Usually somewhat
more than twice the precision you trust in the data and seek in the results is enough."
- W. Kahan, "How Futile are Mindless Assessments of Roundoff in Floating Point Computation", 2006
