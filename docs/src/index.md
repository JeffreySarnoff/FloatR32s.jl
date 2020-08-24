# Robust32s.jl
### Increased accuracy using Float32 data

#### Copyright ©2020 by Jeffrey Sarnoff.  Released under the MIT License.

----
https://jeffreysarnoff.github.io/Robust32s.jl/dev/
----

## How To Use This Package

There are two exported types: `Robust32` and `ComplexR32`.  Use them as you would `Float32` and `ComplexF32`.
It Just Works.

### The Only Requirements

- You are comfortable working with your data as `Float32s`.
- You are comfortable reporting your results as `Float32s`.

### Bringing in data

Your raw data `rawdata` is stored as `Ints` or as `Float32s` or as `Float64s` or other `Real` type.

`data = Robust32.(rawdata);`

### Computing on the data

Just do it.

`processed_data = your_computation(data)`

### Collecting results

Your processed data `processed_data` is stored as `Robust32s`. Report the results as `Float32s`.

`results = Float32.(processed_data);`


