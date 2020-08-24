## How It Works

### The Essential Requirements

- You are comfortable working with all input data as `Float32s`.
- You are comfortable reporting all the results using `Float32s`.
- You respect `Robust32` values (just let them work, no peeking).

### The Pervasive Conventions

- any information introduced to the computation is converted to `Float32s`.
- programmed calculations and computations are performed using `Float64s`.
- all information retrieved from the computation is converted to `Float32s`.

## Why It Works


> An old Rule-of-Thumb offers the simplest way to evade numerical embarrassment.
 Perform your compututations carrying somewhat more than twice the precision
 of your data and somewhat more than twice the precision of you seek in your results.
 This rule has long served statistics, optimization, root-finding, geometry,
 and differential equations. Rare exceptions exist, of course.
 > (from W. Kahan, amended)
 
 This package implements that Rule-of-Thumb in a highly performant manner.
 We work from `Float32` data and provide `Float32` results.
 Internally, `Float64` math is used. 
 `Float64s` have 2x + 5 more significant bits than `Float32s`.
 So, all is well.

