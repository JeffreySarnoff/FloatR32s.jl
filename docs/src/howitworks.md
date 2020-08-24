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
 > (adapted from W. Kahan)
 
 This package implements that Rule-of-Thumb in a highly performant manner.
 To offer the desired performance, this package works with `Float32` data
 and provides `Float32` results.  This is handled automatically; so as long
 as you are committed to working with `Float32` data and obtaining `Float32`
 results, all is well.

