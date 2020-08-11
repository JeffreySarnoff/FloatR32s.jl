# Float32Rs

A more robust Float32 that preserves `float` performance.

> An old Rule-of-Thumb offers the simplest way to evade numerical embarrassment.
 Perform your compututations carrying somewhat more than twice the precision
 of your data and somewhat more than twice the precision of you seek in your results.
 This rule has long served statistics, optimization, root-finding, geometry,
 and differential equations. Rare exceptions exist, of course.
 > (adapted from W. Kahan)
 
 This package implements that Rule-of-Thumb in a highly performant manner.
 To offer the desired performance, this package works with `Float32` data
 and provides `Float32` results.  This is handled automatically; so as long
 as you are committed to working with data and results that are `Float32s`,
 all is well.
 
 
".. the simplest way to evade numerical embarrassment is to perform computation carrying extravagantly
more precision throughout than you think necessary, and pray that it is enough. Usually somewhat
more than twice the precision you trust in the data and seek in the results is enough."
- W. Kahan, "How Futile are Mindless Assessments of Roundoff in Floating Point Computation", 2006
