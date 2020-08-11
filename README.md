# Float32Rs

A more robust Float32 that preserves `float` performance.

"An old Rule-of-Thumb offers the simplest way to evade numerical embarrassment.
 Perform your compututations carrying somewhat more than twice the precision
 of your data and somewhat more than twice the precision of your results.
 This rule has long served statistics, optimization, root-finding, geometry,
 and differential equations. Rare exceptions exist, of course."
 - (edited together for content from multiple writings of W. Kahan)
 
"An old Rule-of-Thumb renders roundoff extremely unlikely to cause embarrassment:
In all intermediate computation, perform arithmetic carrying somewhat more than
twice as many significant digits as are trusted in the data and desired in the final result.
This rule has long served statistics, optimization, root-finding, geometry, structural
analysis and differential equations. Rare exceptions exist, of course. Nothing is perfect."
- W. Kahan, 2016

".. the simplest way to evade numerical embarrassment is to perform computation carrying extravagantly
more precision throughout than you think necessary, and pray that it is enough. Usually somewhat
more than twice the precision you trust in the data and seek in the results is enough."
- W. Kahan, "How Futile are Mindless Assessments of Roundoff in Floating Point Computation", 2006
