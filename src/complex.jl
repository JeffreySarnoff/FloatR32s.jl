negate(x::ComplexR32) = ComplexR32(-x.re, -x.im)
flipsign(x::ComplexR32, y::T) where {T<:Real} = signbit(y) ? negate(x) : x

Bool(x::ComplexR32) = iszero(x.im) ? isone(x.re) : false
ComplexR32(x::Bool) = x ? one(ComplexR32) : zero(ComplexR32)
