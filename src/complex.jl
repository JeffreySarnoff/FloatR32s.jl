negate(x::ComplexR32) = ComplexR32(-x.re, -x.im)
flipsign(x::ComplexR32, y::T) where {T<:Real} = signbit(y) ? negate(x) : x

Bool(x::ComplexR32) = iszero(x.im) ? isone(x.re) : false
ComplexR32(x::Bool) = x ? one(ComplexR32) : zero(ComplexR32)

const ComplexR32_0 = ComplexR32(Robust32_0, Robust32_0)
const ComplexR32_1 = ComplexR32(Robust32_1, Robust32_0)
const ComplexR32_2 = ComplexR32(Robust32_2, Robust32_0)

zero(::Type{ComplexR32}) = ComplexR32_0
one(::Type{ComplexR32})  = ComplexR32_1
two(::Type{ComplexR32})  = ComplexR32_2
zero(x::ComplexR32) = ComplexR32_0
one(x::ComplexR32)  = ComplexR32_1
two(x::ComplexR32)  = ComplexR32_2

iszero(x::ComplexR32) = x === ComplexR32_0
isone(x::ComplexR32)  = x === ComplexR32_1
istwo(x::ComplexR32)  = x === ComplexR32_2


