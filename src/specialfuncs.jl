import SpecialFunctions
const SF = SpecialFunctions

# functions of 1 argument returning a floating point Real value
for F in (:besselj0, :besselj1, :bessely0, :bessely1, 
          :cosint, :dawson, :ellipe, :ellipk, 
          :erf, :erfc, :erfcinv, :erfcx, :erfi, :erfinv, 
          :gamma, :invdigamma, :logabsgamma, 
          :logerfc, :logerfcx, :sinint)
  @eval SF.$F(x::Robust32) = Robust32(SF.$F(Float64(x)))
end

# functions of 1 argument returning a floating point Complex value
for F in (:airyaiprime, :airyaiprimex, :airyaix, :airybi, 
          :airybiprime, :airybiprimex, :airybix, 
          :besselj0, :besselj1, :bessely0, :bessely1, 
          :dawson, :erf, :erfc, :erfcx, :erfi, 
          :gamma, :loggamma)
  @eval SF.$F(x::ComplexR32) = ComplexR32(SF.$F(ComplexF64(x)))
end

# functions of 2 arguments returning a floating point Complex value
for F in (:besseli, :besselix, :besselj, :besseljx,
          :besselk, :besselkx, :bessely, :besselyx)
  @eval begin
    SF.$F(x::Robust32, y::ComplexR32) = ComplexR32(SF.$F(Float64(x), ComplexF64(y)))
    SF.$F(x::Robust32, y::ComplexF64) = ComplexR32(SF.$F(Float64(x), y))
    SF.$F(x::Float64, y::ComplexR32) =  ComplexR32(SF.$F(x, ComplexF64(y)))                
  end                             
end

# functions of 2 arguments returning a floating point 2-tuple value
for F in (:besselj, :bessely)
  @eval begin
    SF.$F(x::Int32, y::Robust32) = Robust32(SF.$F(x, Float64(y)))
  end                  
end

for F in (:gamma_inc_inv,)
  @eval begin
    SF.$F(x::Robust32, y::Robust32, z::Robust32) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
  end                             
end
