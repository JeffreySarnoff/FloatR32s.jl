import SpecialFunctions
const SF = SpecialFunctions

for F in (:besselj0, :besselj1, :bessely0, :bessely1, 
          :cosint, :dawson, :ellipe, :ellipk, 
          :erf, :erfc, :erfcinv, :erfcx, :erfi, :erfinv, 
          :gamma, :invdigamma, :logabsgamma, 
          :logerfc, :logerfcx, :sinint)
  @eval SF.$F(x::Robust32) = Robust32(SF.$F(Float64(x)))
end

for F in (:airyaiprime, :airyaiprimex, :airyaix, :airybi, 
          :airybiprime, :airybiprimex, :airybix, 
          :besselj0, :besselj1, :bessely0, :bessely1, 
          :dawson, :erf, :erfc, :erfcx, :erfi, 
          :gamma, :loggamma)
  @eval SF.$F(x::ComplexR32) = ComplexR32(SF.$F(ComplexF64(x)))
end

for F in (:besseli, :besselix, :besselj, :besselj, :besseljx,
          :besselk, :besselkx, :bessely, :bessely, :besselyx)
  @eval begin
    SF.$F(x::Robust32, y::ComplexR64) = ComplexR32(SF.$F(Float64(x), ComplexF64(y)))
    SF.$F(x::Robust32, y::ComplexF64) = ComplexR32(SF.$F(Float64(x), y))
    SF.$F(x::Float64, y::ComplexR64) =  ComplexR32(SF.$F(x, ComplexF64(y)))                
  end                             
end

for F in (:besselj, :bessely)
  @eval SF.$F(x::Integer, y::Robust32) = Robust32(SF.$F(Int32(x), Float64(y))) 
end

for F in (:gamma_inc_inv)
  @eval begin
    SF.$F(x::Robust32, y::Robust32, z::Robust32) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
  end                             
end
