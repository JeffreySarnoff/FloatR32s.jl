import SpecialFunctions
const SF = SpecialFunctions

# functions of 1 argument returning a value of Float type
for F in (:besselj0, :besselj1, :bessely0, :bessely1, 
          :cosint, :dawson, :ellipe, :ellipk, 
          :erf, :erfc, :erfcinv, :erfcx, :erfi, :erfinv, 
          :gamma, :invdigamma, :logabsgamma, 
          :logerfc, :logerfcx, :sinint)
  @eval SF.$F(x::Robust32) = Robust32(SF.$F(Float64(x)))
end

# functions of 1 argument returning a value of Complex type
for F in (:airyaiprime, :airyaiprimex, :airyaix, :airybi, 
          :airybiprime, :airybiprimex, :airybix, 
          :besselj0, :besselj1, :bessely0, :bessely1, 
          :dawson, :erf, :erfc, :erfcx, :erfi, 
          :gamma, :loggamma)
  @eval SF.$F(x::ComplexR32) = ComplexR32(SF.$F(ComplexF64(x)))
end

# functions of 2 arguments returning a value of Complex type
for F in (:besseli, :besselix, :besselj, :besseljx,
          :besselk, :besselkx, :bessely, :besselyx)
  @eval begin
    SF.$F(x::Robust32, y::ComplexR32) = ComplexR32(SF.$F(Float64(x), ComplexF64(y)))
    SF.$F(x::Robust32, y::ComplexF64) = ComplexR32(SF.$F(Float64(x), y))
    SF.$F(x::Float64, y::ComplexR32) =  ComplexR32(SF.$F(x, ComplexF64(y)))                
    SF.$F(x::AbstractFloat, y::ComplexR32) =  ComplexR32(SF.$F(x, ComplexF64(y)))                
  end                             
end

# functions of 2 arguments returning a value of type typeof(y)
for F in (:besselj, :bessely)
  @eval begin    SF.$F(x::Robust32, y::Robust32, z::Robust32) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
    SF.$F(x::Int32, y::Robust32) = Robust32(SF.$F(x, Float64(y)))
    SF.$F(x::Int32, y::ComplexR32) = ComplexR32(SF.$F(x, ComplexF64(y)))
  end                  
end

# functions of 3 arguments returning a value of Float type
for F in (:gamma_inc_inv,)
  @eval begin
    SF.$F(x::Robust32, y::Robust32, z::Robust32) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
    SF.$F(x::Robust32, y::Robust32, z::AbstractFloat) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
    SF.$F(x::Robust32, y::AbstractFloat, z::Robust32) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
    SF.$F(x::AbstractFloat, y::Robust32, z::Robust32) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
    SF.$F(x::Robust32, y::AbstractFloat, z::AbstractFloat) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
    SF.$F(x::AbstractFloat, y::AbstractFloat, z::Robust32) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))                    
    SF.$F(x::AbstractFloat, y::Robust32, z::AbstractFloat) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))                    
  end                             
end

# functions of 3 arguments returning a value of NTuple{2,Float} type
for F in (:gamma_inc)
  @eval begin
    SF.$F(x::Robust32, y::Robust32, IND::Integer) = Robust32(SF.$F(Float64(x), Float64(y), IND))
    SF.$F(x::Robust32, y::AbstractFloat, IND::Integer) = Robust32(SF.$F(Float64(x), Float64(y), IND))
    SF.$F(x::AbstractFloat, y::Robust32, IND::Integer) = Robust32(SF.$F(Float64(x), Float64(y), IND))
  end
end

# functions of 3 arguments returning a value of NTuple{2,Float} type
for F in (:beta_inc)
  @eval begin
    SF.$F(x::Robust32, y::Robust32, z::Robust32) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
    SF.$F(x::Robust32, y::AbstractFloat, z::Robust32) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
    SF.$F(x::Robust32, y::Robust32, z::AbstractFloat) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
    SF.$F(x::Robust32, y::AbstractFloat, z::AbstractFloat) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
    SF.$F(x::AbstractFloat, y::Robust32, z::AbstractFloat) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
    SF.$F(x::AbstractFloat, y::AbstractFloat, z::Robust32) = Robust32(SF.$F(Float64(x), Float64(y), Float64(z)))
  end
end

# functions of 4 arguments returning a value of Float type
for F in (:ncbeta, :ncF)
  @eval begin
    SF.$F(a::Robust32, b::Robust32, c::Robust32, d::Robust32) =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
  end
end

# functions of 4 arguments returning a value of NTuple{2,Float} type
for F in (:beta_inc_inv,)
  @eval begin
    SF.$F(a::Robust32, b::Robust32, c::Robust32, d::Robust32) =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
  end
end

#=
# functions of 4 arguments returning a value of Float type
for F in (:ncbeta, :ncF)
  @eval begin
    SF.$F(a::Robust32, b::Robust32, c::Robust32, d::Robust32) =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32, b::T, c::T, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T, b::Robust32, c::T, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T, b::T, c::Robust32, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T, b::T, c::T, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32, b::Robust32, c::T, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32 b::T, c::Robust32, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32 b::T, c::T, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T b::Robust32, c::T, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T b::T, c::Robust32, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32 b::Robust32, c::Robust32, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32 b::Robust32, c::T, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32 b::T, c::Robust32, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T, b::Robust32, c::Robust32, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
end

# functions of 4 arguments returning a value of NTuple{2,Float} type
for F in (:beta_inc_inv,)
  @eval begin
    SF.$F(a::Robust32, b::Robust32, c::Robust32, d::Robust32) =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32, b::T, c::T, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T, b::Robust32, c::T, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T, b::T, c::Robust32, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T, b::T, c::T, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32, b::Robust32, c::T, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32 b::T, c::Robust32, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32 b::T, c::T, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T b::Robust32, c::T, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T b::T, c::Robust32, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32 b::Robust32, c::Robust32, d::T) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32 b::Robust32, c::T, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::Robust32 b::T, c::Robust32, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
    SF.$F(a::T, b::Robust32, c::Robust32, d::Robust32) where {T<:AbstractFloat} =
      Robust32(SF.$F(Float64(a), Float64(b), Float64(c), Float64(d)))
  end
end

=#
