#=     
     provide(x) uses reinterpret
     
     provide(x::Vector{Robust32}), provide(x::Matrix{Robust32})
        presents an Array of Robust32s as an Array of Float64s

     provide(x::Vector{Float64}), provide(x::Matrix{Float64})
        presents an Array of Float64s as an Array of Robust32s

     rewrap(x) uses unsafe_wrap
     
     rewrap(x::Vector{Robust32}, rewrap(x::Matrix{Robust32})
        obtains an Array of Robust32s as an Array of Float64s

     rewrap(x::Vector{Float64}, rewrap(x::Matrix{Float64})
        obtains an Array of Float64s  as an Array of Robust32s

=#

for (R,F) in ((:Robust32, :Float64), (:ComplexR32, :ComplexF64))
  for A in (:Vector, :Matrix)
    @eval begin
        provide(x::$A{$R}) = reinterpret($F, x)
        provide(x::$A{$F}) = reinterpret($R, x)

        rewrap(x::$A{$R}) = unsafe_wrap($A{$F}, cvtptr($F, x), size(x))             
        rewrap(x::$A{$F}) = unsafe_wrap($A{$R}, cvtptr($R, x), size(x))
    end
  end
end

@inline cvtptr(::Type{T}, m::Array{T,N}) where {N,T} =
    convert(Ptr{T}, pointer(m,1))
