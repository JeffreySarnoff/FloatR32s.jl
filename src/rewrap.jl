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

provide(m::Vector{Float64}) = reinterpret(Robust32, m)
provide(m::Matrix{Float64}) = reinterpret(Robust32, m)

provide(m::Vector{Robust32}) = reinterpret(Float64, m)
provide(m::Matrix{Robust32}) = reinterpret(Float64, m)

provide(m::Vector{ComplexF64}) = reinterpret(ComplexR32, m)
provide(m::Matrix{ComplexF4})  = reinterpret(ComplexR32, m)

provide(m::Vector{ComplexR32}) = reinterpret(ComplexF64, m)
provide(m::Matrix{ComplexR32}) = reinterpret(ComplexF64, m)


@inline ptr(::Type{T}, m::Array{T,N}) where {N,T} = convert(Ptr{T}, pointer(m,1))

rewrap(m::Vector{Float64}) =
    unsafe_wrap(Array{Robust32,1}, ptr(Robust32, m), length(m))
rewrap(m::Vector{Robust32}) =
    unsafe_wrap(Array{Float64,1}, ptr(Float64, m), length(m))
rewrap(m::Matrix{Float64}) =
    unsafe_wrap(Array{Robust32,2}, ptr(Robust32, m), size(m))
rewrap(m::Matrix{Robust32}) =
    unsafe_wrap(Array{Float64,2}, ptr(Float64, m), size(m))

rewrap(m::Vector{ComplexF64}) =
    unsafe_wrap(Array{ComplexR32,1}, ptr(ComplexR32, m), length(m))
rewrap(m::Vector{ComplexR32}) =
    unsafe_wrap(Array{ComplexF64,1}, ptr(ComplexF64, m), length(m))
rewrap(m::Matrix{ComplexF64}) =
    unsafe_wrap(Array{ComplexR32,2}, ptr(ComplexR32, m), size(m))
rewrap(m::Matrix{ComplexR32}) =
    unsafe_wrap(Array{ComplexF64,2}, ptr(ComplexF64, m), size(m))

#=

rewrap(m::Vector{Float64}) =
    unsafe_wrap(Array{Robust32,1}, Ptr{Robust32}(pointer(m,1)), length(m))
rewrap(m::Vector{Robust32}) =
    unsafe_wrap(Array{Float64,1}, Ptr{Float64}(pointer(m,1)), length(m))
rewrap(m::Matrix{Float64}) =
    unsafe_wrap(Array{Robust32,2}, Ptr{Robust32}(pointer(m,1)), size(m))
rewrap(m::Matrix{Robust32}) =
    unsafe_wrap(Array{Float64,2}, Ptr{Float64}(pointer(m,1)), size(m))

rewrap(m::Vector{ComplexF64}) =
    unsafe_wrap(Array{ComplexR32,1}, Ptr{ComplexR32}(pointer(m,1)), length(m))
rewrap(m::Vector{ComplexR32}) =
    unsafe_wrap(Array{ComplexF64,1}, Ptr{ComplexF64}(pointer(m,1)), length(m))
rewrap(m::Matrix{ComplexF64}) =
    unsafe_wrap(Array{ComplexR32,2}, Ptr{ComplexR32}(pointer(m,1)), size(m))
rewrap(m::Matrix{ComplexR32}) =
    unsafe_wrap(Array{ComplexF64,2}, Ptr{ComplexF64}(pointer(m,1)), size(m))

=#    
