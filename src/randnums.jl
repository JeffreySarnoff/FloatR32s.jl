function rand(rng::AbstractRNG, ::Random.SamplerTrivial{Random.CloseOpen01{Robust32}})
    hi, lo  = rand(rng, T, 2)
    if hi === zero(T)
        if lo === zero(T)
            return zero(DoubleFloat(T))
        end
        hi, lo = lo, hi
    end
    frlo, xplo  = frexp(lo)
    xplo = Base.exponent(hi) - min(1, fld(xplo,4)) - abs(Base.exponent(eps(hi)))
    lo = ldexp(frlo, xplo)
    lo = rand(rng, Bool) ? lo : -lo

    DoubleFloat(hi, lo)
end

function rand(rng::AbstractRNG, ::Random.SamplerTrivial{Random.CloseOpen01{Complex{Robust32}}})
    re = rand(rng, Robust32)
    im = rand(rng, Robust32)
    return Complex{Robust32}(re, im)
end

function randpm(::Type{Robust32})
    r = rand(Robust32)
    r = rand(Bool) ? r : -r
    return r
end

function randpm(rng::MersenneTwister, ::Type{Robust32})
    r = rand(rng, Robust32)
    r = rand(rng, Bool) ? r : -r
    return r
end

function randpm(::Type{Complex{Robust32}})
    re = randpm(Robust32)
    im = randpm(Robust32)
    return Complex{Robust32}(re, im)
end

function randpm(rng::MersenneTwister, ::Type{Complex{Robust32}})
    re = randpm(rng, Robust32)
    im = randpm(rng, Robust32)
    return Complex{Robust32}(re, im)
end

function randpm(::Type{Robust32}, n::Int)
    rs = rand(Robust32, n)
    sgns = ones(n) .- (2 .* rand(Bool, n))
    return rs .* sgns
end

function randpm(::Type{Complex{Robust32}}, n::Int)
    res = rand(Robust32, n)
    ims = rand(Robust32, n)
    o = ones(n)
    rsgns = o .- (2 .* rand(Bool,n))
    isgns = o .- (2 .* rand(Bool,n))
    res = res .* rsgns
    ims = ims .* isgns
    return map((re,im)->Complex{Robust32}(re,im), res, ims)
end

function randpm(rng::MersenneTwister, ::Type{Robust32}, n::Int)
    rs = rand(rng, Robust32, n)
    sgns = ones(n) .- (2 .* rand(rng, Bool, n))
    return rs .* sgns
end

function randpm(rng::MersenneTwister, ::Type{Complex{Robust32}}, n::Int)
    res = rand(rng, Robust32, n)
    ims = rand(rng, Robust32, n)
    o = ones(n)
    rsgns = o .- (2 .* rand(rng, Bool,n))
    isgns = o .- (2 .* rand(rng, Bool,n))
    res = res .* rsgns
    ims = ims .* isgns
    return map((re,im)->Complex{Robust32}(re,im), res, ims)
end

# normal variates

function randn(rng::AbstractRNG, ::Type{Robust32})
    urand1, urand2 = rand(rng, Robust32, 2)
    urand1 = urand1 + urand1 - 1
    urand2 = urand2 + urand2 - 1
    s = urand1*urand1 + urand2*urand2

    while s >= 1 || s === 0
        urand1, urand2 = rand(rng, Robust32, 2)
        urand1 = urand1 + urand1 - 1
        urand2 = urand2 + urand2 - 1
        s = urand1*urand1 + urand2*urand2
    end

    s = sqrt( -log(s) / s )
    return (urand1 + urand2) * s
end
