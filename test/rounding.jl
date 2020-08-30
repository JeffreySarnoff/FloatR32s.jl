macro roundsig(source, sigdigs, target)
  :(streq( round(Robust32($source), sigdigits=$sigdigs), Robust32($target) ))
end
macro roundsig(source, sigdigs, base, target)
  :(streq( round(Robust32($source), sigdigits=$sigdigs, base=$base), Robust32($target) ))
end
macro rounddig(source, digs, target)
  :(streq( round(Robust32($source), digits=$digs), Robust32($target) ))
end
macro rounddig(source, digs, base, target)
  :(streq( round(Robust32($source), digits=$digs, base=$base), Robust32($target) ))
end

@testset "significant digits" begin
    @test @roundsig(123.456, 1, 100.0)
    @test @roundsig(123.456, 3, 123.0)
    @test @roundsig(123.456, 5, 123.46)
    @test @roundsig(123.456, 8, 2, 123.5)
    @test @roundsig(123.456, 2, 4, 128.0)
    
    @test @roundsig(0.0, 1, 0.0)
    @test @roundsig(-0.0, 5, -0.0)
     
    @test @roundsig(1.2, 2, 1.2)
    @test @roundsig(1.2, 3, 1.2)
    @test @roundsig(1.2, 4, 1.2)
    
    @test @roundsig(0.6, 1, 0.6)
    @test @roundsig(0.6, 2, 0.6)
    @test @roundsig(0.6, 3, 0.6)
 
    @test @roundsig(7.262839104539736f0, 2, 7.3)
    @test @roundsig(7.262839104539736f0, 3, 7.26)
    @test @roundsig(7.262839104539736f0, 4, 7.263)
 
    @test isinf(round(Robust32(Inf), sigdigits=3))
    @test isnan(round(Robust32(NaN), sigdigits=3))
 end
