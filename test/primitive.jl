two32  = Float32(2)
twoR32 = Robust32(2)
sqrt_two32 = sqrt(two32)
sqrt_twoR32 = sqrt(twoR32)

@testset "low level functions" begin
  @test hash(sqrt_two32) == hash(sqrt_twoR32)
  @test string(sqrt_two32) == string(sqrt_twoR32)

  @test maxintfloat(Float32) == maxintfloat(Robust32)
  @test floatmin(Float32) == floatmin(Robust32)
  @test typemax(Float32) == typemax(Robust32)
end  
