third32  = Float32(1)/Float32(3)
thirdR32 = Robust32(1)/Robust32(3)
fourthirds = Float32(4)/Float32(3)
two32  = Float32(2)
twoR32 = Robust32(2)

@testset "add/sub" begin
  @test two32 + two32 == twoR32 + twoR32  
  @test two32 - third32 == twoR32 - thirdR32  
  # @test fourthirds - (two32 - 2*third32) > 2.75 * (fourthirds - (twoR32 - 2*thirdR32))
end

@testset "mul/div" begin
  @test abs(2/9 -  (third32 * 2*third32)) > 9 *( abs(2/9 - (thirdR32 * 2*thirdR32)))  
  @test two32 - third32 == twoR32 - thirdR32  
  @test fourthirds - (two32 - 2*third32) > 2.75 * (fourthirds - (twoR32 - 2*thirdR32))
end

@testset "root/power" begin
  @test abs(two32 - sqrt(sqrt(two32))^4) > abs(twoR32 - sqrt(sqrt(twoR32))^4)
  @test abs(two32 - sqrt(cbrt(two32))^6) > abs(twoR32 - sqrt(cbrt(twoR32))^6)
end
