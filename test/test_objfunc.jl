module test_objfunc

using Test, CashParadox

@testset "testing objective functions" begin
    
@testset "testing specific output values for eqn_Regime201610 " begin

  @test isapprox(eqn_Regime201610(1.0,0.0,0.0,0.0,0.0,0),-1.0,atol=1e-5)

  @test isapprox(eqn_Regime201610(2.0,0.5,-3.0,4.0,1.0,0.5),-1.1464466094067267262,atol=1e-5)

end


end


end