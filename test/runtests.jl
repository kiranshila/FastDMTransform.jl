using FastDMTransform, NPZ, Test

@testset "FastDMTransform.jl" begin
    pulse = npzread("../pulse.npz")
    dedispersed = npzread("../dedispersed.npz")
    dedisped = fdmt(pulse,1500,1200,1e-3,0,2000)
    @test isapprox(dedisped,dedispersed; atol=0.01)
end
