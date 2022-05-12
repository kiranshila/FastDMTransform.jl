using FastDMTransform, NPZ, Test, Statistics

function standardize(A; dims=1)
    μ = mean(A; dims=dims)
    σ = std(A; mean=μ, dims=dims)
    return @. (A - μ) / σ
end

@testset "FastDMTransform.jl" begin
    pulse = npzread("../pulse.npz")
    dedispersed = standardize(npzread("../dedispersed.npz"))
    dedisped = standardize(fdmt(pulse, 1500, 1200, 1e-3, 0, 2000))
    @test findmax(dedisped) == findmax(dedispersed)
end
