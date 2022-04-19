module FDMT

const KDM = 4.148808e3 # MHz^2 s pc^-1 cm^3

# Get SNR from a DM result
function standardize(A)
    μ = mean(A; dims=1)
    σ = std(A; mean=μ, dims=1)
    return @. (A - μ) / σ
end

include("cpu.jl")
include("gpu.jl")

export standardize

end