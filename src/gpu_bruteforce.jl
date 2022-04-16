using Plots, CUDA, BenchmarkTools, SIGPROC, LoopVectorization, StaticArrays, Statistics,
      FDMT
using DimensionalData.Dimensions
using DimensionalData

const KDM = 4.1488064239e3

@dim DM YDim "Dispersion Measure"

@inline circmod(x, y) = mod(x - 1, y) + 1

function Δt(f_min, f_max, DM, δt)
    return round(Int, KDM * DM * (f_min^-2 - f_max^-2) / δt) + 1
end


dm_min = 0
dm_max = 500
n_dm = 5795

fb = Filterbank("/home/kiran/Downloads/candidate_ovro_20200428.fil")
pulse = cu(fb.data.data)
n_samp, n_chan = size(pulse)
freqs = cu(collect(fb.data.dims[2]))
dms = cu(collect(range(dm_min, dm_max; length=n_dm)))
δt = step(fb.data.dims[1])
f_min, f_max = extrema(freqs)

function dedisp_kernel!(output, source, plan)
    # Pull out constants
    n_samp, n_chan = size(source)

    samp_idx = (blockIdx().y - 1) * blockDim().x + threadIdx().x
    dm_idx = blockIdx().x

    if samp_idx > n_samp
        return nothing
    end

    # Copy shifts into dynamic memory
    shifts = @cuDynamicSharedMem(UInt32, n_chan)
    @inbounds for chan_idx in (threadIdx().x):(blockDim().x):n_chan
        shifts[chan_idx] = plan[chan_idx, dm_idx]
    end

    sync_threads()

    # Sum
    @inbounds for chan_idx in 1:n_chan
        shifted_samp_idx = circmod(samp_idx + shifts[chan_idx],n_samp)
        output[samp_idx, dm_idx] += source[shifted_samp_idx, chan_idx]
    end

    # Kernel is side-effecting, doesn't return
    return nothing
end

function dedisp(source, plan)
    n_samp, n_chan = size(source)
    _, n_dm = size(plan)
    output = CUDA.zeros(n_samp, n_dm)

    # Compile kernel and grab capabilities
    kernel = @cuda launch = false dedisp_kernel!(output, source, plan)
    config = launch_configuration(kernel.fun)
    threads = config.threads # 1024 on a 2080Ti
    blocks = (n_dm, cld(n_samp, threads))

    # Run kernel
    kernel(output, pulse, plan; threads=threads, blocks=blocks,
           shmem=sizeof(UInt32) * n_chan)
    return output
end

function standardize(A)
    μ = mean(A; dims=1)
    σ = std(A; mean=μ, dims=1)
    return @. (A - μ) / σ
end

plan = Δt.(freqs, f_max, dms', δt)

function pipeline(source)
    return standardize(dedisp(source, plan))
end

heatmap(Array(pipeline(pulse)))
