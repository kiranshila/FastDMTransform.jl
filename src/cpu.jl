using LoopVectorization, FFTW, LinearAlgebra, Base.Threads
import Base.Threads.@spawn

struct InputBlock{T<:Real,
                  M<:AbstractMatrix{T},
                  TF<:Real,
                  TT<:Real,
                  TKD<:Real} <: AbstractMatrix{T}
    data::M
    # Frequency of channel 1 and N
    f_ch1::TF
    f_chn::TF
    # Time step in seconds
    t_samp::TT
    # Dispersion delay across the bands covered by the input block
    Δkdisp::TKD
end

struct OutputBlock{T<:Real,M<:Matrix{T}} <: AbstractMatrix{T}
    data::M
    # Dispersion delay (s) for first DM trial
    y_min::Int
    # Dispersion delay (s) for last DM trial
    y_max::Int
end

# Methods for custom array types
Base.size(ib::InputBlock) = size(ib.data)
Base.size(ob::OutputBlock) = size(ob.data)
Base.getindex(ib::InputBlock, i::Int, j::Int) = ib.data[i, j]
Base.getindex(ob::OutputBlock, i::Int, j::Int) = ob.data[i, j]

# Default constructors
function InputBlock(data, f_ch1, f_chn, t_samp)
    @assert f_ch1 >= f_chn "Channel 1 must be higher than Channel N"
    Δkdisp = KDM * (f_chn^-2 - f_ch1^-2)
    return InputBlock(data, f_ch1, f_chn, t_samp, Δkdisp)
end

function OutputBlock(ib::InputBlock{T}, y_min, y_max) where {T}
    n_trials = y_max - y_min + 1
    n_samp = size(ib)[1]
    data = zeros(T, n_samp, n_trials)
    return OutputBlock(data, y_min, y_max)
end

# Channel spacing
f_step(ib::InputBlock) = (ib.f_ch1 - ib.f_chn) / size(ib)[2]

# Natural DM step of FastDMTransform for this block
dm_step(ib::InputBlock) = ib.t_samp / ib.Δkdisp

function split(ib::InputBlock)
    # The split occurs at the channel where half the dispersion
    # delay has been accounted for
    f = (0.5 * ib.f_chn^-2 + 0.5 * ib.f_ch1^-2)^-0.5
    i = round(Int, (ib.f_ch1 - f) / f_step(ib))

    # New frequency ranges
    step = f_step(ib)
    head_start, head_end = promote(ib.f_ch1, ib.f_ch1 - i * step)
    tail_start, tail_end = promote(ib.f_ch1 - (i + 1) * step, ib.f_chn)

    # New data slices
    ch = size(ib)[2]
    head_slice = view(ib.data, :, 1:i)
    tail_slice = view(ib.data, :, (i + 1):ch)

    # New blocks
    head = InputBlock(head_slice, head_start, head_end, ib.t_samp)
    tail = InputBlock(tail_slice, tail_start, tail_end, ib.t_samp)

    return head, tail
end

circmod(x, y) = mod(x - 1, y) + 1

# We really want go get rid of all these allocations, it's slowing us down a lot
function transform_recursive(block::InputBlock, y_min::Int, y_max::Int)
    out = OutputBlock(block, y_min, y_max)
    # Base Case
    if size(block)[2] == 1
        out.data[:, 1] = block.data[:, 1]
        return out
    end
    # Split
    head, tail = split(block)
    # Transform
    y_min_head = round(Int, y_min * head.Δkdisp / block.Δkdisp + 0.5)
    y_max_head = round(Int, y_max * head.Δkdisp / block.Δkdisp + 0.5)
    transformed_head_task = @spawn transform_recursive(head, y_min_head, y_max_head)

    y_min_tail = round(Int, y_min * tail.Δkdisp / block.Δkdisp + 0.5)
    y_max_tail = round(Int, y_max * tail.Δkdisp / block.Δkdisp + 0.5)
    transformed_tail_task = @spawn transform_recursive(tail, y_min_tail, y_max_tail)

    n_samp = size(block)[1]
    j_range = 1:n_samp

    transformed_head = fetch(transformed_head_task)
    transformed_tail = fetch(transformed_tail_task)

    # Merge
    @tturbo for y in y_min:y_max
        # yh = delay across head band
        yh = round(Int, y * head.Δkdisp / block.Δkdisp + 0.5)
        # yt = delay across tail band
        yt = round(Int, y * tail.Δkdisp / block.Δkdisp + 0.5)
        # yb = delay at interface between head and tail
        yb = y - yh - yt
        ih = yh - transformed_head.y_min + 1
        it = yt - transformed_tail.y_min + 1
        i = y - out.y_min + 1
        # Update
        for j in j_range
            j_shift = circmod(j + yh - yb, n_samp)
            out.data[j, i] = transformed_head.data[j, ih] +
                             transformed_tail.data[j_shift, it]
        end
    end

    return out
end

"""
    fdmt(data,f_ch1,f_chn,t_samp,dm_min,dm_max)

Computes the Fast DM Transform for dynamic spectra `data`. This spectra must have time in
the first axis and cover frequencies from `f_ch1` to `f_chn` in descending order (in MHz).
The sample spacing in time (seconds) is given by `t_samp` and the transform covers DMs `dm_min`
to `dm_max`
"""
function fdmt(data::AbstractMatrix, f_ch1::Real, f_chn::Real, t_samp::Real,
              dm_min::Real, dm_max::Real)
    @assert dm_min >= 0 "Minimum DM must be zero"
    @assert dm_max >= dm_min "Maximum DM must be greater than the minimum"

    # Build input block
    block = InputBlock(data, f_ch1, f_chn, t_samp)
    dm_step = t_samp / block.Δkdisp

    # Convert DMs to delays in sample space
    y_min = trunc(Int, dm_min / dm_step)
    y_max = ceil(Int, dm_max / dm_step)

    # Perform transformation
    output = transform_recursive(block, y_min, y_max)
    return output.data
end

export fdmt