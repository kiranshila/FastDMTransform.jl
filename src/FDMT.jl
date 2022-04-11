module FDMT

using LoopVectorization

const KDM = 4.148808e3 # MHz^2 s pc^-1 cm^3

struct InputBlock{T1} <: AbstractMatrix{T1}
    data::AbstractMatrix{T1}
    # Frequency of channel 1 and N
    f_ch1::Real
    f_chn::Real
    # Time step in seconds
    t_samp::Real
    # Dimension time is in
    t_dim::Integer
    # Dimension freq is in
    f_dim::Integer
    # Dispersion delay across the bands covered by the input block
    Δkdisp::Real
end

struct OutputBlock{T1} <: AbstractMatrix{T1}
    data::AbstractMatrix{T1}
    # Dispersion delay (s) for first DM trial
    y_min::Real
    # Dispersion delay (s) for last DM trial
    y_max::Real
    # Dimension time is in
    t_dim::Integer
    # Dimension DM is in
    dm_dim::Integer
end

# Methods for custom array types
Base.size(ib::InputBlock) = size(ib.data)
Base.size(ob::OutputBlock) = size(ob.data)
Base.getindex(ib::InputBlock, i::Int, j::Int) = ib.data[i, j]
Base.getindex(ob::OutputBlock, i::Int, j::Int) = ob.data[i, j]

# Default constructors
function InputBlock(data, f_ch1, f_chn, t_samp, t_dim)
    @assert f_ch1 >= f_chn "Channel 1 must be higher than Channel N"
    Δkdisp = KDM * (f_chn^-2 - f_ch1^-2)
    f_dim = t_dim == 1 ? 2 : 1
    return InputBlock(data, f_ch1, f_chn, t_samp, t_dim, f_dim, Δkdisp)
end

function OutputBlock(ib::InputBlock{T}, y_min, y_max) where {T}
    n_trials = y_max - y_min + 1
    # We want to keep time in the same axis here
    n_samp = size(ib)[ib.t_dim]
    data = ib.t_dim == 1 ? zeros(T, n_samp, n_trials) : zeros(T, n_trials, n_samp)
    dm_dim = ob.t_dim == 1 ? 2 : 1
    return OutputBlock(data, y_min, y_max, ib.t_dim, ib.f_dim)
end

# Channel spacing
f_step(ib::InputBlock) = (ib.f_ch1 - ib.f_chn) / size(ib)[ib.f_dim]

# Natural DM step of FDMT for this block
dm_step(ib::InputBlock) = ib.t_samp / ib.Δkdisp

function split(ib::InputBlock)
    # The split occurs at the channel where half the dispersion
    # delay has been accounted for
    f = (0.5 * ib.f_chn^-2 + 0.5 * ib.f_ch1^-2)^-0.5
    i = round(Int, (ib.f_ch1 - f) / f_step(ib))

    # New frequency ranges
    step = f_step(ib)
    head_end = ib.f_ch1 - i * step
    tail_start = ib.f_ch1 - (i + 1) * step

    # New data slices
    ch = size(ib)[ib.f_dim]
    head_slice = ib.t_dim == 1 ? view(ib.data, :, 1:i) : view(ib.data, 1:i, :)
    tail_slice = ib.t_dim == 1 ? view(ib.data, :, (i + 1):ch) : view(ib.data, (i + 1):ch, :)

    # New blocks
    head = InputBlock(head_slice, ib.f_ch1, head_end, ib.t_samp, ib.t_dim)
    tail = InputBlock(tail_slice, tail_start, ib.f_chn, ib.t_samp, ib.t_dim)

    return head, tail
end

circmod(x, y) = mod(x - 1, y) + 1

function transform_recursive(block, y_min, y_max)
    out = OutputBlock(block, y_min, y_max)
    # Base Case
    if size(block)[block.f_dim] == 1
        if block.t_dim == 1
            out.data[:, 1] = block.data[:, 1]
        else
            out.data[1, :] = block.data[1, :]
        end
        return out
    end
    # Split
    head, tail = split(block)
    # Transform
    y_min_head = trunc(Int, y_min * head.Δkdisp / block.Δkdisp + 0.5)
    y_max_head = trunc(Int, y_max * head.Δkdisp / block.Δkdisp + 0.5)
    transformed_head = transform_recursive(head, y_min_head, y_max_head)

    y_min_tail = trunc(Int, y_min * tail.Δkdisp / block.Δkdisp + 0.5)
    y_max_tail = trunc(Int, y_max * tail.Δkdisp / block.Δkdisp + 0.5)
    transformed_tail = transform_recursive(tail, y_min_tail, y_max_tail)

    n_samp = size(out)[out.t_dim]

    # Merge
    @turbo for y in y_min:y_max
        # yh = delay across head band
        yh = floor(Int, y * head.Δkdisp / block.Δkdisp + 0.5)
        # yt = delay across tail band
        yt = floor(Int, y * tail.Δkdisp / block.Δkdisp + 0.5)
        # yb = delay at interface between head and tail
        yb = y - yh - yt
        ih = yh - transformed_head.y_min + 1
        it = yt - transformed_tail.y_min + 1
        i = y - out.y_min + 1
        # Update FIXME
        for j in 1:n_samp
            j_shift = circmod(j + yh - yb, n_samp)
            out.data[j, i] = transformed_head.data[j, ih] +
                             transformed_tail.data[j_shift, it]
        end
    end

    return out
end

function transform(data, f_ch1, f_chn, t_samp, dm_min, dm_max, t_dim)
    @assert dm_min >= 0 "Minimum DM must be zero"
    @assert dm_max >= dm_min "Maximum DM must be greater than the minimum"

    # Build input block
    block = InputBlock(data, f_ch1, f_chn, t_samp, t_dim)

    # Convert DMs to delays in sample space
    y_min = trunc(Int, dm_min / dm_step(block))
    y_max = ceil(Int, dm_max / dm_step(block))

    # Perform transformation
    return transform_recursive(block, y_min, y_max)
end

export transform

end
