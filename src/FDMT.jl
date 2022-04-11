module FDMT

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
end

struct OutputBlock{T1} <: AbstractMatrix{T1}
    data::AbstractMatrix{T1}
    # Dispersion delay (s) for first DM trial
    y_min::Real
    # Dispersion delay (s) for last DM trial
    y_max::Real
    # Dimension time is in
    t_dim::Integer
end

# Methods for custom array types
Base.size(ib::InputBlock) = size(ib.data)
Base.size(ob::OutputBlock) = size(ob.data)
Base.getindex(ib::InputBlock, i::Int, j::Int) = ib.data[i, j]
Base.getindex(ob::OutputBlock, i::Int, j::Int) = ob.data[i, j]

# Default constructors
function InputBlock(data, f_ch1, f_chn, t_samp)
    @assert f_ch1 >= f_chn "Channel 1 must be higher than Channel N"
    return InputBlock(data, f_ch1, f_chn, t_samp, 1)
end

function OutputBlock(ib::InputBlock{T}, y_min, y_max) where {T}
    n_trials = y_max - y_min + 1
    # We want to keep time in the same axis here
    data = ib.t_dim == 1 ? zeros(T, n_samp(ib), n_trials) : zeros(T, n_trials, n_samp(ib))
    return OutputBlock(data, y_min, y_max, ib.t_dim)
end

# The dimension that isn't time
f_dim(ib::InputBlock) = ib.t_dim == 1 ? 2 : 1

# The dimension that isn't time
dm_dim(ob::OutputBlock) = ob.t_dim == 1 ? 2 : 1

# Number of frequency channels
n_ch(ib::InputBlock) = size(ib)[f_dim(ib)]

# Number of time samples
n_samp(ib::InputBlock) = size(ib)[ib.t_dim]
n_samp(ob::OutputBlock) = size(ob)[ob.t_dim]

# Number of DM trials
n_trials(ob::OutputBlock) = size(ob)[dm_dim(ob)]

# Frequency range
freqs(ib::InputBlock) = range(; start=ib.f_ch1, stop=ib.f_chn, length=n_ch(ib))

# Channel spacing
f_step(ib::InputBlock) = (ib.f_ch1 - ib.f_chn) / n_ch(ib)

# Dispersion delay across the bands covered by the input block
Δkdisp(ib::InputBlock) = KDM * (ib.f_chn^-2 - ib.f_ch1^-2)

# Natural DM step of FDMT for this block
dm_step(ib::InputBlock) = ib.t_samp / Δkdisp(ib)

function split(ib::InputBlock)
    # The split occurs at the channel where half the dispersion
    # delay has been accounted for
    f = (0.5 * ib.f_chn^-2 + 0.5 * ib.f_ch1^-2)^-0.5
    i = round(Int, (ib.f_ch1 - f) / f_step(ib))

    # New frequency ranges
    head_end = freqs(ib)[i]
    tail_start = freqs(ib)[i + 1]

    # New data slices
    ch = n_ch(ib)
    head_slice = ib.t_dim == 1 ? view(ib.data, :, 1:i) : view(ib.data, 1:i, :)
    tail_slice = ib.t_dim == 1 ? view(ib.data, :, (i + 1):ch) : view(ib.data, (i + 1):ch, :)

    # New blocks
    head = InputBlock(head_slice, ib.f_ch1, head_end, ib.t_samp, ib.t_dim)
    tail = InputBlock(tail_slice, tail_start, ib.f_chn, ib.t_samp, ib.t_dim)

    return head, tail
end

function transform_recursive(block, y_min, y_max)
    out = OutputBlock(block, y_min, y_max)
    # Base Case
    if n_ch(block) == 1
        if block.t_dim == 1
            out.data[:,1] = block.data[:,1]
        else
            out.data[1,:] = block.data[1,:]
        end
        return out
    end
    # Split
    head, tail = split(block)
    # Transform
    y_min_head = trunc(Int, y_min * Δkdisp(head) / Δkdisp(block) + 0.5)
    y_max_head = trunc(Int, y_max * Δkdisp(head) / Δkdisp(block) + 0.5)
    transformed_head = transform_recursive(head, y_min_head, y_max_head)

    y_min_tail = trunc(Int, y_min * Δkdisp(tail) / Δkdisp(block) + 0.5)
    y_max_tail = trunc(Int, y_max * Δkdisp(tail) / Δkdisp(block) + 0.5)
    transformed_tail = transform_recursive(tail, y_min_tail, y_max_tail)

    _n_samp = n_samp(out)

    # Merge
    @turbo for y in y_min:y_max
        # yh = delay across head band
        yh = floor(Int, y * Δkdisp(head) / Δkdisp(block) + 0.5)
        # yt = delay across tail band
        yt = floor(Int, y * Δkdisp(tail) / Δkdisp(block) + 0.5)
        # yb = delay at interface between head and tail
        yb = y - yh - yt
        ih = yh - transformed_head.y_min + 1
        it = yt - transformed_tail.y_min + 1
        i = y - out.y_min + 1
        # Update FIXME
        for j in 1:_n_samp
            j_shift = mod1(j + yh - yb, _n_samp)
            out.data[j, i] = transformed_head[j, ih] + transformed_tail[j_shift, it]
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
