using LoopVectorization, FFTW, LinearAlgebra

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

# Natural DM step of FDMT for this block
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
    transformed_head = transform_recursive(head, y_min_head, y_max_head)

    y_min_tail = round(Int, y_min * tail.Δkdisp / block.Δkdisp + 0.5)
    y_max_tail = round(Int, y_max * tail.Δkdisp / block.Δkdisp + 0.5)
    transformed_tail = transform_recursive(tail, y_min_tail, y_max_tail)

    n_samp = size(block)[1]
    j_range = 1:n_samp

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

function fft_fdmt_iterate(input::AbstractArray{T,3}, Δt_max, n_chan_source, f_min, f_max,
                          iteration) where {T}
    n_samp, n_chan, _ = size(input)
    # Calculate this iterations frequency range
    δf = (2^iteration) * (f_max - f_min)/n_chan_source
    dF = (f_max - f_min)/n_chan
    # Maximum number of time shifts in sample channels for this iteration
    Δt = ceil(Int, Δt_max * (f_min^-2 - (f_min + δf)^-2) / (f_min^-2 - f_max^-2))
    # Initialize new output
    f_jumps = Int(n_chan / 2)
    output = zeros(T, n_samp, f_jumps, Δt + 1)
    # See remark in paper about this random correction
    correction = dF / 2
    # Compute shift vector
    Δt_shift = ceil(Int,
                    Δt_max * (f_min^-2 - (f_min + δf / 2 + correction)^-2) /
                    (f_min^-2 - f_max^-2)) + 2
    shift_row = fft(Matrix{T}(I, n_samp, Δt_shift), 1)
    # Perform the iteration
    for i_F in 1:f_jumps
        f_start = (f_max - f_min) / f_jumps * (i_F - 1) + f_min
        f_end = (f_max - f_min) / f_jumps * i_F + f_min
        f_middle_lower = (f_end - f_start) / 2 + f_start - correction
        f_middle_upper = (f_end - f_start) / 2 + f_start + correction
        Δt_local = ceil(Int, Δt_max * (f_start^-2 - f_end^-2) / (f_min^-2 - f_max^-2))
        for i_Δt in 0:Δt_local
            Δt_middle_lower = round(Int,
                                    i_Δt * (f_middle_lower^-2 - f_start^-2) /
                                    (f_end^-2 - f_start^-2))
            Δt_middle_upper = round(Int,
                                    i_Δt * (f_middle_upper^-2 - f_start^-2) /
                                    (f_end^-2 - f_start^-2))
            Δt_rest = i_Δt - Δt_middle_upper
            if !iszero(Δt_middle_lower)
                output[:, i_F, i_Δt + 1] = @views input[:, 2 * i_F - 1,
                                                        Δt_middle_lower + 1] +
                                                  input[:, 2 * i_F, Δt_rest + 1] .*
                                                  shift_row[:, Δt_middle_upper + 1]
            else
                output[:, i_F, i_Δt + 1] = @views input[:, 2 * i_F - 1,
                                                        Δt_middle_lower + 1] +
                                                  input[:, 2 * i_F, Δt_rest + 1]
            end
        end
    end
    return output
end

function fft_fdmt(data::AbstractMatrix, f_min::Real, f_max::Real, Δt_max=nothing;
                  dtype=Float32)
    ##### Initialization
    @assert f_min < f_max
    n_samp, n_chan = size(data)
    if isnothing(Δt_max)
        Δt_max = n_chan
    end
    δf = (f_max - f_min) / n_chan
    # Maximum number of time shifts in sample channels for the first iteration
    Δt = ceil(Int, Δt_max * (f_min^-2 - (f_min + δf)^-2) / (f_min^-2 - f_max^-2))
    # Initialize output
    # Axes are [samples,channels,delta_t] because column-major
    state = zeros(Complex{dtype}, n_samp, n_chan, Δt + 1)
    state[:, :, 1] .= fft(data, 1)
    shifter = fft(cumsum(Matrix{dtype}(I, n_samp, Δt + 1); dims=2), 1)
    for i_Δt in 2:(Δt + 1)
        state[:, :, i_Δt] = state[:, :, 1] .* shifter[:, i_Δt]
    end
    ##### Iterations
    for i_t in 1:Int(log2(n_chan))
        state = fft_fdmt_iterate(state, Δt_max, n_chan, f_min, f_max, i_t)
    end
    state = real.(ifft!(state, 1))
    _, _, n_dm = size(state)
    return reshape(state, n_samp, n_dm)
end

export fdmt