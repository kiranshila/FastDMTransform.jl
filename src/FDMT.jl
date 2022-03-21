module FDMT

using PaddedViews

const KDM = 4.148808e3 # MHz^2 s pc^-1 cm^3

function disp_shift(DM, f₁, f₂)
    return KDM * DM * (f₁^-2 - f₂^-2)
end

# This code is bad, don't actually use this yet.
# The paper from arxiv is riddled with errors and missing details
# the published paper less so, but still not enough to implement from
# just the pseudocode. This is more or less verbatim from the matlab impl
# and is not idiomatic.

function init(I, f_min, f_max, Δt_max, eltype)
    N_t, N_f = size(I)
    δf = (f_max - f_min) / N_f
    δt = ceil(Int,
              Δt_max * (1 / f_min^2 - 1 / (f_min + δf)^2) / (1 / f_min^2 - 1 / f_max^2))
    A = zeros(eltype, N_t, N_f, δt + 1)
    I_padded = PaddedView(0, I, (N_t + δt + 1, N_f))
    A[:, :, 1] .= I
    for i in 2:(δt + 1)
        A[:, :, i] = @views A[:, :, i - 1] + I_padded[i:(i + N_t - 1), :]
    end
    return A
end

function iter!(A, Δt_max, f_min, f_max, i)
    N_t,N_f,_ = size(A)

    # frequency difference between adjacent frequency sub-bands
    Δf = 2^(i) * (f_max - f_min) / N_f

    # frequency difference between adjacent frequency bins
    dF = (f_max - f_min) / N_f

    # the maximum deltaT needed to calculate at the i'th iteration
    Δt = ceil(Int,
              Δt_max * (1 / f_min^2 - 1 / (f_min + Δf)^2) / (1 / f_min^2 - 1 / f_max^2))

    F_jumps = N_f ÷ 2

    # allowing for a shift to occur between subbands
    correction = dF / 2

    for i_F in 1:F_jumps,
        f_start in (f_max - f_min) / F_jumps * (i_F - 1) + f_min

        f_end = (f_max - f_min) / F_jumps * (i_F) + f_min
        f_middle = (f_end - f_start) / 2 + f_start

        Δt_local = ceil(Int,
                        Δt_max * (1 / f_start^2 - 1 / (f_end)^2) /
                        (1 / f_min^2 - 1 / f_max^2))
        for i_dT in 0:Δt_local

            # determining all the needed shift constants for the iteration
            dT_middle = round(Int,
                              i_dT * (1 / (f_middle - correction)^2 - 1 / f_start^2) /
                              (1 / f_end^2 - 1 / f_start^2))
            dT_middle_index = dT_middle + 1
            dT_middle_larger = round(Int,
                                     i_dT *
                                     (1 / (f_middle + correction)^2 - 1 / f_start^2) /
                                     (1 / f_end^2 - 1 / f_start^2))

            dT_rest = i_dT - dT_middle_larger
            dT_rest_index = dT_rest + 1

            i_T_min = 1
            i_T_max = dT_middle_larger + 1

            # Alternative addition rule!
            A[i_T_min:i_T_max, i_F, i_dT + 1] = A[i_T_min:i_T_max,
                                                        2 * i_F - 1,
                                                        dT_middle_index]

            i_T_min = dT_middle_larger + 2
            i_T_max = N_t

            # Addition rule!
            A[i_T_min:i_T_max, i_F, i_dT + 1] = A[i_T_min:i_T_max,
                                                         2 * i_F - 1,
                                                         dT_middle_index] +
                                                       A[(i_T_min - dT_middle_larger):(i_T_max - dT_middle_larger),
                                                         2 * i_F, dT_rest_index]
        end
    end
end

function fdmt(I,f_min,f_max,Δt_max;eltype=Float64)
    _, N_f = size(I)
    A = init(I,f_min,f_max,Δt_max,eltype)
    for i ∈ 1:log2(N_f)
        iter!(A,Δt_max,f_min,f_max,i)
    end
    return A
end

end
