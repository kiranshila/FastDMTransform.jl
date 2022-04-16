# function dedisp(source, output, plan)
#     for i in 1:n_samp
#         @views output[i, :] = sum(circshift(source, (-(i - 1), 0))[plan]; dims=1)
#     end
# end

#dedisp(pulse, output, plan)

# function dedisp_cpu(source, freqs, dms, δt)
#     n_samp, n_chan = size(source)
#     n_dm = length(dms)
#     f_max = maximum(freqs)
#     output = zeros(Float32, n_samp, n_dm)
#     @tturbo for i in 1:n_dm
#         for j in 1:n_chan
#             for k in 1:n_samp
#                 dm = dms[i]
#                 f = freqs[j]
#                 dt = Δt(f, f_max, dm, δt, n_samp)
#                 source_idx = circmod(dt + k - 1, n_samp)
#                 output[k, i] += source[source_idx, j]
#             end
#         end
#     end
#     return output
# end