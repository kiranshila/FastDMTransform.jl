using SIGPROC,Plots, FDMT, BenchmarkTools, CUDA

δt = 10e-6
n_samp = 100000
n_chan = 4096
f_min = 1280
f_max = 1530

t_total = n_samp * δt
dm_max = t_total / (FDMT.KDM * (f_min^-2 - f_max^-2))
dm_min = 2
n_dm = 2048

fb = SIGPROC.fake_pulse(300, f_max, f_min; samples=n_samp, channels=n_chan,w=5,t_step=δt,dtype=UInt8)
pulse = cu(fb.data.data)
freqs = cu(collect(fb.data.dims[2]))
dms = cu(collect(range(dm_min, dm_max; length=n_dm)))

# fb = Filterbank("/home/kiran/Downloads/candidate_ovro_20200428.fil")
# pulse = cu(fb.data.data)
# n_samp, n_chan = size(pulse)
# freqs = cu(collect(fb.data.dims[2]))
# dms = cu(collect(range(dm_min, dm_max; length=n_dm)))
# δt = step(fb.data.dims[1])
# f_min, f_max = extrema(freqs)

##### Non chunked
plan = plan_dedisp(freqs,f_max,dms,δt)
out = dedisp(pulse,plan)

##### Chunked
#n_chunk = 128
#plan = plan_chunked_dedisp(freqs,f_max,dms,δt,n_chunk)
#tmp = CUDA.zeros(n_samp, n_dm, n_chunk)
#out = dedisp_in_chunks!(tmp,pulse, plan)