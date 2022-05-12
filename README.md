# FastDMTransform

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kiranshila.github.io/FastDMTransform.jl/dev)
[![Build Status](https://github.com/kiranshila/FastDMTransform.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kiranshila/FastDMTransform.jl/actions/workflows/CI.yml?query=branch%3Amain)

A pure-julia implementation of the "Fast DM Transform" from:

[An accurate and efficient algorithm for detection of radio bursts with an unknown dispersion measure, for single dish telescopes and interferometers](https://arxiv.org/abs/1411.5373) - Zackay 2014, et. al.

## Implementation Details

Here, we follow the Python implementation of provided by Vincent Morello, [pyfdmt](https://bitbucket.org/vmorello/pyfdmt/src/master/), as they implemented a very nice recursive version, a big improvement over the confusing nested loops from the source paper. In this code, we make the interface and abstractions more Julian, and heavily optimize for performance. This includes transforming the recursion into parallel sums

In this code, we make heavy use of all available AVX extensions using [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl), so performance may be reasonably dependent on CPU, as Intel CPUs have AVX512.

## Benchmarks
Here we will try to FastDMTransform a synthetic pulse from a DM of 0 to 2000 with 4096 frequency channels, and 3000 1ms time samples.
These benchmarks are run on a machine with 64 GB of RAM and an i9-9900X

### "Canonical" MATLAB
Tested here was the original source from the author.
I only have Octave as I don't want to use nonfree software, and Octave doesn't have `timeit`, so I run it a few times with `tic` and `toc` and averaged around 1.3s

### pyfdmt
```python
%timeit -n 1 -r 10 pyfdmt.transform(pulse.T,1500,1200,1e-3,0,2000)
446 ms ± 2.19 ms per loop (mean ± std. dev. of 10 runs, 1 loop each)
```

### FastDMTransform.jl
```julia
julia> @benchmark fdmt(pulse,1500,1200,1e-3,0,2000)
BenchmarkTools.Trial: 71 samples with 1 evaluation.
 Range (min … max):  50.494 ms … 150.447 ms  ┊ GC (min … max): 0.00% … 28.43%
 Time  (median):     55.299 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   71.116 ms ±  30.091 ms  ┊ GC (mean ± σ):  9.60% ± 10.24%

  ▁█▁                                                           
  ███▄▅▅█▁▄▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▁▃▁▄▄▃▃▁▃▁▃▁▃▃▁▁▃▁▁▁▃▃ ▁
  50.5 ms         Histogram: frequency by time          144 ms <

 Memory estimate: 479.17 MiB, allocs estimate: 281374.
 ```

### [Dedisp.jl](https://github.com/kiranshila/Dedisp.jl)

Here, I've been working on an "optimized" brute force approach

```julia
plan = plan_dedisp(range(;start=1500,stop=1200,length=4096),1500,range(;start=0,stop=2000,length=2076),1e-3)
julia> @benchmark Dedisp.dedisp(pulse,plan)
BenchmarkTools.Trial: 1 sample with 1 evaluation.
 Single result which took 7.142 s (0.00% GC) to evaluate,
 with a memory estimate of 23.76 MiB, over 8 allocations.
```

But, FDMT is still not as fast as the GPU-accelerated brute-force version, but reasonably close
```julia
output = CUDA.zeros(3000,2076)
plan = cu(plan)
pulse = cu(pulse)
julia> @benchmark CUDA.@sync dedisp!(output,pulse,plan)
BenchmarkTools.Trial: 85 samples with 1 evaluation.
 Range (min … max):  55.555 ms … 71.742 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     58.945 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   59.135 ms ±  1.516 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%

                                    ▃ █▆▆▄▄▄▃ ▁    ▃           
  ▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▁▁▁▁▁▆▄█▇███████▁█▇▇▇▆█▁▆▁▆▄▁▁▁▄ ▁
  55.6 ms         Histogram: frequency by time        60.5 ms <

 Memory estimate: 8.72 KiB, allocs estimate: 150.
```

 Our implementation here is the fastest (CPU version) among the lot - almost 10x faster than the reference python implementation, where we can cover every possible DM in a 3 second chunk in 50ms, that's 60x faster than realtime!

