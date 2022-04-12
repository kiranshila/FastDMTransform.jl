# FDMT

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kiranshila.github.io/FDMT.jl/dev)
[![Build Status](https://github.com/kiranshila/FDMT.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kiranshila/FDMT.jl/actions/workflows/CI.yml?query=branch%3Amain)

A pure-julia implementation of the "Fast DM Transform" from:

[An accurate and efficient algorithm for detection of radio bursts with an unknown dispersion measure, for single dish telescopes and interferometers](https://arxiv.org/abs/1411.5373) - Zackay 2014, et. al.

## Implementation Details

Here, we follow the Python implementation of provided by Vincent Morello, [pyfdmt](https://bitbucket.org/vmorello/pyfdmt/src/master/), as they implemented a very nice recursive version, a big improvement over the confusing nested loops from the source paper. In this code, we make the interface and abstractions more Julian, and heavily optimize for performance.

In this code, we make heavy use of all available AVX extensions using [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl), so performance may be reasonably dependent on CPU, as Intel CPUs have AVX512.

## Benchmarks
Here we will try to FDMT a synthetic pulse from a DM of 0 to 2000 with 4096 frequency channels, and 3000 1ms time samples.
These benchmarks are run on a machine with 64 GB of RAM and an i9-9900X

### "Canonical" MATLAB
Tested here was the original source from the author.
I only have Octave as I don't want to use nonfree software, and Octave doesn't have `timeit`, so I run it a few times with `tic` and `toc` and averaged around 1.3s

### pyfdmt
```python
%timeit -n 1 -r 10 pyfdmt.transform(pulse.T,1500,1200,1e-3,0,2000)
446 ms ± 2.19 ms per loop (mean ± std. dev. of 10 runs, 1 loop each)
```

### FDMT.jl
```julia
julia> @benchmark transform(pulse,1500,1200,1e-3,0,2000,1)
BenchmarkTools.Trial: 19 samples with 1 evaluation.
 Range (min … max):  246.667 ms … 425.906 ms  ┊ GC (min … max): 2.48% … 41.62%
 Time  (median):     247.260 ms               ┊ GC (median):    2.48%
 Time  (mean ± σ):   266.404 ms ±  48.233 ms  ┊ GC (mean ± σ):  8.28% ± 11.47%

  █                                                              
  █▁▅▆▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃ ▁
  247 ms           Histogram: frequency by time          426 ms <

 Memory estimate: 407.32 MiB, allocs estimate: 15828.
 ```