# Epidemics.jl

[![License:MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pratikunterwegs.github.io/Epidemics.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pratikunterwegs.github.io/Epidemics.jl/dev/)
[![Build Status](https://github.com/pratikunterwegs/Epidemics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pratikunterwegs/Epidemics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/pratikunterwegs/Epidemics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/pratikunterwegs/Epidemics.jl)

_Epidemics.jl_ is a Julia package that aims to mirror the [R package _epidemics_](https://github.com/epiverse-trace/epidemics), and to provide a robust way to model epidemic and disease outbreak scenarios. _Epidemics.jl_ is a work in progress since it lags the development of _epidemics_.

_Epidemics.jl_ is currently faster than _epidemics_.
Automated benchmarking is in the pipeline.

```julia
BenchmarkTools.Trial: 706 samples with 1 evaluation.
 Range (min … max):  6.431 ms … 16.426 ms  ┊ GC (min … max): 0.00% … 58.82%
 Time  (median):     6.626 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   7.085 ms ±  1.976 ms  ┊ GC (mean ± σ):  6.24% ± 12.32%

  █▅                                                          
  ██▇▃▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▃▂ ▂
  6.43 ms        Histogram: frequency by time        15.8 ms <

 Memory estimate: 8.32 MiB, allocs estimate: 119924.
```
