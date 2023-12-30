# Epidemics.jl

[![License:MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pratikunterwegs.github.io/Epidemics.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pratikunterwegs.github.io/Epidemics.jl/dev/)
[![Build Status](https://github.com/pratikunterwegs/Epidemics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pratikunterwegs/Epidemics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/pratikunterwegs/Epidemics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/pratikunterwegs/Epidemics.jl)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

_Epidemics.jl_ is a Julia package that aims to mirror the [R package {epidemics}](https://github.com/epiverse-trace/epidemics), and to provide a robust way to model epidemic and disease outbreak scenarios.
_Epidemics.jl_ is a work in progress since it lags the development of _epidemics_.

**Note** _Epidemics.jl_ is a personal project where I aim to learn more about the Julia language; it comes with no guarantees of current or future support or maintenance.

_Epidemics.jl_ currently has basic implementations of three models, roughly tracking the R package {epidemics}.

1. `epidemic_default()`: the default model, which is an SEIRV compartmental ODE model allowing for a vaccination regime that confers full immunity with a single dose, as well as (optionally) multiple overlapping interventions to reduce social contacts;

2. `epidemic_vacamole()`: the Vacamole model developed by RIVM for the Covid-19 pandemic, which is a work in progress, but will eventually allow for two-dose leaky vaccination, as well as multiple overlapping interventions;

3. `epidemic_stochastic()`: A simple stochastic compartmental SIR model.

_Epidemics.jl_ is likely to include the Ebola model, as well as features such as time-dependence and rate interventions, from {epidemics} at some point.

## Get started

### Installation

_Epidemics.jl_ can be installed from GitHub using the Julia package manager _Pkg.jl_.

```julia
using Pkg
Pkg.add(url="git@github.com:pratikunterwegs/Epidemics.jl.git")
```

### Running an ODE epidemic model

You can run a simple age-structured epidemic model using the function `epidemic_default()` with its default arguments.

```julia
using Epidemics

# an epidemic of 500 days
time_end = 500.0

# the default model with 3 age groups in the default population
epidemic_default(β=1.3/7, σ=0.5, γ=1/7,
    population = Population(),
    time_end=time_end, increment=1.0)
```

### Running a stochastic epidemic model

You can run a simple stochastic epidemic model using the function `epidemic_stochastic()` with its default arguments.

```julia
using Epidemics

# an epidemic of 500 days
time_end = 500.0

# the default model with 3 age groups in the default population
epidemic_stochastic(population_size = 1010,
    n_infectious = 20, n_recovered = 30,
    β = 9.9, σ = 1.01,
    time_end = sim_time_end,
    time_increment = 0.02)
```

## Benchmarking

_Epidemics.jl_ is currently faster than _epidemics_.
Automated benchmarking is in the pipeline can be found in the [development documentation](https://pratikunterwegs.github.io/Epidemics.jl/dev/).
A static example for the default ODE model is shown below.

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
