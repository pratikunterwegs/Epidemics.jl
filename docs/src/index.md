```@meta
CurrentModule = Epidemics
```

# Epidemics.jl

This is some minimal documentation for [Epidemics.jl](https://github.com/pratikunterwegs/Epidemics.jl).
**Note** that this is a personal project, and comes with no current or future support.
This documentation section is intended as a learning experience (for me) in writing Julia package documentation.

[![License:MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
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

## Index

```@index
```

## Documentation

```@autodocs
Modules = [Epidemics]
```
