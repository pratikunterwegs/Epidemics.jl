```@meta
CurrentModule = Epidemics
```

# Epidemics

Documentation for [Epidemics](https://github.com/pratikunterwegs/Epidemics.jl).

## Benchmarking

This section shows some benchmarking.

```@repl
using BenchmarkTools
using Epidemics

# an epidemic of 500 days
time_end = 500.0

# benchmark the default model with 3 age groups, intervention, and vaccination
@benchmark epidemic_default(time_end=time_end, increment=1.0)
```

```@index
```

```@autodocs
Modules = [Epidemics]
```
