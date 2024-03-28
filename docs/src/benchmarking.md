```@meta
CurrentModule = Epidemics
```

This section shows some benchmarking.

```@example
using BenchmarkTools
using Epidemics

# an epidemic of 500 days
time_end = 500.0

# benchmark the default model with 3 age groups, intervention, and vaccination
@benchmark epidemic_default(β=[1.3/7], σ=[0.5], γ=[1/7], population = Population(), time_end=time_end, increment=1.0)
```
