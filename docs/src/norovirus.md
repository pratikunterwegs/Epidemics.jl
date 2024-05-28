```@meta
CurrentModule = Epidemics
```

This section shows how to run the norovirus model with two levels of vaccination-derived protection against symptomatic disease.

```@example basic_norovirus
using Epidemics
using Plots

# all arguments have appropriate defaults
data = epidemic_norovirus()

# plot exposed in unvaccinated group
plot(data, vars=(0, 5:8))
```

```@example basic_norovirus
# plot exposed in single-vaccinated group
plot(data, vars=(0, 25:28))
```

```@example basic_norovirus
# plot exposed in double-vaccinated group
plot(data, vars=(0, 45:48))
```

## Benchmarking

```@example benchmarking
using Epidemics
using BenchmarkTools

@benchmark epidemic_norovirus()
```
