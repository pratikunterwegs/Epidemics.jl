```@meta
CurrentModule = Epidemics
```

This section shows how to run the Daedalus model.

```@example basic_daedalus
using Epidemics
using Plots

# all arguments have appropriate defaults
data = epidemic_daedalus(time_end=600.0)

# plot exposed group
plot(data, vars=(0, 50:99))
```

```@example basic_daedalus
# plot exposed among vaccinated - should start at vax_time = 200.0
plot(data, vars=(0, 394:442))
```

## Benchmarking

```@example benchmarking
# benchmark for a typical daedalus run of 600 days
using Epidemics
using BenchmarkTools

@benchmark epidemic_daedalus(time_end = 600.0)
```
