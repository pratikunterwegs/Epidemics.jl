```@meta
CurrentModule = Epidemics
```

This section shows how to run the Daedalus model.

```@example basic_daedalus
using Epidemics
using Plots

# all arguments have appropriate defaults
data = epidemic_daedalus()

# plot exposed group
plot(data, vars=(0, 5:8))
```

## Benchmarking

```@example benchmarking
using Epidemics
using BenchmarkTools

@benchmark epidemic_daedalus()
```
