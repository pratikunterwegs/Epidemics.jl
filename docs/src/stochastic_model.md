```@meta
CurrentModule = Epidemics
```

This section shows how to run the simple discrete-time, stochastic SIR compartmental model included in Epidemics.jl.

## Get started with a single run

Here we see how to use `epidemic_stochastic()` to run a single realisation of the stochastic model.

```@example
using Epidemics
using Gadfly
using DataFrames
using Random

# end simulation time at 5.0 time units
sim_time_end = 5.0

# run the model with slightly modified parameters from the defaults
Random.seed!(1) # set a seed

data = epidemic_stochastic(
    population_size = 1010,
    n_infectious = 20, n_recovered = 30,
    β = 9.9, σ = 1.01,
    time_end = sim_time_end,
    time_increment = 0.02
)

# pivot the data to long format and plot
data = stack(data, Not(:time))
rename!(data, :variable => :compartment)

plot(
    data, 
    x = "time",
    y = "value", 
    color = "compartment", 
    Geom.line,
    Guide.xlabel("Time"),
    Guide.ylabel("Individuals"),
    Guide.colorkey("Compartment"),
    Theme(
        key_position=:top
    )
)
```

## Run multiple replicates

Here we see how to use `run_replicates()` to run a multiple realisations of the stochastic model in `epidemic_stochastic()`.

`run_replicates()` takes the function name (such as `epidemic_stochastic`), the number of replicates, and any number of function arguments for the epidemic model function.

Here, we run 100 replicates of the stochastic model, passing the same arguments as in the previous example.

```@example
using Epidemics
using Gadfly
using DataFrames
using Random

# end simulation time at 5.0 time units
sim_time_end = 5.0

# run the model with slightly modified parameters from the defaults
Random.seed!(1) # set a seed

data = run_replicates(
    epidemic_stochastic, 100,
    population_size = 1010,
    n_infectious = 20, n_recovered = 30,
    β = 9.9, σ = 1.01,
    time_end = sim_time_end,
    time_increment = 0.01
)

plot(
    data, 
    x = "time",
    y = "value", 
    color = "compartment",
    group = "replicate",
    Geom.line,
    alpha=[0.1],
    Guide.xlabel("Time"),
    Guide.ylabel("Individuals"),
    Guide.colorkey("Compartment"),
    Coord.cartesian(ymin=0, ymax=1000),
    Theme(
        key_position=:top
    )
)
```
