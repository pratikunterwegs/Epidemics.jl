```@meta
CurrentModule = Epidemics
```

# Epidemics

Documentation for [Epidemics](https://github.com/pratikunterwegs/Epidemics.jl).

## Benchmarking

This section shows some benchmarking.

```@example
using BenchmarkTools
using Epidemics

# an epidemic of 500 days
time_end = 500.0

# benchmark the default model with 3 age groups, intervention, and vaccination
@benchmark epidemic_default(time_end=time_end, increment=1.0)
```

## Get started with an intervention on contacts

```@example
using Epidemics
using Gadfly
using DataFrames

# an epidemic of 500 days
sim_time_end = 500.0

# population of 10 million in three age groups
pop = Population(
    demography_vector = 10e6 .* [0.2, 0.5, 0.3],
    initial_conditions = [1 - 1e-6 0.0 1e-6 0.0 0.0; 1 - 1e-6 0.0 1e-6 0.0 0.0; 1 - 1e-6 0.0 1e-6 0.0 0.0],
    contact_matrix = ones(3, 3) * 5
)

# an infection
pandemic = Infection(
    extra_arguments = (
        preinfectious_period = 2,
    )
)

# an intervention that reduces contacts by 20%
intervention = Npi(
    time_begin = sim_time_end / 4, time_end = sim_time_end / 2, 
    contact_reduction = [0.1, 0.1, 0.1]
)

# a dummy vaccination
no_vax = Vaccination(
    time_begin = [0], time_end = [0], ν = [0.0]
)

# run the default model with 3 age groups, intervention, and vaccination
data = epidemic_default(
    population = pop,
    infection = pandemic,
    intervention = intervention,
    vaccination = no_vax,
    time_end = sim_time_end, increment=1.0
)

# convert to dataframe
data_output = DataFrames.DataFrame(data)

# WIP - function to handle data with correct naming
data_output = prepare_data(ode_solution_df = data_output, n_age_groups = 3)

# filter data for infectious only
data_infectious = filter(:compartment => n -> n == "infectious", data_output)

plot(
    data_infectious, 
    x = "timestamp",
    y = "value", 
    color = "demo_group", 
    Geom.line,
    Guide.xlabel("Time"),
    Guide.ylabel("Individuals infectious"),
    Guide.colorkey("Demographic group"),
    Theme(
        key_position=:top
    )
)
```

## Index

```@index
```

## Documentation

```@autodocs
Modules = [Epidemics]
```
