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
@benchmark epidemic_default(r0 = 1.5, infectious_period = 7, preinfectious_period = 2, population = Population(), time_end=time_end, increment=1.0)
```

## Get started

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

# run the default model with 3 age groups, but no intervention or vaccination
data = epidemic_default(
    r0 = 1.5, infectious_period = 7,
    preinfectious_period = 2,
    population = pop,
    time_end = sim_time_end, increment=1.0
)

# convert to dataframe
data_output = DataFrame(data)

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
    Scale.x_continuous(minvalue=100, maxvalue=500),
    Theme(
        key_position=:top
    )
)
```

## Modelling an intervention on social contacts

```@setup simple_contacts_intervention
using Epidemics
using Gadfly
using DataFrames

# an epidemic of 300 days
sim_time_end = 500.0

# population of 10 million in three age groups
pop = Population(
    demography_vector = 10e6 .* [0.2, 0.5, 0.3],
    initial_conditions = [1 - 1e-6 0.0 1e-6 0.0 0.0; 1 - 1e-6 0.0 1e-6 0.0 0.0; 1 - 1e-6 0.0 1e-6 0.0 0.0],
    contact_matrix = ones(3, 3) * 5
)
```

```@example simple_contacts_intervention
# an intervention that reduces contacts by 50%, 20% and 60% for each age group
# respectively
intervention = Npi(
    time_begin = 200, time_end = 230, 
    contact_reduction = [0.3, 0.1, 0.1]
)

# run a model with no intervention or vaccination
data_baseline = epidemic_default(
    r0 = 1.5, infectious_period = 7,
    preinfectious_period = 2,
    population = pop,
    time_end = sim_time_end, increment = 1.0
)

# run the default model with 3 age groups, intervention, no vaccination
data = epidemic_default(
    r0 = 1.5, infectious_period = 7,
    preinfectious_period = 2,
    population = pop,
    intervention = intervention,
    time_end = sim_time_end, increment = 1.0
)

# convert to dataframe
data_baseline = DataFrame(data_baseline)
data_output = DataFrame(data)

# function to handle data with correct naming
data_baseline = prepare_data(ode_solution_df = data_baseline, n_age_groups = 3)
data_output = prepare_data(ode_solution_df = data_output, n_age_groups = 3)

# assign scenario and combine
data_baseline[!, :scenario] .= "baseline"
data_output[!, :scenario] .= "intervention"

data = vcat(data_baseline, data_output)

# filter data for infectious only
data_infectious = filter(:compartment => n -> n == "infectious", data)

plot(
    data_infectious, 
    x = "timestamp",
    y = "value", 
    color = "demo_group",
    linestyle = "scenario",
    Geom.line,
    Guide.xlabel("Time"),
    Guide.ylabel("Individuals infectious"),
    Guide.colorkey("Demographic group"),
    Scale.x_continuous(minvalue=100, maxvalue=500),
    Theme(
        key_position=:top
    )
)
```

The example here shows how implementing a non-pharmaceutical intervention can reduce the number of infectious individuals in the population.

## Index

```@index
```

## Documentation

```@autodocs
Modules = [Epidemics]
```
