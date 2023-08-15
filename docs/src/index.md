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

## Default example

```@repl
using Epidemics
using Gadfly

# an epidemic of 500 days
time_end = 100.0

# benchmark the default model with 3 age groups, intervention, and vaccination
data = epidemic_default(time_end=time_end, increment=1.0)

# filter data for infectious only
data_infectious = filter(:compartment => n -> n == "infectious", data)

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

## Example with an intervention on contacts

```@repl
using Epidemics
using Gadfly

# an epidemic of 500 days
time_end = 200.0

# define an intervention
intervention = Npi

# benchmark the default model with 3 age groups, intervention, and vaccination
data = epidemic_default(time_end=time_end, increment=1.0)

# filter data for infectious only
data_infectious = filter(:compartment => n -> n == "infectious", data)

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
