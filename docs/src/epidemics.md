```@meta
CurrentModule = Epidemics
```

### Installing Epidemics.jl

_Epidemics.jl_ can be installed from GitHub using the Julia package manager _Pkg.jl_.

```julia
using Pkg
Pkg.add(url="git@github.com:pratikunterwegs/Epidemics.jl.git")
```

## Running the default model

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

# make model parameters using helpers
r0 = 1.5
infectious_period = 7
preinfectious_period = 2

β = r0_to_beta(r0 = r0, infectious_period = infectious_period)
σ = preinfectious_period_to_alpha(preinfectious_period = preinfectious_period)
γ = infectious_period_to_gamma(infectious_period = infectious_period)

# run the default model with 3 age groups, but no intervention or vaccination
data = epidemic_default(
    β=[β], σ=[σ], γ=[γ],
    population = pop,
    time_end = sim_time_end, increment=1.0
)

# convert to dataframe
# NOTE that due to vectorisation, the output is a vector of DataFrames
data_output = DataFrame(data[1])

# WIP - function to handle data with correct naming
data_output = prepare_data(data_output, n_age_groups = 3)

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