```@meta
CurrentModule = Epidemics
```

This section shows how to set up and use the `Npi` struct to model interventions on social contacts.

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

# make model parameters using helpers
r0 = 1.5
infectious_period = 7
preinfectious_period = 2

β = r0_to_beta(r0 = r0, infectious_period = infectious_period)
σ = preinfectious_period_to_alpha(preinfectious_period = preinfectious_period)
γ = infectious_period_to_gamma(infectious_period = infectious_period)

# run a model with no intervention or vaccination
data_baseline = epidemic_default(
    β=[β], σ=[σ], γ=[γ],
    population = pop,
    time_end = sim_time_end, increment = 1.0
)

# run the default model with 3 age groups, intervention, no vaccination
data = epidemic_default(
    β=[β], σ=[σ], γ=[γ],
    population = pop,
    intervention = intervention,
    time_end = sim_time_end, increment = 1.0
)

# convert to dataframe
data_baseline = DataFrame(data_baseline[1])
data_output = DataFrame(data[1])

# function to handle data with correct naming
data_baseline = prepare_data(data_baseline, n_age_groups = 3)
data_output = prepare_data(data_output, n_age_groups = 3)

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
