
using Random
using Distributions

"""
    epidemic_stochastic(; population_size, n_infectious, 
        n_recovered, β, σ, time_end, time_increment
    )

Run a stochastic, discrete-time, compartmental epidemic model with the compartments "susceptible", "infectious", and
    "recovered".

# Named arguments
- `population_size::Number`: The total population size. Defaults to `1000`.
- `n_infectious::Number`: The number of initially infected individuals. Defaults to `10`.
- `n_recovered::Number`: The number of initially recovered individuals. Defaults to `0`.
- `β::Number`: The transmission rate of the infection (denoted ``\\beta``). Defaults to `10.0`.
- `σ::Number`: The recovery rate of the infection (denoted ``\\sigma``, often also denoted ``\\gamma``).
    Defaults to `1.0`.
- `time_end::Number`: The time point at which to end the simulation. Defaults to `5.0`.
- `time_increment::Number`: The increment in model time. Defaults to `0.01`.

# Returns
A `DataFrame` with four columns, "time", "susceptible", "infectious", and "recovered", for the values
    of each compartment at each time point in the simulation. The number of rows should be
    roughly equal to `time_end / time_increment` (and the initial conditions).

# Examples
```julia
# with default arguments
epidemic_stochastic()

# with some user-sepcified named arguments
epidemic_stochastic(population_size=5000, n_infectious=119,
    β=9.9, σ=1.1
)
```
"""
function epidemic_stochastic(;
        population_size::Number = 1000,
        n_infectious::Number = 10,
        n_recovered::Number = 0,
        β::Number = 10.0,
        σ::Number = 1.0,
        time_end::Number = 5.0,
        time_increment::Number = 0.01)

    # check inputs
    @assert (population_size > 0)&&isfinite(population_size) "`population_size` must be a positive number"
    @assert (n_infectious >= 0)&&isfinite(n_infectious) "`n_infectious` must be a positive number or zero"
    @assert (n_recovered >= 0)&&isfinite(n_recovered) "`n_recovered` must be a positive number or zero"
    @assert (β >= 0)&&isfinite(β) "`β` must be a positive number or zero"
    @assert (σ >= 0)&&isfinite(σ) "`σ` must be a positive number or zero"
    @assert (time_end >= 0)&&isfinite(time_end) "`time_end` must be a positive number or zero"
    @assert (time_increment >= 0)&&isfinite(time_increment) "`time_increment` must be a positive number or zero"

    # prepare timesteps for data container
    timesteps = range(0.0, time_end, step = time_increment)
    n_compartments = 3 # fixed to SIR model

    n_susceptibles = population_size - n_infectious - n_recovered

    # prepare data storage, time added later
    # extra row for initial conditions
    data = zeros(Number, length(timesteps), n_compartments)
    data[1, :] = [n_susceptibles, n_infectious, n_recovered]

    # precalculate transmission rate and recovery rate
    β_ = β / population_size * time_increment
    p_recovery = 1.0 - exp(-σ * time_increment)

    # loop over time in increments specified
    for t in 2:length(timesteps)
        # recalculate p(infect) as n_infectious changes
        p_infect = 1.0 - exp(-β_ * n_infectious)
        new_infections = rand(Binomial(n_susceptibles, p_infect))
        new_recoveries = rand(Binomial(n_infectious, p_recovery))

        n_susceptibles = n_susceptibles - new_infections
        n_infectious = n_infectious + new_infections - new_recoveries
        n_recovered = n_recovered + new_recoveries

        # record data
        data[t, :] = [n_susceptibles, n_infectious, n_recovered]
    end

    # convert to dataframe, add time
    data = DataFrame(data, ["susceptible", "infectious", "recovered"])
    data.time = timesteps

    return data
end

"""
run_replicates(model_fn, replicates; args...)

Run any model function multiple times while passing keyword arguments. Especially
    useful for capturing uncertainty due to randomness in stochastic models.

# Arguments
- `model_fn::Function`: A model function. Defaults to `epidemic_stochastic()` but
        may be any model function from this or another package.
- `replicates::Number`: The number of replicates. Defaults to `100`.

# Named arguments
- `args...`: Any extra arguments; these are typically arguments intended to be
    passed to `model_fn`.

# Returns
A `DataFrame` in long format with four columns, "time", "compartment", "value",
    and "replicate", for the values of each compartment at each time point in
    each replicate of the simulation run by `model_fn`.

# Examples
```julia
# pass some arguments to epidemic_stochastic()
run_replicates(epidemic_stochastic, 100,
    population_size=1000, n_infectious=10,
    β=10, σ=1
)
```
"""
function run_replicates(model_fn::Function = epidemic_stochastic,
        replicates::Number = 100;
        args...)
    @assert isfinite(replicates)&&(replicates > 0) "`replicates` must be a finite positive number"

    data = [DataFrame() for i in 1:replicates]

    # good candidate for multi-threading
    for i in 1:replicates
        df = model_fn(; args...)
        # add replicate number
        df[:, :replicate] .= i
        data[i] = df
    end

    # bind data and pivot, rename variable
    data = vcat(data...)
    data = stack(data, Not([:time, :replicate]))
    rename!(data, :variable => :compartment)

    return data
end

export epidemic_stochastic, run_replicates
