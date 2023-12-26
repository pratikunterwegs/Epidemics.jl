
using Random
using Distributions

"""

"""
function epidemic_stochastic(;
        population_size = 1000,
        n_infectious = 10,
        n_recovered = 0,
        β = 10,
        σ = 1,
        time_end = 5.0,
        time_increment = 0.01)

    # check inputs
    @assert isa(population_size, Number)&&(population_size > 0) && isfinite(population_size) "`population_size` must be a positive number"
    @assert isa(n_infectious, Number)&&(n_infectious >= 0) && isfinite(n_infectious) "`n_infectious` must be a positive number or zero"
    @assert isa(n_recovered, Number)&&(n_recovered >= 0) && isfinite(n_recovered) "`n_recovered` must be a positive number or zero"
    @assert isa(β, Number)&&(β >= 0) && isfinite(β) "`β` must be a positive number or zero"
    @assert isa(σ, Number)&&(σ >= 0) && isfinite(σ) "`σ` must be a positive number or zero"
    @assert isa(time_end, Number)&&(time_end >= 0) && isfinite(time_end) "`time_end` must be a positive number or zero"
    @assert isa(time_increment, Number)&&(time_increment >= 0) && isfinite(time_increment) "`time_increment` must be a positive number or zero"

    # prepare timesteps for data container
    timesteps = range(0.0, time_end, step = time_increment)
    n_compartments = 3 # fixed to SIR model

    n_susceptibles = population_size - n_infectious - n_recovered

    # prepare data storage, time added later
    # extra row for initial conditions
    data = zeros(Number, length(timesteps), n_compartments)
    data[1, :] = [n_susceptibles, n_infectious, n_recovered]

    # loop over time in increments specified
    for t in 2:length(timesteps)
        # recalculate p(infect) as n_infectious changes
        p_infect = 1 - exp(-β * n_infectious / population_size * time_increment)
        new_infections = rand(Binomial(n_susceptibles, p_infect))
        new_recoveries = rand(Binomial(n_infectious, σ * time_increment))

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

function run_replicates(replicates = 100, fn = epidemic_stochastic; args...)
    @assert isa(replicates, Number)&&isfinite(replicates) && (replicates > 0) "`replicates` must be a finite positive number"
    @assert isa(fn, Function) "`fn` must be a function"

    data = [DataFrame() for i in 1:replicates]

    for i in 1:replicates
        df = fn(; args...)
        # add replicate number
        df.replicate .= i
        data[i] = df
    end

    # bind data
    data = vcat(data...)

    return data
end

export epidemic_stochastic, run_replicates
