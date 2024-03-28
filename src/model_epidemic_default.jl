
include("helpers.jl")
include("intervention.jl")
include("vaccination.jl")

using OrdinaryDiffEq

"""
    epidemic_default_ode!(du, u, parameters, t)

A simple SEIRV epidemic model function that allows for multiple demographic
groups. This function is intended to be called internally from
[`epidemic_default`](@ref).

The function expects the `parameters` argument to be a four element vector or
tuple with the following elements (which do not have to be named):

- a prepared `contact_matrix`, see [`Population`](@ref);
- ``\\beta``, the transmission rate;
- ``\\sigma``, the rate of conversion from exposed to infectious;
- ``\\gamma``, the rate of recovery;
- a matrix specifying the contacts between demographic groups;
- an `Npi` object specifying the intervention applied to each age
group, see [`Npi`](@ref);
- a `Vaccination` object, see [`Vaccination`](@ref);

    
"""
function epidemic_default_ode!(du, u, parameters, t)
    # assumes that each element of the vector is one of the required params
    contact_matrix, β, α, γ, intervention, vaccination = parameters

    # modify contact matrix if the intervention is active
    cm_temp = contact_matrix .*
              (1 .- cumulative_npi(t = t, npi = intervention))

    # get current ν
    ν_now = current_nu(time = t, vaccination = vaccination)

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, 1]
    E = @view u[:, 2]
    I = cm_temp * @view u[:, 3] # matrix mult for cm * I
    I_ = @view u[:, 3] # unmultiplied I for operations involving only I
    R = @view u[:, 4]
    V = @view u[:, 5]

    # views to the change matrix, dU
    dS = @view du[:, 1]
    dE = @view du[:, 2]
    dI = @view du[:, 3]
    dR = @view du[:, 4]
    dV = @view du[:, 5]

    # calculate change in compartment size and update the change matrix dU
    # note the use of @. for broadcasting, equivalent to .=
    @. dS = (-β * S * I) + (-ν_now * S) # contact matrix cannot be multiplied here
    @. dE = β * S * I - α * E
    @. dI = α * E - (γ * I_) # note use of I_
    @. dR = γ * I_
    @. dV = ν_now * S
end

"""
    epidemic_default(β, σ, γ,
        population, intervention, vaccination,
        time_end, increment
    )

Model the progression of an epidemic, with age- or demographic-group specific
contact patterns and proportions, non-pharmaceutical interventions with group-
specific effects, and group-specific vaccination regimes.

## Arguments

- β: the transmission rate ``\\beta`` of the disease; may be a numeric `Vector`;
- σ: the rate of transition from the exposed to the infectious compartment
    ``\\sigma``; may be a numeric `Vector`;
- γ: the recovery rate ``\\gamma``; may be a numeric `Vector`.

- `population`: A `Population` with population characteristics, importantly
    including a contact matrix describing group-specific social contacts, and
    a dmeography vector describing the number of individuals in each demographic
    group;
- `intervention`: An `Npi` for interventions on social contacts;
- `vaccination`: A `Vaccination` for the vaccination regime applied;
- `time_end`: The time in days at which to end the simulation, defaults to 200;
- `increment`: The increment in simulation time, defaults to 1.0.

    
"""
function epidemic_default(;
        β::Vector{<:Number} = [1.3 / 7.0], σ::Vector{<:Number} = [1.0 / 2.0],
        γ::Vector{<:Number} = [1.0 / 7.0],
        population::Population,
        intervention = nothing,
        vaccination = nothing,
        time_end::Number = 200.0,
        increment::Number = 1.0)

    # input checking
    @assert increment<time_end "`increment` must be less than `time_end`"
    @assert increment>0.0 "`increment` must be a positive number"
    @assert time_end>0.0 "`time_end` must be a positive number"

    # @assert β > 0.0&&β < Inf "β must be a finite positive number"
    # @assert σ > 0.0&&σ < Inf "σ must be a finite positive number"
    # @assert γ > 0.0&&γ < Inf "γ must be a finite positive number"

    # check that population has the right number of compartments
    compartments = ["susceptible", "exposed", "infectious", "recovered", "vaccinated"]
    @assert size(population.initial_conditions)[2]==length(compartments) "`population` must have 5 compartments"

    if isnothing(intervention)
        intervention = no_intervention()
    else
        # check type
        @assert isa(intervention, Npi) "`intervention` must be of type Npi"
        # expect that the intervention is compatible with the population
        @assert (size(intervention.contact_reduction)[1] ==
                 length(population.demography_vector))|(size(intervention.contact_reduction)[1] ==
                                                        1) "`intervention` 'contact_reduction' member rows must match the number of demography groups in `population`"
    end

    if isnothing(vaccination)
        vaccination = no_vaccination()
    else
        # check type
        @assert isa(vaccination, Vaccination) "`vaccination` must be of type Vaccination"
        # check that the vaccination is compatible with the population
        @assert size(vaccination.ν)==length(population.demography_vector) "`Vaccination` 'ν' must match the number of demography groups in `Population`"
    end

    # prepare the initial conditions
    init = prepare_initial_conditions(population = population)

    # prepare population contact matrix
    contact_matrix = prepare_contact_matrix(population = population)

    # prepare parameters
    # NOTE: population, intervention, and vaccination may only be scalars; WIP
    parameters = make_combinations([contact_matrix],
        β, σ, γ,
        [intervention],
        [vaccination])

    # prepare the timespan
    timespan = (0.0, time_end)

    # define the ode problem
    ode_problems = [ODEProblem(epidemic_default_ode!, init, timespan, i)
                    for i in parameters]

    # get the solution
    ode_solutions = [solve(i, AutoTsit5(Rosenbrock23()), saveat = increment)
                     for i in ode_problems]

    return ode_solutions
end

export epidemic_default
