
include("helpers.jl")
include("intervention.jl")
include("vaccination.jl")
include("prepare_data.jl")

using OrdinaryDiffEq

"""
    epidemic_vacamole_ode!(du, u, parameters, t)

An epidemic model adapted from the Vacamole model developed by RIVM, the Dutch
public health agency, to model the effect of leaky, two-dose vaccination on
the Covid-19 pandemic (Ainslie et al. 2022).

This adaptation implements leaky two-dose vaccination, as well as allowing for
"hospitalised" and "dead" compartments.
Individuals who have been doubly vaccinated can be modelled as having lower
susceptibility to infection, and lower hospitalisation and mortality rates.

The model also allows for multiple demographic groups and interventions.
This function is intended to be called internally from
[`epidemic_vacamole`](@ref).

The function expects the `parameters` argument to be a four element vector or
tuple with the following elements (which do not have to be named):

- a `Population` object with a prepared contact matrix, see [`Population`]
(@ref);
- ``\\beta``, the baseline transmission rate;
- ``\\beta_v``, the transmission rate for doubly vaccinated individuals;
- ``\\alpha``, the rate of conversion from exposed to infectious;
- ``\\eta``, the baseline hospitalisation rate;
- ``\\eta_v``, the hospitalisation rate for doubly vaccinated individuals;
- ``\\omega``, the baseline mortality rate;
- ``\\omega_v``, the mortality rate for doubly vaccinated individuals;
- ``\\gamma``, the rate of recovery;
- a matrix specifying the contacts between demographic groups;
- an `Npi` object specifying the intervention applied to each age group, see [`Npi`](@ref);
- a `Vaccination` object, see [`Vaccination`](@ref);

    
"""
function epidemic_vacamole_ode!(du, u, parameters, t)
    # assumes that each element of the vector is one of the required params
    population, β, βv, σ, η, ηv, ω, ωv, γ, intervention, vaccination = parameters

    # modify contact matrix if the intervention is active
    contact_matrix = population.contact_matrix .*
                     (1 .- cumulative_npi(t = t, npi = intervention))

    # get current ν
    ν_now = current_nu(time = t, vaccination = vaccination)

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    # column structure
    # 1| 2| 3|4| 5|6| 7|8| 9|10|11
    # S|V1|V2|E|EV|I|IV|H|HV|D|R
    S = @view u[:, 1]
    V1 = @view u[:, 2]
    V2 = @view u[:, 3]
    E = @view u[:, 4]
    Ev = @view u[:, 5]
    I = contact_matrix * ((@view u[:, 6]) + (@view u[:, 7])) # matrix mult for cm * (I + Iv)
    I_ = @view u[:, 6] # unmultiplied I for operations involving only I
    Iv_ = @view u[:, 7]
    H = @view u[:, 8]
    Hv = @view u[:, 9]

    # views to the change matrix, dU
    dS = @view du[:, 1]
    dV1 = @view du[:, 2]
    dV2 = @view du[:, 3]
    dE = @view du[:, 4]
    dEv = @view du[:, 5]
    dI = @view du[:, 6]
    dIv = @view du[:, 7]
    dH = @view du[:, 8]
    dHv = @view du[:, 9]
    dD = @view du[:, 10]
    dR = @view du[:, 11]

    # calculate change in compartment size and update the change matrix dU
    # note the use of @. for broadcasting, equivalent to .=
    @. dS = -(β * S * I) - (ν_now * S) # contact matrix cannot be multiplied here
    @. dV1 = (ν_now * S) - (ν_now * V1) - (β * V1 * I)
    @. dV2 = (ν_now * V1) - (βv * V2 * I)
    @. dE = (β * S * I) + (β * V1 * I) - (σ * E)
    @. dEv = (βv * V2 * I) - (σ * Ev)
    @. dI = (σ * E) - (η * I_) - (γ * I_) - (ω * I_) # note use of I_
    @. dIv = (σ * Ev) - (ηv * Iv_) - (γ * Iv_) - (ωv * Iv_) # note use of I_
    @. dH = (η * I_) - (γ * H) - (ω * H)
    @. dHv = (ηv * Iv_) - (γ * Hv) - (ωv * Hv)
    @. dR = γ * (I_ + Iv_ + H + Hv)
    @. dD = ω * (I_ + H) + ωv * (Iv_ + Hv)
end

"""
    epidemic_vacamole(β, σ, γ, η, ω, βv, ηv, ωv,
        population, intervention, vaccination,
        time_end, increment
    )

Model the progression of an epidemic, with age- or demographic-group specific
contact patterns and proportions, non-pharmaceutical interventions with group-
specific effects, and group-specific vaccination regimes.

The Vacamole model was developed by RIVM, and is particularly useful for
modelling leaky, two-dose vaccination.

## Arguments

- β: the transmissibility of the disease;
- σ: the rate of transition from the exposed to the infectious compartment;
- γ: the recovery rate;
- η: the hospitalisation rate;
- ω: the mortality rate;
- βv: the transmissibility for fully (doubly) vaccinated individuals;
- ηv: the hospitalisation rate for fully vaccinated individuals;
- ωv: the mortality rate for fully vaccinated individuals;

- `population`: A `Population` with population characteristics, which must have
11 compartments in the initial conditions;
- `intervention`: An `Npi` for interventions on social contacts;
- `vaccination`: A `Vaccination` for the vaccination regime applied;
- `time_end`: The time in days at which to end the simulation, defaults to 200;
- `increment`: The increment in simulation time, defaults to 1.0.

    
"""
function epidemic_vacamole(;
        β::Number, σ::Number, γ::Number,
        η::Number, ω::Number,
        βv::Number = β / 2.0,
        ηv::Number = η / 10.0,
        ωv::Number = ω / 10.0,
        population::Population,
        intervention = nothing,
        vaccination = nothing,
        time_end::Number = 200.0,
        increment::Number = 1.0)

    # input checking
    @assert increment<time_end "`increment` must be less than `time_end`"
    @assert increment>0.0 "`increment` must be a positive number"
    @assert time_end>0.0 "`time_end` must be a positive number"

    @assert β > 0.0&&β < Inf "β must be a finite positive number"
    @assert σ > 0.0&&σ < Inf "σ must be a finite positive number"
    @assert γ > 0.0&&γ < Inf "γ must be a finite positive number"
    @assert η > 0.0&&η < Inf "η must be a finite positive number"
    @assert ω > 0.0&&ω < Inf "ω must be a finite positive number"
    @assert βv > 0.0&&βv < Inf "βv must be a finite positive number"
    @assert ηv > 0.0&&ηv < Inf "ηv must be a finite positive number"
    @assert ωv > 0.0&&ωv < Inf "ωv must be a finite positive number"

    # check that population has the right number of compartments
    compartments = 11
    @assert size(population.initial_conditions)[2]==compartments "`population` must have 11 compartments"

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
    population.contact_matrix = prepare_contact_matrix(population = population)

    # prepare parameters
    parameters = tuple(population,
        β, σ, γ, η, ω, βv, ηv, ωv,
        intervention,
        vaccination)

    # prepare the timespan
    timespan = (0.0, time_end)

    # define the ode problem
    ode_problem = ODEProblem(epidemic_vacamole_ode!, init, timespan, parameters)

    # get the solution
    ode_solution = solve(ode_problem, AutoTsit5(Rosenbrock23()),
        saveat = increment)

    return ode_solution
end

export epidemic_vacamole
