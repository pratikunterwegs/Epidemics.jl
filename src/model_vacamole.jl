
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
- an `Npi` object specifying the intervention applied to each age
group, see [`Npi`](@ref);
- a `Vaccination` object, see [`Vaccination`](@ref);

    
"""
function epidemic_vacamole_ode!(du, u, parameters, t)
    # assumes that each element of the vector is one of the required params
    population, β, βv, α, η, ηv, ω, ωv, γ, intervention, vaccination = parameters

    # modify contact matrix if the intervention is active
    contact_matrix = population.contact_matrix .*
                     (1 .- cumulative_npi(t = t, npi = intervention))

    # get current ν
    ν_now = current_nu(time = t, vaccination = vaccination)

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, 1]
    E = @view u[:, 2]
    I = contact_matrix * @view u[:, 3] # matrix mult for cm * I
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
    epidemic_vacamole(r0, infectious_period, preinfectious_period,
        population, intervention, vaccination,
        time_end, increment
    )

Model the progression of an epidemic, with age- or demographic-group specific
contact patterns and proportions, non-pharmaceutical interventions with group-
specific effects, and group-specific vaccination regimes.
    
"""
function epidemic_vacamole(;
    r0::Number, 
    infectious_period::Number,
    preinfectious_period::Number,
    eta::Number, omega::Number,
    susc_reduction_vax::Number,
    hosp_reduction_vax::Number,
    mort_reduction_vax::Number
    population::Population,
    intervention = nothing,
    vaccination = nothing,
    time_end::Number = 200.0,
    increment::Number = 0.1)

    # input checking
    @assert increment<time_end "`increment` must be less than `time_end`"
    @assert increment>0.0 "`increment` must be a positive number"
    @assert time_end>0.0 "`time_end` must be a positive number"

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
        r0_to_beta(r0 = r0, infectious_period = infectious_period),
        preinfectious_period_to_alpha(preinfectious_period = preinfectious_period),
        infectious_period_to_gamma(infectious_period = infectious_period),
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
