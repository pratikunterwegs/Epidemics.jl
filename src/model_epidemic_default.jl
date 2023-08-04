
include("helpers.jl")
include("intervention.jl")
include("population.jl")
include("pathogen.jl")
include("vaccination.jl")
include("prepare_args_default.jl")
include("prepare_data.jl")

using OrdinaryDiffEq

"""
    epidemic_default_ode!(du, u, parameters, t)

A simple SEIRV epidemic model function that allows for multiple demographic
    groups. This function is intended to be called internally from
    [`epidemic_default`](@ref).
    The function expects the `parameters` argument to be a four element vector
    with the following elements:
    - a `Population` object with a prepared contact matrix, see [`Population`]
    (@ref);
    - ``\\beta``, the transmission rate;
    - ``\\alpha``, the rate of conversion from exposed to infectious;
    - ``\\gamma``, the rate of recovery;
    - a matrix specifying the contacts between demographic groups;
    - an `Npi` object specifying the intervention applied to each age
    group, see [`Npi`](@ref);
    - a `Vaccination` object, see [`Vaccination`](@ref);

    
"""
function epidemic_default_ode!(du, u, parameters, t)
    # assumes that each element of the vector is one of the required params
    population, β, α, γ, intervention, vaccination = parameters

    # modify contact matrix if the intervention is active
    if (t > intervention.time_begin) & (t < intervention.time_end)
        population.contact_matrix = population.contact_matrix .*
                                    (1.0 .- intervention.contact_reduction)
    end

    # get current ν
    ν_now = current_nu(time=t, vaccination=vaccination)

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, 1]
    E = @view u[:, 2]
    I = population.contact_matrix * @view u[:, 3] # matrix mult for cm * I
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
    epidemic_default(population, infection, intervention,
        time_end, increment
    )

Model the progression of an epidemic, with age- or demographic-group specific
contact patterns and proportions, non-pharmaceutical interventions with group-
specific effects, and group-specific vaccination regimes.
    
"""
function epidemic_default(;
    population=Population(),
    infection=Infection(),
    intervention=Npi(),
    vaccination=Vaccination(),
    time_end=200.0,
    increment=0.1)

    # TODO: input checking goes here

    # prepare the initial conditions
    init = prepare_initial_conditions(
        population=population
    )

    # prepare parameters
    parameters = prepare_args_default(
        population=population,
        infection=infection,
        intervention=intervention,
        vaccination=vaccination
    )

    # prepare the timespan
    timespan = (0.0, time_end)

    # define the ode problem
    ode_problem = OrdinaryDiffEq.ODEProblem(
        epidemic_default_ode!, init, timespan, parameters
    )

    # get the solution
    ode_solution = OrdinaryDiffEq.solve(
        ode_problem, OrdinaryDiffEq.AutoTsit5(OrdinaryDiffEq.Rosenbrock23()),
        saveat=increment
    )

    # convert to dataframe
    data_output = DataFrames.DataFrame(ode_solution)

    # WIP - function to handle data with correct naming
    data_output = prepare_data(ode_solution_df=data_output)

    return data_output

end

export epidemic_default
