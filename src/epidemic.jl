
include("helpers.jl")
include("intervention.jl")
include("population.jl")
include("pathogen.jl")
include("vaccination.jl")
include("prepare_args_default.jl")
include("model_epidemic_default.jl")
include("prepare_data.jl")

using OrdinaryDiffEq

"""
    epidemic(model_name, population, infection, intervention,
        time_end, increment
    )

Model the progression of an epidemic, with age- or demographic-group specific
contact patterns and proportions, non-pharmaceutical interventions with group-
specific effects, and group-specific vaccination regimes.
    
"""
function epidemic(;
    model_name="seir", # named arguments begin here
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

    # pick the correct model function - can be done from database later
    if model_name == "seir"
        fn_epidemic = epidemic_default!
    end

    # define the ode problem
    ode_problem = OrdinaryDiffEq.ODEProblem(
        fn_epidemic, init, timespan, parameters
    )

    # get the solution
    ode_solution = OrdinaryDiffEq.solve(
        ode_problem, OrdinaryDiffEq.AutoTsit5(OrdinaryDiffEq.Rosenbrock23()),
        saveat=increment
    )

    # convert to dataframe
    data_output = DataFrames.DataFrame(ode_solution)

    # WIP - function to handle data with correct naming

    return prepare_data(ode_solution_df=data_output)

end

export epidemic
