
include("helpers.jl")
include("intervention.jl")
include("population.jl")
include("pathogen.jl")
include("seir.jl")
include("prepare_data.jl")

using OrdinaryDiffEq

"""
    epidemic(model_name, init, contact_matrix, demography_vector, r0,
        preinfectious_period, infectious_period, interv, time_end, increment
    )

Model the progression of an epidemic, with age- or demographic-group specific
contact patterns and proportions, epidemiological parameters, and interventions.
    
"""
function epidemic(;
    model_name="seir", # named arguments begin here
    population=Population(),
    pathogen=Pathogen(),
    intervention=Npi(),
    time_end=200.0,
    increment=0.1)

    # input checking goes here

    # prepare the contact matrix
    contact_matrix = prepare_contact_matrix(
        population=population
    )

    # prepare the initial conditions
    init = prepare_initial_conditions(
        population=population
    )

    # prepare parameters
    parameters = [
        r0_to_beta(r0=pathogen.r0,
            infectious_period=pathogen.infectious_period
        ),
        preinfectious_period_to_beta2(
            preinfectious_period=pathogen.preinfectious_period
        ),
        infectious_period_to_gamma(
            infectious_period=pathogen.infectious_period
        ),
        contact_matrix,
        intervention
    ]

    # prepare the timespan
    timespan = (0.0, time_end)

    # pick the correct model function - can be done from database later
    if model_name == "seir"
        fn_epidemic = seir!
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
