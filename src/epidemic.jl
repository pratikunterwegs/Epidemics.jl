
include("helpers.jl")
include("classes.jl")
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
    init=[1.0-1e-6 1.0-1e-6 1.0-1e-6
        1e-6/2.0 1e-6/2.0 1e-6/2.0
        1e-6/2.0 1e-6/2.0 1e-6/2.0
        0.0 0.0 0.0
    ]',
    contact_matrix=ones(3, 3) * 5.0,
    demography_vector=67e6 * [0.23, 0.4, 0.37],
    r0=[1.5, 1.5, 1.5],
    preinfectious_period=[3.0, 3.0, 3.0],
    infectious_period=[7.0, 7.0, 7.0],
    interv=npi(70.0, 100.0, [0.3, 0.2, 0.8]),
    time_end=200.0,
    increment=0.1)

    # input checking goes here

    # prepare the contact matrix
    contact_matrix = prepare_contact_matrix(
        contact_matrix=contact_matrix,
        demography_vector=demography_vector
    )

    # prepare the initial conditions
    init = prepare_initial_conditions(
        initial_conditions=init, demography_vector=demography_vector
    )

    # prepare parameters
    parameters = [
        r0_to_beta(r0=r0, infectious_period=infectious_period),
        preinfectious_period_to_beta2(preinfectious_period=preinfectious_period),
        infectious_period_to_gamma(infectious_period=infectious_period),
        contact_matrix,
        interv
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
