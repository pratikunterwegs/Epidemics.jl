
include("helpers.jl")
include("seir.jl")

using OrdinaryDiffEq
using DataFrames
using LinearAlgebra

function epidemic(;
    model_name="seir", # named arguments begin here
    init=[0.9 0.9 0.9; 0.09 0.09 0.09; 0.01 0.01 0.01; 0.0 0.0 0.0]',
    contact_matrix, demography_vector,
    parameters=[[0.14, 0.12, 0.2], [0.5, 0.5, 0.5], [0.05, 0.05, 0.05]],
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

    # prepare the timespan
    timespan = (0.0, time_end)

    # pick the correct model function - can be done from database later
    if model_name == "seir"
        fn_epidemic = seir
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

    return data_output

end

export epidemic
