
include("seir.jl")

using OrdinaryDiffEq
using DataFrames

function epidemic(;
    model_name="seir", # named arguments begin here
    init=[1.0 - 1e-4, 1e-4, 0.0, 0.0],
    parameters=[0.1, 0.5, 0.05, 0.01],
    time_end=200.0,
    increment=0.1)

    # input checking goes here

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

    # rename dataframe columns
    data_output = DataFrames.rename(
        data_output,
        Dict(
            "value1" => "p_susceptible", "value2" => "p_exposed",
            "value3" => "p_infected", "value4" => "p_recovered"
        )
    )

    return data_output

end

export epidemic
