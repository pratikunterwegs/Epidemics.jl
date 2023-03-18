
include("seir.jl")

using OrdinaryDiffEq
using DataFrames
using LinearAlgebra

function prepare_contact_matrix(;contact_matrix, demography_vector)
    # input checking here
    # WIP check that square matrix has as many rows as demographic groups

    # scale the contact matrix by its largest real eigenvalue
    # eigvals is assumed to return only the real eigenvalues
    # may need correction
    contact_matrix = contact_matrix / maximum(LinearAlgebra.eigvals(contact_matrix))
    
    # scale by the demography, divide each row by corresponding demography
    contact_matrix = contact_matrix ./ demography_vector
    
    return contact_matrix
end

function prepare_initial_conditions(;initial_conditions, demography_vector)
    # input checking here; check that initial initial_conditions
    # is a matrix and has as many rows as demography_vector

    # rowwise multiplication of initial conditions by demography_vector
    init = initial_conditions .* demography_vector

    # return linearised values
    return vec(init)
end

function epidemic(;
    model_name="seir", # named arguments begin here
    init=[1.0 - 1e-4, 1e-4, 0.0, 0.0],
    contact_matrix, demography_vector,
    parameters=[0.1, 0.5, 0.05, 0.01],
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
