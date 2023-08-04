using DataFrames

function prepare_data(; ode_solution_df::DataFrame, n_age_groups::Number=3,
    compartment_names::Vector{String}=["susceptible", "exposed", "infectious",
        "recovered", "vaccinated"])

    # input checking here

    # create dummy age group names
    demo_groups = string.("demo_grp_", 1:n_age_groups)

    # prepare vector of age-and-compartment names
    df_names = Vector{String}(undef, n_age_groups * length(compartment_names))
    # prepare iterator
    k = 1

    # generate names
    for i in compartment_names
        for j in demo_groups
            df_names[k] = j * "." * i
            k = k + 1
        end
    end

    # check for length of ouput

    # rename dataframe
    rename!(
        ode_solution_df,
        Dict(
            zip(
                filter(e -> e != "timestamp", names(ode_solution_df)),
                df_names
            )
        )
    )
    # convert to long format, requires reassignment
    ode_solution_df = stack(ode_solution_df, df_names)

    # split the variable column
    transform!(
        ode_solution_df,
        :variable => ByRow(x -> split(x, '.')) =>
            [:demo_group, :compartment]
    )
    select!(ode_solution_df, Not([:variable]))

    # return reshaped and renamed data
    return ode_solution_df
end
