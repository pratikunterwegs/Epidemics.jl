using DataFrames

function prepare_data!(; ode_solution_df, n_age_groups=3,
    compartment_names=["susceptible", "exposed", "infectious", "recovered"])

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
    DataFrames.rename!(
        ode_solution_df,
        Dict(
            zip(
                filter(e -> e != "timestamp", names(ode_solution_df)),
                df_names
            )
        )
    )
    # convert to long format, requires reassignment
    ode_solution_df = DataFrames.stack(ode_solution_df, df_names)
    # split the variable column
    DataFrames.transform!(
        ode_solution_df,
        :variable => DataFrames.ByRow(x -> split(x, '.')) =>
            [:demo_group, :compartment]
    )
    DataFrames.select!(ode_solution_df, Not([:variable]))

    # no return type as the function modifies in place
end
