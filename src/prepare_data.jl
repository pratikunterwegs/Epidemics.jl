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

    # Get the names of columns to rename (excluding the unchanged column)
    columns_to_rename = setdiff(names(ode_solution_df), ["timestamp"])
    # Create an array of Pairs for renaming selected columns
    name_pairs = [old => new for (old, new) in zip(columns_to_rename, df_names)]
    # Rename selected columns using rename!
    rename!(ode_solution_df, name_pairs)

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
