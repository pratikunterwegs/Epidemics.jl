using DataFrames

function prepare_data_names(df::DataFrame,
        n_age_groups::Number = 3, compartments::Vector{String} = ["susceptible", "exposed",
            "infectious",
            "recovered", "vaccinated"])

    # create dummy age group names
    demo_groups = string.("demo_grp_", 1:n_age_groups)

    # prepare vector of age-and-compartment names
    df_names = Vector{String}(undef, n_age_groups * length(compartments))
    # prepare iterator
    k = 1

    # generate names
    for i in compartments
        for j in demo_groups
            df_names[k] = j * "." * i
            k = k + 1
        end
    end

    # return names
    df_names
end

function prepare_data(ode_output::DataFrame; n_age_groups::Number = 3,
        compartments::Vector{String} = ["susceptible", "exposed",
            "infectious",
            "recovered", "vaccinated"])
    # no input checking on this internal function
    # create dummy age group names
    demo_groups = string.("demo_grp_", 1:n_age_groups)

    # prepare vector of age-and-compartment names
    df_names = prepare_data_names(ode_output, n_age_groups, compartments)

    # Get the names of columns to rename (excluding the unchanged column)
    columns_to_rename = setdiff(names(ode_output), ["timestamp"])
    # Create an array of Pairs for renaming selected columns
    name_pairs = [old => new for (old, new) in zip(columns_to_rename, df_names)]

    # Rename selected columns using rename!
    rename!(ode_output, name_pairs)

    # convert to long format, requires reassignment
    data = stack(ode_output, df_names)

    # split the variable column
    transform!(data,
        :variable => ByRow(x -> split(x, '.')) => [:demo_group, :compartment])
    select!(data, Not([:variable]))

    # return reshaped and renamed data
    return data
end

export prepare_data
