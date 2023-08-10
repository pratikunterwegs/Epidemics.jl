"""
    prepare_args_default(parameters)

A function to prepare the parameters passed to the default epidemic model.
This function is intended to be called internally from
[`epidemic_default`](@ref).

## Named arguments

The function expects the parameters argument to be a four element `tuple`
with the following elements:
- `population`, for the `population` affected by the epidemic;
- `infection`, for the `infection` causing the epidemic;
- `intervention`, for any `Npi`s applied to the epidemic;
- `vaccination`, for any `vaccination`s applied to the population;
    
## Return

Returns a `tuple` of the parameters for the default epidemic model.
"""
function prepare_args_default(; population::Population, infection::Infection,
    intervention::Npi, vaccination::Vaccination)
    # prepare the population object by scaling the contact matrix
    population.contact_matrix = prepare_contact_matrix(
        population=population
    )

    # prepare parameters
    parameters = tuple(
        population,
        r0_to_beta(r0=infection.r0,
            infectious_period=infection.infectious_period
        ),
        preinfectious_period_to_alpha(
            preinfectious_period=infection.extra_arguments.preinfectious_period
        ),
        infectious_period_to_gamma(
            infectious_period=infection.infectious_period
        ),
        intervention,
        vaccination
    )

    return parameters
end
