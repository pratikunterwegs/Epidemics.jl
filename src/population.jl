"""
  Population(name, demography_vector, initial_conditions, contact_matrix)

A structure to hold population characteristics, including:
- 'name': A name for the population.
- 'demography_vector': A numeric vector of the number of individuals in each
age or demographic group of the population.
- 'initial_conditions': A numeric matrix representing the proportions of each
age or demographic group that are in one of the epidemiological compartments.
- 'contact_matrix': A matrix giving the contacts between the demographic groups
in the population. Must be a square matrix.

"""
mutable struct Population
    name::String
    demography_vector::Vector
    initial_conditions::Matrix{Number}
    contact_matrix::Matrix{Number}
end

# external constructor allowing passing of some default values
function Population(; demography_vector::Vector=67e6 * [0.23, 0.4, 0.37],
    initial_conditions::Matrix=default_initial_conditions(),
    contact_matrix::Matrix=default_contact_matrix() * 5)
    # input checking goes here

    # use default constructor with unnamed 
    return Population("none",
        demography_vector,
        initial_conditions,
        contact_matrix
    )
end
