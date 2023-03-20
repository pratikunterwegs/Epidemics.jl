
using LinearAlgebra

function prepare_contact_matrix(; contact_matrix, demography_vector)
    # input checking here
    # WIP check that square matrix has as many rows as demographic groups

    # scale the contact matrix by its largest real eigenvalue
    # eigvals is assumed to return only the real eigenvalues
    # may need correction
    contact_matrix = contact_matrix /
                     maximum(LinearAlgebra.eigvals(contact_matrix))

    # scale by the demography, divide each row by corresponding demography
    contact_matrix = contact_matrix ./ demography_vector

    return contact_matrix
end

function prepare_initial_conditions(; initial_conditions, demography_vector)
    # input checking here; check that initial initial_conditions
    # is a matrix and has as many rows as demography_vector

    # rowwise multiplication of initial conditions by demography_vector
    init = initial_conditions .* demography_vector

    # return linearised values
    return init
end

function r0_to_beta(; r0, infectious_period)
    return r0 ./ (infectious_period)
end

function preinfectious_period_to_beta2(; preinfectious_period)
    return 1.0 ./ preinfectious_period
end

function infectious_period_to_gamma(; infectious_period)
    return 1.0 ./ infectious_period
end
