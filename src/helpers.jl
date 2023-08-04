
include("population.jl")

using LinearAlgebra

"""
    prepare_contact_matrix(; n_groups)

Prepare a population contact matrix for an epidemic model.

## Named Arguments
- `population`: A `Population` object with information about the population
    affected by the epidemic. Must include a contact matrix, and
    a demography vector.

## Returns
A matrix of `ones` with the dimensions `n_groups * n_groups`.
"""
function prepare_contact_matrix(; population::Population)
    # input checking here
    # WIP check that square matrix has as many rows as demographic groups

    # scale the contact matrix by its largest real eigenvalue
    # eigvals is assumed to return only the real eigenvalues
    # may need correction
    contact_matrix = population.contact_matrix /
                     maximum(LinearAlgebra.eigvals(population.contact_matrix))

    # scale by the demography, divide each row by corresponding demography
    contact_matrix = contact_matrix ./ population.demography_vector

    return contact_matrix
end

"""
    default_initial_conditions(; n_groups, p_infected=, p_exposed)

Create some useful initial conditions for the default SEIRV epidemic model.

## Named Arguments
- `n_groups`: A `Number` for the number of demographic groups in the population.
- `p_infected`: The proportion of each demographic group that is infected and
    also infectious. This is expected to be a `Vector` of `Numbers` with each
    value between 0.0 and 1.0.
- `p_exposed`: The proportion of each demographic group that is exposed but not
    yet infectious. This is expected to be a `Vector` of `Numbers` with each
    value between 0.0 and 1.0.

## Returns
A matrix with the dimensions `n_groups * 5`, with each row representing a
    demographic groups and each column representing one of the five
    epidemiological compartments of the default model.
"""
function default_initial_conditions(; n_groups::Number=3,
    p_infected::Vector{Number}=[1e-6], p_exposed::Vector{Number}=[0])

    @assert isa(n_groups, Number) "Argument `n_groups` must be a single number."
    @assert all(p_infected .>= 0.0) && all(p_infected .<= 1.0)
    "Argument `p_infected` must have values between 0.0 and 1.0."
    @assert all(p_exposed .>= 0.0) && all(p_exposed .<= 1.0)
    "Argument `p_infected` must have values between 0.0 and 1.0."
    @assert all((p_infected .+ p_exposed) .<= 1.0)
    "Sum of proportion infected and exposed cannot be greater than 1.0"

    # prepare initial conditions for the default model
    default_s = 1.0 .- (p_infected .+ p_exposed)
    default_i = p_infected
    default_e = p_exposed
    default_r = repeat([0.0], n_groups)
    default_v = repeat([0.0], n_groups)

    return [
        default_s default_e default_i default_r default_v
    ] # rows represent age groups
end

"""
    default_contact_matrix(; n_groups)

Create a uniform contact matrix.

## Named Arguments
- `n_groups`: A `Number` for the number of demographic groups in the population.

## Returns
A matrix of `ones` with the dimensions `n_groups * n_groups`.
"""
function default_contact_matrix(; n_groups::Number=3)
    return ones(n_groups, n_groups)
end

"""
    prepare_initial_conditions(; population)

Prepare initial conditions from a `Population` object.

## Named Arguments
- `population`: A `Population` object with information about the population
    affected by the epidemic. Must include a matrix of initial conditions, and
    a demography vector.

## Returns
A matrix with the initial conditions scaled by the demography vector.
"""
function prepare_initial_conditions(; population::Population)
    # input checking here; check that initial initial_conditions
    # is a matrix and has as many rows as demography_vector

    # rowwise multiplication of initial conditions by demography_vector
    init = population.initial_conditions .* population.demography_vector

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
