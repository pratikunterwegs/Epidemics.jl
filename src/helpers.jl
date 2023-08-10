
include("population.jl")

using LinearAlgebra

"""
    prepare_contact_matrix(; population)

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
                     maximum(eigvals(population.contact_matrix))

    # scale by the demography, divide each row by corresponding demography
    contact_matrix = contact_matrix ./ population.demography_vector

    return contact_matrix
end

"""
    default_initial_conditions(; n_groups=3, p_infected=[1e-6, 1e-6, 1e-6],
        p_exposed=[0.0, 0.0, 0.0])

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
    p_infected::Vector=[1e-6, 1e-6, 1e-6], p_exposed::Vector=[0.0, 0.0, 0.0])

    @assert all(p_infected .>= 0.0) && all(p_infected .<= 1.0)
    "Argument `p_infected` must have values between 0.0 and 1.0."
    @assert all(p_exposed .>= 0.0) && all(p_exposed .<= 1.0)
    "Argument `p_exposed` must have values between 0.0 and 1.0."
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
    default_contact_matrix(; n_groups=3)

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
    # return linearised values
    return population.initial_conditions .* population.demography_vector
end

"""
    r0_to_beta(; r0, infectious_period)

Calculate the transmission rate ``\\beta`` from the basic reproductive number
``R_0`` and the mean infectious period.

## Named Arguments

- `r0`: A single number for the basic reproduction number of an infection.

- `infectious_period`: A single number for the mean duration in days that
individuals are infectious.

## Returns
A single number representing the transmission rate of the infection ``\\beta``
"""
function r0_to_beta(; r0::Number, infectious_period::Number)
    @assert r0 > 0.0
    "`r0` must be greater than 0.0"
    @assert infectious_period > 0.0
    "`infectious_period` must be greater than 0.0"
    return r0 / (infectious_period)
end

"""
    preinfectious_period_to_alpha(; r0, preinfectious_period)

Calculate the rate of transition from the 'exposed' to the 'infectious'
compartment, ``\\alpha`` from the mean period between exposure and the initial
occurrence of symptoms.

## Named Arguments

- `preinfectious_period`: A single number for the mean duration in days between
individuals being exposed to infection and becoming infectious.

## Returns
A single number representing the transmission rate of the infectious ``\\alpha``
"""
function preinfectious_period_to_alpha(; preinfectious_period::Number)
    @assert preinfectious_period > 0.0
    "`preinfectious_period` must be greater than 0.0"
    return 1.0 / preinfectious_period
end

"""
    infectious_period_to_gamma(; infectious_period)

Calculate the recovery rate of the infection ``\\gamma`` from the infectious
period.

## Named Arguments

- `infectious_period`: A single number for the mean duration in days that
individuals are infectious.

## Returns
A single number representing the recovery rate of the infectious ``\\gamma``
"""
function infectious_period_to_gamma(; infectious_period::Number)
    @assert infectious_period > 0.0
    "`infectious_period` must be greater than 0.0"
    return 1.0 / infectious_period
end

export default_initial_conditions, r0_to_beta, preinfectious_period_to_alpha,
    infectious_period_to_gamma
