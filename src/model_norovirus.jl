
using OrdinaryDiffEq
using LinearAlgebra

# copied from socialmixr::polymod UK
function noromod_contacts()
    return [[1.9157895 1.2818462 4.999847 0.3658614]
            [0.5994315 7.9460784 6.575480 0.6006796]
            [0.4341422 1.2209536 9.169207 1.0025452]
            [0.1306272 0.4586235 4.122357 1.7142857]]
end

function noromod_demography()
    return [3453670.0, 7385454.0, 39774569.0, 9673058.0]
end

function noromod_aging()
    ages = [4, 14, 64, 80]
    da = diff([0; ages])
    aging = Matrix(Diagonal(-1.0 ./ da))
    aging[2, 1] = 0.25
    aging[3, 2] = 0.1
    aging[4, 3] = 0.02

    return aging / 365.0
end

function seasonal_forcing(t, w1, w2)
    return 1.0 + w1 * cos((2.0 * pi * t / 364.0) + w2)
end

function noromod_initial_state()
    init = [[3857263.0, 8103718.0, 42460865.0, 12374961.0]
            [100.0, 0.0, 0.0, 0.0]
            [10.0, 0.0, 0.0, 0.0]
            [0, 0, 0, 0]
            [0.0, 0, 0, 0]]
    compartments = 5
    age_groups = 4

    init = reshape(init, age_groups, compartments)
    init_vax = zeros(age_groups, compartments)
    init_vax2 = copy(init_vax)

    return cat(init, init_vax, init_vax2, dims = 3)
end

"""
    epidemic_norovirus_ode!(du, u, parameters, t)

The ODE system for the norovirus model. This function is intended to be called
    internally from [`epidemic_norovirus`](@ref).    
"""
function epidemic_norovirus_ode!(du, u, parameters, t)
    # du auto-magically takes the type of u (?)
    # each element of the tuple is one of the required params
    contacts, sigma, rho, delta, epsilon, psi, gamma, b, d,
    aging, w1, w2, q, phi, upsilon = parameters

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, 1, :]
    E = @view u[:, 2, :]
    Is = @view u[:, 3, :]
    Ia = @view u[:, 4, :]
    R = @view u[:, 5, :]

    # calculate seasonal term
    seasonal_term = seasonal_forcing(t, w1, w2)

    # calculate infection potential
    infection_potential = q .*
                          (seasonal_term * (contacts * sum(Is .+ (Ia * rho), dims = 2)))

    # calculate new infections and re-infections
    # NOTE: element-wise multiplication
    new_I = S .* infection_potential
    re_I = R .* infection_potential

    # calculate births
    births = b * sum(@view u[:, 1:5, :])
    births = reshape([births; repeat([0], 11)], 4, 3)

    # vaccination and waning in and out
    S_vax_out = S .* phi
    S_waning_out = S .* upsilon
    R_vax_out = R .* phi
    R_waning_out = R .* upsilon
    R_S_direct_waning = R * delta * gamma
    R_S_direct_waning[:, 1] .= 0.0 # remove extra waning term from R -> S

    # views to the change array slice
    dS = @view du[:, 1, :]
    dE = @view du[:, 2, :]
    dIs = @view du[:, 3, :]
    dIa = @view du[:, 4, :]
    dR = @view du[:, 5, :]

    # aging changes vector
    aging_vec = [aging * @view u[:, i, :] for i in 1:5]

    # calculate change in compartment size and update the change matrix dU
    # note the use of @. for broadcasting, equivalent to .=
    # change in susceptibles
    @. dS = -new_I + (delta * R) + (aging_vec[1] + births) -
            (d * S) -
            (S_vax_out + S_waning_out) +
            (S_vax_out[:, [3, 1, 2]] + S_waning_out[:, [2, 3, 1]]) +
            R_S_direct_waning

    # change in exposed
    @. dE = new_I - (epsilon * E) + (aging_vec[2]) - (d * E)

    # calculate exposed to Is and Ia
    E_sigma = (epsilon * E) * sigma
    E_inv_sigma = (epsilon * E) - E_sigma

    # change in infectious symptomatic
    @. dIs = -(psi * Is) + E_sigma + (aging_vec[3]) - (d * Is)

    # change in infectious asymptomatic
    @. dIa = re_I + (psi * Is) + E_inv_sigma -
             (gamma * Ia) + (aging_vec[4]) - (d * Ia)

    # change in recovered
    @. dR = -re_I + (gamma * Ia) - (delta * R) + (aging_vec[5]) -
            (d * R) -
            (R_vax_out + R_waning_out) +
            (R_vax_out[:, [3, 1, 2]] + R_waning_out[:, [2, 3, 1]]) -
            R_S_direct_waning
end

"""
    epidemic_norovirus(initial_state, contacts, aging, sigma, phi, upsilon,
        rho, w1, w2, delta, q1, q2, b, d, epsilon, psi, gamma, n_age_groups,
        time_end, increment
    )

Model the progression of a norovirus epidemic with multiple optional vaccination
strata.
"""
function epidemic_norovirus(;
        initial_state = noromod_initial_state(),
        contacts = noromod_contacts(),
        demography = noromod_demography(),
        aging = noromod_aging(),
        sigma = Matrix(Diagonal([0.82, 0.41, 0.41])),
        phi = [[1e-4, 0, 0,1e-4] [1e-4, 0, 0,1e-4] [0, 0, 0,0]],
        upsilon = 1 ./ ([[0, 0, 0,0] [4.4, 4.4, 4.4,4.4] [4.4, 4.4, 4.4,4.4]] .* 365.0),
        rho = 0.05,
        w1 = 3.6 / 100.0,
        w2 = 0.5 / 100.0,
        delta = 1.0 / (4.4 * 365.0),
        q1 = 0.195,
        q2 = 0.039,
        b = 11.4e-3 / 365.0,
        d = 0.0,
        epsilon = 1.0,
        psi = 0.5,
        gamma = 0.1,
        time_end::Number = 11100.0,
        increment::Number = 1.0)

    # prepare age-specific transmission probability
    q = [q1, q2, q2, q2]
    # fix division by zero in waning parameter
    upsilon[isinf.(upsilon)] .= 0.0

    # scale contacts by demography matrix
    contacts = transpose(transpose(contacts) ./ demography)

    # no seasonal offsettiing for this scenario model
    parameters = tuple(contacts, sigma, rho, delta,
        epsilon, psi, gamma, b, d, aging, w1, w2, q, phi, upsilon)

    # prepare the timespan
    timespan = (0.0, time_end)

    # define the ode problem
    ode_problem = ODEProblem(epidemic_norovirus_ode!, initial_state, timespan, parameters)

    # get the solution
    ode_solutions = solve(ode_problem,
        AutoTsit5(Rosenbrock23()),
        saveat = increment)

    return ode_solutions
end

export noromod_aging, noromod_contacts, noromod_demography, noromod_initial_state
export epidemic_norovirus
