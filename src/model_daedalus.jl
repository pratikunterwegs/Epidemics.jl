
using OrdinaryDiffEq
using LinearAlgebra

# copied from socialmixr::polymod UK
function daedalus_contacts()
    return [[1.9157895 1.2818462 4.999847 0.3658614]
            [0.5994315 7.9460784 6.575480 0.6006796]
            [0.4341422 1.2209536 9.169207 1.0025452]
            [0.1306272 0.4586235 4.122357 1.7142857]]
end

function australia_demography()
    return [1669707, 4769747, 14926119, 4134308]
end

function australia_contacts()
    return [[3.7187500 2.5982168 5.739112 0.2728101]
            [0.9095369 13.0623374 5.741992 0.5229291]
            [0.6420045 1.8348941 11.256655 1.0003495]
            [0.1347582 0.6540519 3.760931 2.5421895]]
end

function daedalus_demography()
    return [3453670.0, 7385454.0, 39774569.0, 9673058.0]
end

"""
    australia_initial_state()

Initial state copied from the noromod code, with five compartments: S, E, Is,
Ia, R.
"""
function australia_initial_state()
    init = [[1669707, 4769747, 14926119, 4134308]
            [100.0, 0.0, 0.0, 0.0]
            [10.0, 0.0, 0.0, 0.0]
            [0, 0, 0, 0]
            [0.0, 0, 0, 0]]
    compartments = 5
    age_groups = 4

    init = reshape(init, age_groups, compartments)

    return init
end

"""
    condition(u, t, integrator)

A condition function that checks if a root is found.
"""
function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
    # pick a reasonable threshold
    threshold = 5e4
    Is = @view u[:, 3]
    total_infected = sum(Is)

    total_infected - threshold
end

"""
    affect!(integrator)

An event function.
"""
function affect!(integrator)
    # scale contacts by 0.5
    print("modifying contacts!")
    integrator.u[1] = integrator.u[1] .* 0.1
end

"""
    epidemic_daedalus_ode!(du, u, parameters, t)

The ODE system for the DAEDALUS model. This function is intended to be called
    internally from [`epidemic_daedalus`](@ref).    
"""
function epidemic_daedalus_ode!(du, u, parameters, t)
    # du auto-magically takes the type of u (?)
    # each element of the tuple is one of the required params
    contacts, sigma, rho, delta, epsilon, psi, gamma, b, d,
    q = parameters

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, 1]
    E = @view u[:, 2]
    Is = @view u[:, 3]
    Ia = @view u[:, 4]
    R = @view u[:, 5]

    # calculate infection potential
    infection_potential = q .* (contacts * sum(Is .+ (Ia * rho), dims = 2))

    # calculate new infections and re-infections
    # NOTE: element-wise multiplication
    new_I = S .* infection_potential
    re_I = R .* infection_potential

    # calculate births
    births = b * sum(@view u[:, 1:5])
    births = reshape([births; repeat([0], 3)], 4, 1)

    # views to the change array slice
    dS = @view du[:, 1]
    dE = @view du[:, 2]
    dIs = @view du[:, 3]
    dIa = @view du[:, 4]
    dR = @view du[:, 5]

    # calculate change in compartment size and update the change matrix dU
    # note the use of @. for broadcasting, equivalent to .=
    # change in susceptibles
    @. dS = -new_I + (delta * R) + births -
            (d * S)

    # change in exposed
    @. dE = new_I - (epsilon * E) - (d * E)

    # calculate exposed to Is and Ia
    E_sigma = (epsilon * E) * sigma
    E_inv_sigma = (epsilon * E) - E_sigma

    # change in infectious symptomatic
    @. dIs = -(psi * Is) + E_sigma - (d * Is)

    # change in infectious asymptomatic
    @. dIa = re_I + (psi * Is) + E_inv_sigma -
             (gamma * Ia) - (d * Ia)

    # change in recovered
    @. dR = -re_I + (gamma * Ia) - (delta * R) -
            (d * R)
end

"""
    epidemic_daedalus(initial_state, contacts, sigma, phi, upsilon,
        rho, w1, w2, delta, q1, q2, b, d, epsilon, psi, gamma, n_age_groups,
        time_end, increment
    )

Model the progression of a daedalus epidemic with multiple optional vaccination
strata.
"""
function epidemic_daedalus(;
        initial_state = australia_initial_state(),
        contacts = australia_contacts(),
        demography = australia_demography(),
        sigma = 0.82,
        rho = 0.05,
        delta = 1.0 / (4.4 * 365.0),
        q1 = 0.195,
        q2 = 0.039,
        b = 11.4e-3 / 365.0,
        d = 0.0,
        epsilon = 1.0,
        psi = 0.5,
        gamma = 0.1,
        time_end::Number = 300.0,
        increment::Number = 1.0)

    # prepare age-specific transmission probability
    q = [q1, q2, q2, q2]

    # scale contacts by demography matrix
    contacts = transpose(transpose(contacts) ./ demography)

    # no seasonal offsettiing for this scenario model
    parameters = tuple(contacts, sigma, rho, delta,
        epsilon, psi, gamma, b, d, q)

    # prepare the timespan
    timespan = (0.0, time_end)

    # define the ode problem
    ode_problem = ODEProblem(epidemic_daedalus_ode!, initial_state, timespan, parameters)

    # cb = ContinuousCallback(condition, affect!)

    # get the solution
    ode_solution = solve(ode_problem,
        Rosenbrock23(),
        saveat = increment)#,
        # callback = cb)

    return ode_solution
end

export australia_contacts, australia_demography, australia_initial_state
export epidemic_daedalus
