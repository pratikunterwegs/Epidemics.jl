
using OrdinaryDiffEq
using LinearAlgebra

# copied from jameel-institute/daedalus
function australia_demography()
    return [1669707, 4769747, 14926119, 4134308]
end

function australia_contacts()
    return [[3.7187500 2.5982168 5.739112 0.2728101]
            [0.9095369 13.0623374 5.741992 0.5229291]
            [0.6420045 1.8348941 11.256655 1.0003495]
            [0.1347582 0.6540519 3.760931 2.5421895]]
end

"""
    australia_initial_state()

Initial state copied from the noromod code, with five compartments: S, E, Is,
Ia, R.
"""
function australia_initial_state(demography)
    p_infected = 1e-6
    p_susc = 1 - p_infected
    init = [repeat([p_susc], 4)
            repeat([p_infected], 4)
            [0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0]
            [0, 0, 0, 0]
            [0, 0, 0, 0]]
    compartments = 7
    age_groups = 4

    init = reshape(init, age_groups, compartments)

    return init .* demography
end

"""
    condition(u, t, integrator)

A condition function that checks if a root is found.
"""
function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
    # pick a reasonable threshold
    threshold = 1e3
    H = @view u[:, 5]
    total_hosp = sum(H)

    total_hosp - threshold
end

"""
    affect!(integrator)

An event function.
"""
function affect!(integrator)
    # scale beta by 0.5
    # ISSUE: cannot change parameter values as passed in tuple
    setindex!(integrator.p, integrator.p[2] * 0.5, 2) 
end

"""
    epidemic_daedalus_ode!(du, u, p, t)

The ODE system for the DAEDALUS model. This function is intended to be called
    internally from [`epidemic_daedalus`](@ref).    
"""
function epidemic_daedalus_ode!(du, u, p, t)
    # du auto-magically takes the type of u (?)
    # each element of the tuple is one of the required params
    contacts, beta, sigma, p_sigma, epsilon,
    rho, eta, omega, gamma_Ia, gamma_Is, gamma_H = p

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, 1]
    E = @view u[:, 2]
    Is = @view u[:, 3]
    Ia = @view u[:, 4]
    H = @view u[:, 5]
    R = @view u[:, 6]
    D = @view u[:, 7]

    # calculate new infections and re-infections
    foi = beta * contacts * sum(Is .+ (Ia * epsilon), dims = 2)
    # NOTE: element-wise multiplication
    new_I = S .* foi

    # views to the change array slice
    dS = @view du[:, 1]
    dE = @view du[:, 2]
    dIs = @view du[:, 3]
    dIa = @view du[:, 4]
    dH = @view du[:, 5]
    dR = @view du[:, 6]
    dD = @view du[:, 7]

    # calculate change in compartment size and update the change matrix dU
    # note the use of @. for broadcasting, equivalent to .=
    # change in susceptibles
    @. dS = -new_I + (rho * R)

    # change in exposed
    @. dE = new_I - (sigma * E)

    # calculate exposed to Is and Ia
    E_sigma = (sigma * E) * p_sigma
    E_inv_sigma = (sigma * E) * (1 - p_sigma)

    # change in infectious symptomatic
    @. dIs = E_sigma - ((gamma_Is + eta) .* Is)

    # change in infectious asymptomatic
    @. dIa = E_inv_sigma - (gamma_Ia * Ia)

    # change in hospitalised
    @. dH = (eta .* Is) - ((gamma_H + omega) .* H)

    # change in recovered
    @. dR = (gamma_Ia * Ia) + (gamma_Is * Is) +
            (gamma_H .* H) - (rho * R)

    # change in dead
    @. dD = omega .* H
end

"""
    epidemic_daedalus()

Model the progression of a daedalus epidemic with multiple optional vaccination
strata.
"""
function epidemic_daedalus(;
        initial_state = australia_initial_state(australia_demography()),
        contacts = australia_contacts(),
        demography = australia_demography(),
        r0 = 1.3,
        sigma = 0.217,
        p_sigma = 0.867,
        epsilon = 0.58,
        rho = 0.003,
        eta = [0.018, 0.082, 0.018, 0.246],
        omega = [0.012, 0.012, 0.012, 0.012],
        gamma_Ia = 0.476,
        gamma_Is = 0.25,
        gamma_H = [0.034, 0.034, 0.034, 0.034],
        time_end::Number = 300.0,
        increment::Number = 1.0)

    # prepare transmission parameter beta as r0 / max(eigenvalue(contacts))
    beta = r0 / maximum(eigvals(contacts))

    # scale contacts by demography; divide col-wise
    contacts = contacts * diagm(1 ./ demography)

    # no seasonal offsettiing for this scenario model
    parameters = tuple(contacts, beta, sigma, p_sigma, epsilon,
        rho, eta, omega, gamma_Ia, gamma_Is, gamma_H)

    # prepare the timespan
    timespan = (0.0, time_end)

    # define the ode problem
    ode_problem = ODEProblem(epidemic_daedalus_ode!, initial_state, timespan, parameters)

    # cb = ContinuousCallback(condition, affect!)

    # get the solution
    ode_solution = solve(ode_problem,
        Rosenbrock23(),
        saveat = increment)
        # callback = cb)

    return ode_solution
end

export australia_contacts, australia_demography, australia_initial_state
export epidemic_daedalus
