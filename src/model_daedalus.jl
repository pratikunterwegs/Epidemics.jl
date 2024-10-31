
using OrdinaryDiffEq
using LinearAlgebra

# copied from jameel-institute/daedalus
function australia_demography()
    return [1669707, 4769747, 14926119, 4134308]
end

function australia_contacts()
    cm = [[3.7187500 2.5982168 5.739112 0.2728101]
          [0.9095369 13.0623374 5.741992 0.5229291]
          [0.6420045 1.8348941 11.256655 1.0003495]
          [0.1347582 0.6540519 3.760931 2.5421895]]

    return cm
end

function aus_workers()
    return [
        331714, 10730, 110740, 113170, 65189, 224517, 35631, 47879, 45956,
        8604, 31481, 18322, 30266, 30367, 52435, 69643, 38263, 21966,
        59115, 33101, 32696, 195972, 79259, 73591, 1199314, 1783483, 319770,
        66707, 62051, 151331, 89280, 899587, 91912, 90495, 270878, 418059,
        166920, 823060, 567062, 859198, 1059016, 1686004, 277154, 268246, 0
    ] .+ 1
end

function worker_contacts(workers = aus_workers())
    # in proportion to workforce and scaled by workforce
    return (2 .+ workers / sum(workers)) ./ workers
end

# Function to prepare 49x49 community contacts matrix
function prepare_contacts(cm = australia_contacts())
    cm_x = ones(49, 49) .* cm[3, 3]
    cm_x[1:4, 1:4] = cm
    cm_x[1:4, 5:49] .= cm[:, 3]
    cm_x[5:49, 1:4] .= reshape(cm[3, :], 1, 4)

    return cm_x
end

# Function to prepare 49 element demography vector for I/N in FOI calculation
function prepare_demog(demog = australia_demography(), workers = aus_workers())
    demog_x = ones(49) * demog[3]
    demog_x[1:4] = demog

    return demog_x
end

"""
    australia_initial_state()

Initial state with seven compartments: S, E, Is, Ia, R, D, and H, and four
age groups and 49 economic sectors (comprised of working age individuals), and
two vaccination strata for unvaccinated and vaccinated individuals.
"""
function australia_initial_state(
        demography = australia_demography(), workers = aus_workers())
    p_infected = 1e-6
    p_susc = 1 - p_infected
    compartments = 7
    age_groups = 4
    econ_groups = 45
    vaccine_strata = 2

    init = [p_susc, 0.0, p_infected, 0.0, 0.0, 0.0, 0.0]
    init = reshape(init, 1, compartments)
    init = repeat(init, age_groups + econ_groups)

    dummy = zeros(age_groups + econ_groups, compartments, vaccine_strata)

    # assign to first layer
    dummy[:, :, 1] = init

    # process demography to calculate non-working working age, and concat
    # demography and worker counts
    i_working_age = 3
    inactive_workers = demography[i_working_age] - sum(workers)
    demography[i_working_age] = inactive_workers

    demography = [demography; workers]

    # multiply by demography - row i times element i
    init = dummy .* demography

    return init
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
    i_age_groups = 1:4
    i_working_age = 3
    i_econ_groups = 5:49

    # du auto-magically takes the type of u (?)
    # each element of the tuple is one of the required params
    contacts, cw, beta, sigma, p_sigma, epsilon,
    rho, eta, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi,
    t_vax = p

    # time dependent vaccination
    nu = t > t_vax ? nu : 0.0

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, 1, :]
    E = @view u[:, 2, :]
    Is = @view u[:, 3, :]
    Ia = @view u[:, 4, :]
    H = @view u[:, 5, :]
    R = @view u[:, 6, :]
    D = @view u[:, 7, :]

    # calculate new infections and re-infections
    community_infectious = sum(Is .+ Ia * epsilon, dims = 2)
    foi = beta * contacts * community_infectious

    # NOTE: element-wise multiplication
    new_I = S .* foi

    # workplace
    new_I_work = S[i_econ_groups, :] .* cw .*
                 (Is[i_econ_groups, :] .+ (Ia[i_econ_groups, :] * epsilon))

    new_I[i_econ_groups, :] += new_I_work

    # views to the change array slice
    dS = @view du[:, 1, :]
    dE = @view du[:, 2, :]
    dIs = @view du[:, 3, :]
    dIa = @view du[:, 4, :]
    dH = @view du[:, 5, :]
    dR = @view du[:, 6, :]
    dD = @view du[:, 7, :]

    # new vaccinations from susceptibles
    new_Svax = @view S[:, 1]
    new_Swane = @view S[:, 2]

    # new vaccinations from recovered
    new_Rvax = @view R[:, 1]
    new_Rwane = @view R[:, 2]

    # calculate change in compartment size and update the change matrix dU
    # note the use of @. for broadcasting, equivalent to .=
    # change in susceptibles
    @. dS = -new_I + (rho * R)
    @. dS[:, 1] += (-new_Svax * nu + new_Swane * psi)
    @. dS[:, 2] += (new_Svax * nu - new_Swane * psi)

    # change in exposed
    @. dE = new_I - (sigma * E)

    # calculate exposed to Is and Ia
    E_sigma = (sigma * E) * p_sigma
    E_inv_sigma = (sigma * E) * (1 - p_sigma)

    # change in infectious symptomatic
    @. dIs = E_sigma - ((gamma_Is .+ eta) .* Is)

    # change in infectious asymptomatic
    @. dIa = E_inv_sigma - (gamma_Ia * Ia)

    # change in hospitalised
    @. dH = (eta .* Is) - ((gamma_H + omega) .* H)

    # change in recovered
    @. dR = (gamma_Ia * Ia) + (gamma_Is * Is) +
            (gamma_H .* H) - (rho * R)

    @. dR[:, 1] += (-new_Rvax * nu + new_Rwane * psi)
    @. dR[:, 2] += (new_Rvax * nu - new_Rwane * psi)

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
        cw = worker_contacts(),
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
        nu = (2 / 100) / 7,
        psi = 1 / 270,
        t_vax = 200,
        time_end::Number = 300.0,
        increment::Number = 1.0)

    # prepare transmission parameter beta as r0 / max(eigenvalue(contacts))
    beta = r0 / maximum(eigvals(contacts))

    # scale contacts by demography; divide col-wise
    contacts = contacts * diagm(1 ./ demography)

    # prepare parameters to account for economic sector groups
    i_working_age = 3
    econ_sectors = 45

    eta = [eta; repeat([eta[i_working_age]], econ_sectors)]
    omega = [omega; repeat([omega[i_working_age]], econ_sectors)]
    gamma_H = [gamma_H; repeat([gamma_H[i_working_age]], econ_sectors)]

    # combined parameters into a tuple: these cannot be modified during solving
    parameters = tuple(contacts, cw, beta, sigma, p_sigma, epsilon,
        rho, eta, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi, t_vax)

    # prepare the timespan
    timespan = (0.0, time_end)

    # define the ode problem
    ode_problem = ODEProblem(
        epidemic_daedalus_ode!, initial_state, timespan, parameters
    )

    # cb = ContinuousCallback(condition, affect!)

    # get the solution
    # NOTE: not compatible with autodiff!
    ode_solution = solve(ode_problem,
        Rosenbrock23(autodiff = false))
    # callback = cb)

    return ode_solution
end

export australia_contacts, australia_demography, australia_initial_state
export epidemic_daedalus
