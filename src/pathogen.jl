"""
  Pathogen(r0, preinfectious_period, infectious_period)

A structure to hold the age-specific pathogen or infection parameters. These
    currently include:
    - 'r0': the basic reproductive number of the infection ``R_0``;
    - 'preinfectious_period': the average period (in simulation time - taken as days
    ) between individuals being exposed to the pathogen and becoming infectious;
    - 'infectious_period': the average period (in simulation time) for which
    individuals are infectious.

    The default model provided in Epidemics.jl is [`epidemic_default!`](@ref), 
    which calculates:
    - ``\\beta`` the transmission rate, which is the rate at which individuals 
    move from the susceptible to the exposed compartment, as 
    ``\\beta = R_0 / \\text{infectious period}``
    - ``\\alpha`` the rate at which exposed individuals enter the infectious
    compartment, calculated as ``\\alpha = 1 / \\text{preinfectious_period}``
    - ``\\gamma`` the recovery rate, the rate at which indviduals move from the
    infectious to the recovered compartment, calculated as
    ``\\gamma = 1 / \\text{infectious_period}``.

    The default model [`epidemic_default!`](@ref) supports only a single,
    population-wide value for each of the transition rates.

"""

mutable struct Pathogen
    r0::Vector{Number}
    preinfectious_period::Vector{Number}
    infectious_period::Vector{Number}
end

function Pathogen(; r0=1.5, preinfectious_period=3, infectious_period=7)
    # convert contact reduction to a vector if a single number
    if length(r0) == 1
        r0 = [r0]
    end
    if length(preinfectious_period) == 1
        preinfectious_period = [preinfectious_period]
    end
    if length(infectious_period) == 1
        infectious_period = [infectious_period]
    end

    return Pathogen(r0, preinfectious_period, infectious_period)
end
