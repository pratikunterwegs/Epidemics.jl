"""
  Infection(name, r0, infectious_period, extra_arguments)

A structure to hold the infection parameters. These
currently include:
- 'name': a string for the name of the infection;
- 'r0': the basic reproductive number of the infection ``R_0``;
- 'infectious_period': the average period (in simulation time) for which
individuals are infectious;
- 'extra_arguments': a named `tuple` with extra parameters of the infection,
these may include values such as the pre-infectious period, or the
hospitalisation rate.

The default model provided in Epidemics.jl is [`epidemic_default`](@ref), 
which calculates:
- ``\\beta`` the transmission rate, which is the rate at which individuals 
move from the susceptible to the exposed compartment, as 
``\\beta = R_0 / \\text{infectious period}``
- ``\\alpha`` the rate at which exposed individuals enter the infectious
compartment, calculated as ``\\alpha = 1 / \\text{preinfectious_period}``
- ``\\gamma`` the recovery rate, the rate at which indviduals move from the
infectious to the recovered compartment, calculated as
``\\gamma = 1 / \\text{infectious_period}``.

The default model [`epidemic_default`](@ref) supports only a single,
population-wide value for each of the transition rates.

"""

mutable struct Infection
    name::String
    r0::Number
    infectious_period::Number
    extra_arguments::NamedTuple
end

function Infection(; name::String="none", r0::Number=1.5,
    infectious_period::Number=7,
    extra_arguments::NamedTuple=(preinfectious_period=5,))

    return Infection(name, r0, infectious_period, extra_arguments)
end
