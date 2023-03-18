
a = 1
using OrdinaryDiffEq

include("seir.jl")

# Parameter
# a vector of vectors of beta, beta2, and gamma
parameters = [[0.14, 0.12, 0.2], [0.5, 0.5, 0.5], [0.05, 0.05, 0.05]]
init = [0.9 0.9 0.9; 0.1 0.1 0.1; 0.0 0.0 0.0; 0.0 0.0 0.0]'
timespan = (0.0, 100)

function Φ!(du, u, p, t)
    # assumes that each element of the vector is a 
    # vector of rates
    β = @view p[1] # access transmission rates
    β2 = @view p[2] # access transmission rates
    γ = @view p[3] # acess recovery rates

    S = @view u[:, 1]
    E = @view u[:, 2]
    I = @view u[:, 3]
    R = @view u[:, 4]
    
    dS = @view du[:, 1]
    dE = @view du[:, 2]
    dI = @view du[:, 3]
    dR = @view du[:, 4]

    @. du[:, 1] = -β * S * I
    @. du[:, 2] = β * S * I - (β2 * E)
    @. du[:, 3] = β2 * E - γ * I
    @. du[:, 4] = γ * I
end

# Problem Definition
problem = ODEProblem(Φ!, init, timespan, parameters)
problem = ODEProblem(seir, init, timespan, parameters)

# Problem Solution
solution = solve(problem, Tsit5());

using Plots
plot(solution,xlabel="Time",ylabel="Number")

