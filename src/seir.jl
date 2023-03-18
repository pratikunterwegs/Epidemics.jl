
# define function ODEs
function seir(du, u, p, t)
    S, E, I, R = u

    # read beta, beta2, and gamma
    β = @view p[1]
    β2 = @view p[2]  ## WIP: Use NamedArrays to reduce confusion
    γ = @view p[3]

    # access the initial condition matrix
    S = @view u[:, 1]
    E = @view u[:, 2]
    I = @view u[:, 3]
    R = @view u[:, 4]

    # change in compartments
    dS = @view du[:,1]
    dE = @view du[:,2]
    dI = @view du[:,3]
    dR = @view du[:,4]

    @. du[:, 1] = -β * S * I
    @. du[:, 2] = (β * S * I) - (β2 * E)
    @. du[:, 3] = (β2 * E) - (γ * I)
    @. du[:, 4] = γ * I
end
