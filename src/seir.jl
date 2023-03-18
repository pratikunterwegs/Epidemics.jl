"""
    seir(du, u, p, t)

A simple SEIR epidemic model function that allows for multiple demographic
    groups. This function is intended to be called internally from [`epidemic`](@ref).

"""
function seir(du, u, p, t)
    # assumes that each element of the vector is a 
    # vector of rates
    β, β2, γ = p

    # view the values of each compartment per each age group
    S = @view u[:, 1]
    E = @view u[:, 2]
    I = @view u[:, 3]
    R = @view u[:, 4]

    # calculate change in compartment size
    @. du[:, 1] = -β .* S .* I
    @. du[:, 2] = β * S * I - (β2 * E)
    @. du[:, 3] = β2 * E - γ * I
    @. du[:, 4] = γ * I
end
