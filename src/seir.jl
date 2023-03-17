
# define function ODEs
function seir(du, u, p, t)
    S, E, I, R = u
    b, b2, g = p # beta, beta2 (E -> I), gamma, gamma2 (R -> S)
    du[1] = -b * S * I # change in susceptibles
    du[2] = b * S * I - b2 * E # change in exposed
    du[3] = b2 * E - g * I # change in infectious(ed)
    du[4] = g * I # change in recovered
end
