"""
    epidemic_default!(du, u, parameters, t)

A simple SEIR epidemic model function that allows for multiple demographic
    groups. This function is intended to be called internally from
    [`epidemic`](@ref).
    The function expects the parameters argument to be a four element vector
    with the following elements:
    - ``\\beta``, the transmission rate;
    - ``\\alpha``, the rate of conversion from exposed to infectious;
    - ``\\gamma``, the rate of recovery;
    - a matrix specifying the contacts between demographic groups;
    - an `Npi` object specifying the intervention applied to each age
    group, see [`Npi`](@ref)
    
"""
function epidemic_default!(du, u, parameters, t)
    # assumes that each element of the vector is one of the required params
    β, α, γ, contact_matrix, intervention, vaccination = parameters

    # modify contact matrix if the intervention is active
    if (t > intervention.time_begin) & (t < intervention.time_end)
        contact_matrix = contact_matrix .*
                         (1.0 .- intervention.contact_reduction)
    end

    # get current ν
    ν_now = current_nu(time=t, vaccination=vaccination)

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, 1]
    E = @view u[:, 2]
    I = contact_matrix * @view u[:, 3] # matrix mult for cm * I
    I_ = @view u[:, 3] # unmultiplied I for operations involving only I
    R = @view u[:, 4]
    V = @view u[:, 5]

    # views to the change matrix, dU
    dS = @view du[:, 1]
    dE = @view du[:, 2]
    dI = @view du[:, 3]
    dR = @view du[:, 4]
    dV = @view du[:, 5]

    # calculate change in compartment size and update the change matrix dU
    # note the use of @. for broadcasting, equivalent to .=
    @. dS = (-β * S * I) + (-ν_now * S) # contact matrix cannot be multiplied here
    @. dE = β * S * I - α * E
    @. dI = α * E - (γ * I_) # note use of I_
    @. dR = γ * I_
    @. dV = ν_now * S
end
