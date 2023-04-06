"""
    seir!(du, u, parameters, t)

A simple SEIR epidemic model function that allows for multiple demographic
groups. This function is intended to be called internally from
[`epidemic`](@ref). The function expects the parameters argument to be a four 
element vector with the following elements: a vector of ``\\beta``, the 
transmission rate, where each element of the vector represents the transmission
rate of the pathogen within a specific demographic group; a vector of 
``\\alpha``, the group-specific rate of conversion from exposed  to infectious; 
a vector of ``\\gamma``, the group-specific rate of recovery; a matrix 
specifying the contacts between demographic groups; a matrix of the intervention
applied to each age group, see [`Npi`](@ref).
    
"""
function seir!(du, u, parameters, t)
    # assumes that each element of the vector is a 
    # vector of rates
    β, α, γ, contact_matrix, intervention = parameters

    # modify contact matrix if the intervention is active
    if (t > intervention.time_begin) & (t < intervention.time_end)
        contact_matrix = contact_matrix .*
                         (1.0 .- intervention.contact_reduction)
    end

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, 1]
    E = @view u[:, 2]
    I = contact_matrix * @view u[:, 3] # matrix mult for cm * I
    I_ = @view u[:, 3] # unmultiplied I for operations involving only I
    R = @view u[:, 4]

    # views to the change matrix, dU
    dS = @view du[:, 1]
    dE = @view du[:, 2]
    dI = @view du[:, 3]
    dR = @view du[:, 4]

    # calculate change in compartment size and update the change matrix dU
    # note the use of @. for broadcasting, equivalent to .=
    @. dS = -β * S * I # contact matrix cannot be multiplied here due to @.?
    @. dE = β * S * I - α * E
    @. dI = α * E - (γ * I_) # note use of I_
    @. dR = γ * I_
end
