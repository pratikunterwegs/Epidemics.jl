"""
  Vaccination(time_begin, time_end, ν)

A structure to hold the end points and strength of a non-pharmaceutical
intervention.

"""
mutable struct Vaccination
    time_begin::Vector{Number}
    time_end::Vector{Number}
    ν::Vector{Number}
end

function Vaccination(; time_begin=50, time_end=80, contact_reduction=0.25)
    # convert parameters to vectors if a single number
    if length(time_begin) == 1
        time_begin = [time_begin]
    end

    if length(time_end) == 1
        time_end = [time_end]
    end

    if length(contact_reduction) == 1
        contact_reduction = [contact_reduction]
    end

    return Vaccination(time_begin, time_end, contact_reduction)
end

# helper function for Vaccination object
function current_nu(; time::Number, vaccination::Vaccination)
    current_nu = vaccination.ν .* (
        (time .> vaccination.time_begin) .& (time .< vaccination.time_end)
    ) # broadcast bitwise & as broadcast logical operator gives bitwise vector
    return current_nu
end
