"""
  Vaccination(time_begin, time_end, ν)

A structure to hold the end points and strength of a non-pharmaceutical
intervention.

"""
mutable struct Vaccination
    time_begin::Vector
    time_end::Vector
    ν::Vector
end

function Vaccination(; time_begin::Vector = [50],
                     time_end::Vector = [80],
                     ν::Vector = [0.25])

    # check for inputs
    @assert (length(time_begin) == length(time_end))&&(length(time_begin) == length(ν)) "All arguments to Vaccination() must be Number Vectors of the same length"

    return Vaccination(time_begin, time_end, ν)
end

# helper function for Vaccination object
function current_nu(; time::Number, vaccination::Vaccination)
    current_nu = vaccination.ν .*
                 ((time .> vaccination.time_begin) .& (time .< vaccination.time_end)) # broadcast bitwise & as broadcast logical operator gives bitwise vector
    return current_nu
end
