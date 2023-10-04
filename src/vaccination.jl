"""
  Vaccination(time_begin, time_end, ν)

A structure to hold the end points and strength of a non-pharmaceutical
intervention.

"""
mutable struct Vaccination
    time_begin::Vector
    time_end::Vector
    ν::Vector

    # internal constructor
    function Vaccination(; time_begin::Vector = [50],
        time_end::Vector = [80],
        ν::Vector = [0.25])

        # check for inputs
        @assert (length(time_begin) == length(time_end))&&(length(time_begin) == length(ν)) "All arguments to Vaccination() must be Number Vectors of the same length"

        new(time_begin, time_end, ν)
    end
end

function no_vaccination(; doses::Number = 1)
    return Vaccination(time_begin = [0], time_end = [0], ν = [0])
end

# helper function for Vaccination object
function current_nu(; time::Number, vaccination::Vaccination)
    # broadcast bitwise & as broadcast logical operator gives bitwise vector
    current_nu = vaccination.ν .*
                 ((time .> vaccination.time_begin) .& (time .< vaccination.time_end))
    return current_nu
end

export Vaccination, current_nu, no_vaccination
