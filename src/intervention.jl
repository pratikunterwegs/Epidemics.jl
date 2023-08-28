"""
  Npi(time_begin, time_end, contact_reduction)

A structure to hold the end points and strength of a non-pharmaceutical
intervention.
    
"""
mutable struct Npi
    time_begin::Vector
    time_end::Vector
    contact_reduction::Matrix

    ## internal constructor
    function Npi(time_begin, time_end, contact_reduction)
        # basic input checking
        @assert length(time_begin)==length(time_end) "`time_begin` and `time_end` must be the same length"
        if any(time_begin .< 0.0)
            error("`time_begin` must be > 0.0")
        end
        if any(time_end .< 0.0)
            error("`time_end` must be > 0.0")
        end
        if any(time_end .< time_begin)
            error("All `time_end` must be > corresponding `time_begin`")
        end
        if any((contact_reduction .> 1.0) .| (contact_reduction .< 0.0))
            error("`contact_reduction` values must be in the range [0.0, 1.0]")
        end
        new(time_begin, time_end, contact_reduction)
    end
end

function Npi(; time_begin::Number = 50, time_end::Number = 80,
             contact_reduction::Vector = [0.25])
    # convert contact reduction to a matrix with one column
    contact_reduction = reshape(contact_reduction, length(contact_reduction), 1)
    return Npi([time_begin], [time_end], contact_reduction)
end

function no_intervention()::Npi
    return Npi(time_begin = 0.0, time_end = 0.0, contact_reduction = [0.0])
end

# Define a constructor method to combine Npis
function c(args::Npi...)
    # collect Npi elements
    time_begin = vcat([a.time_begin for a in args]...)
    time_end = vcat([a.time_end for a in args]...)

    contact_reduction = [a.contact_reduction for a in args]
    contact_reduction = hcat(contact_reduction...)

    return Npi(time_begin, time_end, contact_reduction)
end

# function for a cumulative intervention
function cumulative_npi(; t, npi::Npi)
    # which Npis are active
    active = ((t .>= npi.time_begin) .& (t .<= npi.time_end))

    cr = fill(0, size(npi.contact_reduction)[1])
    for i in eachindex(active)
        if active[i]
            cr = cr .+ npi.contact_reduction[:, i]
        end
    end

    # handle values greater than 1.0
    cr[cr .> 1.0] .= 1.0

    return cr
end

export Npi, c, cumulative_npi, no_intervention
