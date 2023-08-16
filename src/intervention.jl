"""
  Npi(time_begin, time_end, contact_reduction)

A structure to hold the end points and strength of a non-pharmaceutical
intervention.
    
"""
mutable struct Npi
    time_begin::Vector
    time_end::Vector
    contact_reduction::Vector
end

function Npi(; time_begin::Number = 50, time_end::Number = 80,
             contact_reduction::Vector = [0.25])
    # convert contact reduction to matrix
    return Npi([time_begin], [time_end], contact_reduction)
end

function Npi(; time_begin::Vector = [50], time_end::Vector = [80],
             contact_reduction::Vector{Vector})
    # convert contact reduction to matrix
    return Npi(time_begin, time_end, contact_reduction)
end

# Define a constructor method to combine Npis
function c(args::Npi...)
    time_begin = vcat([a.time_begin for a in args]...)
    time_end = vcat([a.time_end for a in args]...)
    contact_reduction = [a.contact_reduction for a in args]
    return Npi(time_begin, time_end, contact_reduction)
end

# function for a cumulative intervention
function cumulative_npi(; t, npi::Npi)
  # which Npis are active
  active = ((t .> npi.time_begin) .& (t .< npi.time_end))
  cr = npi.contact_reduction[active]
  cr = reduce((x, y) -> (1 .- x) .* (1 .- y), cr)
  return cr
end

export Npi, c, cumulative_npi
