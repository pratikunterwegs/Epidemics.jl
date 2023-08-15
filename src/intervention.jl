"""
  Npi(time_begin, time_end, contact_reduction)

A structure to hold the end points and strength of a non-pharmaceutical
intervention.
    
"""
mutable struct Npi
  time_begin::Vector
  time_end::Vector
  contact_reduction::Matrix
end

function Npi(; time_begin::Number = 50, time_end::Number = 80,
             contact_reduction::Vector = [0.25])
    # convert contact reduction to matrix
    cr_matrix = reshape(contact_reduction, length(contact_reduction), 1)
    return Npi(time_begin, time_end, cr_matrix)
end

# function for a cumulative intervention
function cumulative_npi(;t, npi::Npi)
  return npi
end

export Npi
