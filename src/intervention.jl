"""
  Npi(time_begin, time_end, contact_reduction)

A structure to hold the end points and strength of a non-pharmaceutical
intervention.
    
"""
mutable struct Npi
  time_begin::Number
  time_end::Number
  contact_reduction::Vector{Number}
end

function Npi(; time_begin=50, time_end=80, contact_reduction=0.25)
  # convert contact reduction to a vector if a single number
  if length(contact_reduction) == 1
    contact_reduction = [contact_reduction]
  end

  return Npi(time_begin, time_end, contact_reduction)
end
