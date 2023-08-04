"""
  Npi(time_begin, time_end, contact_reduction)

A structure to hold the end points and strength of a non-pharmaceutical
intervention.
    
"""
mutable struct Npi
  time_begin::Number
  time_end::Number
  contact_reduction::Vector
end

function Npi(; time_begin::Number=50, time_end::Number=80,
  contact_reduction::Vector=[0.25])
  
  return Npi(time_begin, time_end, contact_reduction)
end
