"""
  npi(time_begin, time_end, contact_reduction)

A structure to hold the end points and strength of a non-pharmaceutical
intervention.
    
"""
abstract type intervention end
mutable struct npi <: intervention
  time_begin::Number
  time_end::Number
  contact_reduction::Vector{Number}
end
