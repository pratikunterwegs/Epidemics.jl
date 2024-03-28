module Epidemics

export epidemic_default  # This line exports the function

# Write your package code here
include("model_epidemic_default.jl")
include("model_vacamole.jl")
include("model_stochastic.jl")
include("tools.jl")

end
