module Epidemics

export epidemic_default  # This line exports the function

# Write your package code here
include("model_epidemic_default.jl")
include("model_vacamole.jl")
include("model_stochastic.jl")
include("model_norovirus.jl")
include("tools.jl")
include("prepare_data.jl")

end
