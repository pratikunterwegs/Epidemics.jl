
using BenchmarkTools
using Epidemics

time_end = 1000.0

# benchmark the default model with 3 age groups, intervention, and vaccination
@benchmark Epidemics.epidemic_default(time_end=time_end, increment=1.0)
