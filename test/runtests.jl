using Epidemics
using Test

@testset "Epidemics.jl" begin
    x = 2
    y = 2
    @test Epidemics.sum_values(x,y) == 4
end
