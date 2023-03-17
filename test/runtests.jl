using Epidemics
using Test

using DataFrames

@testset "Epidemics.jl" begin
    x = 2
    y = 2
    @test Epidemics.sum_values(x,y) == 4
end

@testset "SEIR model" begin
    data = Epidemics.epidemic()
    @test typeof(data) == DataFrames.DataFrame
end