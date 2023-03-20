using Epidemics
using Test

using DataFrames

@testset "Epidemics.jl" begin
    x = 2
    y = 2
    @test Epidemics.sum_values(x,y) == 4
end

@testset "SEIR model return type" begin
    time_end = 100.0
    n_age_groups = 3.0
    n_compartments = 4.0
    
    data = Epidemics.epidemic(time_end = time_end, increment = 1.0)
    @test typeof(data) == DataFrames.DataFrame
    @test size(data, 2) == 4 # test for four cols
    # count initial rows for t = 0.0
    initial_rows = n_age_groups * n_compartments
    # test for expected number of columns
    @test size(data, 1) == Int((time_end * n_age_groups * n_compartments) + (initial_rows))

    # test that final values sum to the same as initial values
    initial_pop = sum(filter(:timestamp => x -> x == 0.0, data).value)
    final_pop = sum(filter(:timestamp => x -> x == time_end, data).value)
    @test initial_pop ≈ final_pop atol=1e-6

end