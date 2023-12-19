using Epidemics
using Test

using DataFrames

@testset "Default model return type" begin
    time_end = 100.0
    n_age_groups = 3.0
    n_compartments = 5 #SEIRV

    # create a dummy population and infection
    population = Population()

    # generate parameters using helpers
    r0 = 1.3
    infectious_period = 5
    preinfectious_period = 2

    β = r0_to_beta(r0 = r0, infectious_period = infectious_period)
    σ = preinfectious_period_to_alpha(preinfectious_period = preinfectious_period)
    γ = infectious_period_to_gamma(infectious_period = infectious_period)

    # run the model
    data = epidemic_default(β = β, σ = σ, γ = γ,
        population = population,
        time_end = time_end,
        increment = 1.0)

    # convert to data.frame and apply `prepare_data`
    data = prepare_data(ode_solution_df = DataFrame(data), n_age_groups = 3)

    @test typeof(data) == DataFrames.DataFrame
    @test size(data, 2) == 4 # test for four cols
    # count initial rows for t = 0.0
    initial_rows = n_age_groups * n_compartments
    # test for expected number of columns
    @test size(data, 1) == Int((time_end * n_age_groups * n_compartments) + (initial_rows))

    # test that final values sum to the same as initial values
    initial_pop = sum(filter(:timestamp => x -> x == 0.0, data).value)
    final_pop = sum(filter(:timestamp => x -> x == time_end, data).value)
    @test initial_pop≈final_pop atol=1e-6
end

# tests for helpers relating to initial conditions
@testset "Helper functions for initial conditions" begin

    # test for default initial conditions
    # returns a matrix
    @test isa(default_initial_conditions(), AbstractMatrix{<:Number})
    # breaks when inputs are bad
    # p_infected outside limits
    @test_throws AssertionError default_initial_conditions(p_infected = [1.1])
    @test_throws AssertionError default_initial_conditions(p_infected = [-1e-6])
    # p_exposed outside limits
    @test_throws AssertionError default_initial_conditions(p_exposed = [1.1])
    @test_throws AssertionError default_initial_conditions(p_exposed = [-1e-6])
    # p_infected and p_exposed outside limits
    @test_throws AssertionError default_initial_conditions(p_infected = [0.5],
        p_exposed = [0.6])
end

# tests for helpers to convert parameters
@testset "Helper functions for parameter conversion" begin

    # tests for default parameter conversion
    # functions throw assertion errors
    # r0_to_beta
    @test_throws AssertionError r0_to_beta(r0 = 0.0, infectious_period = 5)
    @test_throws AssertionError r0_to_beta(r0 = 1.3, infectious_period = -5)
    # preinfectious period to alpha
    @test_throws AssertionError preinfectious_period_to_alpha(preinfectious_period = 0.0)
    # infectious period to gamma
    @test_throws AssertionError infectious_period_to_gamma(infectious_period = 0)
end

# tests for the intervention class
@testset "Intervention class" begin
    @test isa(no_intervention(), Npi)
    @test_throws ErrorException Npi(time_begin = 1, time_end = 0)
    @test_throws ErrorException Npi(time_begin = 1, time_end = 0)
    @test_throws ErrorException Npi(contact_reduction = [1.1, 0.1])
end

@testset "Vacamole model return type" begin
    time_end = 100.0
    n_age_groups = 3.0
    n_compartments = 5 #SEIRV

    # create a dummy population and infection
    population = Population(demography_vector = 10e6 .* [0.2, 0.5, 0.3],
        initial_conditions = [1-1e-6 0.0 0.0 0.0 0.0 1e-6 0.0 0.0 0.0 0.0 0.0;
            1-1e-6 0.0 0.0 0.0 0.0 1e-6 0.0 0.0 0.0 0.0 0.0;
            1-1e-6 0.0 0.0 0.0 0.0 1e-6 0.0 0.0 0.0 0.0 0.0],
        contact_matrix = ones(3, 3) * 5)

    # generate parameters using helpers
    r0 = 1.3
    infectious_period = 5
    preinfectious_period = 2

    β = r0_to_beta(r0 = 1.3, infectious_period = infectious_period)
    σ = preinfectious_period_to_alpha(preinfectious_period = preinfectious_period)
    γ = infectious_period_to_gamma(infectious_period = infectious_period)

    # run the model
    data = epidemic_vacamole(β = β, σ = σ, γ = γ, η = 1 / 100, ω = 1 / 1000,
        population = population,
        time_end = time_end,
        increment = 1.0)

    # convert to data.frame and apply `prepare_data`
    data = prepare_data(ode_solution_df = DataFrame(data), n_age_groups = 3)

    @test typeof(data) == DataFrames.DataFrame
    @test size(data, 2) == 4 # test for four cols
    # count initial rows for t = 0.0
    initial_rows = n_age_groups * n_compartments
    # test for expected number of columns
    @test size(data, 1) == Int((time_end * n_age_groups * n_compartments) + (initial_rows))

    # test that final values sum to the same as initial values
    initial_pop = sum(filter(:timestamp => x -> x == 0.0, data).value)
    final_pop = sum(filter(:timestamp => x -> x == time_end, data).value)
    @test initial_pop≈final_pop atol=1e-6
end