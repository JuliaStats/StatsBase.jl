using Test
using StatsBase
using Aqua

@testset "Aqua tests (performance)" begin
    # This tests that we don't accidentally run into
    # https://github.com/JuliaLang/julia/issues/29393
    # Aqua.test_unbound_args(StatsBase)
    ua = Aqua.detect_unbound_args_recursively(StatsBase)
    @test length(ua) == 0

    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
    # Test that we're not introducing method ambiguities across deps
    ambs = Aqua.detect_ambiguities(StatsBase; recursive = true)
    # Uncomment for debugging:
    # for method_ambiguity in ambs
    #     @show method_ambiguity
    # end
    @test length(ambs) ≤ 5
    potentially_pirated_methods = Aqua.Piracy.hunt(StatsBase)
    # Uncomment for debugging:
    # for pirate in potentially_pirated_methods
    #     @show pirate
    # end
    @test length(potentially_pirated_methods) ≤ 27
end

@testset "Aqua tests (additional)" begin
    Aqua.test_undefined_exports(StatsBase)
    Aqua.test_stale_deps(StatsBase)
    Aqua.test_deps_compat(StatsBase)
    Aqua.test_project_extras(StatsBase)
    Aqua.test_project_toml_formatting(StatsBase)
end

nothing
