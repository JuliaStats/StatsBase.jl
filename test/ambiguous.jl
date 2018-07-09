using StatsBase
using Test
using Statistics

@testset "Ambiguities" begin
    # Ambiguities with Base, Core, and stdlib Statistics introduced by this package
    @test_broken isempty(detect_ambiguities(StatsBase, Base, Core, Statistics))
end
