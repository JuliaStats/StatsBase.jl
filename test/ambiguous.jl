using StatsBase
using Test
if VERSION >= v"0.7.0-beta.85"
    using Statistics
end

@testset "Ambiguities" begin
    # Ambiguities with Base and Core introduced by this package
    tocheck = [StatsBase, Base, Core]
    VERSION >= v"0.7.0-beta.85" && push!(tocheck, Statistics)
    @test_broken isempty(detect_ambiguities(tocheck...))
end
