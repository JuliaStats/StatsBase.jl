using Aqua
using StatsBase
using Test

@testset "Aqua" begin
    Aqua.test_all(
        StatsBase;
        # Ignore an ambiguity introduced by SortingAlgorithms
        # Disable test on Julia < 1.6 due to ambiguities in ChainRulesCore
        # (fixed in more recent versions not supported by Julia < 1.6)
        ambiguities=VERSION < v"1.6" ? false : (exclude = [sort!],),
        # `describe`, `pairwise`, `pairwise!`, and `weights` are defined in StatsAPI but owned by StasBase
        # `quantile(::AbstractVector{<:Real})` is a deprecated method that should be removed
        piracy=(treat_as_own = [describe, pairwise, pairwise!, weights, quantile],),
    ) 
end
