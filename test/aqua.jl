using Aqua
using StatsBase
using Test

@testset "Aqua" begin
    Aqua.test_all(
        StatsBase;
        # Ambiguity checks error due to a bug in Aqua on older Julia versions:
        # https://github.com/JuliaTesting/Aqua.jl/issues/141
        # We also ignore an ambiguity introduced by SortingAlgorithms
        ambiguities=VERSION >= v"1.6" ? (exclude = [sort!],) : false,
        # We should ignore `describe`, `pairwise`, `pairwise!` etc. which are owned by StatsBase once supported by Aqua:
        # https://github.com/JuliaTesting/Aqua.jl/pull/140
        piracy=(broken = true,),
    ) 
end
