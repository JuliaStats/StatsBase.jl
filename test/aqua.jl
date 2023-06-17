using Aqua
using StatsBase
using Test

@testset "Aqua" begin
    Aqua.test_all(
        StatsBase;
        ambiguities=(exclude = [sort!],), # ambiguity in SortingAlgorithms (ignore)
        piracy=(broken = true,), # should ignore `describe`, `pairwise`, `pairwise!` etc. which are owned by StatsBase once supported by Aqua
    ) 
end
