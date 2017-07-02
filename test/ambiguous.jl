using StatsBase, Base.Test

@testset "Ambiguities" begin
    # Ambiguities with Base and Core introduced by this package
    @test_broken isempty(detect_ambiguities(StatsBase, Base, Core))
end
