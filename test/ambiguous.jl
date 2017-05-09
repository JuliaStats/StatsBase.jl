using Base.Test

@testset "Ambiguities" begin
    # Preexisting ambiguities between Base and Core -- we want to ignore these
    base_core_ambig = detect_ambiguities(Base, Core)

    # Ambiguities introduced by this package
    our_ambig = detect_ambiguities(StatsBase, Base, Core)

    @test isempty(setdiff(our_ambig, base_core_ambig))
end
