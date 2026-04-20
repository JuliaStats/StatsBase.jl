using StatsBase, Test
import JET

@testset "JET" begin
    # Check that there are no undefined global references and undefined field accesses
    JET.test_package(StatsBase; target_modules = (StatsBase,), mode = :typo)

    # Default error analysis
    # Note: This analysis is not enough strict to guarantee that there are no runtime errors!
    kwargs = if isdefined(JET, :LastFrameModuleExact) # JET@0.11
        (; target_modules = (StatsBase,))
    else
        (; target_defined_modules = true)
    end
    # The (deprecated) `model_response(::StatisticalModel)` calls the interface
    # function `response(::StatisticalModel)` for which no method exists yet
    # Ref https://github.com/aviatesk/JET.jl/issues/495
    res = JET.report_package(StatsBase; kwargs..., mode = :basic, toplevel_logger = nothing)
    println(res)
    reports = JET.get_reports(res)
    @test_broken isempty(reports)
    @test length(reports) <= 1
end
