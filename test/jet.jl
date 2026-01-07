using StatsBase, Test
using JET: JET

# JET has only experimental support for Julia 1.12 currently
# It throws an internal `AssertionError` in the tests below
if VERSION < v"1.12-"
    @testset "JET" begin
        # Check that there are no undefined global references and undefined field accesses
        JET.test_package("StatsBase"; target_defined_modules=true, mode=:typo)

        # Default error analysis for common problem fails since JET errors on interface definitions
        # The (deprecated) `model_response(::StatisticalModel)` calls the interface
        # function `response(::StatisticalModel)` for which no method exists yet
        # Note: This analysis is not enough strict to guarantee that there are no runtime errors!
        # Ref https://github.com/aviatesk/JET.jl/issues/495
        res = JET.report_package("StatsBase"; target_defined_modules=true, mode=:basic,
                                 toplevel_logger=nothing)
        println(res)
        reports = JET.get_reports(res)
        @test_broken isempty(reports)
        @test length(reports) <= 1
    end
end
