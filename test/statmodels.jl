using StatsBase
using Base.Test

## format_pvc: Formatting of p-values
@test StatsBase.format_pvc(1.0) == "1.0000"
@test StatsBase.format_pvc(1e-1) == "0.1000"
@test StatsBase.format_pvc(1e-5) == "<1e-4"
@test StatsBase.format_pvc(NaN) == "NaN"
@test_throws ErrorException StatsBase.format_pvc(-0.1)
@test_throws ErrorException StatsBase.format_pvc(1.1)
