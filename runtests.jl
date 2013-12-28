using Stats

tests = ["means", "scalarstats", "counts", "ranking", "corr"]

println("Running tests:")

for t in tests
	tfile = joinpath("test", "$(t).jl")
	println(" * $(tfile) ...")
	include(tfile)
end

