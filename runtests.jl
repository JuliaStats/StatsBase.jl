using Stats

tests = ["01", "means", "scalarstats", "intstats"]

println("Running tests:")

for t in tests
	tfile = joinpath("test", "$(t).jl")
	println(" * $(tfile) ...")
	include(tfile)
end

