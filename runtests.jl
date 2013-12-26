using Stats

tests = ["warrays", "means", "scalarstats", "counts"]

println("Running tests:")

for t in tests
	tfile = joinpath("test", "$(t).jl")
	println(" * $(tfile) ...")
	include(tfile)
end

