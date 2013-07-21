#
# Correctness Tests
#

# require("test.jl")

using Stats

my_tests = ["test/01.jl"]
            #"test/statquiz.jl"]

println("Running tests:")

for my_test in my_tests
    println(" * $(my_test)")
    include(my_test)
end
