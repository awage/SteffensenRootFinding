using Profile
using ProfileView
using BenchmarkTools
using DrWatson
@quickactivate
# using CodecZlib
# using LaTeXStrings
# using Statistics
# using Profile
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))

# Profile.init()
ε = 1.e-14;  max_it = 50; 
setprecision(BigFloat, 50; base = 10)

# Iterations
x0 = big(rand())
N = stephenson_map_real(F_list[1], fam_list_real[2])
ds = FunIterator(N, [x0])

function funfun(k)
    for _ in 1:k
        x0 = big(4*(rand()-0.5))
        set_state!(ds, [x0])
        n = _get_iterations!(ds, ε, max_it)
        # @show n
    end
end


funfun(1)

# @profview  funfun(10000)
@btime  funfun(10000)
# funfun(10000)

# Profile.print()
