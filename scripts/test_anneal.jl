using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))

function iterate(ds, x, ε, max_it)
    d = length(x) 
    set_state!(ds, d == 1 ? x[1] : x) 
    n = @timed get_iterations!(ds, ε, max_it)
    return n.value[1], n.time
end


ε = 1e-25;  max_it = 100; 
# ε = 1e-8;  max_it = 100; 
setprecision(BigFloat, 50; base = 10)


# ds = setup_iterator(F_list[19], g_list[1], X0[19]; algtype = :Steffensen)
# @show n, xf, q = compute_figure(ds, ε, max_it)

g(x) = g_list[1](x,ε)
ds = setup_iterator(F_list[1], g, X0[1]; algtype = :normal)
@show n, xf = iterate(ds, X0[1], ε, max_it)

g(x) = g_list[1](x,ε)
ds = setup_iterator(F_list[7], g, big(X0[7]); algtype = :normal)
@show n, xf = iterate(ds, big(rand()), ε, max_it)

# ds = setup_iterator(F_list[26], g_list[1], X0[26]; algtype = :Steffensen)
# @show n, xf, q = compute_figure(ds, ε, max_it)

# ds = setup_iterator(F_list[19], g_list[1], X0[19]; algtype = :accelerated)
# @show n, xf, q = compute_figure(ds, ε, max_it)

# ds = setup_iterator(F_list[1], g_list[1], X0[1]; algtype = :accelerated)
# @show n, xf, q = compute_figure(ds, ε, max_it)

# ds = setup_iterator(F_list[1], g_list[1], X0[1]; algtype = :accelerated_secant)
# @show n, xf, q = compute_figure(ds, ε, max_it)

g(x) = g_list[1](x,ε)

ds = setup_iterator(optim_f[2], g, X0_optim[1]; algtype = :accelerated)
# ds = setup_iterator(F_list[23], g, X0[23]; algtype = :Steffensen)
@show n, t = iterate(ds, X0_optim[1], 1e-8, 10000)
