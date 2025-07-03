using DrWatson
@quickactivate
using CairoMakie
using CodecZlib
using LaTeXStrings
using Statistics:mean
include("../src/function_stuff.jl")
include("../src/function_list.jl")
include("../src/color_stuff.jl")
include("../src/basins_compute.jl")

function plot_basins(ds_it, grid, res, ε = 1e-8, max_it =50; force = false, shaded = true, show_attractors = false, prefix = "stephenson")
    data = _get_basins(ds_it, grid, res, ε, max_it; force, prefix)
    @unpack basins, iterations, roots, grid = data
    @show bas_num = unique(basins)
    if  1 < length(bas_num) < 100
        fig = plot_heatmap(grid, basins, iterations, roots; ukeys = bas_num, shaded, show_attractors, xticksvisible = false, yticksvisible = false, xticklabelsvisible = false, yticklabelsvisible = false)
        s = plotsdir(savename(string(prefix), @dict(res,ε),"png"))
        save(s, fig)
    end
end


# Plot all basins 
res = 200
xg = yg = range(-10, 10; length = res)
grid = (xg, yg)
force = true
ε = 1e-8
max_it = 100
for i in [1, 6, 9, 11]
    for (k,g) in enumerate(g_list)
        F = [ x -> real(F_list[i](x[1]+im*x[2])),  x -> imag(F_list[i](x[1]+im*x[2]))]
        gg(x) = g(x,ε/2)
        alg = :accelerated
        ds_it = setup_iterator(F, gg, rand(2); algtype = alg)
        plot_basins(ds_it, grid, res, ε, max_it;  prefix = string("steffenson_", alg, "_f",i, "_g", k), force , shaded = true)
    end
    for (k,g) in enumerate(g_list)
        F = [ x -> real(F_list[i](x[1]+im*x[2])),  x -> imag(F_list[i](x[1]+im*x[2]))]
        gg(x) = g(x,ε/2)
        alg = :normal
        ds_it = setup_iterator(F, gg, rand(2); algtype = alg)
        plot_basins(ds_it, grid, res, ε, max_it;  prefix = string("steffenson_", alg, "_f",i, "_g", k), force , shaded = true)
    end
end


for i in [15, 18] 
    for (k,g) in enumerate(g_list)
        gg(x) = g(x,ε/2)
        alg = :accelerated
        ds_it = setup_iterator(F_list[i], gg, rand(2); algtype = alg)
        plot_basins(ds_it, grid, res, ε, max_it;  prefix = string("steffenson_", alg, "_f",i, "_g", k), force , shaded = true)
    end
    for (k,g) in enumerate(g_list)
        gg(x) = g(x,ε/2)
        alg = :normal
        ds_it = setup_iterator(F_list[i], gg, rand(2); algtype = alg)
        plot_basins(ds_it, grid, res, ε, max_it;  prefix = string("steffenson_", alg, "_f",i, "_g", k), force , shaded = true)
    end
end
