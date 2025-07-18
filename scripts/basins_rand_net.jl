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
    if  1 < length(bas_num) < 100
        fig = plot_heatmap(grid, basins, iterations, roots; ukeys = bas_num, shaded, show_attractors, xticksvisible = false, yticksvisible = false, xticklabelsvisible = false, yticklabelsvisible = false)
        s = plotsdir(savename(string(prefix), @dict(res,ε),"png"))
        save(s, fig)
    end
end


function plot_basins(d)
    @unpack res, ε, N, g_ind, max_it = d
    xg = yg = range(-5, 5; length = res)
    grid = (xg, yg)
    ε = 1e-8
    max_it = 150
    s = 4.5; gg = 2.5; N = 8
    rng = MersenneTwister(123);

    v = rand(rng, N)
    P1 = (v .> 0.5)
    # v = rand(rng, N)
    P2 = 1 .- P1
    P1 = P1 ./ norm(P1)
    P2 = P2 ./ norm(P2)

    # P1 = [ 1; zeros(N-1)]
    # P2 = [ 0; 1; zeros(N-2)]

    W = randn(rng,N,N)
    for k in 1:N; W[k,k] = 0.0; end;
    # F(u) = -u .+ s*tanh.(u) .+ gg/sqrt(N)*W*tanh.(u) 
    F = Vector{Function}(undef,N)
    for n in 1:N
        F[n] = u -> (-u[n] .+ s*tanh.(u[n]) .+ gg/sqrt(N)*W[n,:]'*tanh.(u))
    end


    g(x) = g_list[g_ind](x,ε/2)
    alg = :accelerated
    ds_it = setup_iterator(F, g, rand(N); algtype = alg)

    x1, _ =  get_state(ds_it)
    roots =  typeof(x1)[] 
    basins = zeros(Int32,res,res); iterations = zeros(Int16,res,res)

    @showprogress for (i,x) in enumerate(xg), (j,y) in enumerate(yg) 
        set_state!(ds_it, x*P1 + y*P2)
        n = @timed get_iterations!(ds_it, ε, max_it)
        it = n.value[1]
        if it ≥ max_it
            # the alg. did not converge
            basins[i,j] = -1
        else
            # We identify the root with the mapper.
            xf, _ = get_state(ds_it)
            basins[i,j] = custom_mapper(xf, roots, 0.01)
        end
        iterations[i,j] = it
    end

    return @strdict(grid, basins, iterations, roots)
end


# Plot all basins 
res = 1000
ε = 1e-8
max_it = 100 
N = 12
force = false

for g_ind in 1:3
    d = @dict(res, ε, max_it, N, g_ind) # parametros
    data, file = produce_or_load(
        datadir(""), # path
        d, # container for parameter
        plot_basins, # function
        prefix = "basins_rand_net", # prefix for savename
        force = force, # true for forcing sims
        wsave_kwargs = (;compress = true)
    )
    @unpack basins, grid, roots = data
    xg, yg = grid
    s = plotsdir(savename("fig_rand_net",d,"png"))
    f = Figure(size = (400,400))
    ax = Axis(f[1,1], xticksvisible = false, yticksvisible = false, xticklabelsvisible = false, yticklabelsvisible = false) #0xlabel = L"x_1", ylabel = L"x_2") #, yscale = log10);

    using Colors, ColorSchemes
    # cmap = ColorScheme([RGB(1,1,1), RGB(0,1,0), RGB(0.1,0.1,0.1), RGB(1,0.46, 0.46), RGB(0.34,0.34,1)] )
    # cmap = ColorScheme([RGB(1,0,0), RGB(0,0,0), RGB(0.2,0.0,0.0), RGB(1,1,1)])
    cmap = :flag
    heatmap!(ax, xg, yg, basins; colormap = cmap)
    save(s,f)
end
