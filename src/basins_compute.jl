using LinearAlgebra:norm
using ProgressMeter
using Random

include(srcdir("function_stuff.jl"))


""" 
    function _get_basins(N,i,res,ε,max_it) -> data

Convenience function to compute and store the basins
and attractors of the funcion i. with the proximity algorithm

"""
function _get_basins(ds_it, grid, res, ε, max_it; prefix = "basins_", force = false)
    d = @dict(ds_it, grid, res, ε, max_it) # parametros
    data, file = produce_or_load(
        datadir(""), # path
        d, # container for parameter
        compute_basins, # function
        prefix = prefix, # prefix for savename
        force = force, # true for forcing sims
        wsave_kwargs = (;compress = true)
    )
    return data
end


function get_stats(ds, Nsamples, grid, ε, max_it; seed = 123, prefix = "stats_", force = false)
    d = @dict(ds, Nsamples, ε, max_it, grid, seed) # parametros
    data, file = produce_or_load(
        datadir(""), # path
        d, # container for parameter
        compute_stats, # function
        prefix = prefix, # prefix for savename
        force = force, # true for forcing sims
        wsave_kwargs = (;compress = true)
    )
    return data
end



"""
    compute_basins(d) -> Dict

Compute the basins using first AttractorsViaRecurrences to 
locate the attractors (roots) of the function N_β. These 
attractors are passed to AttractorsViaProximity with a 
convergence criterion ε such that the algorithm stops 
when |f(x) - r| < ε. 
The basins, the iteration matrix, the metrics and the attractors
are returned into a name dictionnary. 
"""
function compute_basins(d)
    @unpack ds_it, grid, res, ε, max_it = d
    x1, _ =  get_state(ds_it)
    roots =  typeof(x1)[] 
    xg = yg = range(-10, 10; length = res)
    grid = (xg, yg)

    basins = zeros(Int32,res,res); iterations = zeros(Int16,res,res)
    exec_time = zeros(res,res)

@showprogress for (i,x) in enumerate(xg), (j,y) in enumerate(yg) 
        set_state!(ds_it, [x,y])
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
        exec_time[i,j] = n.time
    end

    return @strdict(grid, basins, iterations, exec_time, roots)
end


function custom_mapper(xf, roots, ε)

    if isempty(roots) 
        push!(roots, xf) 
        return 1
    end

    for (k,r) in enumerate(roots)
        if norm(r .- xf) < ε
            # @show norm(r .- xf) 
            return k
        end
    end
    push!(roots, xf) 
    return length(roots) 
end

"""
Compute stats!
"""
function compute_stats(d)
    @unpack ds, seed, Nsamples, ε, max_it, grid = d
    dim = length(grid)
    iterations = 0.0
    exec_time = 0.0
    nc = 0

    samp = sampler(grid, seed)

    for k in 1:Nsamples
        set_state!(ds, samp())
        # set_state!(ds, big(samp()))
        n = @timed get_iterations!(ds, ε, max_it)
        it = n.value[1]
        if it ≥ max_it
            # the alg. did not converge
            nc += 1
        else
            iterations += it
            exec_time += n.time
        end
    end
    exec_time = exec_time/iterations
    iterations = iterations/(Nsamples - nc)
    nc = nc/Nsamples
    return @strdict(grid, Nsamples, iterations, exec_time, nc)
end


# Estimate order
function  estimate_ACOC!(T, yy)
    qn = 0.
    for k in 3:T-2
        num = log(norm(yy[k+1] - yy[k])) - log(norm(yy[k] - yy[k-1])) 
        den = log(norm(yy[k] - yy[k-1]))- log(norm(yy[k-1] - yy[k-2]))
        qn = num/den
    end
    return qn
end



function sampler(grid, seed::Int)::Function

  rng = MersenneTwister(seed)  

  if length(grid) == 1
      generate_point = () -> rand(rng) * (grid[1][end] - grid[1][1]) + grid[1][1]
  else
      generate_point = () -> [rand(rng) * (grid[i][end] - grid[i][1]) + grid[i][1] for i in 1:length(grid)]
  end

  return generate_point
end
