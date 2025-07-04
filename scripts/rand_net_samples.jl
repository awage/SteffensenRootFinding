using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))

using Printf
macro display_float_2_digits(arr)
    quote
    _display_float_2_digits_impl($(esc(arr)))
    end
end

function _display_float_2_digits_impl(arr::AbstractArray{<:AbstractFloat})
    formatted_elements = [@sprintf("%.2f", x) for x in arr]
    println("[", join(formatted_elements, ", "), "]")
end

    
function compute_figure(ds, ε, max_it)
    # xf, fx = get_state(ds) 
    n, yy = get_iterations!(ds, ε, max_it)
    xf, fx = get_state(ds) 
    # @show Float64.(xf), Float64.(fx)
    # @display_float_2_digits xf
    # @display_float_2_digits fx

    if 5 ≤ n < max_it
        q = estimate_ACOC!(n, yy)
    else
        q = 0
    end
    return n, xf, q
end

function get_roots_number(N, Nsamples, max_it, rng) 
    ε = 1e-8;  
    s = 4.5; gg = 2.5; 
    W = randn(rng,N,N)
    for k in 1:N; W[k,k] = 0.0; end;
    # F(u) = -u .+ s*tanh.(u) .+ gg/sqrt(N)*W*tanh.(u) 
    F = Vector{Function}(undef,N)
    for n in 1:N
        F[n] = u -> (-u[n] .+ s*tanh.(u[n]) .+ gg/sqrt(N)*W[n,:]'*tanh.(u))
    end

    # F = F_list[18]
    # Nsamples = 50000
    alg = :accelerated
    r_num = zeros(Int,3)
    for (j,g) in enumerate(g_list)
        roots =  typeof(rand(N))[] 
        g_eps(x) = g(x,ε/2)
        for k in 1:Nsamples
            X0 = 5*(rand(N) .- 0.5)*2 
            ds = setup_iterator(F, g_eps, X0; algtype = alg)
            n, xf, q = compute_figure(ds, ε, max_it)
            if n < max_it
                custom_mapper(xf, roots, 0.01)
            end
        end
        r_num[j] = length(roots)
    end
    @show r_num
    return r_num
end



function _roots_number(d) 
    @unpack dims, Nsamples, max_it = d 
    roots_N = zeros(Int, length(dims), length(g_list), length(Nsamples))
    rng = MersenneTwister(123);
    for (j,N) in enumerate(dims) 
        rr = zeros(Int,length(Nsamples),3)
        for (h,Ns) in enumerate(Nsamples)
            roots_N[j,:,h] = get_roots_number(N, Ns, max_it, rng)
        end
        @show roots_N[j,:,:]
    end
    return @strdict(dims, Nsamples, roots_N)
end
 # using JLD2
 # @save "tmp_Nsamples_roots.jld2" dims Nsamples roots_N

a = 0.9313508638295191
b = -0.12406996404465495
nroots_fit(x) = exp(a*x + b)

max_it = 200; dims = 3:12
Nsamples = round.(Int, logrange(1000,100000, length = 5))
force = false
d = @dict(dims, Nsamples, max_it) # parametros
data, file = produce_or_load(
    datadir(""), # path
    d, # container for parameter
    _roots_number, # function
    prefix = "roots_rand_net", # prefix for savename
    force = force, # true for forcing sims
    wsave_kwargs = (;compress = true)
)
@unpack roots_N = data

using CairoMakie
f = Figure(); 
ax = Axis(f[1,1], xlabel = L"N_{dim}", ylabel = L"\Delta N_{roots}") #, yscale = log10);

ind = 3
plot!(ax, dims, roots_N[:,1,ind] .- roots_N[:,3,ind], label = L"g_1")
plot!(ax, dims, roots_N[:,2,ind] .- roots_N[:,3,ind], label = L"g_2")
axislegend(ax; position = :rb) 
# plot!(ax, dims, nroots.(dims))

f
