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

ε = 1e-8;  max_it = 1000; 
setprecision(BigFloat, 100; base = 10)

s = 4.5; gg = 2.5; N = 10
rng = MersenneTwister(123);
W = randn(rng,N,N)
for k in 1:N; W[k,k] = 0.0; end;
# F(u) = -u .+ s*tanh.(u) .+ gg/sqrt(N)*W*tanh.(u) 
F = Vector{Function}(undef,N)
for n in 1:N
    F[n] = u -> (-u[n] .+ s*tanh.(u[n]) .+ gg/sqrt(N)*W[n,:]'*tanh.(u))
end

# F = F_list[18]

alg = :accelerated
roots =  typeof(rand(N))[] 

for k in 1:40000
    X0 = 5*randn(N)

    # for g in g_list
    g_eps(x) = g_list[1](x,ε/2)
        ds = setup_iterator(F, g_eps, X0; algtype = alg)
        n, xf, q = compute_figure(ds, ε, max_it)
        custom_mapper(xf, roots, 0.01)
        # @show n
        # println(" ---------")
    # end

    # println(" ---------")
    # println(" ---------")
end

@show length(roots)

roots2 =  typeof(rand(N))[] 

for k in 1:40000
    X0 = 5*randn(N)

    # for g in g_list
    g_eps(x) = g_list[3](x,ε/2)
        ds = setup_iterator(F, g_eps, X0; algtype = alg)
        n, xf, q = compute_figure(ds, ε, max_it)
        custom_mapper(xf, roots2, 0.01)
        # @show n
        # println(" ---------")
    # end

    # println(" ---------")
    # println(" ---------")
end

@show length(roots2)
