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

function get_roots_number(N, Nsamples) 
    ε = 1e-8;  max_it = 1000; 
    s = 4.5; gg = 2.5; 
    rng = MersenneTwister(123);
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
            X0 = 5*randn(N)
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


dims = 5:15
Nsamples = round.(Int, logrange(1000,50000, length = 10))
roots_N = []
for N in dims 
    rr = zeros(Int,length(Nsamples),3)
    for (h,Ns) in enumerate(Nsamples)
        rr[h,:] = get_roots_number(N, Ns)
    end
    @show rr
    push!(roots_N, rr)
end

using JLD2
@save "tmp_Nsamples_roots.jld2" dims Nsamples roots_N


