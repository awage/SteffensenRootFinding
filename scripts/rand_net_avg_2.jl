using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
# Add NonlinearSolve and LinearAlgebra for the comparison
using NonlinearSolve
using LinearAlgebra 

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
    n, yy = get_iterations!(ds, ε, max_it)
    xf, fx = get_state(ds) 
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
    
    # 1. Define Function for Custom Solver (Vector of functions)
    F = Vector{Function}(undef,N)
    for n in 1:N
        F[n] = u -> (-u[n] .+ s*tanh.(u[n]) .+ gg/sqrt(N)*W[n,:]'*tanh.(u))
    end

    # 2. Define Function for NonlinearSolve.jl (In-place f!(du, u, p))
    # p = (s, gg, N, W)
    function f_lib!(du, u, p)
        s, gg, N, W = p
        u_tanh = tanh.(u)
        # du = -u + s*tanh(u) + (gg/sqrt(N)) * W * tanh(u)
        # Using mul! for W * u_tanh is more efficient, but direct * works too
        du .= -u .+ s .* u_tanh .+ (gg/sqrt(N)) .* (W * u_tanh)
    end
    p_lib = (s, gg, N, W)

    alg = :accelerated
    
    # Increase size by 1 to hold the library result at the end
    num_methods = length(g_list)
    r_num = zeros(Int, num_methods + 1) 

    # --- Run Custom Methods ---
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

    # --- Run NonlinearSolve.jl ---
    roots_lib = typeof(rand(N))[]
    for k in 1:Nsamples
        X0 = 5*(rand(N) .- 0.5)*2 
        prob = NonlinearProblem(f_lib!, X0, p_lib)
        # NewtonRaphson is a robust choice comparable to accelerated fixed-point
        # sol = solve(prob, NewtonRaphson(), abstol=ε, reltol=ε, maxiters=max_it)
        sol = solve(prob, FastShortcutNonlinearPolyalg(), abstol=ε, reltol=ε, maxiters=max_it)
        
        # Check convergence (SciMLBase.successful_retcode(sol) or specific codes)
        if sol.retcode == ReturnCode.Success
            custom_mapper(sol.u, roots_lib, 0.01)
        end
    end
    r_num[end] = length(roots_lib)

    @show r_num
    return r_num
end

function _roots_number(d) 
    @unpack dims, Nsamples, Navg, max_it = d 
    
    # Resize roots_N to hold custom methods + 1 library method
    # Assuming g_list is defined globally or imported from function_list.jl
    n_methods = length(g_list) + 1
    
    roots_N = zeros(Int, length(dims), n_methods, Navg)
    rng = MersenneTwister(123);
    
    for (j,N) in enumerate(dims) 
        # rr is temporary storage
        for h in 1:Navg
            roots_N[j,:,h] = get_roots_number(N, Nsamples, max_it, rng)
        end
        @show roots_N[j,:,:]
    end
    return @strdict(dims, Nsamples, Navg, roots_N)
end

a = 0.9313508638295191
b = -0.12406996404465495
nroots_fit(x) = exp(a*x + b)

max_it = 200; dims = 3:12
Nsamples = 20000
Navg = 10
force = true
d = @dict(dims, Navg, Nsamples, max_it)

data, file = produce_or_load(
    datadir(""),
    d,
    _roots_number,
    prefix = "roots_rand_net_avg_comp", # Changed prefix to avoid overwriting old data
    force = force,
    wsave_kwargs = (;compress = true)
)
@unpack roots_N = data

using CairoMakie
f = Figure(size = (800, 600)); 
ax = Axis(f[1,1]; xlabel = L"N_{dim}", ylabel = L"\Delta N_{roots}", xlabelsize = 25, ylabelsize = 25)

ind = 3 # Index of the baseline method in g_list
rr = mean(roots_N, dims = 3)

# Calculate standard deviations for error bars
rs1 = std(roots_N[:,1,:] .- roots_N[:,3,:], dims = 2)
rs2 = std(roots_N[:,2,:] .- roots_N[:,3,:], dims = 2)
# Library comparison (Last index vs Baseline)
rs_lib = std(roots_N[:,end,:] .- roots_N[:,3,:], dims = 2)

# Plot Differences
errorbars!(ax, dims, rr[:,1] .- rr[:,3], vec(rs1))
errorbars!(ax, dims, rr[:,2] .- rr[:,3], vec(rs2))
errorbars!(ax, dims, rr[:,end] .- rr[:,3], vec(rs_lib), color=:green) # Library errors

plot!(ax, dims, rr[:,1] .- rr[:,3], label = L"g_1")
plot!(ax, dims, rr[:,2] .- rr[:,3], label = L"g_2")
plot!(ax, dims, rr[:,end] .- rr[:,3], label = "NonlinearSolve", color=:green)

axislegend(ax; position = :rb) 

# --- Inset Plot (Total Roots) ---
ax_inset = Axis(f[1, 1],
    width=Relative(0.4),
    height=Relative(0.4),
    halign=0.2,
    valign=0.9,
    title="Total Roots found", 
    yscale = log10)

# Plot g_1
plot!(ax_inset, dims, rr[:,1], label=L"g_1")
# Plot Library Solver
plot!(ax_inset, dims, rr[:,end], label="NLS", color=:green)

translate!(ax_inset.blockscene, 0, 0, 150)
ax_inset.xticks = [3,5, 7 ,9, 11,13]
ax_inset.yticks = [10, 100, 1000]
ax.xticks = dims
axislegend(ax_inset; position = :rb, labelsize=10)

save("fig_roots_comparison.pdf",f)
