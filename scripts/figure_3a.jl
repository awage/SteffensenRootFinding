using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))
# Add NonlinearSolve
using NonlinearSolve

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

# ---------------------------------------------------------
# NonlinearSolve.jl Helper Structures and Functions
# ---------------------------------------------------------
struct NonlinearSolveParams
    s::Float64
    gg::Float64
    W::Matrix{Float64}
    N::Int
end

# Out-of-place system function for NonlinearSolve
function F_system(u, p::NonlinearSolveParams)
    N = p.N
    s = p.s
    gg = p.gg
    W_mat = p.W
    
    tanh_u = tanh.(u) 
    # Calculation: -u + s*tanh(u) + (gg/sqrt(N)) * W * tanh(u)
    # Using matrix multiplication for the W term is usually faster and cleaner
    # than the element-wise loop if W is standard Matrix
    term_W = (gg / sqrt(N)) .* (W_mat * tanh_u)
    
    return -u .+ s .* tanh_u .+ term_W
end
# ---------------------------------------------------------

function compute_figure(ds, ε, max_it)
    # xf, fx = get_state(ds) 
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

    # 1. Setup Custom Function (Vector of functions)
    F_custom = Vector{Function}(undef,N)
    for n in 1:N
        F_custom[n] = u -> (-u[n] .+ s*tanh.(u[n]) .+ gg/sqrt(N)*W[n,:]'*tanh.(u))
    end

    # 2. Setup NonlinearSolve Function and Params
    p_nl = NonlinearSolveParams(s, gg, W, N)
    
    alg_custom = :accelerated
    
    # Initialize storage for unique roots found by each method
    # Structure: Vector of Vectors. 
    # Indices: 1..length(g_list) are custom, +1 is Newton, +2 is FastShortcut
    num_methods = length(g_list) + 2
    roots_storage = [typeof(rand(N))[] for _ in 1:num_methods]

    # --- MAIN LOOP ---
    for k in 1:Nsamples
        # A. Generate ONE Initial Condition for this sample
        X0_base = 5*(rand(N) .- 0.5)*2 

        # B. Run Custom Solvers (g_list)
        for (j, g) in enumerate(g_list)
            # Use copy(X0_base) to ensure no solver mutates the starting point for others
            X0 = copy(X0_base) 
            g_eps(x) = g(x, ε/2)
            ds = setup_iterator(F_custom, g_eps, X0; algtype = alg_custom)
            n, xf, q = compute_figure(ds, ε, max_it)
            
            if n < max_it
                custom_mapper(xf, roots_storage[j], 0.01)
            end
        end

        # C. Run NonlinearSolve: NewtonRaphson
        # Index: length(g_list) + 1
        X0_newton = copy(X0_base)
        prob_newton = NonlinearProblem(F_system, X0_newton, p_nl)
        sol_newton = solve(prob_newton, NewtonRaphson(), abstol = ε)
        
        if sol_newton.retcode == ReturnCode.Success
            custom_mapper(sol_newton.u, roots_storage[length(g_list) + 1], 0.01)
        end

        # D. Run NonlinearSolve: FastShortcutNonlinearPolyalg
        # Index: length(g_list) + 2
        X0_fast = copy(X0_base)
        prob_fast = NonlinearProblem(F_system, X0_fast, p_nl)
        sol_fast = solve(prob_fast, FastShortcutNonlinearPolyalg(), abstol = ε)
        
        if sol_fast.retcode == ReturnCode.Success
            custom_mapper(sol_fast.u, roots_storage[length(g_list) + 2], 0.01)
        end
    end

    # Convert vector of vectors to vector of counts
    r_num = length.(roots_storage)
    @show r_num
    return r_num
end

function _roots_number(d) 
    @unpack dims, Nsamples, Navg, max_it = d 
    # Adjusted size: length(g_list) + 2 to accommodate the 2 new solvers
    roots_N = zeros(Int, length(dims), length(g_list) + 2, Navg)
    rng = MersenneTwister(123);
    for (j,N) in enumerate(dims) 
        for h in 1:Navg
            roots_N[j,:,h] = get_roots_number(N, Nsamples, max_it, rng)
        end
        @show roots_N[j,:,:]
    end
    return @strdict(dims, Nsamples, Navg, roots_N)
end

# --- Plotting Section ---

max_it = 200; dims = 3:12
Nsamples = 10000
Navg = 10
force = true
d = @dict(dims, Navg, Nsamples, max_it) 

data, file = produce_or_load(
    datadir(""), 
    d, 
    _roots_number, 
    prefix = "roots_rand_net_avg_4", 
    force = force, 
    wsave_kwargs = (;compress = true)
)
@unpack roots_N = data

using CairoMakie
f = Figure(); 
ax = Axis(f[1,1]; xlabel = L"N_{dim}", ylabel = L"\Delta N_{roots}", xlabelsize = 25, ylabelsize = 25)

# Indices definition
idx_ref = 4            # Reference method (e.g., one of the g_list methods)
idx_newton = length(g_list) + 1
idx_fast = length(g_list) + 2

rr = mean(roots_N, dims = 3)

# Calculate Std Dev relative to reference
rs1 = std(roots_N[:,1,:] .- roots_N[:,idx_ref,:], dims = 2)
rs2 = std(roots_N[:,2,:] .- roots_N[:,idx_ref,:], dims = 2)
rs3 = std(roots_N[:,3,:] .- roots_N[:,idx_ref,:], dims = 2)
rs_newton = std(roots_N[:,idx_newton,:] .- roots_N[:,idx_ref,:], dims = 2)
rs_fast = std(roots_N[:,idx_fast,:] .- roots_N[:,idx_ref,:], dims = 2)

# Plot Errorbars
errorbars!(ax, dims, rr[:,1] .- rr[:,idx_ref], vec(rs1), color = :blue)
errorbars!(ax, dims, rr[:,2] .- rr[:,idx_ref], vec(rs2), color = :red)
errorbars!(ax, dims, rr[:,3] .- rr[:,idx_ref], vec(rs3), color = :black)
errorbars!(ax, dims, rr[:,idx_newton] .- rr[:,idx_ref], vec(rs_newton), color = :green)
errorbars!(ax, dims, rr[:,idx_fast] .- rr[:,idx_ref], vec(rs_fast), color = :orange)

# Plot Lines
plot!(ax, dims, rr[:,1] .- rr[:,idx_ref], label = L"g_1", color = :blue)
plot!(ax, dims, rr[:,2] .- rr[:,idx_ref], label = L"g_2", color = :red)
plot!(ax, dims, rr[:,3] .- rr[:,idx_ref], label = L"g_3", color = :black)
plot!(ax, dims, rr[:,idx_newton] .- rr[:,idx_ref], label = "NewtonRaphson", color = :green)
plot!(ax, dims, rr[:,idx_fast] .- rr[:,idx_ref], label = "FastShortcut", color = :orange)

axislegend(ax; position = :lt) 

# Inset plot for absolute number of roots
ax_inset = Axis(f[1, 1],
    width=Relative(0.35),
    height=Relative(0.35),
    halign=0.15,
    valign=0.15,
    title="Roots number (SM)", 
    yscale = log10)

plot!(ax_inset, dims, rr[:,1])
translate!(ax_inset.blockscene, 0, 0, 150)
ax_inset.xticks = [3,5, 7 ,9, 11,13]
ax_inset.yticks = [10, 100, 1000]
ax.xticks = dims

save("fig_roots_comparison.pdf",f)
