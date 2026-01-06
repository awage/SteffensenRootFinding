using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
using LinearAlgebra
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))
using NonlinearSolve
using Printf

struct ExactKuramotoParams
    N::Int
    K::Float64
    omega::Vector{Float64}
end

function F_kuramoto_system(u, p::ExactKuramotoParams)
    N = p.N
    K = p.K
    omega = p.omega
    F = similar(u)
    
    for i in 1:N
        sum_sin = 0.0
        theta_i = u[i]
        for j in 1:N
            sum_sin += sin(u[j] - theta_i)
        end
        F[i] = omega[i] + (K / N) * sum_sin
    end
    return F
end

function compute_performance(ds, ε, max_it, theta_exact)
    n, yy = get_iterations!(ds, ε, max_it)
    xf, fx = get_state(ds) 
    
    if n < max_it && !any(isnan, xf)
        q = 1
    else
        q = 0
    end
    return n, q
end

function get_exact_success_rate(N, Nsamples, max_it, rng) 
    ε = 1e-8
    K = 5.0 # Strong coupling to ensure the chosen cluster is stable
    
    # Pre-allocate output
    num_methods = length(g_list) + 2
    success_counts = zeros(Int, num_methods)

    # --- MAIN LOOP ---
    for k in 1:Nsamples
        
        theta_exact = (rand(rng, N) .- 0.5) .* (π/2)
        
        #  CALCULATE OMEGA 
        omega_fixed = zeros(N)
        for i in 1:N
            interaction = 0.0
            for j in 1:N
                interaction += sin(theta_exact[j] - theta_exact[i])
            end
            # For theta_exact to be a root, omega + (K/N)*interaction = 0
            omega_fixed[i] = -(K / N) * interaction
        end

        # 3. Setup Problem
        p_nl = ExactKuramotoParams(N, K, omega_fixed)
        
        # Setup Custom Function Wrapper
        F_custom = Vector{Function}(undef, N)
        for i in 1:N
            F_custom[i] = u -> begin
                sum_sin = 0.0
                theta_i = u[i]
                for j in 1:N
                    sum_sin += sin(u[j] - theta_i)
                end
                return omega_fixed[i] + (K / N) * sum_sin
            end
        end

        # 4. Define Initial Condition
        perturbation = (rand(rng, N) .- 0.5) .* pi/2 # Random noise up to +/- 1.0 rad
        X0 = theta_exact .+ perturbation

        alg_custom = :accelerated

        for (j, g) in enumerate(g_list)
            X0_c = copy(X0) 
            g_eps(x) = g(x, ε/2)
            ds = setup_iterator(F_custom, g_eps, X0_c; algtype = alg_custom)
            
            n, q = compute_performance(ds, ε, max_it, theta_exact)
            success_counts[j] += q
        end

        # Run NonlinearSolve: NewtonRaphson
        idx_newton = length(g_list) + 1
        X0_newton = copy(X0)
        prob_newton = NonlinearProblem(F_kuramoto_system, X0_newton, p_nl)
        
        try
            sol = solve(prob_newton, NewtonRaphson(), abstol = ε)
            if sol.retcode == ReturnCode.Success
                success_counts[idx_newton] += 1
            end
        catch
        end

        # Run NonlinearSolve: FastShortcutNonlinearPolyalg
        idx_fast = length(g_list) + 2
        X0_fast = copy(X0)
        prob_fast = NonlinearProblem(F_kuramoto_system, X0_fast, p_nl)
        
        try
            sol = solve(prob_fast, FastShortcutNonlinearPolyalg(), abstol = ε)
            if sol.retcode == ReturnCode.Success
                success_counts[idx_fast] += 1
            end
        catch
        end
    end

    rates = success_counts ./ Nsamples
    @show rates
    return rates
end

function _success_rates_data(d) 
    @unpack dims, Nsamples, Navg, max_it = d 
    
    success_rates_N = zeros(Float64, length(dims), length(g_list) + 2, Navg)
    rng = MersenneTwister(123);
    
    for (j,N) in enumerate(dims) 
        for h in 1:Navg
            success_rates_N[j,:,h] = get_exact_success_rate(N, Nsamples, max_it, rng)
        end
        avg_for_dim = mean(success_rates_N[j,:,:], dims=2)
        println("Dim $N finished. Avg rates: $avg_for_dim")
    end
    
    return @strdict(dims, Nsamples, Navg, success_rates_N)
end

# --- Plotting / Execution Section ---

max_it = 200
dims = 10:10:50
Nsamples = 1000 
Navg = 3
force = true
d = @dict(dims, Navg, Nsamples, max_it) 

data, file = produce_or_load(
    datadir(""), 
    d, 
    _success_rates_data, 
    prefix = "success_rates_kuramoto_exact", 
    force = force, 
    wsave_kwargs = (;compress = true)
)

@unpack success_rates_N = data

mean_success = mean(success_rates_N, dims=3) 
println("Mean Success Rates (Exact Kuramoto Recovery):")
display(mean_success)

let 
m = mean_success[:,:]; ind = 1
for r in eachrow(m)
    print(dims[ind]); ind += 1
    for c in r
    print(" & ",  round(Float64((1-c)), digits =2))
end
println("\\\\")
end
end
