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

function print_table_all()
    ε = 1.e-8;  max_it = 100; force = true; Nsamples = Int(5e4)
    setprecision(BigFloat, 100; base = 10)

    open("table3_dat.txt","w") do io
    # for i in [1:15 ; 17:21]
    for i in 1:21
        println(io,L"{\footnotesize $f_{", i, L"}$}" )

        grid = ntuple(i -> range(-10, 10, length = 10), length(X0[i]))
        it = zeros(length(g_list))
        ex = zeros(length(g_list))
        cv = zeros(length(g_list))
        for alg in [:normal :accelerated]
            if alg == :normal
                println(io,"& {\\footnotesize (norm.)}" )
            else 
                println(io,"& {\\footnotesize (accel.)}" )
            end

            ds = [setup_iterator(F_list[i], x -> g(x,ε), X0[i]; algtype = alg) for g in g_list]

            for k in eachindex(g_list)
                d = get_stats(ds[k], Nsamples, grid, ε, max_it; seed = 123, prefix = string("stats_", alg, "_f", i, "_g",k ), force = force)
                @unpack nc, iterations, exec_time = d
                it[k] = iterations; ex[k] = exec_time; cv[k] = nc
                print(io," & ",  round(Float64((cv[k])*100), digits =1))
            end

            # MEASURE TIMING and iterations FOR CONVERGING IC. We select only ic that converges for all functions. 
            nmb = 0; cnt = 0;
            tt = zeros(length(g_list))
            it = zeros(length(g_list))
            samp = sampler(grid,123)
            while nmb < 500  && cnt < Int(1e4)
                x0 = samp()
                n, tm, ex_code = get_timing(ds, x0, ε, max_it) 
                if ex_code
                    tt .= tt .+ tm
                    it .= it .+ n
                    nmb = nmb + 1
                end
                cnt = cnt  + 1
            end
            tt = tt./it

            println(io," ")
            for k in 1:length(g_list)
                print(io," & ",  round(it[k]/nmb, digits = 1))
            end

            println(io," ")
            for k in 1:length(g_list)
                print(io," & ",  round(Float64((tt[k]/tt[length(g_list)])), digits = 2))
            end
            
            setprecision(BigFloat, 100; base = 10)
            ds = [setup_iterator(F_list[i], x -> g(x,1e-25), big.(X0[i]); algtype = alg) for g in g_list]
            nmb = 0; cnt = 0;
            q = zeros(length(g_list))
            samp = sampler(grid,123)
            for k in eachindex(ds)
                while cnt < Int(1e4)
                    x0 = samp()
                    d = length(x0) 
                    set_state!(ds[k], d == 1 ? x0[1] : x0) 
                    n, yy = get_iterations!(ds[k], 1e-25, max_it)
                    if 6 ≤ n < max_it
                        q[k] = estimate_ACOC!(n, yy)
                        @show q 
                        break
                    end
                    cnt = cnt  + 1
                end
            end

            println(io," ")
            for k in 1:length(g_list)
                print(io," & ",  round(Float64(q[k]), digits = 2))
            end
            println(io," \\\\")
        end
    end
end
end

function get_timing(ds_v, x0, ε, max_it)
    tt = zeros(length(ds_v))
    n = zeros(Int, length(ds_v))
    for (k,ds) in enumerate(ds_v) 
        n[k], tt[k] = iterate(ds, x0, ε, max_it)
        if n[k] ≥ max_it
            return n, tt, false
        end
    end

    return n, tt, true
end

print_table_all()


