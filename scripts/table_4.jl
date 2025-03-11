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
    ε = 1.e-8;  max_it = 200; force = false; Nsamples = Int(5e4)

    open("table4_dat.txt","w") do io
    for i in 1:5
        println(io,"{\\footnotesize ", optim_flist[i], "}" )

        grid = ntuple(i -> range(-10, 10, length = 10), length(X0_optim[i]))
        it = zeros(length(g_list))
        ex = zeros(length(g_list))
        cv = zeros(length(g_list))
        for alg in [:normal :accelerated]
                if alg == :normal
                    println(io,"& {\\footnotesize (norm.)}" )
                else 
                    println(io,"& {\\footnotesize (accel.)}" )
                end

        ds = [setup_iterator(optim_f[i], x -> g(x,ε), X0_optim[i]; algtype = alg) for g in g_list]

        for k in eachindex(g_list)
            d = get_stats(ds[k], Nsamples, grid, ε, max_it; seed = 123, prefix = string("stats_", alg, "_f", i, "_g",k ), force = force)
            @unpack nc, iterations, exec_time = d
            it[k] = iterations; ex[k] = exec_time; cv[k] = nc
            print(io," & ",  round(Float64((cv[k])*100), digits =1))
        end

        println(io," ")
        for k in 1:length(g_list)
            print(io," & ",  round(Float64((it[k])), digits =1))
        end

        # MEASURE TIMING FOR CONVERGING IC. We select only ic that converges for all functions. 
        k = 0; cnt = 0;
        tt = zeros(length(g_list))
        sampler, = statespace_sampler(grid)
        while k < 500  && cnt < Int(1e4)
            x0 = sampler()
            tm, ex_code = get_timing(ds, x0, ε, max_it) 
            if ex_code
                tt .= tt .+ tm
                k = k  + 1
            end
            cnt = cnt  + 1
        end

        println(io," ")
        for k in 1:length(g_list)
            print(io," & ",  round(Float64((tt[k]/tt[length(g_list)])), digits = 2))
        end
        println(io," \\\\")
        println(io," \\hline")
    end
    end
end
end

function get_timing(ds_v, x0, ε, max_it)
    tt = zeros(length(ds_v))
    for (k,ds) in enumerate(ds_v) 
        n, tt[k] = iterate(ds, x0, ε, max_it)
        if n ≥ max_it
            return tt, false
        end
    end

    return tt, true
end

print_table_all()


