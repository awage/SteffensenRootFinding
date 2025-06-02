using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))

function compute_figure(ds, ε, max_it)
@show     xf, fx = get_state(ds) 
    n, yy = get_iterations!(ds, ε, max_it)
@show     xf, fx = get_state(ds) 
    if 5 ≤ n < max_it
        q = estimate_ACOC!(n, yy)
    else
        q = 0
    end
    return n, xf, q
end

function print_table_all()
    ε = 1e-25;  max_it = 1000; 
    setprecision(BigFloat, 100; base = 10)

    open("table2_dat.txt","w") do io
    for i in 1:6
        println(io,L"{\footnotesize $f_{", i, L"}$}" )

        for alg in [:normal :accelerated]
                if alg == :normal
                    println(io,"& {\\footnotesize (norm.)}" )
                else 
                    println(io,"& {\\footnotesize (accel.)}" )
                end
            # Iterations
            xf_v = []
            q_v = []

            for (k,g) in enumerate(g_list)
                g_eps(x) = g(x,ε/2)
                ds = setup_iterator(F_list[i], g_eps, big.(X0[i]); algtype = alg)
                n, xf, q = compute_figure(ds, ε, max_it)
                @show xf,n
                println(" ---------")
                push!(xf_v, xf)
                push!(q_v, q)
                print(io," & ", n)
            end
            
            println(io," ")
            #  Final point 
             for k in 1:length(g_list)
                 if length(xf_v[k]) > 1
                     print(io," & (")
                     for x in xf_v[k]; print(io, round(Float64(x), digits =2), ", "); end
                     print(io,")")
                 else
                  print(io, " & ", round(Float64(xf_v[k]), digits =2)," ");
                 end
             end

            println(io," ")
             # convergence order. 
             for k in 1:length(g_list)
                 print(io," & ", round(Float64(q_v[k]), digits =1))
             end

            println(io," \\\\")

        end # algtype
    end

    end
end

print_table_all()



