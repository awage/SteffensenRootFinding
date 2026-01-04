using DrWatson
@quickactivate
using LaTeXStrings
include(srcdir("function_list.jl"))
using CairoMakie
f = Figure(); 
ax = Axis(f[1,1]; ylabel = L"g_i(x)", xlabel = L"x", xlabelsize = 25, ylabelsize = 25)

x = range(-4,4, length = 100)
# Plot Lines
lines!(ax, x, g_list[1].(x,1e-8), label = L"g_1(x) = tanh(x)", color = :blue)
lines!(ax, x, g_list[2].(x,1e-8), label = L"g_2(x) = \text{sign($x$)max($x,\varepsilon)$}", color = :red)
lines!(ax, x, 1e6*g_list[3].(1e-3*x,1e-8), label = L"g_3(x) = max(min(z^2, \varepsilon)  \text{ (scaled)}", color = :green)

axislegend(ax; position = :lt) 
save("fig1.pdf",f)

