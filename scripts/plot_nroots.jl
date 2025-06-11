using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
using CairoMakie
using Printf
using JLD2
@load "tmp_Nsamples_roots.jld2" dims Nsamples roots_N

f = Figure(); 
ax = Axis(f[1,1], xlabel = L"N_{dim}", ylabel = L"\Delta N_{roots}") #, yscale = log10);

a = 0.9313508638295191
b = -0.12406996404465495
nroots(x) = exp(a*x + b)

# plot!(ax, dims, roots_N[:,1,end]) 
ind = 10
plot!(ax, dims, roots_N[:,1,ind] .- roots_N[:,3,ind])
plot!(ax, dims, roots_N[:,2,ind] .- roots_N[:,3,ind])
# plot!(ax, dims, nroots.(dims))

f

# M = zeros(Int, length(dims), 3)
# for (k,m) in enumerate(roots_N)
#     @show m[end,:]
#     M[k,:] = m[end,:]
# end

# plot!(ax, dims, M[:,1] - M[:,3]);

# f
