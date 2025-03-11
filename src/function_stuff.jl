using LinearAlgebra
using Test  # Import the Test module for assertions

"""
Represents the state of an iterative solver.

Fields:
- `x`: Current solution vector/scalar.
- `fx`: Function evaluation at `x`.
- `dfx`: Derivative/Jacobian/approximation thereof at `x`.
"""
mutable struct State{N <: Number}
    x::Union{Vector{N}, N}
    fx::Union{Vector{N}, N}
    dfx::Union{Matrix{N}, Vector{N}, N}
end


"""
Iterator for root-finding algorithms, encapsulating the state and update map.

Fields:
- `update_map!`:  In-place update function that modifies the state.  Corresponds to the chosen Steffensen variant.
- `f`:  The function(s) to find the root(s) of.  Can be a single function or a vector of functions (for systems).
- `state`:  The current state of the solver.
"""
mutable struct FunIterator
    update_map!::Function
    f::Union{Function, Vector{Function}}
    state::State
end


"""
Sets up the `FunIterator` based on the chosen algorithm.

Args:
- `f`: The function or vector of functions to solve.
- `g`: A function used in the Steffensen-type update.
- `x`: Initial guess for the solution.
- `algtype`:  Symbol indicating the Steffensen variant to use: `:normal` or `:accelerated`.

Returns:
A `FunIterator` object, ready for iteration.
"""
function setup_iterator(f::Union{Function,Vector{Function}}, g::Function, x; algtype::Symbol = :normal)
    update_map! = case_steffenson_map(f, g, length(x), algtype)
    fi = FunIterator(update_map!, f, State(zero(x), zero(x), zero(x)))
    set_state!(fi, x)  # Initialize the state
    return fi
end

"""
Helper function to select and create the appropriate update map.
"""
function case_steffenson_map(f::Union{Function,Vector{Function}}, g::Function, d::Int, algtype::Symbol)
    if algtype == :normal
        gg = (fx, dfx) -> g(fx) 
    elseif algtype == :accelerated
        gg = (fx, dfx) -> g(-fx/dfx) 
    else
        throw(ArgumentError("Invalid algtype: $algtype. Choose :normal or :accelerated.")) 
    end

    return steffenson_map(f, gg, d) 
end



"""
Constructor for `FunIterator` that initializes the `State` struct.
"""
function FunIterator(update_map!::Function, f, x::Union{Vector{T}, T} where T <: Number; n::Int = length(x))
    fx = isa(f, Vector) ? map(h -> h(x), f) : f(x)
    dfx = if isa(f, Vector)
        ones(eltype(x), n, length(f))
    elseif n == 1
        one(eltype(x))
    else
        ones(eltype(x), n)
    end

    state = State(x, fx, dfx)
    return FunIterator(update_map!, f, state)
end


"""
Performs a single iteration step using the update map.
"""
function step!(fi::FunIterator)
    if norm(fi.state.x) > 1e12
        throw(DomainError(fi.state.x, "x is too large"))
    end
    fi.update_map!(fi.state)  
end


"""
Returns the current state (x, f(x)).
"""
get_state(fi::FunIterator) = fi.state.x, fi.state.fx


"""
Sets the state variables x, fx, and dfx. If fx or dfx are not provided, they are computed.
"""
function set_state!(fi::FunIterator, x; fx = nothing, dfx = nothing)
    fi.state.x = x

    if isnothing(fx)
        fi.state.fx = isa(fi.f, Vector) ? map(h -> h(x), fi.f) : fi.f(x)
    else
        fi.state.fx = fx
    end

    if isnothing(dfx)
        n = length(x)
        if isa(fi.f,Vector)
            fi.state.dfx =  ones(eltype(x), n, length(fi.f))
        elseif n > 1
            fi.state.dfx = ones(eltype(x), size(fi.state.dfx))
        else
            fi.state.dfx = one(eltype(x))
        end
    else
        fi.state.dfx = dfx
    end

    return nothing
end


"""
Performs iterations until a stopping criterion is met.

Args:
- `fi`: The `FunIterator` object.
- `ε`: Tolerance for convergence (norm(f(x)) < ε).
- `max_it`: Maximum number of iterations.

Returns:
- `k`: Number of iterations performed.
- `yy`: Vector of solution states (`x` values) at each iteration.  Includes the initial guess.
"""
function get_iterations!(fi::FunIterator, ε::Real, max_it::Int)
    xn, fx = get_state(fi)
    yy = Vector{typeof(xn)}(undef, max_it + 1)
    yy[1] = xn
    k = 1

    try
        while norm(fx) > ε && k < max_it
            step!(fi)
            xn, fx = get_state(fi)
            k += 1
            yy[k] = xn
        @show xn, fx
        end
    catch e
        @warn "Iteration failed at step $k: $e"
        @show xn
        return max_it, yy
    end
    return k, yy
end



# ----------------------- Steffensen Map Implementations -----------------------

# Barrier function.
function steffenson_map(f::Function, g::Function, d::Int)
    if d == 1
        return _steffenson_map(f, g)
    else
        return _steffenson_map(f, g, d)
    end
end


# Generalized Steffensen method real values
function _steffenson_map(f::Function, g::Function)
    function N!(S::State)
        x, fx, dfx = S.x, S.fx, S.dfx
        gx = g(fx, dfx)
        fx_h = f(x + gx)
        dfx_new = (fx_h - fx)/gx
        x = x - fx/dfx_new
        S.x = x; S.fx = f(x); S.dfx = dfx_new
        return nothing 
    end
    return N!
end


# Generalized Steffensen method for R^d → R
function _steffenson_map(f::Function, g::Function, d)
    function construct_gradient(f, g, d, S::State)
        fx = S.fx; Jx = S.dfx; x = S.x
        J = zeros(eltype(x), d)
        gx = [g(fx, Jx[k]) for k in 1:d] 
        G = zeros(eltype(x), d)
        for  k in 1:d
            G .= 0.0
            G[k] = gx[k]
            J[k] = (f(x .+ G) - fx)/gx[k]
        end
        return J
    end
    function N!(S::State)
        x = S.x; fx = S.fx
        Jx = construct_gradient(f, g, d, S)
        nJ = norm(Jx)
        if nJ > 0
            x_new = x - fx*Jx/nJ^2
        else
            x_new = x
        end
        S.x = x_new; S.fx = f(x_new); S.dfx = Jx  # Store the gradient in dfx
        return nothing 
    end
    return N!
end




# Generate a estimated jacobian function using the same technique for functions
# from R^k -> R^k._
function steffenson_map(f::Array{Function}, g::Function, d::Int)
    function construct_jacobian(f, g, d, S::State)
        Jx = S.dfx; fx = S.fx; x = S.x 
        J = zeros(eltype(x), d, d)

        gx = [ g(fx[n], Jx[n,k]) for n in 1:d, k in 1:d ] 

        G = zeros(eltype(x), d)
        for n in 1:d, k in 1:d
            G .= 0.0
            G[k] = gx[n,k]
            J[n,k] = (f[n](x .+ G) - fx[n])/gx[n,k]
        end
        return J
    end

    function N!(S::State)
        x = S.x; fx = S.fx
        Jx = construct_jacobian(f, g, d, S)
        if 0 < abs(det(Jx)) < Inf
            x_new =  x - inv(Jx)*fx
            S.x = x_new
        else
            S.x = x
        end
        S.fx = map(h -> h(S.x), f)
        S.dfx = Jx #Store the Jacobian here
        return nothing 
    end
    return N!
end
