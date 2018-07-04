"""
$(TYPEDEF)

Represents an expression that may be evaluated at a later time, by storing
both a function, `f`, and a tuple of function arguments, `args`.

`LazyExpression`s are typically not manually constructed by a user, and hence
are not exported. Instead, `LazyExpression`s should be created using the [`@expression`]
macro.

A `LazyExpression` may be evaluated by simply calling it with no arguments.

# Example

```julia
julia> a = ones(2); b = ones(2);

julia> expr = SimpleQP.LazyExpression(+, a, b)
LazyExpression{Base.#+, …}(…)

julia> expr()
2-element Array{Float64,1}:
 2.0
 2.0

julia> b .= 2
2-element Array{Float64,1}:
 2.0
 2.0

julia> expr()
2-element Array{Float64,1}:
 3.0
 3.0
```
"""
struct LazyExpression{F, A}
    f::F
    args::A
    LazyExpression(f::F, args...) where {F} = new{F, typeof(args)}(f, args)
end

Base.show(io::IO, expr::LazyExpression{F}) where {F} = print(io, "LazyExpression{$F, …}(…)")
function Base.show(io::IO, expr::LazyExpression{F}) where {F <: FunctionWrapper}
    print(io, "LazyExpression{FunctionWrapper{…}($(expr.f.obj[]))}(…)")
end


# Evaluation
@inline evalarg(x) = x
@inline evalarg(x::Parameter) = x()
@inline evalarg(x::LazyExpression) = x()
@inline (expr::LazyExpression{F, A})() where {F, A} = evalexpr(expr.f, expr.args)
@generated function evalexpr(f::F, args::Tuple{Vararg{Any, N}}) where {F, N}
    # equivalent to f(map(evalarg, args)...), minus the inference issues
    argexprs = [:(evalarg(args[$i])) for i = 1 : N]
    quote
        Base.@_inline_meta
        f($(argexprs...))
    end
end


# Wrapping
const WrappedExpression{T} = LazyExpression{FunctionWrapper{T, Tuple{}}, Tuple{}}

Base.convert(::Type{WrappedExpression{T}}, expr::LazyExpression) where {T} =
    LazyExpression(FunctionWrapper{T, Tuple{}}(expr))
Base.convert(::Type{WrappedExpression{T}}, expr::WrappedExpression{T}) where {T} = expr
Base.convert(::Type{WrappedExpression{T}}, value::T) where {T} =
    convert(WrappedExpression{T}, LazyExpression(identity, value))

"""
$(SIGNATURES)

Wrap a `LazyExpression` in a `FunctionWrappers.FunctionWrapper` and return a new
`LazyExpression` with the `FunctionWrapper` as the function `f` and an empty tuple
as the arguments `arg`.

The type parameters of the returned `LazyFunction` depend only on the type of the value
returned by `expr`. This is useful when a common interface is needed for different
`LazyExpression`s that share the same return value type.
"""
function wrap(expr::LazyExpression)
    T = typeof(expr())
    convert(WrappedExpression{T}, expr)
end


# expression macro
"""
$(SIGNATURES)

Create a new [`LazyExpression`](@ref) and apply optimizations to it to reduce allocations
and improve performance.

Expressions that do not depend on [`Parameter`](@ref)s or other [`LazyExpression`]s
are simply evaluated straight away.

# Examples

Creating an expression that represents `p * x1`, where `p` is a parameter that always evaluates to 2:

```julia
julia> model = SimpleQP.MockModel(); # a 'mock model' used only for demonstrations and tests

julia> x1 = Variable(model)
SimpleQP.Functions.Variable(1)

julia> p = Parameter{Int}(() -> 2, model)
Parameter{Int64, …}(…)

julia> expr = @expression p * x1
LazyExpression{FunctionWrapper{…}(LazyExpression{Base.#*, …}(…))}(…)

julia> expr()
2 * x1
```

Creating an expression that represents `p ⋅ x`, where `p` is a parameter that evaluates to [1, 2] and
`x` is a vector of two variables:

```julia
julia> model = SimpleQP.MockModel();

julia> x = Variable.(1 : 2);

julia> p = Parameter(identity, [1, 2], model)
Parameter{Array{Int64,1}, …}(…)

julia> expr = @expression p ⋅ x
LazyExpression{FunctionWrapper{…}(LazyExpression{SimpleQP.Functions.#vecdot!, …}(…))}(…)

julia> expr()
1 * x1 + 2 * x2 + 0

julia> @allocated expr()
0
```

Note that evaluating the expression does not allocate, because the ⋅ operation is
optimized and transformed into a call to the in-place `Functions.vecdot!` function.
"""
macro expression(expr)
    # We call expand(expr) to lower the code and make some useful conversions like:
    #   * [x; y] turns from Expr(:vcat, :x, :y) into Expr(:call, :vcat, :x, :y)
    #   * x .+ y turns from Expr(:call, :.+, :x, :y) into Expr(:call, :broadcast, :+, :x, :y) (on v0.6)
    postwalk(expand(expr)) do x
        if @capture(x, f_(args__))
            :(SimpleQP.optimize_toplevel(SimpleQP.LazyExpression($f, $(args...))))
        else
            if x isa Expr && x.head ∉ [:block, :line]
                buf = IOBuffer()
                dump(buf, expr)
                msg =
                    """
                    Unhandled expression head: $(x.head). expr:
                    $(String(take!(buf)))
                    """
                return :(throw(ArgumentError($msg)))
            end
            esc(x)
        end
    end
end


# Optimizations
function optimize_toplevel(@nospecialize expr::LazyExpression)
    expr isa WrappedExpression && throw(ArgumentError(("Cannot optimize wrapped expressions")))
    if any(arg -> arg isa Parameter || arg isa LazyExpression, expr.args)
        argtypes = map(arg -> typeof(evalarg(arg)), expr.args) # TODO: use Core.typeof in 0.7 to improve handling of Types
        return wrap(optimize(expr, argtypes...))
    else
        # simply evaluate
        return expr()
    end
end

optimizearg(arg) = arg
optimizearg(expr::LazyExpression{typeof(identity)}) = expr.args[1]

optimize(expr::LazyExpression, argtypes...) = expr

function optimize(expr::LazyExpression{typeof(*)}, ::Type{<:AbstractMatrix{T}}, ::Type{<:AbstractVector{<:Union{Variable, AffineFunction}}}) where T
    A, x = expr.args
    dest = deepcopy(expr())
    LazyExpression(Functions.matvecmul!, dest, A, x)
end

function optimize(expr::LazyExpression{typeof(*)},
        ::Type{<:TransposeVector{Variable}},
        ::Type{<:AbstractMatrix{T}},
        ::Type{<:AbstractVector{Variable}}) where {T}
    x, Q, y = expr.args
    dest = zero(QuadraticFunction{T})
    LazyExpression(Functions.bilinearmul!, dest, Q, x, y)
end

function optimize(expr::LazyExpression{<:Union{typeof(vecdot), typeof(dot)}}, ::Type{<:AbstractVector}, ::Type{<:AbstractVector})
    x, y = expr.args
    dest = deepcopy(expr())
    LazyExpression(Functions.vecdot!, dest, x, y)
end

function optimize(expr::LazyExpression{typeof(+)}, ::Type, ::Type, ::Type, ::Type...)
    optimize_toplevel(LazyExpression(+, optimize_toplevel(LazyExpression(+, expr.args[1], expr.args[2])), expr.args[3 : end]...))
end

function optimize(expr::LazyExpression{typeof(+)}, ::Type, ::Type)
    ret = expr()
    if ret isa AffineFunction || ret isa QuadraticFunction
        LazyExpression(Functions.add!, zero(typeof(ret)), expr.args...)
    elseif ret isa AbstractVector{<:AffineFunction}
        LazyExpression(Functions.vecadd!, deepcopy(ret), expr.args...)
    else
        expr
    end
end

function optimize(expr::LazyExpression{typeof(-)}, ::Type, ::Type)
    ret = expr()
    if ret isa AffineFunction || ret isa QuadraticFunction
        LazyExpression(Functions.subtract!, zero(typeof(ret)), expr.args...)
    elseif ret isa AbstractVector{<:AffineFunction}
        LazyExpression(Functions.vecsubtract!, deepcopy(ret), expr.args...)
    else
        expr
    end
end

function optimize(expr::LazyExpression{typeof(*)}, ::Type{<:AffineFunction}, ::Type{<:Union{Number, Variable, LinearTerm}})
    LazyExpression(Functions.mul!, deepcopy(expr()), expr.args...)
end

function optimize(expr::LazyExpression{typeof(*)}, ::Type{<:QuadraticFunction}, ::Type{<:Union{Number, Variable}})
    LazyExpression(Functions.mul!, deepcopy(expr()), expr.args...)
end

function optimize(expr::LazyExpression{typeof(*)}, ::Type{<:Union{Number, Variable, LinearTerm}}, ::Type{<:AffineFunction})
    LazyExpression(Functions.mul!, deepcopy(expr()), expr.args...)
end

function optimize(expr::LazyExpression{typeof(*)}, ::Type{<:Union{Number, Variable}}, ::Type{<:QuadraticFunction})
    LazyExpression(Functions.mul!, deepcopy(expr()), expr.args...)
end

function optimize(expr::LazyExpression{typeof(vcat)}, ::Type{<:AbstractVector{<:AffineFunction}}...)
    LazyExpression(Functions.vcat!, deepcopy(expr()), expr.args...)
end

function optimize(expr::LazyExpression{typeof(convert)}, ::Type, ::Type{<:AbstractVector})
    LazyExpression(Compat.copyto!, deepcopy(expr()), expr.args[2])
end

function optimize(expr::LazyExpression{typeof(*)}, ::Type{<:Number}, ::Type{<:AbstractVector{<:Union{Variable, AffineFunction}}})
    LazyExpression(Functions.scale!, deepcopy(expr()), expr.args...)
end

function optimize(expr::LazyExpression{typeof(*)}, ::Type{<:AbstractVector{<:Union{Variable, AffineFunction}}}, ::Type{<:Number})
    LazyExpression(Functions.scale!, deepcopy(expr()), expr.args...)
end

function optimize(expr::LazyExpression{typeof(Base.vect)}, ::Type{<:Union{Variable, LinearTerm, AffineFunction}})
    LazyExpression(deepcopy(expr()), expr.args...) do dest, x
        @boundscheck size(dest) == (1,) || throw(DimensionMismatch())
        dest[1] = x
        dest
    end
end
