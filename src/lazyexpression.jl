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

lazy_wrap(x) = esc(x)

function lazy_wrap(x::Expr)
    if x.head == :line
        return esc(x)
    elseif x.head == :$ && length(x.args) == 1
        return esc(x.args[1])
    elseif x.head == :block
        return Expr(x.head, lazy_wrap.(x.args)...)
    elseif x.head == :call
        if x.args[1] == GlobalRef(Core, :getfield)  # TODO: probably different on 0.7
            return esc(x)
        else
            return :(wrap(SimpleQP.optimize_toplevel(SimpleQP.LazyExpression($(lazy_wrap.(x.args)...)))))
        end
    else
        buf = IOBuffer()
        dump(buf, x)
        msg =
            """
            Unhandled expression head: $(x.head). expr:
            $(String(take!(buf)))
            """
        return :(throw(ArgumentError($msg)))
    end
end


macro expression(expr)
    lazy_wrap(expand(expr))
end


# Optimizations
function optimize_toplevel(@nospecialize expr::LazyExpression)
    if any(arg -> arg isa Parameter || arg isa LazyExpression, expr.args)
        expr′ = LazyExpression(expr.f, map(optimizearg, expr.args)...)
        argtypes = map(arg -> typeof(evalarg(arg)), expr.args) # TODO: use Core.typeof in 0.7 to improve handling of Types
        return optimize(expr′, argtypes...)
    else
        # simply evaluate
        return LazyExpression(identity, expr())
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

# Wrapping
const WrappedExpression{T} = LazyExpression{FunctionWrapper{T, Tuple{}}, Tuple{}}

Base.convert(::Type{WrappedExpression{T}}, expr::WrappedExpression{T}) where {T} = expr
Base.convert(::Type{WrappedExpression{T}}, expr::LazyExpression) where {T} =
    LazyExpression(FunctionWrapper{T, Tuple{}}(expr))

function wrap(expr::LazyExpression)
    T = typeof(expr())
    convert(WrappedExpression{T}, expr)
end

