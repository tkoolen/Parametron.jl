struct LazyExpression{F, A}
    f::F
    args::A
    LazyExpression(f::F, args...) where {F} = new{F, typeof(args)}(f, args)
end

Base.show(io::IO, expr::LazyExpression) = print(io, "LazyExpression{…}($(string(expr.f)), …)")


# Evaluation
evalarg(x) = x
evalarg(x::Parameter) = x()
evalarg(x::LazyExpression) = x()
(expr::LazyExpression)() = evalexpr(expr.f, expr.args...)
@generated function evalexpr(f, args::Vararg{Any, N}) where N
    # equivalent to f(map(evalarg, args)...), minus the inference issues
    argexprs = [:(evalarg(args[$i])) for i = 1 : N]
    :(f($(argexprs...)))
end


# expression macro
macro expression(expr)
    postwalk(expr) do x
        if @capture(x, f_(args__))
            :(SimpleQP.optimize_toplevel(SimpleQP.LazyExpression($f, $(args...))))
        else
            if x isa Expr && x.head ∉ [:block, :line, :globalref]
                dump(expr)
                error("Unhandled expression head: $(x.head)")
            end
            x
        end
    end |> esc
end


# Optimizatons
function optimize_toplevel(@nospecialize expr::LazyExpression)
    if any(arg -> arg isa Parameter || arg isa LazyExpression, expr.args)
        expr′ = LazyExpression(expr.f, map(optimizearg, expr.args)...)
        argtypes = map(arg -> typeof(evalarg(arg)), expr.args)
        return optimize(expr′, argtypes...)
    else
        # simply evaluate
        return LazyExpression(identity, expr())
    end
end

optimizearg(arg) = arg
optimizearg(expr::LazyExpression{typeof(identity)}) = expr.args[1]

optimize(expr::LazyExpression, argtypes...) = expr

function optimize(expr::LazyExpression{typeof(*)}, ::Type{<:DenseMatrix{T}}, ::Type{<:DenseVector{Variable}}) where T <: LinearAlgebra.BlasFloat
    A, x = expr.args
    rows, _ = size(A())
    y = [zero(AffineFunction{Float64}) for _ = 1 : rows]
    LazyExpression(Functions.matvecmul!, y, A, x)
end

function optimize(expr::LazyExpression{typeof(*)}, ::Type{<:TransposeVector{Variable}}, ::Type{<:AbstractMatrix{T}}, ::Type{<:DenseVector{Variable}}) where {T <: Number}
    x, Q, y = expr.args
    dest = zero(QuadraticFunction{T})
    LazyExpression(Functions.bilinearmul!, dest, Q, x, y)
end

function optimize(expr::LazyExpression{<:Union{typeof(vecdot), typeof(dot)}}, ::Type{<:DenseVector}, ::Type{<:DenseVector})
    x, y = expr.args
    dest = zero(expr())
    LazyExpression(Functions.vecdot!, dest, x, y)
end

function optimize(expr::LazyExpression{typeof(+)}, ::Type, ::Type, ::Type, ::Type...)
    optimize_toplevel(LazyExpression(+, optimize_toplevel(LazyExpression(+, expr.args[1], expr.args[2])), expr.args[3 : end]...))
end

function optimize(expr::LazyExpression{typeof(-)}, ::Type, ::Type, ::Type, ::Type...)
    optimize_toplevel(LazyExpression(-, optimize_toplevel(LazyExpression(-, expr.args[1], expr.args[2])), expr.args[3 : end]...))
end

function optimize(expr::LazyExpression{typeof(+)}, ::Type, ::Type)
    ret = expr()
    x, y = expr.args
    if ret isa AffineFunction || ret isa QuadraticFunction
        LazyExpression(Functions.add!, zero(typeof(ret)), x, y)
    elseif ret isa Vector{<:AffineFunction}
        LazyExpression(Functions.vecadd!, similar(ret), x, y)
    else
        expr
    end
end

function optimize(expr::LazyExpression{typeof(-)}, ::Type, ::Type)
    ret = expr()
    x, y = expr.args
    if ret isa AffineFunction || ret isa QuadraticFunction
        LazyExpression(Functions.subtract!, zero(typeof(ret)), x, y)
    elseif ret isa DenseVector{<:AffineFunction}
        LazyExpression(Functions.vecsubtract!, similar(ret), x, y)
    else
        expr
    end
end


# Wrapping
struct WrappedExpression{T}
    f::FunctionWrapper{T, Tuple{}}
end
(expr::WrappedExpression)() = expr.f()
Base.convert(::Type{WrappedExpression{T}}, expr::LazyExpression) where {T} = WrappedExpression{T}(expr)
