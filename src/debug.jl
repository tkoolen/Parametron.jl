_debug_args(expr::LazyExpression) = expr.args
_debug_args(expr::LazyExpression{<:FunctionWrapper}) = expr.f.obj[].args

function findallocs(io::IO, x, depth = 0, argnum = nothing)
    depth > 0 && print(io, "  "^depth)
    argnum != nothing && print(io, "[$argnum]: ")
    if x isa LazyExpression || x isa Parameter
        x isa Parameter && setdirty!(x)
        allocs = @allocated x()
        print(io, "$x: ")
        color = allocs > 0 ? :light_red : :green
        printstyled(io, allocs, color=color)
        print(io, " bytes")
        println(io)
    else
        println(io, typeof(x))
    end
    if x isa LazyExpression
        for (argnum, arg) in enumerate(_debug_args(x))
            findallocs(io, arg, depth + 1, argnum)
        end
    end
end

"""
$(SIGNATURES)

Utility function that can be used to track down allocations in [`LazyExpression`](@ref)s.

# Examples

The following session shows the output of `findallocs` if the expression doesn't allocate:

```julia
julia> model = Parametron.MockModel(); x = [Variable(model) for i in 1 : 2];

julia> param = Parameter{Int}(() -> 2, model)
Parameter{Int64, …}(…)

julia> expr = @expression param * x
LazyExpression{FunctionWrapper{…}(LazyExpression{Parametron.Functions.#scale!, …}(…))}(…)

julia> Parametron.findallocs(expr)
LazyExpression{FunctionWrapper{…}(LazyExpression{Parametron.Functions.#scale!, …}(…))}(…): 0 bytes
  [1]: Array{Parametron.Functions.LinearTerm{Int64},1}
  [2]: Parameter{Int64, …}(…): 0 bytes
  [3]: Array{Parametron.Functions.Variable,1}
```

In this session, `param` allocates, and `findallocs` reports the allocation:

```julia
julia> model = Parametron.MockModel(); x = [Variable(model) for i in 1 : 2];

julia> param = Parameter(() -> zeros(2), model)
Parameter{Array{Float64,1}, …}(…)

julia> expr = @expression param ⋅ x
LazyExpression{FunctionWrapper{…}(LazyExpression{Parametron.Functions.#vecdot!, …}(…))}(…)

julia> Parametron.findallocs(expr)
LazyExpression{FunctionWrapper{…}(LazyExpression{Parametron.Functions.#vecdot!, …}(…))}(…): 0 bytes
  [1]: Parametron.Functions.AffineFunction{Float64}
  [2]: Parameter{Array{Float64,1}, …}(…): 96 bytes
  [3]: Array{Parametron.Functions.Variable,1}
```
"""
findallocs(x) = findallocs(stdout, x)
