function addparameter! end

"""
$(TYPEDEF)

Represents a 'placeholder' for data; a value that may be filled in later.

`Parameter`s can be evaluated by simply calling them with no arguments.

`Parameter`s keep track of whether they've already been evaluated using `dirty` flag.
To reevalute a parameter, the dirty flag must first be set using `setdirty!(parameter)`.
The update function will then be called when the parameter itself is called.

# Examples

```julia
julia> model = Parametron.MockModel() # a 'mock model' used only for demonstrations and tests
Parametron.MockModel(Parametron.Parameter[], Base.RefValue{Int64}(1))

julia> value = Ref(1)
Base.RefValue{Int64}(1)

julia> p = Parameter{Int}(() -> value[], model)
Parameter{Int64, …}(…)

julia> p()
1

julia> value[] = 2
2

julia> p()
1

julia> Parametron.setdirty!(p); p()
2
```
"""
struct Parameter{T, F, InPlace}
    dirty::Base.RefValue{Bool}
    f::F
    val::Base.RefValue{T}

    # out-of-place
    """
    $(SIGNATURES)

    Create a new `Parameter` with an update function `f` that takes
    no arguments and returns a value of type `T`.
    """
    Parameter{T}(f::F, model) where {T, F} = addparameter!(model, new{T, F, false}(Ref(true), f, Base.RefValue{T}()))

    # in-place
    """
    $(SIGNATURES)

    Create a new 'in-place' `Parameter` with an update function `f` that takes
    `val` as its argument and updates it in place.
    """
    Parameter(f::F, val::T, model) where {T, F} = addparameter!(model, new{T, F, true}(Ref(true), f, Base.RefValue(val)))
end
Parameter(f, model) = Parameter{typeof(f())}(f, model)

isinplace(::Type{Parameter{T, F, InPlace}}) where {T, F, InPlace} = InPlace
Base.show(io::IO, param::Parameter{T, F, InPlace}) where {T, F, InPlace} = print(io, "Parameter{$T, …}(…)")

@inline function (parameter::Parameter{T})() where {T}
    if parameter.dirty[]
        update!(parameter)
        parameter.dirty[] = false
    end
    parameter.val[]::T
end

update!(parameter::Parameter{T, F, true}) where {T, F} = (parameter.f(parameter.val[]); nothing)
update!(parameter::Parameter{T, F, false}) where {T, F} = (parameter.val[] = parameter.f()::T; nothing)

setdirty!(parameter::Parameter) = (parameter.dirty[] = true; nothing)
