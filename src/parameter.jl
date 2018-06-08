function addparameter! end

struct Parameter{T, F, InPlace}
    dirty::Base.RefValue{Bool}
    f::F
    val::Base.RefValue{T}

    # out-of-place
    Parameter{T}(f::F, model) where {T, F} = addparameter!(model, new{T, F, false}(Ref(true), f, Base.RefValue{T}()))

    # in-place
    Parameter(f::F, val::T, model) where {T, F} = addparameter!(model, new{T, F, true}(Ref(true), f, Base.RefValue(val)))
end
Parameter(f, model) = Parameter{typeof(f())}(f, model)

isinplace(::Type{Parameter{T, F, InPlace}}) where {T, F, InPlace} = InPlace
Base.show(io::IO, param::Parameter{T, F, InPlace}) where {T, F, InPlace} = print(io, "Parameter{$T, …}(…)")

function (parameter::Parameter)()
    if parameter.dirty[]
        if isinplace(typeof(parameter))
            parameter.f(parameter.val[])
        else
            parameter.val[] = parameter.f()
        end
        parameter.dirty[] = false
    end
    parameter.val[]
end

setdirty!(parameter::Parameter) = (parameter.dirty[] = true; nothing)
