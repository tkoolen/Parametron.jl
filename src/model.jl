mutable struct Model{T, O<:MOI.AbstractOptimizer}
    params::Vector{Parameter}
    backend::SimpleQPMOIModel{T}
    optimizer::O
    initialized::Bool
    objective::Objective # abstract type to support different types of objectives. RTTI used in update!.
    constraints::Constraints{T}
    model_var_to_optimizer::Vector{MOI.VariableIndex}

    function Model{T}(optimizer::O) where {T, O}
        VI = MOI.VariableIndex
        params = Parameter[]
        backend = SimpleQPMOIModel{T}()
        initialized = false
        objective = Objective(T, @expression zero(AffineFunction{T}))
        MOI.set!(backend, MOI.ObjectiveSense(), MOI.OptimizationSense(Minimize))
        MOI.set!(backend, MOI.ObjectiveFunction{typeof(objective.f)}(), objective.f)
        constraints = Constraints{T}()
        model_var_to_optimizer = Vector{MOI.VariableIndex}()
        new{T, O}(params, backend, optimizer, initialized, objective, constraints, model_var_to_optimizer)
    end
end

"""
$(SIGNATURES)

Create a new `Model`, representing an optimization problem to be solved
by the optimizer `optimizer` (a `MathOptInterface.AbstractOptimizer`).
"""
Model(optimizer::MOI.AbstractOptimizer) = Model{Float64}(optimizer)

Base.show(io::IO, m::Model{T, O}) where {T, O} = print(io, "Model{$T, $O}(…)")

"""
$(SIGNATURES)

Mark all parameters associated with the model as 'dirty' (out of date),
meaning they must be updated upon their next evaluation.
"""
setdirty!(model::Model) = foreach(setdirty!, model.params)

addparameter!(model::Model, param::Parameter) = (push!(model.params, param); param)

"""
$(SIGNATURES)

Create a new decision variable (`Variable`) associated with the model.
"""
function Variable(m::Model)
    m.initialized && error("Model has already been initialized.")
    index = MOI.addvariable!(m.backend)
    Variable(index)
end

"""
$(SIGNATURES)

Set the objective function and optimization sense (`Minimize` or `Maximize`).
"""
function setobjective!(m::Model{T}, sense::Sense, expr) where T
    m.initialized && error("Model was already initialized. setobjective! can only be called before initialization.")
    m.objective = Objective(T, expr)
    MOI.set!(m.backend, MOI.ObjectiveSense(), MOI.OptimizationSense(sense))
    MOI.set!(m.backend, MOI.ObjectiveFunction{typeof(m.objective.f)}(), m.objective.f)
    nothing
end

function addconstraint!(m::Model, constraint::Constraint)
    m.initialized && error("Model was already initialized. addconstraint! can only be called before initialization.")
    constraint.modelindex = MOI.addconstraint!(m.backend, constraint.f, constraint.set)
    push!(m.constraints, constraint)
    nothing
end

_add_nonnegative_constraint!(m::Model{T}, expr, val::AffineFunction) where {T} =
    addconstraint!(m, Constraint(T, expr, MOI.GreaterThan(zero(T))))
_add_nonpositive_constraint!(m::Model{T}, expr, val::AffineFunction) where {T} =
    addconstraint!(m, Constraint(T, expr, MOI.LessThan(zero(T))))
_add_zero_constraint!(m::Model{T}, expr, val::AffineFunction) where {T} =
    addconstraint!(m, Constraint(T, expr, MOI.EqualTo(zero(T))))

_add_nonnegative_constraint!(m::Model{T}, expr, val::AbstractVector{<:AffineFunction}) where {T} =
    addconstraint!(m, Constraint(T, expr, MOI.Nonnegatives(length(val))))
_add_nonpositive_constraint!(m::Model{T}, expr, val::AbstractVector{<:AffineFunction}) where {T} =
    addconstraint!(m, Constraint(T, expr, MOI.Nonpositives(length(val))))
_add_zero_constraint!(m::Model{T}, expr, val::AbstractVector{<:AffineFunction}) where {T} =
    addconstraint!(m, Constraint(T, expr, MOI.Zeros(length(val))))

add_nonnegative_constraint!(m::Model, expr) = _add_nonnegative_constraint!(m, expr, evalarg(expr))
add_nonpositive_constraint!(m::Model, expr) = _add_nonpositive_constraint!(m, expr, evalarg(expr))
add_zero_constraint!(m::Model, expr) = _add_zero_constraint!(m, expr, evalarg(expr))

add_integer_constraint!(m::Model, x::Variable) = addconstraint!(m, Constraint(MOI.SingleVariable(MOI.VariableIndex(x)), MOI.Integer()))
add_binary_constraint!(m::Model, x::Variable) = addconstraint!(m, Constraint(MOI.SingleVariable(MOI.VariableIndex(x)), MOI.ZeroOne()))

function mapindices!(m::Model, idxmap)
    mapindices!(m.constraints, idxmap)
    backend_var_indices = MOI.get(m.backend, MOI.ListOfVariableIndices())
    resize!(m.model_var_to_optimizer, length(backend_var_indices))
    for index in backend_var_indices
        m.model_var_to_optimizer[index.value] = idxmap[index]
    end
end

"""
$(SIGNATURES)

Copy the problem to be solved to the optimizer.

Users should generally not need to call this function directly, as it is automatically
called the first time [`solve!`](@ref) is called on a `Model`.
"""
@noinline function initialize!(m::Model)
    result = MOI.copy!(m.optimizer, m.backend)
    if result.status == MOI.CopySuccess
        mapindices!(m, result.indexmap)
    else
        error("Copy failed with status ", result.status, ". Message: ", result.message)
    end
    m.initialized = true
    nothing
end

"""
$(SIGNATURES)

Re-evaluate the expressions used to build the constraints and objective function of `Model` `m`.

Users should generally not need to call this function directly, as it is automatically
called in [`solve!`](@ref).
"""
function update!(m::Model{T}) where T
    setdirty!(m)
    if m.objective isa Objective{AffineFunction{T}, MOI.ScalarAffineFunction{T}}
        update!(m.objective::Objective{AffineFunction{T}, MOI.ScalarAffineFunction{T}}, m.optimizer, m.model_var_to_optimizer)
    elseif m.objective isa Objective{QuadraticFunction{T}, MOI.ScalarQuadraticFunction{T}}
        update!(m.objective::Objective{QuadraticFunction{T}, MOI.ScalarQuadraticFunction{T}}, m.optimizer, m.model_var_to_optimizer)
    else
        error("Unhandled objective type")
    end
    update!(m.constraints, m.optimizer, m.model_var_to_optimizer)
    nothing
end

"""
$(SIGNATURES)

Solve the model `m`. (Re-)evaluate constraint and objective expressions, update the
optimizer's internal representation of the problem, and start the optimization procedure.
"""
function solve!(m::Model)
    if !m.initialized
        initialize!(m)
        m.initialized = true
    end
    update!(m)
    MOI.optimize!(m.optimizer)
    nothing
end

"""
$(SIGNATURES)

Return the value of variable `x` as determined by the optimizer.
"""
value(m::Model, x::Variable) = MOI.get(m.optimizer, MOI.VariablePrimal(), m.model_var_to_optimizer[x.index])

"""
$(SIGNATURES)

Return the value of the objective function at the solution found by the optimizer.
"""
objectivevalue(m::Model) = MOI.get(m.optimizer, MOI.ObjectiveValue())

"""
$(SIGNATURES)

Return the termination status of the solver.
"""
terminationstatus(m::Model) = MOI.get(m.optimizer, MOI.TerminationStatus())

"""
$(SIGNATURES)

Return information regarding the primal of the problem.
"""
primalstatus(m::Model) = MOI.get(m.optimizer, MOI.PrimalStatus())

"""
$(SIGNATURES)

Return information regarding the dual of the problem.
"""
dualstatus(m::Model) = MOI.get(m.optimizer, MOI.DualStatus())

"""
$(SIGNATURES)

Add a constraint to the model using operators `==`, `<=`, `>=`, or `in`/`∈`.

`in`/`∈` may only be used for single variables with a right hand side that is
one of:

* ℤ or Integers
* {0, 1} or ZeroOne

# Examples

Let `model` be a `Model` instance. The constraint `x >= zeros(2)` can
be added as follows:

```julia
julia> x = [Variable(model) for i = 1 : 2];

julia> @constraint(model, x >= zeros(2))
```

The constraint that variable `x[1]` should be an integer can be
expressed using:

```julia
julia> @constraint(model, x ∈ ℤ)
```
"""
macro constraint(model, expr)
    if !(expr isa Expr) || expr.head != :call || length(expr.args) != 3
        return :(throw(ArgumentError("Expected expression of the form `a relation b`")))
    end
    relation = expr.args[1]
    lhs = expr.args[2]
    rhs = expr.args[3]
    if relation == :(>=)
        ret = :(SimpleQP.add_nonnegative_constraint!($model, SimpleQP.@expression $lhs - $rhs))
    elseif relation == :(<=)
        ret = :(SimpleQP.add_nonpositive_constraint!($model, SimpleQP.@expression $lhs - $rhs))
    elseif relation == :(==)
        ret = :(SimpleQP.add_zero_constraint!($model, SimpleQP.@expression $lhs - $rhs))
    elseif relation ∈ [:∈, :in]
        if rhs ∈ [:ℤ, :Integers]
            ret = :(SimpleQP.add_integer_constraint!($model, $lhs))
        elseif rhs ∈ [:({0, 1}), :ZeroOne]
            ret = :(SimpleQP.add_binary_constraint!($model, $lhs))
        else
            return :(throw(ArgumentError("'in' only supports ℤ/Integers and {0, 1}/ZeroOne")))
        end
    else
        return :(throw(ArgumentError("Relation not recognized")))
    end
    esc(ret)
end

"""
$(SIGNATURES)

Set the objective function of the model.

# Examples

Let `model` be a `Model` instance. The objective 'minimize x ⋅ x' can
be added as follows:

```julia
julia> x = [Variable(model) for i = 1 : 2];

julia> @objective model Minimize x ⋅ x
"""
macro objective(model, sense, expr)
    quote
        SimpleQP.setobjective!($model, $sense, SimpleQP.@expression $expr)
    end |> esc
end
