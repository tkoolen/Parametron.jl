mutable struct Objective{T}
    sense::Sense
    expr::WrappedExpression{QuadraticFunction{T}}
    moi_f::MOI.ScalarQuadraticFunction{T}

    function Objective{T}(sense::Sense, expr) where {T}
        converted = @expression convert(QuadraticFunction{T}, expr)
        wrapped = convert(WrappedExpression{QuadraticFunction{T}}, converted)
        new{T}(sense, wrapped, MOI.ScalarQuadraticFunction(wrapped()))
    end
end

# TODO: ScalarAffineFunction constraints
mutable struct Constraint{T, S<:MOI.AbstractSet}
    expr::WrappedExpression{Vector{AffineFunction{T}}}
    moi_f::MOI.VectorAffineFunction{T}
    set::S
    modelindex::MOI.ConstraintIndex{MOI.VectorAffineFunction{T}, S}
    optimizerindex::MOI.ConstraintIndex{MOI.VectorAffineFunction{T}, S}

    function Constraint{T}(expr, set::S) where {T, S<:MOI.AbstractSet}
        converted = @expression convert(Vector{AffineFunction{T}}, expr)
        wrapped = convert(WrappedExpression{Vector{AffineFunction{T}}}, converted)
        new{T, S}(wrapped, MOI.VectorAffineFunction(wrapped()), set)
    end
end

mutable struct Model{T, O<:MOI.AbstractOptimizer}
    params::Vector{Parameter}
    backend::SimpleQPMOIModel{T}
    optimizer::O
    initialized::Bool
    objective::Objective{T}
    nonnegconstraints::Vector{Constraint{T, MOI.Nonnegatives}}
    nonposconstraints::Vector{Constraint{T, MOI.Nonpositives}}
    zeroconstraints::Vector{Constraint{T, MOI.Zeros}}
    user_var_to_optimizer::Vector{MOI.VariableIndex}

    function Model{T}(optimizer::O) where {T, O}
        VI = MOI.VariableIndex
        params = Parameter[]
        backend = SimpleQPMOIModel{T}()
        initialized = false
        objective = Objective{T}(Minimize, @expression zero(QuadraticFunction{T}))
        nonnegconstraints = Constraint{T, MOI.Nonnegatives}[]
        nonposconstraints = Constraint{T, MOI.Nonpositives}[]
        zeroconstraints = Constraint{T, MOI.Zeros}[]
        user_var_to_optimizer = Vector{MOI.VariableIndex}()
        new{T, O}(params, backend, optimizer, initialized, objective, nonnegconstraints, nonposconstraints, zeroconstraints, user_var_to_optimizer)
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
    m.objective = Objective{T}(sense, expr)
    MOI.set!(m.backend, MOI.ObjectiveSense(), MOI.OptimizationSense(sense))
    MOI.set!(m.backend, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), m.objective.moi_f)
    nothing
end

addconstraint!(m::Model{T}, c::Constraint{T, MOI.Nonnegatives}) where {T} = push!(m.nonnegconstraints, c)
addconstraint!(m::Model{T}, c::Constraint{T, MOI.Nonpositives}) where {T} = push!(m.nonposconstraints, c)
addconstraint!(m::Model{T}, c::Constraint{T, MOI.Zeros}) where {T} = push!(m.zeroconstraints, c)

function addconstraint!(m::Model{T}, f, set::MOI.AbstractVectorSet) where T
    m.initialized && error()
    constraint = Constraint{T}(f, set)
    constraint.modelindex = MOI.addconstraint!(m.backend, MOI.VectorAffineFunction(constraint.expr()), set)
    addconstraint!(m, constraint)
    nothing
end

constraintdim(expr::LazyExpression) = length(expr())
constraintdim(val) = length(val)

add_nonnegative_constraint!(m::Model, f) = addconstraint!(m, f, MOI.Nonnegatives(constraintdim(f)))
add_nonpositive_constraint!(m::Model, f) = addconstraint!(m, f, MOI.Nonpositives(constraintdim(f)))
add_zero_constraint!(m::Model, f) = addconstraint!(m, f, MOI.Zeros(constraintdim(f)))

function mapindices(m::Model, idxmap)
    for c in m.nonnegconstraints
        c.optimizerindex = idxmap[c.modelindex]
    end
    for c in m.nonposconstraints
        c.optimizerindex = idxmap[c.modelindex]
    end
    for c in m.zeroconstraints
        c.optimizerindex = idxmap[c.modelindex]
    end
    backend_var_indices = MOI.get(m.backend, MOI.ListOfVariableIndices())
    resize!(m.user_var_to_optimizer, length(backend_var_indices))
    for index in backend_var_indices
        m.user_var_to_optimizer[index.value] = idxmap[index]
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
        mapindices(m, result.indexmap)
    else
        error("Copy failed with status ", result.status, ". Message: ", result.message)
    end
    m.initialized = true
    nothing
end

function update!(obj::Objective, m::Model)
    update!(obj.moi_f, obj.expr(), m.user_var_to_optimizer)
    MOI.set!(m.optimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), obj.moi_f)
    nothing
end

function update!(constraint::Constraint, m::Model)
    update!(constraint.moi_f, constraint.expr(), m.user_var_to_optimizer)
    MOI.set!(m.optimizer, MOI.ConstraintFunction(), constraint.optimizerindex, constraint.moi_f)
    nothing
end

"""
$(SIGNATURES)

Re-evaluate the expressions used to build the constraints and objective function of `Model` `m`.

Users should generally not need to call this function directly, as it is automatically
called in [`solve!`](@ref).
"""
function update!(m::Model)
    setdirty!(m)
    update!(m.objective, m)
    for c in m.nonnegconstraints
        update!(c, m)
    end
    for c in m.nonposconstraints
        update!(c, m)
    end
    for c in m.zeroconstraints
        update!(c, m)
    end
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
value(m::Model, x::Variable) = MOI.get(m.optimizer, MOI.VariablePrimal(), m.user_var_to_optimizer[x.index])

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

Add a constraint to the model using operators `==`, `<=`, or `>=`.

# Examples

Let `model` be a `Model` instance. The constraint `x >= zeros(2)` can
be added as follows:

```julia
julia> x = [Variable(model) for i = 1 : 2];

julia> @constraint(model, x >= zeros(2))
```
"""
macro constraint(model, expr)
    addcon = if @capture(expr, >=(lhs_, rhs_))
        :(SimpleQP.add_nonnegative_constraint!)
    elseif @capture(expr, ==(lhs_, rhs_))
        :(SimpleQP.add_zero_constraint!)
    elseif @capture(expr, <=(lhs_, rhs_))
        :(SimpleQP.add_nonpositive_constraint!)
    else
        return :(throw(ArgumentError("Relation not recognized")))
    end
    quote
        $addcon($model, @expression $lhs - $rhs)
    end |> esc
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
        setobjective!($model, $sense, @expression $expr)
    end |> esc
end
