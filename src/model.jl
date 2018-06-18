mutable struct Objective{T}
    sense::Sense
    expr::WrappedExpression{QuadraticFunction{T}}
    moi_f::MOI.ScalarQuadraticFunction{T}

    function Objective{T}(sense::Sense, expr::WrappedExpression{QuadraticFunction{T}}) where {T}
        new{T}(sense, expr, MOI.ScalarQuadraticFunction(expr()))
    end

    function Objective{T}(sense:: Sense, expr::LazyExpression) where {T}
        Objective(sense, WrappedExpression{QuadraticFunction{T}}(expr))
    end
end
Objective(sense::Sense, expr::WrappedExpression{QuadraticFunction{T}}) where {T} = Objective{T}(sense, expr)

# TODO: ScalarAffineFunction constraints
mutable struct Constraint{T, S<:MOI.AbstractSet}
    expr::WrappedExpression{Vector{AffineFunction{T}}}
    moi_f::MOI.VectorAffineFunction{T}
    set::S
    modelindex::MOI.ConstraintIndex{MOI.VectorAffineFunction{T}, S}
    optimizerindex::MOI.ConstraintIndex{MOI.VectorAffineFunction{T}, S}

    function Constraint{T}(expr::WrappedExpression{Vector{AffineFunction{T}}}, set::S) where {T, S<:MOI.AbstractSet}
        new{T, S}(expr, MOI.VectorAffineFunction(expr()), set)
    end

    function Constraint{T}(expr::LazyExpression, set::S) where {T, S<:MOI.AbstractSet}
        Constraint(WrappedExpression{Vector{AffineFunction{T}}}(expr), set)
    end
end
Constraint(expr::WrappedExpression{Vector{AffineFunction{T}}}, set::S) where {T, S<:MOI.AbstractSet} = Constraint{T}(expr, set)

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
        objfun = zero(QuadraticFunction{T})
        objective = Objective{T}(Minimize, @expression identity(objfun))
        nonnegconstraints = Constraint{T, MOI.Nonnegatives}[]
        nonposconstraints = Constraint{T, MOI.Nonpositives}[]
        zeroconstraints = Constraint{T, MOI.Zeros}[]
        user_var_to_optimizer = Vector{MOI.VariableIndex}()
        new{T, O}(params, backend, optimizer, initialized, objective, nonnegconstraints, nonposconstraints, zeroconstraints, user_var_to_optimizer)
    end
end
Model(optimizer::MOI.AbstractOptimizer) = Model{Float64}(optimizer)

Base.show(io::IO, m::Model{T, O}) where {T, O} = print(io, "Model{$T, $O}(…)")

setdirty!(model::Model) = foreach(setdirty!, model.params)
addparameter!(model::Model, param::Parameter) = (push!(model.params, param); param)

function Variable(m::Model)
    m.initialized && error()
    index = MOI.addvariable!(m.backend)
    Variable(index)
end

function setobjective!(m::Model, sense::Sense, expr)
    m.initialized && error("Model was already initialized. setobjective! can only be called before initialization.")
    m.objective.sense = sense
    m.objective.expr = expr
    m.objective.moi_f = MOI.ScalarQuadraticFunction(expr())
    MOI.set!(m.backend, MOI.ObjectiveSense(), MOI.OptimizationSense(sense))
    MOI.set!(m.backend, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), m.objective.moi_f)
    nothing
end

addconstraint!(m::Model{T}, c::Constraint{T, MOI.Nonnegatives}) where {T} = push!(m.nonnegconstraints, c)
addconstraint!(m::Model{T}, c::Constraint{T, MOI.Nonpositives}) where {T} = push!(m.nonposconstraints, c)
addconstraint!(m::Model{T}, c::Constraint{T, MOI.Zeros}) where {T} = push!(m.zeroconstraints, c)

function addconstraint!(m::Model{T}, f, set::MOI.AbstractVectorSet) where T
    f′ = @expression convert(Vector{AffineFunction{T}}, f) # to handle e.g. constraint functions that return SVectors
    m.initialized && error()
    constraint = Constraint{T}(f′, set)
    constraint.modelindex = MOI.addconstraint!(m.backend, MOI.VectorAffineFunction(f′()), set)
    addconstraint!(m, constraint)
    nothing
end

add_nonnegative_constraint!(m::Model, f) = addconstraint!(m, f, MOI.Nonnegatives(length(f())))
add_nonpositive_constraint!(m::Model, f) = addconstraint!(m, f, MOI.Nonpositives(length(f())))
add_zero_constraint!(m::Model, f) = addconstraint!(m, f, MOI.Zeros(length(f())))

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
    MOI.modifyconstraint!(m.optimizer, constraint.optimizerindex, constraint.moi_f)
    nothing
end

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

function solve!(m::Model)
    if !m.initialized
        initialize!(m)
        m.initialized = true
    end
    update!(m)
    MOI.optimize!(m.optimizer)
    nothing
end

value(m::Model, x::Variable) = MOI.get(m.optimizer, MOI.VariablePrimal(), m.user_var_to_optimizer[x.index])
objectivevalue(m::Model) = MOI.get(m.optimizer, MOI.ObjectiveValue())
terminationstatus(m::Model) = MOI.get(m.optimizer, MOI.TerminationStatus())
primalstatus(m::Model) = MOI.get(m.optimizer, MOI.PrimalStatus())

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
        $addcon($(esc(model)), @expression $(esc(:($lhs - $rhs))))
    end
end

macro objective(model, sense, expr)
    quote
        setobjective!($(esc(model)), $(esc(sense)), @expression $(esc(expr)))
    end
end
