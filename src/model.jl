struct Constraint{S<:MOI.AbstractSet}
    index::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64}, S}
    fun::DataPair{AffineFunction, MOI.VectorAffineFunction{Float64}}
    set::S
end

function Constraint(index::MOI.ConstraintIndex, f::AffineFunction, set::MOI.AbstractSet)
    Constraint(index, DataPair(f, MOI.VectorAffineFunction(f)), set)
end

mutable struct Model{O<:MOI.AbstractOptimizer}
    backend::SimpleQPModel{Float64}
    optimizer::O
    initialized::Bool
    objective::DataPair{QuadraticForm, MOI.ScalarQuadraticFunction{Float64}}
    nonnegconstraints::Vector{Constraint{MOI.Nonnegatives}}
    nonposconstraints::Vector{Constraint{MOI.Nonpositives}}
    zeroconstraints::Vector{Constraint{MOI.Zeros}}
    optimizer_to_backend::MOIU.IndexMap # FIXME: make type stable

    function Model(optimizer::O) where O
        VI = MOI.VariableIndex
        backend = SimpleQPModel{Float64}()
        initialized = false
        objective = DataPair(QuadraticForm(), MOI.ScalarQuadraticFunction(VI[], Float64[], VI[], VI[], Float64[], 0.0))
        nonnegconstraints = Constraint{MOI.Nonnegatives}[]
        nonposconstraints = Constraint{MOI.Nonpositives}[]
        zeroconstraints = Constraint{MOI.Zeros}[]
        optimizer_to_backend = MOIU.IndexMap()
        new{O}(backend, optimizer, initialized, objective, nonnegconstraints, nonposconstraints, zeroconstraints, optimizer_to_backend)
    end
end

function Variable(m::Model)
    m.initialized && error()
    index = MOI.addvariable!(m.backend)
    Variable(index)
end

function setobjective!(m::Model, sense::Senses.Sense, f::QuadraticForm)
    m.initialized && error()
    MOI.set!(m.backend, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), MOI.ScalarQuadraticFunction(f))
    MOI.set!(m.backend, MOI.ObjectiveSense(), MOI.OptimizationSense(sense))
    copyto!(m.objective.native, f)
    nothing
end

Base.push!(m::Model, c::Constraint{MOI.Nonnegatives}) = push!(m.nonnegconstraints, c)
Base.push!(m::Model, c::Constraint{MOI.Nonpositives}) = push!(m.nonposconstraints, c)
Base.push!(m::Model, c::Constraint{MOI.Zeros}) = push!(m.zeroconstraints, c)

function addconstraint!(m::Model, f::AffineFunction, set::MOI.AbstractVectorSet)
    m.initialized && error()
    index = MOI.addconstraint!(m.backend, MOI.VectorAffineFunction(f), set)
    push!(m, Constraint(index, f, set))
    nothing
end

add_nonnegative_constraint!(m::Model, f::AffineFunction) = addconstraint!(m, f, MOI.Nonnegatives(outputdim(f)))
add_nonpositive_constraint!(m::Model, f::AffineFunction) = addconstraint!(m, f, MOI.Nonpositives(outputdim(f)))
add_zero_constraint!(m::Model, f::AffineFunction) = addconstraint!(m, f, MOI.Zeros(outputdim(f)))

function initialize!(m::Model)
    # Copy
    result = MOI.copy!(m.optimizer, m.backend)
    if result.status == MOI.CopySuccess
        m.optimizer_to_backend = result.indexmap
    else
        error("Copy failed with status ", result.status, ". Message: ", result.message)
    end

    # Map variables
    # FIXME

    m.initialized = true
    nothing
end

function updateobjective!(m::Model)
    update!(m.objective)
    MOI.set!(m.optimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), m.objective.moi)
end

function updateconstraint!(m::Model, constraint::Constraint)
    update!(constraint.fun)
    MOI.set!(m.optimizer, MOI.ConstraintFunction(), constraint.index, constraint.fun.moi)
end

function update!(m::Model)
    # FIXME: map variables
    updateobjective!(m)
    updateconstraint!.(m, m.nonnegconstraints)
    updateconstraint!.(m, m.nonposconstraints)
    updateconstraint!.(m, m.zeroconstraints)
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
