mutable struct Model{O<:MOI.AbstractOptimizer}
    backend::SimpleQPMOIModel{Float64}
    optimizer::O
    initialized::Bool
    objective::DataPair{QuadraticFunction, MOI.ScalarQuadraticFunction{Float64}}
    nonnegconstraints::Vector{AffineConstraint{MOI.Nonnegatives}}
    nonposconstraints::Vector{AffineConstraint{MOI.Nonpositives}}
    zeroconstraints::Vector{AffineConstraint{MOI.Zeros}}
    varboundindices::Dict{Variable, VarBoundIndexPair} # TODO: CardinalDict?
    user_var_to_optimizer::Vector{MOI.VariableIndex}

    function Model(optimizer::O) where O
        VI = MOI.VariableIndex
        backend = SimpleQPMOIModel{Float64}()
        initialized = false
        objective = DataPair(QuadraticFunction(), MOI.ScalarQuadraticFunction(VI[], Float64[], VI[], VI[], Float64[], 0.0))
        nonnegconstraints = AffineConstraint{MOI.Nonnegatives}[]
        nonposconstraints = AffineConstraint{MOI.Nonpositives}[]
        zeroconstraints = AffineConstraint{MOI.Zeros}[]
        varboundindices = Dict{Variable, VarBoundIndexPair}()
        user_var_to_optimizer = Vector{MOI.VariableIndex}()
        new{O}(backend, optimizer, initialized, objective, nonnegconstraints, nonposconstraints, zeroconstraints, varboundindices, user_var_to_optimizer)
    end
end

function Variable(m::Model)
    m.initialized && error()
    index = MOI.addvariable!(m.backend)
    Variable(index)
end

function setobjective!(m::Model, sense::Senses.Sense, f)
    m.initialized && error()
    MOI.set!(m.backend, MOI.ObjectiveSense(), MOI.OptimizationSense(sense)) # TODO: consider putting in a DataPair as well
    m.objective.native = convert(QuadraticFunction, f)
    m.objective.moi = MOI.ScalarQuadraticFunction(m.objective.native)
    # update!(m.objective)
    MOI.set!(m.backend, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), m.objective.moi)
    nothing
end

Base.push!(m::Model, c::AffineConstraint{MOI.Nonnegatives}) = push!(m.nonnegconstraints, c)
Base.push!(m::Model, c::AffineConstraint{MOI.Nonpositives}) = push!(m.nonposconstraints, c)
Base.push!(m::Model, c::AffineConstraint{MOI.Zeros}) = push!(m.zeroconstraints, c)

function addconstraint!(m::Model, f, set::MOI.AbstractVectorSet)
    m.initialized && error()
    f_affine = convert(AffineFunction, f)
    constraint = AffineConstraint(f_affine, set)
    constraint.indexpair.modelindex = MOI.addconstraint!(m.backend, MOI.VectorAffineFunction(f_affine), set)
    push!(m, constraint)
    nothing
end

add_nonnegative_constraint!(m::Model, f) = addconstraint!(m, f, MOI.Nonnegatives(outputdim(f)))
add_nonpositive_constraint!(m::Model, f) = addconstraint!(m, f, MOI.Nonpositives(outputdim(f)))
add_zero_constraint!(m::Model, f) = addconstraint!(m, f, MOI.Zeros(outputdim(f)))

function setbounds!(m::Model, x::Variable, lower::Float64, upper::Float64)
    interval = MOI.Interval(lower, upper)
    if m.initialized
        MOI.modifyconstraint!(m.optimizer, m.varboundindices[x].optimizerindex, interval)
    else
        if haskey(m.varboundindices, x)
            index = m.varboundindices[x].modelindex
            MOI.modifyconstraint!(m.backend, index, interval)
        else
            indexpair = VarBoundIndexPair()
            indexpair.modelindex = MOI.addconstraint!(m.backend, MOI.SingleVariable(MOI.VariableIndex(x)), interval)
            m.varboundindices[x] = indexpair
        end
    end
    nothing
end

function mapindices!(indexpair::ConstraintIndexPair, idxmap)
    indexpair.optimizerindex = idxmap[indexpair.modelindex]
end

function mapindices!(m::Model, idxmap)
    for c in m.nonnegconstraints
        mapindices!(c.indexpair, idxmap)
    end
    for c in m.nonposconstraints
        mapindices!(c.indexpair, idxmap)
    end
    for c in m.zeroconstraints
        mapindices!(c.indexpair, idxmap)
    end
    for indexpair in values(m.varboundindices)
        mapindices!(indexpair, idxmap)
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
        mapindices!(m, result.indexmap)
    else
        error("Copy failed with status ", result.status, ". Message: ", result.message)
    end
    m.initialized = true
    nothing
end

function updateobjective!(m::Model)
    update!(m.objective, m.user_var_to_optimizer)
    MOI.set!(m.optimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), m.objective.moi)
end

function updateconstraint!(m::Model, constraint::AffineConstraint)
    update!(constraint.fun, m.user_var_to_optimizer)
    MOI.modifyconstraint!(m.optimizer, constraint.indexpair.optimizerindex, constraint.fun.moi)
end

function update!(m::Model)
    updateobjective!(m)
    for c in m.nonnegconstraints
        updateconstraint!(m, c)
    end
    for c in m.nonposconstraints
        updateconstraint!(m, c)
    end
    for c in m.zeroconstraints
        updateconstraint!(m, c)
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

function value(m::Model, x::Variable)
    MOI.get(m.optimizer, MOI.VariablePrimal(), m.user_var_to_optimizer[x.index])
end

function objectivevalue(m::Model)
    MOI.get(m.optimizer, MOI.ObjectiveValue())
end
