mutable struct Model{O<:MOI.AbstractOptimizer}
    backend::SimpleQPModel{Float64}
    optimizer::O
    initialized::Bool
    objective::QuadraticForm
    moiobjective::MOI.ScalarQuadraticFunction{Float64}
    optimizer_to_backend::MOIU.IndexMap # FIXME: make type stable
    # constraintfuns::Vector{AffineFunction}

    function Model(optimizer::O) where O
        VI = MOI.VariableIndex
        backend = SimpleQPModel{Float64}()
        initialized = false
        objective = QuadraticForm()
        moiobjective = MOI.ScalarQuadraticFunction(VI[], Float64[], VI[], VI[], Float64[], 0.0)
        optimizer_to_backend = MOIU.IndexMap()
        new{O}(backend, optimizer, initialized, objective, moiobjective, optimizer_to_backend)
    end
end

function Variable(m::Model)
    m.initialized && error()
    index = MOI.addvariable!(m.backend)
    Variable(index)
end

function setobjective!(m::Model, sense::Senses.Sense, f::QuadraticForm)
    m.initialized && error()
    m.objective = f
    MOI.set!(m.backend, MOI.ObjectiveSense(), MOI.OptimizationSense(sense))
    nothing
end

function initialize!(m::Model)
    # Objective
    set!(m.moiobjective, m.objective)
    MOI.set!(m.backend, MOI.ObjectiveFunction{typeof(m.moiobjective)}(), m.moiobjective)

    # Constraints

    # Copy
    result = MOI.copy!(m.optimizer, m.backend)
    if result.status == MOI.CopySuccess
        m.optimizer_to_backend = result.indexmap
    else
        error("Copy failed with status ", result.status, ". Message: ", result.message)
    end
    m.initialized = true
    nothing
end

function update!(m::Model)
    m.initialized || error()

    # Objective
    set!(m.moiobjective, m.objective)
    MOI.set!(m.optimizer, MOI.ObjectiveFunction{typeof(m.moiobjective)}(), m.moiobjective)

    # Constraints

    nothing
end
