struct DataPair{A, B}
    native::A
    moi::B
end

mutable struct Model{O<:MOI.AbstractOptimizer}
    backend::SimpleQPModel{Float64}
    optimizer::O
    initialized::Bool
    objective::DataPair{QuadraticForm, MOI.ScalarQuadraticFunction{Float64}}
    optimizer_to_backend::MOIU.IndexMap # FIXME: make type stable
    # constraintfuns::Vector{AffineFunction}

    function Model(optimizer::O) where O
        VI = MOI.VariableIndex
        backend = SimpleQPModel{Float64}()
        initialized = false
        objective = DataPair(QuadraticForm(), MOI.ScalarQuadraticFunction(VI[], Float64[], VI[], VI[], Float64[], 0.0))
        optimizer_to_backend = MOIU.IndexMap()
        new{O}(backend, optimizer, initialized, objective, optimizer_to_backend)
    end
end

function Variable(m::Model)
    m.initialized && error()
    index = MOI.addvariable!(m.backend)
    Variable(index)
end

function setobjective!(m::Model, sense::Senses.Sense, f::QuadraticForm)
    m.initialized && error()
    copyto!(m.objective.native, f)
    MOI.set!(m.backend, MOI.ObjectiveSense(), MOI.OptimizationSense(sense))
    nothing
end

function initialize!(m::Model)
    # Objective
    set!(m.objective.moi, m.objective.native)
    MOI.set!(m.backend, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), m.objective.moi)

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
    set!(m.objective.moi, m.objective.native)
    MOI.set!(m.optimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), m.objective.moi)

    # Constraints

    nothing
end

function solve!(m::Model)
    if m.initialized
        update!(m)
    else
        initialize!(m)
    end
    MOI.optimize!(m.optimizer)
    nothing
end
