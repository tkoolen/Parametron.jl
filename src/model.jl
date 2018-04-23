mutable struct Model{O<:MOI.AbstractOptimizer}
    backend::SimpleQPModel{Float64}
    optimizer::O
    initialized::Bool
    objective::QuadraticForm
    optimizer_to_backend::MOIU.IndexMap # FIXME: make type stable
    # constraintfuns::Vector{AffineFunction}

    function Model(optimizer::O) where O
        backend = SimpleQPModel{Float64}()
        initialized = false
        objective = QuadraticForm()
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
    m.objective = f
    MOI.set!(m.backend, MOI.ObjectiveSense(), MOI.OptimizationSense(sense))
    nothing
end

function initialize!(m::Model)
    moiobjective = MOI.ScalarQuadraticFunction(m.objective)
    MOI.set!(m.backend, MOI.ObjectiveFunction{typeof(moiobjective)}(), moiobjective)
    result = MOI.copy!(m.optimizer, m.backend)
    if result.status == MOI.CopySuccess
        m.optimizer_to_backend = result.indexmap
    else
        error("Copy failed with status ", result.status, ". Message: ", result.message)
    end
    m.initialized = true
    nothing
end
