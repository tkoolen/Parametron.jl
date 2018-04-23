MOIU.@model(SimpleQPModel, # modelname
    (), # scalarsets
    (), # typedscalarsets
    (Zeros, Nonnegatives, Nonpositives), # vectorsets
    (), # typedvectorsets
    (), # scalarfunctions
    (ScalarQuadraticFunction,), # typedscalarfunctions
    (VectorAffineFunction,), # vectorfunctions
    () # typedvectorfunctions
)

mutable struct Model{O<:MOI.AbstractOptimizer}
    backend::SimpleQPModel{Float64}
    optimizer::O
    initialized::Bool
    objective::QuadraticForm
    # constraintfuns::Vector{AffineFunction}

    function Model(optimizer::O) where O
        backend = SimpleQPModel{Float64}()
        initialized = false
        objective = QuadraticForm()
        new{O}(backend, optimizer, initialized, objective)
    end
end

function Variable(m::Model)
    m.initialized && error()
    index = MOI.addvariable!(m.backend)
    Variable(index)
end

module Senses
import MathOptInterface
const MOI = MathOptInterface

@enum Sense Min Max

function MOI.OptimizationSense(sense::Sense)
    if sense == Min
        MOI.MinSense
    elseif sense == Max
        MOI.MaxSense
    else
        error()
    end
end
end

function setobjective!(m::Model, sense::Senses.Sense, f::QuadraticForm)
    m.initialized && error()
    m.objective = f
    MOI.set!(m.backend, MOI.ObjectiveSense(), MOI.OptimizationSense(sense))
    nothing
end
