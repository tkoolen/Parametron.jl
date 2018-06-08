__precompile__()

module SimpleQP

# types
export
    Variable,
    LinearTerm,
    QuadraticTerm,
    AffineFunction,
    QuadraticFunction,
    Model,
    Parameter

# enum values
export
    Minimize,
    Maximize

# functions
export
    setobjective!,
    add_nonnegative_constraint!,
    add_nonpositive_constraint!,
    add_zero_constraint!,
    solve!,
    value,
    objectivevalue

# macros
export
    @expression,
    @constraint,
    @objective

using Compat
using Compat.LinearAlgebra
import FunctionWrappers: FunctionWrapper
import MathOptInterface
import MacroTools: @capture, postwalk

const LinearAlgebra = Compat.LinearAlgebra
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

@static if VERSION >= v"0.7-"
    const TransposeVector{T, V<:AbstractVector{T}} = Transpose{T, V}
    const AdjointVector{T, V<:AbstractVector{T}} = Adjoint{T, V}
else
    const TransposeVector{T, V<:AbstractVector{T}} = RowVector{T, V}
    const AdjointVector{T, V<:AbstractVector{T}} = RowVector{T, ConjVector{T, V}}
end

@enum Sense Minimize Maximize

# include("util.jl")
include("functions.jl")
include("parameter.jl")

using .Functions

include("lazyexpression.jl")
include("moi_interop.jl")
include("model.jl")
include("mockmodel.jl")
include("debug.jl")

end # module
