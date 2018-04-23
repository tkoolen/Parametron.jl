__precompile__()

module SimpleQP

export
    Variable,
    LinearFunction,
    AffineFunction,
    QuadraticForm,
    Model

export
    Senses

export
    quad,
    setobjective!,
    add_nonnegative_constraint!,
    add_nonpositive_constraint!,
    add_zero_constraint!,
    solve!

using Compat
import MathOptInterface
# import CardinalDicts: CardinalDict

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
const SparseSymmetric64 = Symmetric{Float64,SparseMatrixCSC{Float64,Int}}

include("util.jl")
include("functions.jl")
include("moi_interop.jl")
include("model.jl")

end # module
