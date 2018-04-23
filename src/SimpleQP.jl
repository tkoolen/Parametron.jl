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
    setobjective!

import MathOptInterface

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
const SparseSymmetric64 = Symmetric{Float64,SparseMatrixCSC{Float64,Int}}

include("util.jl")
include("functions.jl")
include("model.jl")

end # module
