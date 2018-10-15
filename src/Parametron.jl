module Parametron

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
    objectivevalue,
    terminationstatus,
    primalstatus,
    dualstatus,
    findallocs

# macros
export
    @expression,
    @constraint,
    @objective

using LinearAlgebra
using DocStringExtensions

include("FunctionWrappersQuickFix.jl")
using .FunctionWrappersQuickFix: FunctionWrapper

import MathOptInterface
import MacroTools: @capture, postwalk

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

@enum Sense Minimize Maximize

include("util.jl")
include("functions.jl")
include("parameter.jl")

using .Functions

include("lazyexpression.jl")
include("moi_interop.jl")
include("model.jl")
include("mockmodel.jl")
include("debug.jl")

end # module
