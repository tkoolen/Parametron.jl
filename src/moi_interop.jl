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

function set!(sqf::MOI.ScalarQuadraticFunction, qf::QuadraticForm)
    sqf.constant == 0.0 || error()
    empty!(sqf.affine_variables)
    empty!(sqf.affine_coefficients)
    resize!(sqf.quadratic_rowvariables, 0)
    resize!(sqf.quadratic_colvariables, 0)
    resize!(sqf.quadratic_coefficients, 0)
    for i = 1 : numterms(qf)
        scale = prod(s -> s[], qf.scales[i])
        A = qf.As[i]
        Adata = A.data
        x = qf.xs[i]
        @inbounds for col = 1 : Adata.n, k = Adata.colptr[col] : (Adata.colptr[col + 1] - 1) # from sparse findn
            row = Adata.rowval[k]
            if row <= col # upper triangle
                push!(sqf.quadratic_rowvariables, x[row].index)
                push!(sqf.quadratic_colvariables, x[col].index)
                push!(sqf.quadratic_coefficients, scale * Adata.nzval[k])
            end
        end
    end
    sqf
end

function MOI.ScalarQuadraticFunction(f::QuadraticForm)
    VI = MOI.VariableIndex
    ret = MOI.ScalarQuadraticFunction(VI[], Float64[], VI[], VI[], Float64[], 0.0)
    set!(ret, f)
    ret
end
