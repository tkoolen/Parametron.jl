MOIU.@model(SimpleQPMOIModel, # modelname
    (), # scalarsets
    (), # typedscalarsets
    (Zeros, Nonnegatives, Nonpositives), # vectorsets
    (), # typedvectorsets
    (), # scalarfunctions
    (ScalarQuadraticFunction,), # typedscalarfunctions
    (), # vectorfunctions
    (VectorAffineFunction,) # typedvectorfunctions
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
end # module

function set!(moi_f::MOI.ScalarQuadraticFunction, f::QuadraticForm)
    moi_f.constant == 0.0 || error()
    empty!(moi_f.affine_variables)
    empty!(moi_f.affine_coefficients)
    empty!(moi_f.quadratic_rowvariables)
    empty!(moi_f.quadratic_colvariables)
    empty!(moi_f.quadratic_coefficients)
    for i = 1 : numterms(f)
        s = scale(f, i)
        A = f.As[i]
        Adata = A.data
        x = f.xs[i]
        @inbounds for col = 1 : Adata.n, k = Adata.colptr[col] : (Adata.colptr[col + 1] - 1) # from sparse findn
            row = Adata.rowval[k]
            if row <= col # upper triangle
                push!(moi_f.quadratic_rowvariables, x[row].index)
                push!(moi_f.quadratic_colvariables, x[col].index)
                push!(moi_f.quadratic_coefficients, s * Adata.nzval[k])
            end
        end
    end
    moi_f
end

function MOI.ScalarQuadraticFunction(f::QuadraticForm)
    VI = MOI.VariableIndex
    ret = MOI.ScalarQuadraticFunction(VI[], Float64[], VI[], VI[], Float64[], 0.0)
    set!(ret, f)
    ret
end

function set!(moi_f::MOI.VectorAffineFunction, f::AffineFunction)
    linear = f.linear
    empty!(moi_f.outputindex)
    empty!(moi_f.variables)
    empty!(moi_f.coefficients)
    @inbounds for i = 1 : numterms(linear)
        s = scale(linear, i)
        A = linear.As[i]
        x = linear.xs[i]
        indices = CartesianIndices(A)
        for i in eachindex(A)
            index = indices[i]
            row = index[1]
            col = index[2]
            push!(moi_f.outputindex, row)
            push!(moi_f.variables, x[col].index)
            push!(moi_f.coefficients, A[i])
        end
    end
    resize!(moi_f.constant, outputdim(f))
    copyto!(moi_f.constant, f.constant)
    moi_f
end

function MOI.VectorAffineFunction(f::AffineFunction)
    VI = MOI.VariableIndex
    ret = MOI.VectorAffineFunction(Int[], VI[], Float64[], Float64[])
    set!(ret, f)
    ret
end

struct DataPair{A, B}
    native::A
    moi::B
end

function update!(pair::DataPair)
    # TODO: hash?
    set!(pair.moi, pair.native)
end
