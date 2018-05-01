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

Functions.Variable(v::MOI.VariableIndex) = Variable(v.value)
MOI.VariableIndex(v::Variable) = MOI.VariableIndex(v.index)

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

struct IdentityVarMap
end
Base.getindex(::IdentityVarMap, i) = MOI.VariableIndex(i)

function update!(moi_f::MOI.ScalarQuadraticFunction, f::QuadraticFunction, varmap = IdentityVarMap())
    affine = f.affine
    constant = 0.0
    @inbounds for scaled in affine.constant.terms
        constant += scaled.scalar * scaled.val.v[1]
    end
    moi_f.constant == constant || throw(ArgumentError("MOI ScalarQuadraticFunction constant can't be modified."))

    empty!(moi_f.affine_variables)
    empty!(moi_f.affine_coefficients)
    for scaled in affine.linear.terms
        s = scaled.scalar
        linearterm = scaled.val
        A = linearterm.A
        x = linearterm.x
        for col in 1 : size(A, 2)
            push!(moi_f.affine_variables, varmap[x[col].index])
            push!(moi_f.affine_coefficients, A[col])
        end
    end

    empty!(moi_f.quadratic_rowvariables)
    empty!(moi_f.quadratic_colvariables)
    empty!(moi_f.quadratic_coefficients)
    quadratic = f.quadratic
    for scaled in quadratic.terms
        s = scaled.scalar
        quadraticterm = scaled.val
        Q = quadraticterm.Q
        x = quadraticterm.x
        Q.uplo == 'U' || error()
        Qdata = Q.data
        @inbounds for col = 1 : Qdata.n, k = Qdata.colptr[col] : (Qdata.colptr[col + 1] - 1) # from sparse findn
            row = Qdata.rowval[k]
            if row <= col # upper triangle
                push!(moi_f.quadratic_rowvariables, varmap[x[row].index])
                push!(moi_f.quadratic_colvariables, varmap[x[col].index])
                push!(moi_f.quadratic_coefficients, 2 * s * Qdata.nzval[k])
            end
        end
    end
    moi_f
end

function MOI.ScalarQuadraticFunction(f::QuadraticFunction)
    VI = MOI.VariableIndex
    affine = f.affine
    constant = 0.0
    @inbounds for scaled in affine.constant.terms
        constant += scaled.scalar * scaled.val.v[1]
    end
    ret = MOI.ScalarQuadraticFunction(VI[], Float64[], VI[], VI[], Float64[], constant)
    update!(ret, f)
    ret
end

function update!(moi_f::MOI.VectorAffineFunction, f::AffineFunction, varmap = IdentityVarMap())
    empty!(moi_f.outputindex)
    empty!(moi_f.variables)
    empty!(moi_f.coefficients)
    @inbounds for scaled in f.linear.terms
        s = scaled.scalar
        linearterm = scaled.val
        A = linearterm.A
        x = linearterm.x
        indices = CartesianIndices(A)
        for i in eachindex(A)
            index = indices[i]
            row = index[1]
            col = index[2]
            push!(moi_f.outputindex, row)
            push!(moi_f.variables, varmap[x[col].index])
            push!(moi_f.coefficients, A[i])
        end
    end
    resize!(moi_f.constant, outputdim(f))
    moi_f.constant[:] = 0
    @inbounds for scaled in f.constant.terms
        s = scaled.scalar
        c = scaled.val
        moi_f.constant .+= s .* c.v
    end
    moi_f
end

function MOI.VectorAffineFunction(f::AffineFunction)
    VI = MOI.VariableIndex
    ret = MOI.VectorAffineFunction(Int[], VI[], Float64[], Float64[])
    update!(ret, f)
    ret
end

mutable struct DataPair{A, B}
    native::A
    moi::B
end

function update!(pair::DataPair, varmap::Vector{MOI.VariableIndex})
    # TODO: hash?
    update!(pair.moi, pair.native, varmap)
end
