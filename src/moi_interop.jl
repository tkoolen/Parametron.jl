# TODO: Scalar constraints
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

# Variable conversion
Functions.Variable(v::MOI.VariableIndex) = Variable(v.value)
MOI.VariableIndex(v::Variable) = MOI.VariableIndex(v.index)


# Sense conversion
function MOI.OptimizationSense(sense::Sense)
    if sense == Minimize
        MOI.MinSense
    elseif sense == Maximize
        MOI.MaxSense
    else
        error()
    end
end


# update!
struct IdentityVarMap end
Base.getindex(::IdentityVarMap, i) = MOI.VariableIndex(i)

function update!(moi_f::MOI.ScalarAffineFunction, f::AffineFunction, varmap = IdentityVarMap())
    moi_f.constant == f.constant[] || throw(ArgumentError("MOI.ScalarAffineFunction constant can't be modified."))
    resize!(moi_f.terms, length(f.linear))
    @inbounds for i in eachindex(f.linear)
        term = f.linear[i]
        moi_f.terms[i] = MOI.ScalarAffineTerm(term.coeff, varmap[term.var.index])
    end
    moi_f
end

function update!(moi_f::MOI.ScalarQuadraticFunction, f::QuadraticFunction, varmap = IdentityVarMap())
    affine = f.affine
    moi_f.constant == affine.constant[] || throw(ArgumentError("MOI.ScalarQuadraticFunction constant can't be modified."))
    resize!(moi_f.affine_terms, length(affine.linear))
    @inbounds for i in eachindex(affine.linear)
        term = affine.linear[i]
        moi_f.affine_terms[i] = MOI.ScalarAffineTerm(term.coeff, varmap[term.var.index])
    end
    resize!(moi_f.quadratic_terms, length(f.quadratic))
    @inbounds for i in eachindex(f.quadratic)
        term = f.quadratic[i]
        var_index_1 = varmap[term.rowvar.index]
        var_index_2 = varmap[term.colvar.index]
        coeff = ifelse(term.rowvar == term.colvar, 2 * term.coeff, term.coeff)
        moi_f.quadratic_terms[i] = MOI.ScalarQuadraticTerm(coeff, var_index_1, var_index_2)
    end
    moi_f
end

function update!(moi_f::MOI.VectorAffineFunction, fs::Vector{<:AffineFunction}, varmap = IdentityVarMap())
    resize!(moi_f.constants, length(fs))
    num_linear_terms = 0
    for f in fs
        num_linear_terms += length(f.linear)
    end
    resize!(moi_f.terms, num_linear_terms)
    i = 1
    @inbounds for row in eachindex(fs)
        f = fs[row]
        for term in f.linear
            moi_f.terms[i] = MOI.VectorAffineTerm(Int64(row), MOI.ScalarAffineTerm(term.coeff, varmap[term.var.index]))
            i += 1
        end
        moi_f.constants[row] = f.constant[]
    end
    moi_f
end

# MOI function construction
function MOI.ScalarAffineFunction(f::AffineFunction{T}) where T
    ret = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{T}[], f.constant)
    update!(ret, f)
end

function MOI.ScalarQuadraticFunction(f::QuadraticFunction{T}) where T
    ret = MOI.ScalarQuadraticFunction(MOI.ScalarAffineTerm{T}[], MOI.ScalarQuadraticTerm{T}[], f.affine.constant[])
    update!(ret, f)
end

function MOI.VectorAffineFunction(fs::Vector{AffineFunction{T}}) where T
    ret = MOI.VectorAffineFunction(MOI.VectorAffineTerm{T}[], T[])
    update!(ret, fs)
end
