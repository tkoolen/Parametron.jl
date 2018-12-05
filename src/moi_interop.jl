# MOI backend AbstractModel
MOIU.@model(ParametronMOIModel, # modelname
    (MOI.ZeroOne, MOI.Integer), # scalarsets
    (MOI.LessThan, MOI.GreaterThan, MOI.EqualTo), # typedscalarsets
    (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives), # vectorsets
    (), # typedvectorsets
    (MOI.SingleVariable,), # scalarfunctions
    (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction), # typedscalarfunctions
    (), # vectorfunctions
    (MOI.VectorAffineFunction,) # typedvectorfunctions
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
    moi_f.constant = f.constant[]
    resize!(moi_f.terms, length(f.linear))
    @inbounds for i in eachindex(f.linear)
        term = f.linear[i]
        moi_f.terms[i] = MOI.ScalarAffineTerm(term.coeff, varmap[term.var.index])
    end
    moi_f
end

function update!(moi_f::MOI.ScalarQuadraticFunction, f::QuadraticFunction, varmap = IdentityVarMap())
    affine = f.affine
    moi_f.constant = affine.constant[]
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

update!(moi_f::MOI.AbstractFunction, f::Nothing, varmap = IdentityVarMap()) = moi_f


# make_moi_equivalent
make_moi_equivalent(::Type{AffineFunction{T}}) where {T} =
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{T}[], zero(T))
make_moi_equivalent(::Type{QuadraticFunction{T}}) where {T} =
    MOI.ScalarQuadraticFunction(MOI.ScalarAffineTerm{T}[], MOI.ScalarQuadraticTerm{T}[], zero(T))
make_moi_equivalent(::Type{<:AbstractVector{AffineFunction{T}}}) where {T} =
    MOI.VectorAffineFunction(MOI.VectorAffineTerm{T}[], T[])


# canonical_function_type
canonical_function_type(::Type{Variable}, ::Type{T}) where {T} = AffineFunction{T}
canonical_function_type(::Type{<:LinearTerm}, ::Type{T}) where {T} = AffineFunction{T}
canonical_function_type(::Type{<:AffineFunction}, ::Type{T}) where {T} = AffineFunction{T}
canonical_function_type(::Type{<:AbstractVector{<:AffineFunction}}, ::Type{T}) where {T} = Vector{AffineFunction{T}}
canonical_function_type(::Type{<:QuadraticTerm}, ::Type{T}) where {T} = QuadraticFunction{T}
canonical_function_type(::Type{<:QuadraticFunction}, ::Type{T}) where {T} = QuadraticFunction{T}


# moi_to_native_type
moi_to_native_type(::Type{MOI.ScalarAffineFunction}) = AffineFunction{T} where T
moi_to_native_type(::Type{MOI.ScalarQuadraticFunction}) = QuadraticFunction{T} where T
moi_to_native_type(::Type{MOI.VectorAffineFunction}) = Vector{AffineFunction{T}} where T
moi_to_native_type(::Type{MOI.VectorQuadraticFunction}) = Vector{QuadraticFunction{T}} where T
moi_to_native_type(::Type{MOI.SingleVariable}) = Nothing


# Objective
struct Objective{E, F}
    expr::WrappedExpression{E}
    f::F
    isconstant::Bool
end

function Objective(::Type{T}, expr) where T
    val = evalarg(expr)
    E = canonical_function_type(typeof(val), T)
    converted = @expression convert(E, expr)
    isconstant = converted isa E # i.e., it's just a value; not a LazyExpression or a Parameter
    wrapped = convert(WrappedExpression{E}, converted)
    f = make_moi_equivalent(E)
    F = typeof(f)
    update!(f, wrapped())
    Objective{E, F}(wrapped, f, isconstant)
end

function update!(objective::Objective{E, F}, optimizer::MOI.AbstractOptimizer, varmap) where {E, F}
    if !objective.isconstant
        update!(objective.f, objective.expr(), varmap)
        MOI.set(optimizer, MOI.ObjectiveFunction{F}(), objective.f)
    end
    nothing
end


# Constraint
mutable struct Constraint{E, F<:MOI.AbstractFunction, S<:MOI.AbstractSet}
    expr::WrappedExpression{E}
    f::F
    set::S
    isconstant::Bool
    modelindex::MOI.ConstraintIndex{F, S}
    optimizerindex::MOI.ConstraintIndex{F, S}

    function Constraint(::Type{T}, expr, set::S) where {T, S<:MOI.AbstractSet}
        val = evalarg(expr)
        E = canonical_function_type(typeof(val), T)
        converted = @expression convert(E, expr)
        isconstant = converted isa E # i.e., it's just a value; not a LazyExpression or a Parameter
        wrapped = convert(WrappedExpression{E}, converted)
        f = make_moi_equivalent(E)
        F = typeof(f)
        update!(f, wrapped())
        new{E, F, S}(wrapped, f, set, isconstant)
    end

    function Constraint(f::F, set::S) where {F<:MOI.AbstractFunction, S<:MOI.AbstractSet}
        E = Nothing
        wrapped = convert(WrappedExpression{E}, @expression nothing)
        new{E, F, S}(wrapped, f, set, true)
    end
end

function update!(constraint::Constraint, optimizer::MOI.AbstractOptimizer, varmap)
    if !constraint.isconstant
        update!(constraint.f, constraint.expr(), varmap)
        MOI.set(optimizer, MOI.ConstraintFunction(), constraint.optimizerindex, constraint.f)
    end
    nothing
end
update!(constraint::Constraint{Nothing}, optimizer::MOI.AbstractOptimizer, varmap) = nothing


# Constraints
let
    constraint_specs = [(MOI.ScalarAffineFunction{T} where T, MOI.GreaterThan{T} where T),
                        (MOI.ScalarAffineFunction{T} where T, MOI.LessThan{T} where T),
                        (MOI.ScalarAffineFunction{T} where T, MOI.EqualTo{T} where T),
                        (MOI.VectorAffineFunction{T} where T, MOI.Nonnegatives),
                        (MOI.VectorAffineFunction{T} where T, MOI.Nonpositives),
                        (MOI.VectorAffineFunction{T} where T, MOI.Zeros),
                        (MOI.ScalarQuadraticFunction{T} where T, MOI.GreaterThan{T} where T),
                        (MOI.ScalarQuadraticFunction{T} where T, MOI.LessThan{T} where T),
                        (MOI.ScalarQuadraticFunction{T} where T, MOI.EqualTo{T} where T),
                        (MOI.VectorQuadraticFunction{T} where T, MOI.Nonnegatives),
                        (MOI.VectorQuadraticFunction{T} where T, MOI.Nonpositives),
                        (MOI.VectorQuadraticFunction{T} where T, MOI.Zeros),
                        (MOI.SingleVariable, MOI.Integer),
                        (MOI.SingleVariable, MOI.ZeroOne)]
    fieldnames = Symbol[]
    constrainttypes = []
    for (F, S) in constraint_specs
        E = moi_to_native_type(F)
        funcsym = F isa UnionAll ? F.body.name.name : F.name.name
        setsym = S isa UnionAll ? S.body.name.name : S.name.name
        ctype = if F isa UnionAll && S isa UnionAll
            Constraint{E{T}, F{T}, S{T}} where T
        elseif F isa UnionAll
            Constraint{E{T}, F{T}, S} where T
        elseif S isa UnionAll
            Constraint{E, F, S{T}} where T
        else
            Constraint{E, F, S}
        end
        push!(fieldnames, Symbol(lowercase(String(funcsym)) * "_in_" * lowercase(String(setsym))))
        push!(constrainttypes, ctype)
    end

    # Constraints struct
    @eval begin
        struct Constraints{T}
            $([:($name::Vector{$(ctype isa UnionAll ? :($ctype{T}) : ctype)}) for (name, ctype) in zip(fieldnames, constrainttypes)]...)
        end

        function Constraints{T}() where {T}
            Constraints{T}($([:(Vector{$(ctype isa UnionAll ? :($ctype{T}) : ctype)}()) for ctype in constrainttypes]...))
        end
    end

    # Base.push! methods
    for ((F, S), name, ctype) in zip(constraint_specs, fieldnames, constrainttypes)
        if ctype isa UnionAll
            @eval function Base.push!(constraints::Constraints{T}, constraint::$(:($ctype{T}))) where T
                push!(constraints.$name, constraint)
            end
        else
            @eval Base.push!(constraints::Constraints, constraint::$(ctype)) = push!(constraints.$name, constraint)
        end
    end

    # update!
    @eval begin
        function update!(constraints::Constraints, optimizer::MOI.AbstractOptimizer, varmap)
            $(map(fieldnames) do fieldname
                quote
                    for constraint in constraints.$fieldname
                        update!(constraint, optimizer, varmap)
                    end
                end
            end...)
            nothing
        end
    end

    # mapindices!
    @eval begin
        function mapindices!(constraints::Constraints, idxmap)
            $(map(fieldnames) do fieldname
                quote
                    for constraint in constraints.$fieldname
                        constraint.optimizerindex = idxmap[constraint.modelindex]
                    end
                end
            end...)
            nothing
        end
    end
end
