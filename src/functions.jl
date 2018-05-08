module Functions

export
    Variable,
    # ScalarLinearTerm,
    # VectorLinearTerm,
    # ScalarQuadraticTerm,
    Parameter,
    ScalarLinearFunction,
    VectorLinearFunction,
    QuadraticForm,
    ScalarAffineFunction,
    VectorAffineFunction,
    QuadraticFunction

export
    wrap
    # outputdim

using Compat
import FunctionWrappers: FunctionWrapper


# Variable
struct Variable
    index::Int
end
Base.transpose(v::Variable) = v


# ScalarLinearTerm, ScalarQuadraticTerm, VectorLinearTerm
struct ScalarLinearTerm{T}
    coeff::T
    var::Variable
end
getcoeff(term::ScalarLinearTerm) = term.coeff
setcoeff(term::ScalarLinearTerm, coeff) = ScalarLinearTerm(coeff, term.var)

struct VectorLinearTerm{T}
    row::Int
    scalarterm::ScalarLinearTerm{T}
end
getcoeff(term::VectorLinearTerm) = getcoeff(term.scalarterm)
setcoeff(term::VectorLinearTerm, coeff) = VectorLinearTerm(term.row, setcoeff(term.scalarterm, coeff))

struct ScalarQuadraticTerm{T}
    coeff::T
    rowvar::Variable
    colvar::Variable
end
getcoeff(term::ScalarQuadraticTerm) = term.coeff
setcoeff(term::ScalarQuadraticTerm, coeff) = ScalarQuadraticTerm(coeff, term.rowvar, term.colvar)

for Term in [:ScalarLinearTerm, :VectorLinearTerm, :ScalarQuadraticTerm]
    @eval begin
        Base.convert(::Type{$Term{T}}, term::$Term{T}) where {T} = term
        Base.convert(::Type{$Term{T}}, term::$Term) where {T} = setcoeff(term, convert(T, getcoeff(term)))
        Base.:-(term::$Term) = setcoeff(term, -getcoeff(term))
        Base.:*(c::Real, term::$Term) = setcoeff(term, c * getcoeff(term))
        Base.:*(term::$Term, c::Real) = c * term
    end
end


# Parameter
struct Parameter{T, N, F}
    size::NTuple{N, Int}
    f::F
end

Parameter{T}(s::NTuple{N, Int}, f::F) where {T, N, F} = Parameter{T, N, F}(s, f)
Parameter(val::T) where {T} = (f = () -> val; Parameter{T}(size(val), f))
(p::Parameter{T})() where {T} = p.f()::T
wrap(p::Parameter{T}) where {T} = Parameter{T}(p.size, FunctionWrapper{T, Tuple{}}(p.f))
Base.size(p::Parameter) = p.size

Base.:+(p::Parameter) = p

function Base.:-(p::Parameter{T}) where T <: Real
    f = let inner = p.f
        () -> -inner()
    end
    Parameter{T}(size(p), f)
end

function Base.:-(p::Parameter{T}) where T <: Vector
    ret = T()
    f = let inner = p.f, ret = ret
        function ()
            vec = inner()
            resize!(ret, length(vec))
            @inbounds map!(-, ret, vec)
            ret
        end
    end
    Parameter{T}(size(p), f)
end
# TODO: handle matrix case?

for op in [:+, :-, :*]
    @eval begin
        function Base.$op(p1::Parameter{T1}, p2::Parameter{T2}) where {T1 <: Real, T2 <: Real}
            T = promote_type(T1, T2)
            f = let f1 = p1.f, f2 = p2.f
                () -> $op(f1(), f2())
            end
            Parameter{T}((), f)
        end
    end
end

function Base.:*(pscalar::Parameter{T1}, pvector::Parameter{Vector{T2}}) where {T1 <: Real, T2 <: Real}
    T = promote_type(T1, T2)
    ret = Vector{T}()
    f = let fscalar = pscalar.f, fvector = pvector.f, ret = ret
        function ()
            scalar = fscalar()
            vector = fvector()
            resize!(ret, length(vector))
            @inbounds map!(x -> scalar * x, ret, vector)
            ret
        end
    end
    Parameter{Vector{T}}(size(pvector), f)
end
Base.:*(pvector::Parameter{Vector{T2}}, pscalar::Parameter{T1}) where {T1 <: Real, T2 <: Real} = pscalar * pvector


# ScalarLinearFunction, VectorLinearFunction, QuadraticForm
for (Fun, Term) in [
        (:ScalarLinearFunction, :ScalarLinearTerm),
        (:QuadraticForm, :ScalarQuadraticTerm)]

    @eval begin
        struct $Fun{T, F}
            terms::F
        end
        $Fun{T}(terms::F) where {T, F} = $Fun{T, F}(terms)
        Base.size(::$Fun) = ()
        Base.similar(::$Fun, ::Type{T}, terms::F) where {T, F} = $Fun{T, F}(terms)
        Base.similar(::$Fun{T}, terms::F) where {T, F} = $Fun{T, F}(terms)
        function Base.zero(::Type{$Fun{T}}) where T
            ret = Vector{$Term{T}}()
            terms = let ret = ret
                () -> ret
            end
            $Fun{T}(terms)
        end
    end
end

struct VectorLinearFunction{T, F}
    length::Int
    terms::F
end
VectorLinearFunction{T}(l::Int, terms::F) where {T, F} = VectorLinearFunction{T, F}(l, terms)
Base.length(fun::VectorLinearFunction) = fun.length
Base.size(fun::VectorLinearFunction) = (fun.length,)
Base.similar(fun::VectorLinearFunction, ::Type{T}, terms::F) where {T, F} = VectorLinearFunction{T, F}(length(fun), terms)
Base.similar(fun::VectorLinearFunction{T}, terms::F) where {T, F} = VectorLinearFunction{T, F}(length(fun), terms)

for (Fun, Term) in [
        (:ScalarLinearFunction, :ScalarLinearTerm),
        (:VectorLinearFunction, :VectorLinearTerm),
        (:QuadraticForm, :ScalarQuadraticTerm)]

    @eval begin
        wrap(l::$Fun{T}) where {T} = similar(l, FunctionWrapper{Vector{$Term{T}}, Tuple{}}(l.terms))

        Base.:+(l::$Fun) = l

        function Base.:-(l::$Fun{T}) where T
            ret = Vector{$Term{T}}()
            terms = let ltermsfun = l.terms, ret = ret
                function ()
                    lterms = ltermsfun()
                    resize!(ret, length(lterms))
                    @inbounds map!(-, ret, lterms)
                    ret
                end
            end
            similar(l, terms)
        end

        function Base.:+(l1::$Fun{T1}, l2::$Fun{T2}) where {T1, T2}
            @boundscheck size(l1) == size(l2) || throw(DimensionMismatch())
            T = promote_type(T1, T2)
            ret = Vector{$Term{T}}()
            terms = let l1termsfun = l1.terms, l2termsfun = l2.terms, ret = ret
                function ()
                    l1terms = l1termsfun()
                    l2terms = l2termsfun()
                    resize!(ret, length(l1terms) + length(l2terms))
                    i = 1
                    @inbounds for term in l1terms
                        ret[i] = term
                        i += 1
                    end
                    @inbounds for term in l2terms
                        ret[i] = term
                        i+= 1
                    end
                    ret
                end
            end
            similar(l1, T, terms)
        end

        function Base.:-(l1::$Fun{T1}, l2::$Fun{T2}) where {T1, T2}
            @boundscheck size(l1) == size(l2) || throw(DimensionMismatch())
            T = promote_type(T1, T2)
            ret = Vector{$Term{T}}()
            terms = let l1termsfun = l1.terms, l2termsfun = l2.terms, ret = ret
                function ()
                    l1terms = l1termsfun()
                    l2terms = l2termsfun()
                    resize!(ret, length(l1terms) + length(l2terms))
                    i = 1
                    @inbounds for term in l1terms
                        ret[i] = term
                        i += 1
                    end
                    @inbounds for term in l2terms
                        ret[i] = -term
                        i+= 1
                    end
                    ret
                end
            end
            similar(l1, T, terms)
        end

        function Base.:*(c::Parameter{T1}, l::$Fun{T2}) where {T1 <: Real, T2}
            T = promote_type(T1, T2)
            ret = Vector{$Term{T}}()
            terms = let c = c, lterms = l.terms, ret = ret
                function ()
                    cval = c()
                    ltermsval = l.terms()
                    resize!(ret, length(ltermsval))
                    @inbounds for i in eachindex(ltermsval)
                        ret[i] = cval * ltermsval[i]
                    end
                    ret
                end
            end
            similar(l, T, terms)
        end

        Base.:*(l::$Fun, c::Parameter{<:Real}) = c * l
    end
end

# VectorLinearFunction-specific
function Base.:*(A::Parameter{Matrix{T}}, x::Vector{Variable}) where T
    @boundscheck size(A)[2] == length(x) || throw(DimensionMismatch())
    ret = Vector{VectorLinearTerm{T}}()
    terms = let A = A, x = x, ret = ret
        function ()
            Aval = A()
            @boundscheck size(Aval, 2) == length(x) || throw(DimensionMismatch())
            n = length(Aval)
            resize!(ret, n)
            i = 1
            @inbounds for col = 1 : size(Aval, 2)
                for row = 1 : size(Aval, 1)
                    ret[i] = VectorLinearTerm(row, ScalarLinearTerm(Aval[i], x[col]))
                    i += 1
                end
            end
            ret
        end
    end
    VectorLinearFunction{T}(size(A)[1], terms)
end

function Base.:*(c::Parameter{T}, x::Vector{Variable}) where T <: Real
    ret = Vector{VectorLinearTerm{T}}()
    terms = let c = c, x = x, ret = ret
        function ()
            cval = c()
            resize!(ret, length(x))
            @inbounds for i in eachindex(x)
                ret[i] = VectorLinearTerm(i, ScalarLinearTerm(cval, x[i]))
            end
            ret
        end
    end
    VectorLinearFunction{T}(length(x), terms)
end
Base.:*(x::Vector{Variable}, c::Parameter{<:Real}) = c * x

# QuadraticForm-specific
function Base.dot(x::Vector{Variable}, l::VectorLinearFunction{T}) where T
    @boundscheck length(x) == length(l) || throw(DimensionMismatch())
    ret = Vector{ScalarQuadraticTerm{T}}()
    terms = let x = x, ltermsfun = l.terms, ret = ret
        function ()
            linearterms = ltermsfun()
            resize!(ret, length(linearterms))
            for i in eachindex(linearterms)
                @inbounds linearterm = linearterms[i]
                quadraticterm = ScalarQuadraticTerm(linearterm.scalarterm.coeff, x[linearterm.row], linearterm.scalarterm.var)
                @inbounds ret[i] = quadraticterm
            end
            ret
        end
    end
    QuadraticForm{T}(terms)
end


# ScalarAffineFunction, VectorAffineFunction
for (AffineFunction, Linear, Constant) in [
        (:ScalarAffineFunction, :ScalarLinearFunction, :(Parameter{T} where T<:Real)),
        (:VectorAffineFunction, :VectorLinearFunction, :(Parameter{Vector{T}} where T<:Real))]
    @eval begin
        struct $AffineFunction{L <: $Linear, C <: $Constant}
            linear::L
            constant::C

            function $AffineFunction(linear::L, constant::C) where {L <: $Linear, C <: $Constant}
                @boundscheck size(linear) == size(constant) || throw(DimensionMismatch())
                new{L, C}(linear, constant)
            end
        end

        wrap(a::$AffineFunction) = $AffineFunction(wrap(a.linear), wrap(a.constant))
    end

    for op in [:+, :-]
        @eval begin
            Base.$op(l::$Linear, c::$Constant) = $AffineFunction(l, $op(c))
            Base.$op(c::$Constant, l::$Linear) = $AffineFunction($op(l), c)
            Base.$op(a::$AffineFunction, c::$Constant) = $AffineFunction(a.linear, $op(a.constant, c))
            Base.$op(c::$Constant, a::$AffineFunction) = $AffineFunction(a.linear, $op(c, a.constant))
            Base.$op(a::$AffineFunction, l::$Linear) = $AffineFunction($op(a.linear, l), a.constant)
            Base.$op(l::$Linear, a::$AffineFunction) = $AffineFunction($op(l, a.linear), a.constant)
        end
    end
end

# ScalarAffineFunction-specific
ScalarAffineFunction(linear::ScalarLinearFunction{T}) where {T} = ScalarAffineFunction(linear, Parameter(zero(T)))
ScalarAffineFunction(constant::Parameter{T}) where T<:Real = ScalarAffineFunction(zero(ScalarLinearFunction{T}), constant)

# VectorAffineFunction-specific
VectorAffineFunction(linear::VectorLinearFunction{T}) where {T} = VectorAffineFunction(linear, Parameter(zeros(size(linear))))
function VectorAffineFunction(constant::Parameter{Vector{T}}) where {T}
    ret = Vector{VectorLinearTerm{T}}()
    terms = let ret = ret
        () -> ret
    end
    linear = VectorLinearFunction{T}(size(constant)[1], terms)
    VectorAffineFunction(linear, constant)
end

# QuadraticFunction
struct QuadraticFunction{Q <: QuadraticForm, A <: ScalarAffineFunction}
    quadratic::Q
    affine::A
end

wrap(q::QuadraticFunction) = QuadraticFunction(wrap(q.quadratic), wrap(q.affine))

for op in [:+, :-]
    @eval begin
        Base.$op(qform::QuadraticForm, c::Parameter{T} where T<:Real) = QuadraticFunction(qform, ScalarAffineFunction(c))
        Base.$op(c::Parameter{T} where T<:Real, qform::QuadraticForm) = qform + c
        # Base.$op(c::)
    end
end

end
