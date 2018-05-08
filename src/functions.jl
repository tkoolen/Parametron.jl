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
struct Parameter{T, F}
    f::F
end

Parameter{T}(f::F) where {T, F} = Parameter{T, F}(f)
Parameter(val::T) where {T} = (f = () -> val; Parameter{T, typeof(f)}(f))
(p::Parameter{T})() where {T} = p.f()::T
wrap(p::Parameter{T}) where {T} = Parameter{T}(FunctionWrapper{T, Tuple{}}(p.f))

Base.:+(p::Parameter) = p

function Base.:-(p::Parameter{T}) where T <: Real
    f = let inner = p.f
        () -> -inner()
    end
    Parameter{T}(f)
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
    Parameter{T}(f)
end
# TODO: handle matrix case?

for op in [:+, :-, :*]
    @eval begin
        function Base.$op(p1::Parameter{T1}, p2::Parameter{T2}) where {T1 <: Real, T2 <: Real}
            T = promote_type(T1, T2)
            f = let f1 = p1.f, f2 = p2.f
                () -> $op(f1(), f2())
            end
            Parameter{T}(f)
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
    Parameter{Vector{T}}(f)
end
Base.:*(pvector::Parameter{Vector{T2}}, pscalar::Parameter{T1}) where {T1 <: Real, T2 <: Real} = pscalar * pvector


# ScalarLinearFunction, VectorLinearFunction, QuadraticForm
for (Fun, Term) in [
        (:ScalarLinearFunction, :ScalarLinearTerm),
        (:VectorLinearFunction, :VectorLinearTerm),
        (:QuadraticForm, :ScalarQuadraticTerm)]

    @eval begin
        struct $Fun{T, F}
            terms::F
        end

        $Fun{T}(terms::F) where {T, F} = $Fun{T, F}(terms)
        wrap(l::$Fun{T}) where {T} = $Fun{T}(FunctionWrapper{Vector{$Term{T}}, Tuple{}}(l.terms))

        Base.:+(l::$Fun) = l

        function Base.:-(l::$Fun{T}) where T
            ret = Vector{$Term{T}}()
            terms = let inner = l.terms, ret = ret
                function ()
                    innerterms = inner()
                    resize!(ret, length(innerterms))
                    @inbounds map!(-, ret, innerterms)
                    ret
                end
            end
            $Fun{T}(terms)
        end

        function Base.:+(l1::$Fun{T1}, l2::$Fun{T2}) where {T1, T2}
            T = promote_type(T1, T2)
            ret = Vector{$Term{T}}()
            terms = let l1 = l1, l2 = l2, ret = ret
                function ()
                    empty!(ret)
                    append!(ret, l1.terms())
                    append!(ret, l2.terms())
                    ret
                end
            end
            $Fun{T}(terms)
        end

        function Base.:-(l1::$Fun{T1}, l2::$Fun{T2}) where {T1, T2}
            T = promote_type(T1, T2)
            ret = Vector{$Term{T}}()
            terms = let l1 = l1, l2 = l2, ret = ret
                function ()
                    empty!(ret)
                    append!(ret, l1.terms())
                    i = length(ret) + 1
                    l2terms = l2.terms()
                    n = length(ret) + length(l2terms)
                    resize!(ret, n)
                    @inbounds for term in l2.terms()
                        ret[i] = term
                        i += 1
                    end
                    ret
                end
            end
            $Fun{T}(terms)
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
            $Fun{T}(terms)
        end

        Base.:*(l::$Fun, c::Parameter{<:Real}) = c * l
    end
end

# VectorLinearFunction-specific
function Base.:*(A::Parameter{Matrix{T}}, x::Vector{Variable}) where T
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
    VectorLinearFunction{T}(terms)
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
    VectorLinearFunction{T}(terms)
end
Base.:*(x::Vector{Variable}, c::Parameter{<:Real}) = c * x

# QuadraticForm-specific
function Base.dot(x::Vector{Variable}, l::VectorLinearFunction{T}) where T
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
        end

        $AffineFunction(l::$Linear{S}) where {S} = $AffineFunction(l, zero($Constant{S}))

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
