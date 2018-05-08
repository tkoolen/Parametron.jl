module Functions

export
    Variable,
    ScalarLinearTerm,
    VectorLinearTerm,
    ScalarQuadraticTerm

export
    Parameter,
    ScalarLinearFunction,
    VectorLinearFunction,
    QuadraticForm,
    ScalarAffineFunction,
    VectorAffineFunction,
    QuadraticFunction

export
    WrappedParameter,
    WrappedScalarLinearFunction,
    WrappedVectorLinearFunction,
    WrappedQuadraticForm,
    WrappedScalarAffineFunction,
    WrappedVectorAffineFunction,
    WrappedQuadraticFunction

export
    wrap

using Compat
import FunctionWrappers: FunctionWrapper


# Variable
struct Variable
    index::Int
end
Base.hash(v::Variable, h::UInt) = hash(v.index, h)
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
        Base.:*(c::Number, term::$Term) = setcoeff(term, c * getcoeff(term))
        Base.:*(term::$Term, c::Number) = c * term
    end
end


# Const
struct Const{T}
    val::T
end
(c::Const)() = c.val


# Parameter
struct Parameter{T, N, F}
    size::NTuple{N, Int}
    f::F
end

Parameter{T}(size::NTuple{N, Int}, f::F) where {T, N, F} = Parameter{T, N, F}(size, f)
Parameter(val::T) where {T} = (Parameter{T}(size(val), Const(copy(val))))
(p::Parameter{T})() where {T} = p.f()::T
Base.size(p::Parameter) = p.size
wrap(p::Parameter{T}) where {T} = Parameter{T}(p.size, FunctionWrapper{T, Tuple{}}(p.f))
const WrappedParameter{T, N} = Parameter{T, N, FunctionWrapper{T, Tuple{}}}

function Base.show(io::IO, p::Parameter{T}) where {T}
    s = size(p)
    dimstring = isempty(s) ? "Scalar" : length(s) == 1 ? "$(s[1])-element" : join(map(string, s), '×')
    print(io, dimstring, " Parameter{$T, …}(…)")
end

Base.:+(p::Parameter) = p

function Base.:-(p::Parameter{T}) where T <: Number
    f = let inner = p.f
        () -> -inner()
    end
    Parameter{T}(size(p), f)
end

function Base.:-(p::Parameter{T}) where T <: Array
    ret = T(undef, size(p))
    f = let inner = p.f, ret = ret
        function ()
            val = inner()
            @boundscheck size(val) == size(ret) || throw(DimensionMismatch())
            @inbounds map!(-, ret, val)
            ret
        end
    end
    Parameter{T}(size(p), f)
end

for op in [:+, :-, :*]
    @eval begin
        function Base.$op(p1::Parameter{T1}, p2::Parameter{T2}) where {T1 <: Number, T2 <: Number}
            T = promote_type(T1, T2)
            f = let f1 = p1.f, f2 = p2.f
                () -> $op(f1(), f2())
            end
            Parameter{T}((), f)
        end
    end
end

for op in [:+, :-]
    @eval begin
        function Base.$op(p1::Parameter{Array{T1, N}}, p2::Parameter{Array{T2, N}}) where {N, T1 <: Number, T2 <: Number}
    T = promote_type(T1, T2)
            @boundscheck size(p1) == size(p2) || throw(DimensionMismatch())
            ret = Array{T, N}(size(p1))
            f = let f1 = p1.f, f2 = p2.f
        function ()
                    map!($op, ret, f1(), f2())
                    ret
                end
            end
            Parameter{Array{T, N}}(size(p1), f)
        end
    end
end

function Base.:*(pscalar::Parameter{T1}, parray::Parameter{A}) where {T1 <: Number, T2 <: Number, N, A <: Array{T2, N}}
    T = promote_type(T1, T2)
    ret = Array{T, N}(size(parray))
    f = let fscalar = pscalar.f, farray = parray.f, ret = ret
        function ()
            scalar = fscalar()
            array = farray()
            @boundscheck size(array) == size(ret) || throw(DimensionMismatch())
            @inbounds map!(x -> scalar * x, ret, array)
            ret
        end
    end
    Parameter{Array{T, N}}(size(parray), f)
end
Base.:*(pvector::Parameter{Vector{T2}}, pscalar::Parameter{T1}) where {T1 <: Number, T2 <: Number} = pscalar * pvector


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
        Base.show(io::IO, ::$Fun{T}) where {T} = print(io, "$(string($Fun)){$T, …}(…)")
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

function (f::ScalarLinearFunction{T1})(vals::Dict{Variable, T2}) where {T1, T2}
    T = promote_type(T1, T2)
    ret = zero(T)
    for term in f.terms()
        ret += term.coeff * vals[term.var]
    end
    ret
end

function (f::QuadraticForm{T1})(vals::Dict{Variable, T2}) where {T1, T2}
    T = promote_type(T1, T2)
    ret = zero(T)
    for term in f.terms()
        ret += term.coeff * vals[term.rowvar] * vals[term.colvar]
    end
    ret
end

struct VectorLinearFunction{T, F}
    length::Int
    terms::F
end
VectorLinearFunction{T}(l::Int, terms::F) where {T, F} = VectorLinearFunction{T, F}(l, terms)
Base.length(fun::VectorLinearFunction) = fun.length
Base.size(fun::VectorLinearFunction) = (fun.length,)
Base.show(io::IO, f::VectorLinearFunction{T}) where {T} = print(io, "VectorLinearFunction{$T, …}(…) with output dimension $(length(f))")
Base.similar(fun::VectorLinearFunction, ::Type{T}, terms::F) where {T, F} = VectorLinearFunction{T, F}(length(fun), terms)
Base.similar(fun::VectorLinearFunction{T}, terms::F) where {T, F} = VectorLinearFunction{T, F}(length(fun), terms)

function (f::VectorLinearFunction{T1})(vals::Dict{Variable, T2}) where {T1, T2}
    T = promote_type(T1, T2)
    ret = zeros(T, length(f))
    for term in f.terms()
        ret[term.row] += term.scalarterm.coeff * vals[term.scalarterm.var]
    end
    ret
end

for (Fun, Term) in [
        (:ScalarLinearFunction, :ScalarLinearTerm),
        (:VectorLinearFunction, :VectorLinearTerm),
        (:QuadraticForm, :ScalarQuadraticTerm)]

    WrappedFun = Symbol(:Wrapped, Fun)
    @eval begin
        wrap(l::$Fun{T}) where {T} = similar(l, FunctionWrapper{Vector{$Term{T}}, Tuple{}}(l.terms))
        const $WrappedFun{T} = $Fun{T, FunctionWrapper{Vector{$Term{T}}, Tuple{}}}

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

        function Base.:*(c::Parameter{T1}, l::$Fun{T2}) where {T1 <: Number, T2}
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

        Base.:*(l::$Fun, c::Parameter{<:Number}) = c * l
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

function Base.:*(c::Parameter{T}, x::Vector{Variable}) where T <: Number
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
Base.:*(x::Vector{Variable}, c::Parameter{<:Number}) = c * x

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
for (AffineFunction, Linear, Constant, N) in [
        (:ScalarAffineFunction, :ScalarLinearFunction, :(Parameter{R} where R<:Real), 0),
        (:VectorAffineFunction, :VectorLinearFunction, :(Parameter{Vector{R}} where R<:Real), 1)]
    WrappedFun = Symbol(:Wrapped, AffineFunction)
    WrappedLinear = Symbol(:Wrapped, Linear)

    @eval begin
        struct $AffineFunction{T, L <: $Linear{T}, C <: $Constant{T}}
            linear::L
            constant::C

            function $AffineFunction(linear::L, constant::C) where {T, L <: $Linear{T}, C <: $Constant{T}}
                @boundscheck size(linear) == size(constant) || throw(DimensionMismatch())
                new{T, L, C}(linear, constant)
            end
        end

        (f::$AffineFunction)(vals::Dict{Variable}) = f.linear(vals) + f.constant()
        wrap(a::$AffineFunction) = $AffineFunction(wrap(a.linear), wrap(a.constant))
        const $WrappedFun{T} = $AffineFunction{T, $WrappedLinear{T}, $Constant{T, $N, FunctionWrapper{T, Tuple{}}}}

        Base.:+(a::$AffineFunction) = a
        Base.:-(a::$AffineFunction) = $AffineFunction(-a.linear, -a.constant)
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
ScalarAffineFunction(constant::Parameter{T}) where {T <: Number} = ScalarAffineFunction(zero(ScalarLinearFunction{T}), constant)
ScalarAffineFunction(constant::Number) = ScalarAffineFunction(Parameter(constant))
Base.zero(::Type{ScalarAffineFunction{T}}) where {T} = ScalarAffineFunction(Parameter(zero(T)))
Base.show(io::IO, f::ScalarAffineFunction{T}) where {T} = print(io, "ScalarAffineFunction{$T, …}(…)")

# VectorAffineFunction-specific
VectorAffineFunction(linear::VectorLinearFunction{T}) where {T} = VectorAffineFunction(linear, Parameter(zeros(size(linear))))
function VectorAffineFunction(constant::Parameter{Vector{T}}) where {T}
    ret = Vector{VectorLinearTerm{T}}()
    terms = () -> ret
    linear = VectorLinearFunction{T}(size(constant)[1], terms)
    VectorAffineFunction(linear, constant)
end
VectorAffineFunction(constant::Vector{<:Number}) = VectorAffineFunction(Parameter(constant))
Base.show(io::IO, f::VectorAffineFunction{T}) where {T} = print(io, "VectorAffineFunction{$T, …}(…) with output dimension $(length(f.linear))")


# QuadraticFunction
struct QuadraticFunction{T, Q <: QuadraticForm{T}, A <: ScalarAffineFunction{T}}
    quadratic::Q
    affine::A
end

QuadraticFunction(quadratic::Q, affine::A) where {T, Q<:QuadraticForm{T}, A<:ScalarAffineFunction{T}} = QuadraticFunction{T, Q, A}(quadratic, affine)
QuadraticFunction(quadratic::QuadraticForm{T}) where {T} = QuadraticFunction(quadratic, zero(ScalarAffineFunction{T}))
QuadraticFunction(affine::ScalarAffineFunction{T}) where {T} = QuadraticFunction(zero(QuadraticForm{T}), affine)
QuadraticFunction(linear::ScalarLinearFunction) = QuadraticFunction(ScalarAffineFunction(linear))
QuadraticFunction(constant::Parameter{<:Number}) = QuadraticFunction(ScalarAffineFunction(constant))
QuadraticFunction(constant::Number) = QuadraticFunction(Parameter(constant))

(f::QuadraticFunction)(vals::Dict{Variable}) = f.quadratic(vals) + f.affine(vals)
Base.zero(::Type{QuadraticFunction{T}}) where {T} = QuadraticFunction(Parameter(zero(T)))
Base.show(io::IO, f::QuadraticFunction{T}) where {T} = print(io, "QuadraticFunction{$T, …}(…)")
wrap(q::QuadraticFunction) = QuadraticFunction(wrap(q.quadratic), wrap(q.affine))
const WrappedQuadraticFunction{T} = QuadraticFunction{T, WrappedQuadraticForm{T}, WrappedScalarAffineFunction{T}}

for op in [:+, :-]
    @eval begin
        Base.$op(qform::QuadraticForm, c::Parameter{<:Number}) = QuadraticFunction(qform, ScalarAffineFunction($op(c)))
        Base.$op(c::Parameter{T} where T<:Number, qform::QuadraticForm) = QuadraticFunction($op(qform), c)
        Base.$op(qform::QuadraticForm, l::ScalarLinearFunction) = QuadraticForm(qform, ScalarAffineFunction($op(l)))
        Base.$op(l::ScalarLinearFunction, qform::QuadraticForm) = QuadraticForm(qform, ScalarAffineFunction($op(l)))
    end
end


# Convenience type aliases
const ScalarFunction{T} = Union{ScalarLinearFunction{T}, QuadraticForm{T}, ScalarAffineFunction{T}, QuadraticFunction{T}}
const VectorFunction{T} = Union{VectorLinearFunction{T}, VectorAffineFunction{T}}
const AnyFunction{T} = Union{ScalarFunction{T}, VectorFunction{T}}


# Interaction with values
Base.:*(A::Matrix, x::Vector{Variable}) = Parameter(A) * x
Base.:*(c::Number, x::Vector{Variable}) = Parameter(c) * x

for T in [Number, Array{<:Number}]
    for op in [:+, :-]
        @eval begin
            Base.$op(p::Parameter{<:$T}, c::$T) = $op(p, Parameter(c))
            Base.$op(c::$T, p::Parameter{<:$T}) = $op(Parameter(c), p)
        end
    end
    @eval begin
        Base.:*(p::Parameter{<:$T}, c::Number) = p * Parameter(c)
        Base.:*(c::Number, p::Parameter{<:$T}) = Parameter(c) * p
    end
end

for op in [:+, :-]
    @eval begin
        Base.$op(f::ScalarFunction, c::Number) = $op(f, Parameter(c))
        Base.$op(c::Number, f::ScalarFunction) = $op(Parameter(c), f)
        Base.$op(f::VectorFunction, c::Vector{<:Number}) = $op(f, Parameter(c))
        Base.$op(c::Vector{<:Number}, f::VectorFunction) = $op(Parameter(c), f)
    end
end

Base.:*(f::AnyFunction, c::Number) = f * Parameter(c)
Base.:*(c::Number, f::AnyFunction) = Parameter(c) * f

end
