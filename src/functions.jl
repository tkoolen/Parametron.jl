module Functions

export
    Variable,
    LinearTerm,
    QuadraticTerm,
    NumberParameter,
    ArrayParameter,
    LinearFunction,
    QuadraticForm,
    AffineFunction


# export
#     LinearFunction,
#     VectorLinearFunction,
#     QuadraticForm,
#     ScalarAffineFunction,
#     VectorAffineFunction,
#     QuadraticFunction

# export
#     WrappedParameter,
#     WrappedScalarLinearFunction,
#     WrappedVectorLinearFunction,
#     WrappedQuadraticForm,
#     WrappedScalarAffineFunction,
#     WrappedVectorAffineFunction,
#     WrappedQuadraticFunction

# export
#     wrap

using Compat
import FunctionWrappers: FunctionWrapper

dimstring(d::Tuple{Vararg{Int}}) = isempty(d) ? "Scalar" : length(d) == 1 ? "$(d[1])-element" : join(map(string, d), '×')

# Const
struct Const{T}
    val::T
end
(c::Const)() = c.val


# LazyMap
struct LazyMap{T, F, A}
    f::F
    dest::Vector{T}
    argfun::A
end

function (m::LazyMap)()
    arg = m.argfun()
    dest = m.dest
    resize!(dest, length(arg))
    @inbounds map!(m.f, dest, arg)
    dest
end

lazymap(f, dest, argfun) = LazyMap(f, dest, argfun)
lazymap(f, dest, argfun::Const) = Const(LazyMap(f, dest, argfun)())


# LazyOpCat
struct LazyOpCat{T, O, F1, F2}
    dest::Vector{T}
    op::O
    arg1fun::F1
    arg2fun::F2
end

function (opcat::LazyOpCat)()
    arg1 = opcat.arg1fun()
    arg2 = opcat.arg2fun()
    dest = opcat.dest
    op = opcat.op
    resize!(dest, length(arg1) + length(arg2))
    i = 1
    @inbounds for term in arg1
        dest[i] = term
        i += 1
    end
    @inbounds for term in arg2
        dest[i] = op(term)
        i+= 1
    end
    dest
end

lazyopcat(dest, op, arg1fun, arg2fun) = LazyOpCat(dest, op, arg1fun, arg2fun)
lazyopcat(dest, op, arg1fun::Const, arg2fun::Const) = Const(LazyOpCat(dest, op, arg1fun, arg2fun)())


# LazyUnaryOp
struct LazyUnaryOp{O, F}
    op::O
    argfun::F
end

(unaryop::LazyUnaryOp)() = unaryop.op(unaryop.argfun())

lazyscalarop(op, argfun) = LazyUnaryOp(op, argfun)
lazyscalarop(op, argfun::Const) = Const(LazyUnaryOp(op, argfun)())


# LazyBinaryOp
struct LazyBinaryOp{O, F1, F2}
    op::O
    argfun1::F1
    argfun2::F2
end

(binaryop::LazyBinaryOp)() = binaryop.op(binaryop.argfun1(), binaryop.argfun2())

lazyscalarop(op, argfun1, argfun2) = LazyBinaryOp(op, argfun1, argfun2)
lazyscalarop(op, argfun1::Const, argfun2::Const) = Const(LazyBinaryOp(op, argfun1, argfun2)())

# Variable
struct Variable <: Number
    index::Int
end
Base.hash(v::Variable, h::UInt) = hash(v.index, h)
Base.one(::Type{Variable}) = 1


# LinearTerm, QuadraticTerm
struct LinearTerm{T <: Number} <: Number
    coeff::T
    var::Variable
end
LinearTerm{T}(var::Variable) where {T <: Number} = LinearTerm{T}(one(T), var)
LinearTerm(var::Variable) = LinearTerm(1, var)
getcoeff(term::LinearTerm) = term.coeff
setcoeff(term::LinearTerm, coeff) = LinearTerm(coeff, term.var)
Base.:*(coeff::Number, var::Variable) = LinearTerm(coeff, var)
Base.:*(var::Variable, coeff::Number) = LinearTerm(coeff, var)
Base.promote_rule(::Type{LinearTerm{T}}, ::Type{Variable}) where {T} = LinearTerm{T}
Base.convert(::Type{LinearTerm{T}}, var::Variable) where {T} = LinearTerm{T}(var)

struct QuadraticTerm{T <: Number}
    coeff::T
    rowvar::Variable
    colvar::Variable
end
getcoeff(term::QuadraticTerm) = term.coeff
setcoeff(term::QuadraticTerm, coeff) = QuadraticTerm(coeff, term.rowvar, term.colvar)
Base.:*(x::Variable, linearterm::LinearTerm) = QuadraticTerm(linearterm.coeff, x, linearterm.var)
Base.:*(linearterm::LinearTerm, x::Variable) = QuadraticTerm(linearterm.coeff, linearterm.var, x)
Base.:*(x::Variable, y::Variable) = QuadraticTerm(1, x, y)

for Term in [:LinearTerm, :QuadraticTerm]
    @eval begin
        Base.promote_rule(::Type{$Term{T1}}, ::Type{$Term{T2}}) where {T1 <: Number, T2 <: Number} = $Term{promote_type(T1, T2)}
        Base.convert(::Type{$Term{T}}, term::$Term{T}) where {T} = term
        Base.convert(::Type{$Term{T}}, term::$Term) where {T} = setcoeff(term, convert(T, getcoeff(term)))
        Base.:+(term::$Term) = term
        Base.:-(term::$Term) = setcoeff(term, -getcoeff(term))
        Base.:*(c::Number, term::$Term) = setcoeff(term, c * getcoeff(term))
        Base.:*(term::$Term, c::Number) = c * term
    end
end

# Base.dot(c::Number, var::Variable) = c * var
# Base.dot(var::Variable, c::Number) = var * c
# Base.dot(x::Variable, y::Variable) = x * y
# Base.dot(l::LinearTerm, x::Variable) = l * x
# Base.dot(x::Variable, l::LinearTerm) = x * l

# Idea: coloring
# Idea: NumberParam, ArrayParam
# Idea: terms functions for matrix param * var vector: distinguish between first and last index
# Idea: LazyMatMul
# Idea: LazyDot

# LinearFunction, QuadraticForm
for (Fun, Term) in [(:LinearFunction, :LinearTerm), (:QuadraticForm, :QuadraticTerm)]
    @eval begin
        struct $Fun{T <: Number, F} <: Number
            terms::F
        end

        $Fun{T}(terms::F) where {T<:Number, F} = $Fun{T, F}(terms)
        $Fun(term::$Term{T}) where {T} = $Fun{T}(Const([term]))

        Base.zero(::Type{<:$Fun{T}}) where {T} = $Fun{T}(Const($Term{T}[]))
        Base.zero(::$Fun{T}) where {T} = zero($Fun{T})
        Base.zero(::Type{$Term{T}}) where {T} = zero($Fun{T})
        Base.zero(::$Term{T}) where {T} = zero($Fun{T})

        Base.show(io::IO, f::$Fun{T, N}) where {T, N} = print(io, "$(string($Fun)){$T,…}(…)")

        Base.promote_rule(::Type{F1}, ::Type{F2}) where {T1, T2, F1<:$Fun{T1}, F2<:$Fun{T2}} = $Fun{promote_type(T1, T2)} # TODO: check if F1, F2 necessary
        Base.convert(::Type{$Fun{T}}, f::$Fun{T}) where {T<:Number} = f
        Base.convert(::Type{$Fun{T}}, f::$Fun) where {T<:Number} = $Fun{T}(lazymap(identity, $Term{T}[], f.terms))

        Base.promote_rule(::Type{F}, ::Type{T}) where {T1<:Number, T2<:Number, F<:$Fun{T1}, T<:$Term{T2}} = $Fun{promote_type(T1, T2)} # TODO: check if F, T necessary
        Base.convert(::Type{$Fun{T}}, term::$Term) where {T<:Number} = $Fun($Term{T}(term))

        Base.:+(f::$Fun) = f
        Base.:-(f::$Fun{T}) where {T} = $Fun{T}(lazymap(-, $Term{T}[], f.terms))
        Base.:*(f::$Fun{T}, s::Number) where {T} = $Fun{T}(lazymap(term -> term * s, $Term{T}[], f.terms))
        Base.:*(s::Number, f::$Fun) = f * s
    end

    for op in [:+, :-]
        @eval begin
            Base.$op(x::$Term{T}, y::$Term{T}) where {T<:Number} = $Fun{T}(Const([x, y]))
            Base.$op(x::$Fun{T}, y::$Fun{T}) where {T<:Number} = $Fun{T}(lazyopcat($Term{T}[], $op, x.terms, y.terms))
            Base.$op(x::$Term, y::$Term) = $op(promote(x, y)...)
            Base.$op(x::$Fun, y::$Fun) = $op(promote(x, y)...)
            Base.$op(f::$Fun, term::$Term) = $op(promote(f, term)...)
            Base.$op(term::$Term, f::$Fun) = $op(promote(term, f)...)
            Base.$op(f::$Fun{T}, x::Variable) where {T<:Number} = $op(f, LinearTerm{T}(x))
            Base.$op(x::Variable, f::$Fun{T}) where {T<:Number} = $op(LinearTerm{T}(x), f)
            Base.$op(f::$Term{T}, x::Variable) where {T<:Number} = $op(f, LinearTerm{T}(x))
            Base.$op(x::Variable, f::$Term{T}) where {T<:Number} = $op(LinearTerm{T}(x), f)
        end
    end
end

Base.:+(x::Variable, y::Variable) = LinearTerm(x) + LinearTerm(y)
Base.:-(x::Variable, y::Variable) = LinearTerm(x) - LinearTerm(y)
Base.:*(l::LinearFunction{T}, x::Variable) where {T} = QuadraticForm{T}(lazymap(term -> term * x, QuadraticTerm{T}[], l.terms))
Base.:*(x::Variable, l::LinearFunction) = l * x
Base.zero(::Type{Variable}) = zero(LinearFunction{Int})

# TODO: necessary for generic matmul; consider removing:
Base.promote_rule(::Type{LinearFunction{T, Const{Vector{LinearTerm{T}}}}}, ::Type{Variable}) where {T<:Number} = LinearFunction{T, Const{Vector{LinearTerm{T}}}}
Base.convert(::Type{LinearFunction{T, Const{Vector{LinearTerm{T}}}}}, var::Variable) where {T<:Number} = LinearFunction(LinearTerm{T}(var))
Base.promote_rule(::Type{LinearFunction{T, Const{Vector{LinearTerm{T}}}}}, ::Type{LinearTerm{T}}) where {T<:Number} = LinearFunction{T, Const{Vector{LinearTerm{T}}}}
Base.convert(::Type{LinearFunction{T, Const{Vector{LinearTerm{T}}}}}, term::LinearTerm{T}) where {T<:Number} = LinearFunction(term)


function (f::LinearFunction{T1})(vals::Dict{Variable, T2}) where {T1, T2}
    T = promote_type(T1, T2)
    ret = zero(T)
    for term in f.terms()
        if term.coeff != 0
            ret += term.coeff * vals[term.var]
        end
    end
    ret
end

function (f::QuadraticForm{T1})(vals::Dict{Variable, T2}) where {T1, T2}
    T = promote_type(T1, T2)
    ret = zero(T)
    for term in f.terms()
        if term.coeff != 0
            ret += term.coeff * vals[term.rowvar] * vals[term.colvar]
        end
    end
    ret
end


# NumberParameter
struct NumberParameter{T<:Number, F} <: Number
    f::F
end
NumberParameter{T}(f::F) where {T, F} = NumberParameter{T, F}(f)
NumberParameter(val::T) where {T<:Number} = NumberParameter{T}(Const(val))
Base.show(io::IO, p::NumberParameter{T}) where {T} = print(io, "NumberParameter{$T,…}(…)")
Base.zero(::Type{<:NumberParameter{T}}) where {T} = NumberParameter(zero(T))
(p::NumberParameter{T})() where {T} = p.f()::T

Base.promote_rule(::Type{<:NumberParameter{T1}}, ::Type{<:NumberParameter{T2}}) where {T1<:Number, T2<:Number} = NumberParameter{promote_type(T1, T2)}
Base.convert(::Type{NumberParameter{T}}, p::NumberParameter{T}) where {T<:Number} = p
Base.convert(::Type{NumberParameter{T}}, p::NumberParameter) where {T<:Number} = NumberParameter{T}(lazyscalarop(x -> convert(T, x), p.f))

Base.promote_rule(::Type{<:NumberParameter{T1}}, ::Type{T2}) where {T1<:Number, T2<:Number} = NumberParameter{promote_type(T1, T2)}
Base.convert(::Type{NumberParameter{T}}, c::Number) where {T} = NumberParameter(convert(T, c))

Base.:+(p::NumberParameter) = p
Base.:-(p::NumberParameter{T}) where {T} = NumberParameter{T}(lazyscalarop(-, p.f))

for op in [:+, :-, :*]
    @eval begin
        Base.$op(p1::NumberParameter{T}, p2::NumberParameter{T}) where {T <: Number} = NumberParameter{T}(lazyscalarop($op, p1.f, p2.f))
        Base.$op(p1::NumberParameter, p2::NumberParameter) = $op(promote(p1, p2)...)
        Base.$op(p1::NumberParameter, p2::Number) = $op(promote(p1, p2)...)
        Base.$op(p1::Number, p2::NumberParameter) = $op(promote(p1, p2)...)
    end
end

Base.:*(x::Variable, p::NumberParameter) = LinearTerm(p, x)
Base.:*(p::NumberParameter, x::Variable) = LinearTerm(p, x)
for Term in [LinearTerm, QuadraticTerm]
    @eval begin
        Base.:*(term::$Term, p::NumberParameter) = setcoeff(term, p * getcoeff(term))
        Base.:*(p::NumberParameter, term::$Term) = setcoeff(term, p * getcoeff(term))
    end
end


# AffineFunction
struct AffineFunction{T, L <: LinearFunction{T}, C <: NumberParameter{T}}
    linear::L
    constant::C
end

AffineFunction(l::L, c::C) where {T, L <: LinearFunction{T}, C <: NumberParameter{T}} = AffineFunction{T, L, C}(l, c)

function AffineFunction(l::LinearFunction{T1}, c::NumberParameter{T2}) where {T1, T2}
    T = promote_type(T1, T2)
    AffineFunction(convert(LinearFunction{T}, l), convert(NumberParameter{T}, c))
end

AffineFunction(l::LinearFunction{T}) where {T} = AffineFunction(l, NumberParameter(zero(T)))
AffineFunction(c::NumberParameter{T}) where {T <: Number} = AffineFunction(zero(LinearFunction{T}), c)
AffineFunction(c::Number) = AffineFunction(NumberParameter(c))
Base.zero(::Type{AffineFunction{T}}) where {T} = AffineFunction(zero(T))
Base.show(io::IO, f::AffineFunction{T}) where {T} = print(io, "AffineFunction{$T, …}(…)")

(f::AffineFunction)(vals::Dict{Variable}) = f.linear(vals) + f.constant()

Base.:+(a::AffineFunction) = a
Base.:-(a::AffineFunction) = AffineFunction(-a.linear, -a.constant)

for op in [:+, :-]
    @eval begin
        Base.$op(l::LinearFunction, c::NumberParameter) = AffineFunction(l, $op(c))
        Base.$op(c::NumberParameter, l::LinearFunction) = AffineFunction($op(l), c)
        Base.$op(a::AffineFunction, c::NumberParameter) = AffineFunction(a.linear, $op(a.constant, c))
        Base.$op(c::NumberParameter, a::AffineFunction) = AffineFunction(a.linear, $op(c, a.constant))
        Base.$op(a::AffineFunction, l::LinearFunction) = AffineFunction($op(a.linear, l), a.constant)
        Base.$op(l::LinearFunction, a::AffineFunction) = AffineFunction($op(l, a.linear), a.constant)
        Base.$op(a1::AffineFunction, a2::AffineFunction) = AffineFunction($op(a1.linear, a2.linear), $op(a1.constant, a2.constant))
        Base.$op(a::AffineFunction, c::Number) = $op(a, NumberParameter(c))
        Base.$op(c::Number, a::AffineFunction) = $op(NumberParameter(c), a)
    end
end


# QuadraticFunction
# struct QuadraticFunction{T, Q <: QuadraticForm{T}, A <: AffineFunction{T}}
#     quadratic::Q
#     affine::A
# end



# ArrayParameter
struct ArrayParameter{T, N, F} <: AbstractArray{T, N}
    f::F
    size::NTuple{N, Int}
end
ArrayParameter{T}(size::NTuple{N, Int}, f::F) where {T, N, F} = ArrayParameter{T, N, F}(f, size)
ArrayParameter(val::Array{T}) where {T} = ArrayParameter{T}(Const(copy(val))), size(val)
Base.size(p::ArrayParameter) = p.size
Base.show(io::IO, p::ArrayParameter{T, N}) where {T, N} = print(io, dimstring(size(p)), " ArrayParameter{$T,$N,…}(…)")
(p::ArrayParameter{T, N})() where {T, N} = p.f()::Array{T, N}


# wrap(p::ArrayParameter{T}) where {T} = ArrayParameter{T}(p.size, FunctionWrapper{T, Tuple{}}(p.f))
# const WrappedParameter{T, N} = ArrayParameter{T, N, FunctionWrapper{T, Tuple{}}}

# Base.:+(p::Parameter) = p

# function Base.:-(p::Parameter{T}) where T <: Number
#     f = let inner = p.f
#         () -> -inner()
#     end
#     Parameter{T}(size(p), f)
# end

# function Base.:-(p::Parameter{T}) where T <: Array
#     ret = T(undef, size(p))
#     f = let inner = p.f, ret = ret
#         function ()
#             val = inner()
#             @boundscheck size(val) == size(ret) || throw(DimensionMismatch())
#             @inbounds map!(-, ret, val)
#             ret
#         end
#     end
#     Parameter{T}(size(p), f)
# end

# for op in [:+, :-, :*]
#     @eval begin
#         function Base.$op(p1::Parameter{T1}, p2::Parameter{T2}) where {T1 <: Number, T2 <: Number}
#             T = promote_type(T1, T2)
#             f = let f1 = p1.f, f2 = p2.f
#                 () -> $op(f1(), f2())
#             end
#             Parameter{T}((), f)
#         end
#     end
# end

# for op in [:+, :-]
#     @eval begin
#         function Base.$op(p1::Parameter{Array{T1, N}}, p2::Parameter{Array{T2, N}}) where {N, T1 <: Number, T2 <: Number}
#     T = promote_type(T1, T2)
#             @boundscheck size(p1) == size(p2) || throw(DimensionMismatch())
#             ret = Array{T, N}(size(p1))
#             f = let f1 = p1.f, f2 = p2.f
#                 function ()
#                     map!($op, ret, f1(), f2())
#                     ret
#                 end
#             end
#             Parameter{Array{T, N}}(size(p1), f)
#         end
#     end
# end

# function Base.:*(pscalar::Parameter{T1}, parray::Parameter{A}) where {T1 <: Number, T2 <: Number, N, A <: Array{T2, N}}
#     T = promote_type(T1, T2)
#     ret = Array{T, N}(size(parray))
#     f = let fscalar = pscalar.f, farray = parray.f, ret = ret
#         function ()
#             scalar = fscalar()
#             array = farray()
#             @boundscheck size(array) == size(ret) || throw(DimensionMismatch())
#             @inbounds map!(x -> scalar * x, ret, array)
#             ret
#         end
#     end
#     Parameter{Array{T, N}}(size(parray), f)
# end
# Base.:*(pvector::Parameter{Vector{T2}}, pscalar::Parameter{T1}) where {T1 <: Number, T2 <: Number} = pscalar * pvector


# struct VectorLinearFunction{T, F}
#     length::Int
#     terms::F
# end
# VectorLinearFunction{T}(l::Int, terms::F) where {T, F} = VectorLinearFunction{T, F}(l, terms)
# Base.length(fun::VectorLinearFunction) = fun.length
# Base.size(fun::VectorLinearFunction) = (fun.length,)
# Base.show(io::IO, f::VectorLinearFunction{T}) where {T} = print(io, "VectorLinearFunction{$T, …}(…) with output dimension $(length(f))")
# Base.similar(fun::VectorLinearFunction, ::Type{T}, terms::F) where {T, F} = VectorLinearFunction{T, F}(length(fun), terms)
# Base.similar(fun::VectorLinearFunction{T}, terms::F) where {T, F} = VectorLinearFunction{T, F}(length(fun), terms)

# function (f::VectorLinearFunction{T1})(vals::Dict{Variable, T2}) where {T1, T2}
#     T = promote_type(T1, T2)
#     ret = zeros(T, length(f))
#     for term in f.terms()
#         ret[term.row] += term.scalarterm.coeff * vals[term.scalarterm.var]
#     end
#     ret
# end

# for (Fun, Term) in [
#         (:LinearFunction, :LinearTerm),
#         (:VectorLinearFunction, :VectorLinearTerm),
#         (:QuadraticForm, :QuadraticTerm)]

#     WrappedFun = Symbol(:Wrapped, Fun)
#     @eval begin
#         wrap(l::$Fun{T}) where {T} = similar(l, FunctionWrapper{Vector{$Term{T}}, Tuple{}}(l.terms))
#         const $WrappedFun{T} = $Fun{T, FunctionWrapper{Vector{$Term{T}}, Tuple{}}}

#         Base.:+(l::$Fun) = l

#         function Base.:-(l::$Fun{T}) where T
#             ret = Vector{$Term{T}}()
#             terms = let ltermsfun = l.terms, ret = ret
#                 function ()
#                     lterms = ltermsfun()
#                     resize!(ret, length(lterms))
#                     @inbounds map!(-, ret, lterms)
#                     ret
#                 end
#             end
#             similar(l, terms)
#         end

#         function Base.:+(l1::$Fun{T1}, l2::$Fun{T2}) where {T1, T2}
#             @boundscheck size(l1) == size(l2) || throw(DimensionMismatch())
#             T = promote_type(T1, T2)
#             ret = Vector{$Term{T}}()
#             terms = let l1termsfun = l1.terms, l2termsfun = l2.terms, ret = ret
#                 function ()
#                     l1terms = l1termsfun()
#                     l2terms = l2termsfun()
#                     resize!(ret, length(l1terms) + length(l2terms))
#                     i = 1
#                     @inbounds for term in l1terms
#                         ret[i] = term
#                         i += 1
#                     end
#                     @inbounds for term in l2terms
#                         ret[i] = term
#                         i+= 1
#                     end
#                     ret
#                 end
#             end
#             similar(l1, T, terms)
#         end

#         function Base.:-(l1::$Fun{T1}, l2::$Fun{T2}) where {T1, T2}
#             @boundscheck size(l1) == size(l2) || throw(DimensionMismatch())
#             T = promote_type(T1, T2)
#             ret = Vector{$Term{T}}()
#             terms = let l1termsfun = l1.terms, l2termsfun = l2.terms, ret = ret
#                 function ()
#                     l1terms = l1termsfun()
#                     l2terms = l2termsfun()
#                     resize!(ret, length(l1terms) + length(l2terms))
#                     i = 1
#                     @inbounds for term in l1terms
#                         ret[i] = term
#                         i += 1
#                     end
#                     @inbounds for term in l2terms
#                         ret[i] = -term
#                         i+= 1
#                     end
#                     ret
#                 end
#             end
#             similar(l1, T, terms)
#         end

#         function Base.:*(c::Parameter{T1}, l::$Fun{T2}) where {T1 <: Number, T2}
#             T = promote_type(T1, T2)
#             ret = Vector{$Term{T}}()
#             terms = let c = c, lterms = l.terms, ret = ret
#                 function ()
#                     cval = c()
#                     ltermsval = l.terms()
#                     resize!(ret, length(ltermsval))
#                     @inbounds for i in eachindex(ltermsval)
#                         ret[i] = cval * ltermsval[i]
#                     end
#                     ret
#                 end
#             end
#             similar(l, T, terms)
#         end

#         Base.:*(l::$Fun, c::Parameter{<:Number}) = c * l
#     end
# end

# # VectorLinearFunction-specific
# function Base.:*(A::Parameter{Matrix{T}}, x::Vector{Variable}) where T
#     @boundscheck size(A)[2] == length(x) || throw(DimensionMismatch())
#     ret = Vector{VectorLinearTerm{T}}()
#     terms = let A = A, x = x, ret = ret
#         function ()
#             Aval = A()
#             @boundscheck size(Aval, 2) == length(x) || throw(DimensionMismatch())
#             n = length(Aval)
#             resize!(ret, n)
#             i = 1
#             @inbounds for col = 1 : size(Aval, 2)
#                 for row = 1 : size(Aval, 1)
#                     ret[i] = VectorLinearTerm(row, LinearTerm(Aval[i], x[col]))
#                     i += 1
#                 end
#             end
#             ret
#         end
#     end
#     VectorLinearFunction{T}(size(A)[1], terms)
# end

# function Base.:*(c::Parameter{T}, x::Vector{Variable}) where T <: Number
#     ret = Vector{VectorLinearTerm{T}}()
#     terms = let c = c, x = x, ret = ret
#         function ()
#             cval = c()
#             resize!(ret, length(x))
#             @inbounds for i in eachindex(x)
#                 ret[i] = VectorLinearTerm(i, LinearTerm(cval, x[i]))
#             end
#             ret
#         end
#     end
#     VectorLinearFunction{T}(length(x), terms)
# end
# Base.:*(x::Vector{Variable}, c::Parameter{<:Number}) = c * x

# # QuadraticForm-specific
# function Base.dot(x::Vector{Variable}, l::VectorLinearFunction{T}) where T
#     @boundscheck length(x) == length(l) || throw(DimensionMismatch())
#     ret = Vector{QuadraticTerm{T}}()
#     terms = let x = x, ltermsfun = l.terms, ret = ret
#         function ()
#             linearterms = ltermsfun()
#             resize!(ret, length(linearterms))
#             for i in eachindex(linearterms)
#                 @inbounds linearterm = linearterms[i]
#                 quadraticterm = QuadraticTerm(linearterm.scalarterm.coeff, x[linearterm.row], linearterm.scalarterm.var)
#                 @inbounds ret[i] = quadraticterm
#             end
#             ret
#         end
#     end
#     QuadraticForm{T}(terms)
# end


# # ScalarAffineFunction, VectorAffineFunction
# for (AffineFunction, Linear, Constant, N) in [
#         (:ScalarAffineFunction, :LinearFunction, :(Parameter{R} where R<:Real), 0),
#         (:VectorAffineFunction, :VectorLinearFunction, :(Parameter{Vector{R}} where R<:Real), 1)]
#     WrappedFun = Symbol(:Wrapped, AffineFunction)
#     WrappedLinear = Symbol(:Wrapped, Linear)

#     @eval begin
#         struct $AffineFunction{T, L <: $Linear{T}, C <: $Constant{T}}
#             linear::L
#             constant::C

#             function $AffineFunction(linear::L, constant::C) where {T, L <: $Linear{T}, C <: $Constant{T}}
#                 @boundscheck size(linear) == size(constant) || throw(DimensionMismatch())
#                 new{T, L, C}(linear, constant)
#             end
#         end

#         (f::$AffineFunction)(vals::Dict{Variable}) = f.linear(vals) + f.constant()
#         wrap(a::$AffineFunction) = $AffineFunction(wrap(a.linear), wrap(a.constant))
#         const $WrappedFun{T} = $AffineFunction{T, $WrappedLinear{T}, $Constant{T, $N, FunctionWrapper{T, Tuple{}}}}

#         Base.:+(a::$AffineFunction) = a
#         Base.:-(a::$AffineFunction) = $AffineFunction(-a.linear, -a.constant)
#     end

#     for op in [:+, :-]
#         @eval begin
#             Base.$op(l::$Linear, c::$Constant) = $AffineFunction(l, $op(c))
#             Base.$op(c::$Constant, l::$Linear) = $AffineFunction($op(l), c)
#             Base.$op(a::$AffineFunction, c::$Constant) = $AffineFunction(a.linear, $op(a.constant, c))
#             Base.$op(c::$Constant, a::$AffineFunction) = $AffineFunction(a.linear, $op(c, a.constant))
#             Base.$op(a::$AffineFunction, l::$Linear) = $AffineFunction($op(a.linear, l), a.constant)
#             Base.$op(l::$Linear, a::$AffineFunction) = $AffineFunction($op(l, a.linear), a.constant)
#         end
#     end
# end

# # ScalarAffineFunction-specific
# ScalarAffineFunction(linear::LinearFunction{T}) where {T} = ScalarAffineFunction(linear, Parameter(zero(T)))
# ScalarAffineFunction(constant::Parameter{T}) where {T <: Number} = ScalarAffineFunction(zero(LinearFunction{T}), constant)
# ScalarAffineFunction(constant::Number) = ScalarAffineFunction(Parameter(constant))
# Base.zero(::Type{ScalarAffineFunction{T}}) where {T} = ScalarAffineFunction(Parameter(zero(T)))
# Base.show(io::IO, f::ScalarAffineFunction{T}) where {T} = print(io, "ScalarAffineFunction{$T, …}(…)")

# # VectorAffineFunction-specific
# VectorAffineFunction(linear::VectorLinearFunction{T}) where {T} = VectorAffineFunction(linear, Parameter(zeros(size(linear))))
# function VectorAffineFunction(constant::Parameter{Vector{T}}) where {T}
#     ret = Vector{VectorLinearTerm{T}}()
#     terms = () -> ret
#     linear = VectorLinearFunction{T}(size(constant)[1], terms)
#     VectorAffineFunction(linear, constant)
# end
# VectorAffineFunction(constant::Vector{<:Number}) = VectorAffineFunction(Parameter(constant))
# Base.show(io::IO, f::VectorAffineFunction{T}) where {T} = print(io, "VectorAffineFunction{$T, …}(…) with output dimension $(length(f.linear))")


# # QuadraticFunction
# struct QuadraticFunction{T, Q <: QuadraticForm{T}, A <: ScalarAffineFunction{T}}
#     quadratic::Q
#     affine::A
# end

# QuadraticFunction(quadratic::Q, affine::A) where {T, Q<:QuadraticForm{T}, A<:ScalarAffineFunction{T}} = QuadraticFunction{T, Q, A}(quadratic, affine)
# QuadraticFunction(quadratic::QuadraticForm{T}) where {T} = QuadraticFunction(quadratic, zero(ScalarAffineFunction{T}))
# QuadraticFunction(affine::ScalarAffineFunction{T}) where {T} = QuadraticFunction(zero(QuadraticForm{T}), affine)
# QuadraticFunction(linear::LinearFunction) = QuadraticFunction(ScalarAffineFunction(linear))
# QuadraticFunction(constant::Parameter{<:Number}) = QuadraticFunction(ScalarAffineFunction(constant))
# QuadraticFunction(constant::Number) = QuadraticFunction(Parameter(constant))

# (f::QuadraticFunction)(vals::Dict{Variable}) = f.quadratic(vals) + f.affine(vals)
# Base.zero(::Type{QuadraticFunction{T}}) where {T} = QuadraticFunction(Parameter(zero(T)))
# Base.show(io::IO, f::QuadraticFunction{T}) where {T} = print(io, "QuadraticFunction{$T, …}(…)")
# wrap(q::QuadraticFunction) = QuadraticFunction(wrap(q.quadratic), wrap(q.affine))
# const WrappedQuadraticFunction{T} = QuadraticFunction{T, WrappedQuadraticForm{T}, WrappedScalarAffineFunction{T}}

# for op in [:+, :-]
#     @eval begin
#         Base.$op(qform::QuadraticForm, c::Parameter{<:Number}) = QuadraticFunction(qform, ScalarAffineFunction($op(c)))
#         Base.$op(c::Parameter{T} where T<:Number, qform::QuadraticForm) = QuadraticFunction($op(qform), c)
#         Base.$op(qform::QuadraticForm, l::LinearFunction) = QuadraticForm(qform, ScalarAffineFunction($op(l)))
#         Base.$op(l::LinearFunction, qform::QuadraticForm) = QuadraticForm(qform, ScalarAffineFunction($op(l)))
#     end
# end


# # Convenience type aliases
# const ScalarFunction{T} = Union{LinearFunction{T}, QuadraticForm{T}, ScalarAffineFunction{T}, QuadraticFunction{T}}
# const VectorFunction{T} = Union{VectorLinearFunction{T}, VectorAffineFunction{T}}
# const AnyFunction{T} = Union{ScalarFunction{T}, VectorFunction{T}}


# # Interaction with values
# Base.:*(A::Matrix, x::Vector{Variable}) = Parameter(A) * x
# Base.:*(c::Number, x::Vector{Variable}) = Parameter(c) * x

# for T in [Number, Array{<:Number}]
#     for op in [:+, :-]
#         @eval begin
#             Base.$op(p::Parameter{<:$T}, c::$T) = $op(p, Parameter(c))
#             Base.$op(c::$T, p::Parameter{<:$T}) = $op(Parameter(c), p)
#         end
#     end
#     @eval begin
#         Base.:*(p::Parameter{<:$T}, c::Number) = p * Parameter(c)
#         Base.:*(c::Number, p::Parameter{<:$T}) = Parameter(c) * p
#     end
# end

# for op in [:+, :-]
#     @eval begin
#         Base.$op(f::ScalarFunction, c::Number) = $op(f, Parameter(c))
#         Base.$op(c::Number, f::ScalarFunction) = $op(Parameter(c), f)
#         Base.$op(f::VectorFunction, c::Vector{<:Number}) = $op(f, Parameter(c))
#         Base.$op(c::Vector{<:Number}, f::VectorFunction) = $op(Parameter(c), f)
#     end
# end

# Base.:*(f::AnyFunction, c::Number) = f * Parameter(c)
# Base.:*(c::Number, f::AnyFunction) = Parameter(c) * f

end
