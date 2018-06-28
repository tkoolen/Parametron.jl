module Functions

export
    Variable,
    LinearTerm,
    QuadraticTerm,
    AffineFunction,
    QuadraticFunction

using Compat
using Compat.LinearAlgebra

@static if VERSION >= v"0.7-"
    const TransposeVector{T, V<:AbstractVector{T}} = Transpose{T, V}
    const AdjointVector{T, V<:AbstractVector{T}} = Adjoint{T, V}
else
    const TransposeVector{T, V<:AbstractVector{T}} = RowVector{T, V}
    const AdjointVector{T, V<:AbstractVector{T}} = RowVector{T, ConjVector{T, V}}
end

const LinearAlgebra = Compat.LinearAlgebra

coefftype(::Type{T}) where {T<:Number} = T

# Variable
struct Variable
    index::Int
end
Base.hash(v::Variable, h::UInt) = hash(v.index, h)
coefftype(::Type{Variable}) = Int


# LinearTerm, QuadraticTerm
struct LinearTerm{T}
    coeff::T
    var::Variable
end
LinearTerm{T}(var::Variable) where {T} = LinearTerm{T}(one(T), var)
Base.show(io::IO, term::LinearTerm) = print(io, term.coeff, " * ", "x", term.var.index)
getcoeff(term::LinearTerm) = term.coeff
setcoeff(term::LinearTerm, coeff) = LinearTerm(coeff, term.var)
Base.:*(coeff::Number, var::Variable) = LinearTerm(coeff, var)
Base.:*(var::Variable, coeff::Number) = LinearTerm(coeff, var)
Base.promote_rule(::Type{LinearTerm{T}}, ::Type{Variable}) where {T} = LinearTerm{T}
Base.convert(::Type{LinearTerm{T}}, var::Variable) where {T} = LinearTerm{T}(var)

struct QuadraticTerm{T}
    coeff::T
    rowvar::Variable
    colvar::Variable
end
Base.show(io::IO, term::QuadraticTerm) = print(io, term.coeff, " * x", term.rowvar.index, " * x", term.colvar.index)
getcoeff(term::QuadraticTerm) = term.coeff
setcoeff(term::QuadraticTerm, coeff) = QuadraticTerm(coeff, term.rowvar, term.colvar)
Base.:*(x::Variable, y::LinearTerm) = QuadraticTerm(y.coeff, x, y.var)
Base.:*(y::LinearTerm, x::Variable) = QuadraticTerm(y.coeff, y.var, x)
Base.:*(x::Variable, y::Variable) = QuadraticTerm(1, x, y)
Base.:*(x::LinearTerm, y::LinearTerm) = QuadraticTerm(x.coeff * y.coeff, x.var, y.var)

for Term in [:LinearTerm, :QuadraticTerm]
    @eval begin
        Base.promote_rule(::Type{$Term{T1}}, ::Type{$Term{T2}}) where {T1, T2} = $Term{promote_type(T1, T2)}
        Base.convert(::Type{$Term{T}}, term::$Term{T}) where {T} = term
        Base.convert(::Type{$Term{T}}, term::$Term) where {T} = setcoeff(term, convert(T, getcoeff(term)))
        @inline coefftype(::Type{$Term{T}}) where {T} = T
        Base.:+(term::$Term) = term
        Base.:-(term::$Term) = setcoeff(term, -getcoeff(term))
        Base.:*(c::Number, term::$Term) = setcoeff(term, c * getcoeff(term))
        Base.:*(term::$Term, c::Number) = c * term
    end
end


# AffineFunction
struct AffineFunction{T}
    linear::Vector{LinearTerm{T}}
    constant::Base.RefValue{T}
end

AffineFunction{T}(linear::Vector{<:LinearTerm}, constant::Number) where {T} =
    AffineFunction{T}(copyto!(Vector{LinearTerm{T}}(undef, length(linear)), linear), Ref(T(constant)))
AffineFunction(linear::Vector{LinearTerm{T}}, constant::S) where {T, S<:Number} =
    AffineFunction{promote_type(T, S)}(linear, constant)

AffineFunction{T}(x::AffineFunction) where {T} = AffineFunction{T}(x.linear, x.constant[])
AffineFunction{T}(x::LinearTerm) where {T} = AffineFunction{T}([LinearTerm{T}(x)], Ref(zero(T)))
AffineFunction{T}(x::Variable) where {T} = AffineFunction{T}(LinearTerm{T}(x))
AffineFunction{T}(x::Number) where {T} = AffineFunction{T}(LinearTerm{T}[], Ref(T(x)))

AffineFunction(x::AffineFunction{T}) where {T} = AffineFunction{T}(x)
AffineFunction(x::LinearTerm{T}) where {T} = AffineFunction{T}(x)
AffineFunction(x::T) where {T<:Number} = AffineFunction{T}(x)

coefftype(::Type{AffineFunction{T}}) where {T} = T

Base.:(==)(x::AffineFunction, y::AffineFunction) = x.linear == y.linear && x.constant[] == y.constant[]
Base.isequal(x::AffineFunction, y::AffineFunction) = isequal(x.linear, y.linear) && isequal(x.constant[], y.constant[])
Base.hash(x::AffineFunction, h::UInt) = (h = hash(x.linear, h); hash(x.constant, h))

Base.zero(::Type{AffineFunction{T}}) where {T} = AffineFunction(LinearTerm{T}[], Ref(zero(T)))
zero!(f::AffineFunction) = (empty!(f.linear); f.constant[] = 0; f)

Base.r_promote_type(::typeof(+), ::Type{LinearTerm{T}}) where {T} = AffineFunction{T}
Base.convert(::Type{AffineFunction{T}}, x::Number) where {T} = AffineFunction(LinearTerm{T}[], convert(T, x))
Base.convert(::Type{AffineFunction{T}}, x::LinearTerm) where {T} = AffineFunction([convert(LinearTerm{T}, x)], zero(T))
Base.convert(::Type{AffineFunction{T}}, x::Variable) where {T} = AffineFunction{T}(x)

function Base.show(io::IO, f::AffineFunction)
    for term in f.linear
        print(io, term, " + ")
    end
    print(io, f.constant[])
end

function (f::AffineFunction{T})(vals::Associative{Variable, S}) where {T, S}
    R′ = Base.promote_op(*, T, S)
    R = Base.promote_op(+, R′, R′)
    ret = convert(R, f.constant[])
    for term in f.linear
        ret += term.coeff * vals[term.var]
    end
    ret
end


# QuadraticFunction
struct QuadraticFunction{T}
    quadratic::Vector{QuadraticTerm{T}}
    affine::AffineFunction{T}
end

QuadraticFunction(quadratic::Vector{QuadraticTerm{T}}, affine::AffineFunction{S}) where {T, S} =
    QuadraticFunction{promote_type(T, S)}(quadratic, affine)

QuadraticFunction{T}(x::QuadraticFunction) where {T} =
    QuadraticFunction{T}(copyto!(Vector{QuadraticTerm{T}}(undef, length(x.quadratic)), x.quadratic), AffineFunction{T}(x.affine))
QuadraticFunction{T}(x::AffineFunction) where {T} = QuadraticFunction{T}(QuadraticTerm{T}[], AffineFunction{T}(x))
QuadraticFunction{T}(x::QuadraticTerm) where {T} = QuadraticFunction([QuadraticTerm{T}(x)], zero(AffineFunction{T}))
QuadraticFunction{T}(x::LinearTerm) where {T} = QuadraticFunction{T}(QuadraticTerm{T}[], AffineFunction{T}(x))
QuadraticFunction{T}(x::Variable) where {T} = QuadraticFunction(AffineFunction{T}(x))
QuadraticFunction{T}(x::Number) where {T} = QuadraticFunction{T}(AffineFunction{T}(x))

QuadraticFunction(x::QuadraticFunction{T}) where {T} = QuadraticFunction{T}(x)
QuadraticFunction(x::AffineFunction{T}) where {T} = QuadraticFunction(QuadraticTerm{T}[], x)
QuadraticFunction(x::QuadraticTerm{T}) where {T} = QuadraticFunction([x], zero(AffineFunction{T}))
QuadraticFunction(x::LinearTerm{T}) where {T<:Number} = QuadraticFunction(QuadraticTerm{T}[], AffineFunction{T}(x))
QuadraticFunction(x::T) where {T<:Number} = QuadraticFunction(QuadraticTerm{T}[], AffineFunction{T}(x))

coefftype(::Type{QuadraticFunction{T}}) where {T} = T

Base.:(==)(x::QuadraticFunction, y::QuadraticFunction) = x.quadratic == y.quadratic && x.affine == y.affine
Base.isequal(x::QuadraticFunction, y::QuadraticFunction) = isequal(x.quadratic, y.quadratic) && isequal(x.affine, y.affine)
Base.hash(x::QuadraticFunction, h::UInt) = (h = hash(x.quadratic, h); hash(x.affine, h))

Base.zero(::Type{QuadraticFunction{T}}) where {T} = QuadraticFunction(QuadraticTerm{T}[], zero(AffineFunction{T}))
zero!(f::QuadraticFunction) = (empty!(f.quadratic); zero!(f.affine); f)

Base.convert(::Type{QuadraticFunction{T}}, x::QuadraticTerm) where {T} =
    QuadraticFunction([convert(QuadraticTerm{T}, x)], zero(AffineFunction{T}))


function Base.show(io::IO, f::QuadraticFunction)
    for term in f.quadratic
        print(io, term, " + ")
    end
    print(io, f.affine)
end

function (f::QuadraticFunction{T})(vals::Associative{Variable, S}) where {T, S}
    ret = f.affine(vals)
    for term in f.quadratic
        ret += term.coeff * vals[term.rowvar] * vals[term.colvar]
    end
    ret
end


# copyto!
Compat.copyto!(f::AffineFunction, x::Number) = (zero!(f); f.constant[] = x; f)
Compat.copyto!(f::AffineFunction, x::LinearTerm) = (zero!(f); push!(f.linear, x); f)
Compat.copyto!(f::AffineFunction{T}, x::Variable) where {T} = copyto!(f, LinearTerm{T}(x))
function Compat.copyto!(f::AffineFunction, x::AffineFunction)
    resize!(f.linear, length(x.linear))
    copyto!(f.linear, x.linear)
    f.constant[] = x.constant[]
    f
end
function Compat.copyto!(f::QuadraticFunction, x::Union{<:Number, <:LinearTerm, Variable, <:AffineFunction})
     empty!(f.quadratic)
     copyto!(f.affine, x)
     f
end
Compat.copyto!(f::QuadraticFunction, x::QuadraticTerm) = (zero!(f); push!(f.quadratic, x); f)
function Compat.copyto!(f::QuadraticFunction, x::QuadraticFunction)
    resize!(f.quadratic, length(x.quadratic))
    copyto!(f.quadratic, x.quadratic)
    copyto!(f.affine, x.affine)
    f
end


# add!
add!(f::AffineFunction, x::Number) = (f.constant[] += x; f)
add!(f::AffineFunction{T}, x::Variable) where {T} = add!(f, LinearTerm{T}(x))
add!(f::AffineFunction, x::LinearTerm) = (push!(f.linear, x); f)
add!(f::AffineFunction, x::AffineFunction) = (append!(f.linear, x.linear); f.constant[] += x.constant[]; f)

add!(f::QuadraticFunction, x::Union{<:Number, <:LinearTerm, Variable, <:AffineFunction}) = (add!(f.affine, x); f)
add!(f::QuadraticFunction, x::QuadraticTerm) = (push!(f.quadratic, x); f)
add!(f::QuadraticFunction, x::QuadraticFunction) = (append!(f.quadratic, x.quadratic); add!(f.affine, x.affine); f)

add!(dest, x, y) = (copyto!(dest, x); add!(dest, y); dest)


# subtract!
subtract!(f::AffineFunction, x::Number) = (f.constant[] -= x; f)
subtract!(f::AffineFunction{T}, x::Variable) where {T} = subtract!(f, LinearTerm{T}(x))
subtract!(f::AffineFunction, x::LinearTerm) = add!(f, -x)
function subtract!(f::AffineFunction, x::AffineFunction)
    offset = length(f.linear)
    resize!(f.linear, offset + length(x.linear))
    @inbounds for i in eachindex(x.linear)
        f.linear[offset + i] = -x.linear[i]
    end
    f.constant[] -= x.constant[]
    f
end

subtract!(f::QuadraticFunction, x::Number) = (subtract!(f.affine, x); f)
subtract!(f::QuadraticFunction, x::LinearTerm) = (subtract!(f.affine, x); f)
subtract!(f::QuadraticFunction, x::Variable) = (subtract!(f.affine, x); f)
subtract!(f::QuadraticFunction, x::AffineFunction) = (subtract!(f.affine, x); f)
subtract!(f::QuadraticFunction, x::QuadraticTerm) = add!(f, -x)
function subtract!(f::QuadraticFunction, x::QuadraticFunction)
    offset = length(f.quadratic)
    resize!(f.quadratic, offset + length(x.quadratic))
    @inbounds for i in eachindex(x.quadratic)
        f.quadratic[offset + i] = -x.quadratic[i]
    end
    subtract!(f.affine, x.affine)
    f
end

subtract!(dest, x, y) = (copyto!(dest, x); subtract!(dest, y); dest)


# muladd!
function muladd!(dest::AffineFunction, x::AffineFunction, y::Number)
    offset = length(dest.linear)
    resize!(dest.linear, offset + length(x.linear))
    @inbounds for i in eachindex(x.linear)
        dest.linear[offset + i] = x.linear[i] * y
    end
    dest.constant[] += x.constant[] * y
    dest
end
muladd!(dest::AffineFunction, x::Number, y::AffineFunction) = muladd!(dest, y, x)

function muladd!(dest::QuadraticFunction, x::QuadraticFunction, y::Number)
    offset = length(dest.quadratic)
    resize!(dest.quadratic, offset + length(x.quadratic))
    @inbounds for i in eachindex(x.quadratic)
        dest.quadratic[offset + i] = x.quadratic[i] * y
    end
    muladd!(dest.affine, x.affine, y)
    dest
end
muladd!(dest::QuadraticFunction, x::Number, y::QuadraticFunction) = muladd!(dest, y, x)

function muladd!(dest::QuadraticFunction, x::AffineFunction, y::Union{Variable, <:LinearTerm})
    offset = length(dest.quadratic)
    resize!(dest.quadratic, offset + length(x.linear))
    @inbounds for i in eachindex(x.linear)
        dest.quadratic[offset + i] = x.linear[i] * y
    end
    add!(dest.affine, x.constant[] * y)
    dest
end
muladd!(dest::QuadraticFunction, x::Union{Variable, <:LinearTerm}, y::AffineFunction) = muladd!(dest, y, x)

function muladd!(dest::QuadraticFunction, x::AffineFunction, y::AffineFunction)
    xlinear = x.linear
    ylinear = y.linear
    quadoffset = length(dest.quadratic)
    resize!(dest.quadratic, quadoffset + length(xlinear) * length(ylinear))
    k = 1
    @inbounds for i in eachindex(xlinear)
        for j in eachindex(ylinear)
            dest.quadratic[quadoffset + k] = xlinear[i] * ylinear[j]
            k += 1
        end
    end
    destaffine = dest.affine
    xconst = x.constant[]
    yconst = y.constant[]
    linoffset = length(destaffine.linear)
    resize!(destaffine.linear, linoffset + length(xlinear) + length(ylinear))
    k = linoffset + 1
    @inbounds for i in eachindex(xlinear)
        destaffine.linear[k] = xlinear[i] * yconst
        k += 1
    end
    @inbounds for i in eachindex(ylinear)
        destaffine.linear[k] = ylinear[i] * xconst
        k += 1
    end
    destaffine.constant[] += xconst * yconst
    dest
end

mul!(dest::Union{<:AffineFunction, <:QuadraticFunction}, x, y) = (zero!(dest); muladd!(dest, x, y))

# Operators
for (op, fun!) in [(:+, add!), (:-, subtract!)]
    @eval begin
        Base.$op(x::LinearTerm{T}, y::LinearTerm{T}) where {T} = AffineFunction{T}([x, $op(y)], zero(T))
        Base.$op(x::LinearTerm, y::LinearTerm) = +(promote(x, y)...)
        Base.$op(x::Variable, y::Variable) = $op(LinearTerm{Int}(x), LinearTerm{Int}(y))
        Base.$op(x::LinearTerm, y::Number) = AffineFunction([x], $op(y))
        Base.$op(x::Number, y::LinearTerm) = AffineFunction([$op(y)], x)
        Base.$op(x::Variable, y::T) where {T<:Number} = $op(LinearTerm{T}(x), y)
        Base.$op(x::T, y::Variable) where {T<:Number} = $op(x, LinearTerm{T}(y))
        Base.$op(x::AffineFunction{T}, y::AffineFunction{S}) where {T, S} =
            $fun!(AffineFunction{promote_type(T, S)}(x), y)
        Base.$op(x::AffineFunction{T}, y::S) where {T, S<:Union{Number, Variable, LinearTerm}} =
            $fun!(AffineFunction{promote_type(T, coefftype(S))}(x), y)
        Base.$op(x::T, y::AffineFunction{S}) where {T<:Union{Number, Variable, LinearTerm}, S} =
            $fun!(AffineFunction{promote_type(coefftype(T), S)}(x), y)

        Base.$op(x::QuadraticTerm{T}, y::QuadraticTerm{T}) where {T} = QuadraticFunction([x, $op(y)], zero(AffineFunction{T}))
        Base.$op(x::QuadraticTerm, y::QuadraticTerm) = +(promote(x, y)...)
        Base.$op(x::QuadraticTerm{T}, y::LinearTerm{T}) where {T} = $fun!(QuadraticFunction(x), y)
        Base.$op(x::LinearTerm{T}, y::QuadraticTerm{T}) where {T} = $fun!(QuadraticFunction(x), y)
        Base.$op(x::LinearTerm{T}, y::QuadraticTerm{S}) where {T, S} = (R = promote_type(T, S); $op(LinearTerm{R}(x), QuadraticTerm{R}(y)))
        Base.$op(x::QuadraticTerm{T}, y::LinearTerm{S}) where {T, S} = (R = promote_type(T, S); $op(QuadraticTerm{R}(x), LinearTerm{R}(y)))
        Base.$op(x::QuadraticTerm{T}, y::Variable) where {T} = $op(x, LinearTerm{T}(y))
        Base.$op(x::Variable, y::QuadraticTerm{T}) where {T} = $op(LinearTerm{T}(x), y)
        Base.$op(x::QuadraticFunction{T}, y::QuadraticFunction{S}) where {T, S} =
            $fun!(QuadraticFunction{promote_type(T, S)}(x), y)
        Base.$op(x::QuadraticFunction{T}, y::S) where {T, S<:Union{Number, Variable, LinearTerm, QuadraticTerm, AffineFunction}} =
            $fun!(QuadraticFunction{promote_type(T, coefftype(S))}(x), y)
        Base.$op(x::T, y::QuadraticFunction{S}) where {T<:Union{Number, Variable, LinearTerm, QuadraticTerm, AffineFunction}, S} =
            $fun!(QuadraticFunction{promote_type(coefftype(T), S)}(x), y)
    end
end

Base.:*(x::AffineFunction{T}, y::Variable) where {T} = muladd!(zero(QuadraticFunction{T}), x, y)
Base.:*(x::Variable, y::AffineFunction{T}) where {T} = muladd!(zero(QuadraticFunction{T}), x, y)
Base.:*(x::AffineFunction{T}, y::Union{LinearTerm{S}, AffineFunction{S}}) where {T, S} = muladd!(zero(QuadraticFunction{promote_type(T, S)}), x, y)
Base.:*(x::Union{LinearTerm{S}, AffineFunction{S}}, y::AffineFunction{T}) where {T, S} = muladd!(zero(QuadraticFunction{promote_type(T, S)}), x, y)
Base.:*(x::AffineFunction{T}, y::AffineFunction{S}) where {T, S} = muladd!(zero(QuadraticFunction{promote_type(T, S)}), x, y)
Base.:*(x::AffineFunction{T}, y::S) where {T, S<:Number} = muladd!(zero(AffineFunction{promote_type(T, S)}), x, y)
Base.:*(x::S, y::AffineFunction{T}) where {T, S<:Number} = muladd!(zero(AffineFunction{promote_type(T, S)}), x, y)
Base.:*(x::QuadraticFunction{T}, y::S) where {T, S<:Number} = muladd!(zero(QuadraticFunction{promote_type(T, S)}), x, y)
Base.:*(x::S, y::QuadraticFunction{T}) where {T, S<:Number} = muladd!(zero(QuadraticFunction{promote_type(T, S)}), x, y)


# Number-like interface
const SimpleQPFunctions = Union{Variable, <:LinearTerm, <:QuadraticTerm, <:AffineFunction, <:QuadraticFunction}
Base.transpose(x::SimpleQPFunctions) = x
Base.dot(x::SimpleQPFunctions, y::SimpleQPFunctions) = x * y
Base.zero(::T) where {T<:SimpleQPFunctions} = zero(T)
Base.one(::T) where {T<:SimpleQPFunctions} = one(T)
if VERSION >= v"0.7-"
    Base.adjoint(x::SimpleQPFunctions) = x
end


# Array operations
# TODO: reduce code duplication
function vecdot!(dest, x::AbstractVector, y::AbstractVector)
    # fallback
    vecdot(x, y)
end

function vecdot!(dest::AffineFunction,
        x::AbstractVector{<:Union{<:Number, Variable, <:LinearTerm}},
        y::AbstractVector{<:Union{<:Number, Variable, <:LinearTerm}})
    zero!(dest)
    @boundscheck length(x) == length(y) || throw(DimensionMismatch())
    linear = dest.linear
    resize!(linear, length(x))
    @inbounds for i in eachindex(x)
        linear[i] = x[i] * y[i]
    end
    dest
end

function vecdot!(dest::QuadraticFunction,
        x::AbstractVector{<:Union{<:Number, Variable, <:LinearTerm}},
        y::AbstractVector{<:Union{<:Number, Variable, <:LinearTerm}})
    zero!(dest)
    @boundscheck length(x) == length(y) || throw(DimensionMismatch())
    quadratic = dest.quadratic
    resize!(quadratic, length(x))
    @inbounds for i in eachindex(x)
        quadratic[i] = x[i] * y[i]
    end
    dest
end

function quadvecdot!(dest::QuadraticFunction, x::AbstractVector, y::AbstractVector)
    zero!(dest)
    @boundscheck length(x) == length(y) || throw(DimensionMismatch())
    for i in eachindex(x)
        muladd!(dest, x[i], y[i])
    end
    dest
end

vecdot!(dest::QuadraticFunction, x::AbstractVector{<:AffineFunction}, y::AbstractVector{<:Union{Variable, <:LinearTerm}}) =
    quadvecdot!(dest, x, y)
vecdot!(dest::QuadraticFunction, x::AbstractVector{<:Union{Variable, <:LinearTerm}}, y::AbstractVector{<:AffineFunction}) =
    quadvecdot!(dest, x, y)
vecdot!(dest::QuadraticFunction, x::AbstractVector{<:AffineFunction}, y::AbstractVector{<:AffineFunction}) =
    quadvecdot!(dest, x, y)

for (vecfun!, scalarfun!) in [(:vecadd!, :add!), (:vecsubtract!, :subtract!)]
    @eval begin
        function $vecfun!(dest::AbstractVector{AffineFunction{T}}, x::AbstractVector, y::AbstractVector) where T
            n = length(x)
            @boundscheck n == length(y) || throw(DimensionMismatch())
            length(dest) == n || resize!(dest, n)
            @inbounds for i in eachindex(dest)
                zero!(dest[i])
                $scalarfun!(dest[i], x[i], y[i])
            end
            dest
        end
    end
end

function matvecmul!(
        y::AbstractVector{AffineFunction{T}},
        A::AbstractMatrix{T},
        x::AbstractVector{Variable}) where T
    rows, cols = size(A)
    @boundscheck length(y) == rows || throw(DimensionMismatch())
    @boundscheck length(x) == cols || throw(DimensionMismatch())
    @inbounds for row in eachindex(y)
        if isassigned(y, row)
            zero!(y[row])
        else
            y[row] = zero(AffineFunction{T})
        end
        resize!(y[row].linear, cols)
    end
    i = 1
    @inbounds for col in Base.OneTo(cols)
        for row in Base.OneTo(rows)
            y[row].linear[col] = A[i] * x[col]
            i += 1
        end
    end
    y
end

function matvecmul!(
    y::AbstractVector{AffineFunction{T}},
    A::AbstractMatrix{<:Number},
    x::AbstractVector{<:AffineFunction}) where T
    rows, cols = size(A)
    @boundscheck length(y) == rows || throw(DimensionMismatch())
    @boundscheck length(x) == cols || throw(DimensionMismatch())
    @inbounds for row in eachindex(y)
        if isassigned(y, row)
            zero!(y[row])
        else
            y[row] = zero(AffineFunction{T})
        end
    end
    i = 1
    @inbounds for col in Base.OneTo(cols)
        for row in Base.OneTo(rows)
            muladd!(y[row], A[i], x[col])
            i += 1
        end
    end
    y
end


function bilinearmul!(
        dest::QuadraticFunction,
        Q::AbstractMatrix,
        x::Union{TransposeVector{Variable, <:AbstractVector{Variable}}, AdjointVector{Variable, <:AbstractVector{Variable}}},
        y::AbstractVector{Variable})
    @boundscheck size(Q) == (length(x), length(y)) || throw(DimensionMismatch())
    zero!(dest)
    quadratic = dest.quadratic
    resize!(quadratic, length(Q))
    k = 1
    @inbounds for row in eachindex(x)
        xrow = x[row]
        for col in eachindex(y)
            quadratic[k] = QuadraticTerm(Q[k], xrow, y[col])
            k += 1
        end
    end
    dest
end

function Base.:*(A::StridedMatrix{T}, x::StridedVector{Variable}) where {T<:LinearAlgebra.BlasFloat}
    matvecmul!(similar(x, AffineFunction{T}, size(A, 1)) , A, x)
end

function scale!(
    dest::AbstractVector{<:LinearTerm},
    x::Number,
    y::AbstractVector{Variable})
    @boundscheck indices(dest) == indices(y) || throw(DimensionMismatch())
    @inbounds for i in eachindex(dest)
        dest[i] = x * y[i]
    end
    dest
end

function scale!(
    dest::AbstractVector{<:LinearTerm},
    x::AbstractVector{Variable},
    y::Number)
    @boundscheck indices(dest) == indices(x) || throw(DimensionMismatch())
    @inbounds for i in eachindex(dest)
        dest[i] = x[i] * y
    end
    dest
end

function scale!(
    dest::AbstractVector{<:AffineFunction},
    x::Number,
    y::AbstractVector{<:AffineFunction})
    @boundscheck indices(dest) == indices(y) || throw(DimensionMismatch())
    @inbounds for i in eachindex(dest)
        mul!(dest[i], x, y[i])
    end
    dest
end

function scale!(
    dest::AbstractVector{<:AffineFunction},
    x::AbstractVector{<:AffineFunction},
    y::Number)
    @boundscheck indices(dest) == indices(x) || throw(DimensionMismatch())
    @inbounds for i in eachindex(dest)
        mul!(dest[i], x[i], y)
    end
    dest
end

if VERSION >= v"0.7-"
    error("TODO: implement mul! for 0.7")
else
    LinearAlgebra.At_mul_B(A::StridedMatrix{T}, x::StridedVector{Variable}) where {T <: LinearAlgebra.BlasFloat} =
        At_mul_B!(Vector{AffineFunction{T}}(undef, size(A, 2)), transpose(A), x)
    LinearAlgebra.Ac_mul_B(A::StridedMatrix{T}, x::StridedVector{Variable}) where {T <: LinearAlgebra.BlasFloat} =
        Ac_mul_B!(Vector{AffineFunction{T}}(undef, size(A, 2)), adjoint(A), x)

    LinearAlgebra.A_mul_B!(y::StridedVector{AffineFunction{T}}, A::StridedMatrix{T}, x::StridedVector{Variable}) where {T <: LinearAlgebra.BlasFloat} =
        matvecmul!(y, A, x)
    LinearAlgebra.At_mul_B!(y::StridedVector{AffineFunction{T}}, A::StridedMatrix{T}, x::StridedVector{Variable}) where {T <: LinearAlgebra.BlasFloat} =
        matvecmul!(y, transpose(A), x)
    LinearAlgebra.Ac_mul_B!(y::StridedVector{AffineFunction{T}}, A::StridedMatrix{T}, x::StridedVector{Variable}) where {T <: LinearAlgebra.BlasFloat} =
        matvecmul!(y, adjoint(A), x)
end

function LinearAlgebra.vecdot(x::AbstractArray{T}, y::AbstractArray{Variable}) where {T<:Number}
    vecdot!(zero(AffineFunction{T}), x, y)
end

function LinearAlgebra.vecdot(x::AbstractArray{Variable}, y::AbstractArray{T}) where {T<:Number}
    vecdot!(zero(AffineFunction{T}), x, y)
end

function LinearAlgebra.vecdot(
        x::AbstractArray{T},
        y::AbstractArray{S}) where {T <: Union{Variable, <:LinearTerm, <:AffineFunction}, S <: Union{Variable, <:LinearTerm, <:AffineFunction}}
    R = T <: Variable ? (S <: Variable ? Int : coefftype(S)) : coefftype(T)
    vecdot!(zero(QuadraticFunction{R}), x, y)
end

function LinearAlgebra.vecdot(x::AbstractArray{LinearTerm{T}}, y::AbstractArray{Variable}) where T
    vecdot!(zero(QuadraticFunction{T}), x, y)
end

function LinearAlgebra.vecdot(x::AbstractArray{Variable}, y::AbstractArray{LinearTerm{T}}) where T
    vecdot!(zero(QuadraticFunction{T}), x, y)
end

# vcat
function _vcat!(dest::AbstractVector{<:AffineFunction},
                i::Integer)
    @boundscheck i == lastindex(dest) + 1 || throw(DimensionMismatch())
    dest
end

function _vcat!(dest::AbstractVector{<:AffineFunction},
                i::Integer,
                source::AbstractVector{<:AffineFunction},
                remaining::Vararg{<:AbstractVector{<:AffineFunction}, N}) where {N}
    @boundscheck i >= firstindex(dest) && (i + length(source) - 1) <= lastindex(dest) || throw(DimensionMismatch())
    @inbounds for s in source
        copyto!(dest[i], s)
        i += 1
    end
    _vcat!(dest, i, remaining...)
end

function vcat!(y::AbstractVector{<:AffineFunction},
               args::Vararg{<:AbstractVector{<:AffineFunction}, N}) where N
    @inbounds for i in eachindex(y)
        zero!(y[i])
    end
    _vcat!(y, 1, args...)
    y
end


end
