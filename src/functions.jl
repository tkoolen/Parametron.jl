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


# Variable
struct Variable
    index::Int
end
Base.hash(v::Variable, h::UInt) = hash(v.index, h)
Base.one(::Type{Variable}) = 1


# LinearTerm, QuadraticTerm
struct LinearTerm{T}
    coeff::T
    var::Variable
end
LinearTerm{T}(var::Variable) where {T} = LinearTerm{T}(one(T), var)
LinearTerm(var::Variable) = LinearTerm(1, var)
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
AffineFunction{T}(f::AffineFunction) where {T} = AffineFunction(copyto!(Vector{LinearTerm{T}}(undef, length(f.linear)), f.linear), T(f.constant[]))
AffineFunction(f::AffineFunction{T}) where {T} = AffineFunction{T}(f)
AffineFunction(linear::Vector{LinearTerm{T}}, constant::Base.RefValue{S}) where {T, S} = AffineFunction{promote_type(T, S)}(linear, constant)
AffineFunction(linear::Vector{LinearTerm{T}}, constant) where {T} = AffineFunction(linear, Ref(constant))
coefftype(::Type{AffineFunction{T}}) where {T} = T

Base.:(==)(x::AffineFunction, y::AffineFunction) = x.linear == y.linear && x.constant[] == y.constant[]
Base.isequal(x::AffineFunction, y::AffineFunction) = isequal(x.linear, y.linear) && isequal(x.constant[], y.constant[])
Base.hash(x::AffineFunction, h::UInt) = (h = hash(x.linear, h); hash(x.constant, h))

Base.zero(::Type{AffineFunction{T}}) where {T} = AffineFunction(LinearTerm{T}[], Ref(zero(T)))
zero!(f::AffineFunction) = (empty!(f.linear); f.constant[] = 0; f)

Base.r_promote_type(::typeof(+), ::Type{LinearTerm{T}}) where {T} = AffineFunction{T}
Base.convert(::Type{AffineFunction{T}}, x::Number) where {T} = AffineFunction(LinearTerm{T}[], convert(T, x))
Base.convert(::Type{AffineFunction{T}}, x::LinearTerm) where {T} = AffineFunction([convert(LinearTerm{T}, x)], zero(T))

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
QuadraticFunction{T}(f::QuadraticFunction) where {T} = QuadraticFunction(copyto!(Vector{QuadraticTerm{T}}(undef, length(f.quadratic)), f.quadratic), AffineFunction{T}(f.affine))
QuadraticFunction(f::QuadraticFunction{T}) where {T} = QuadraticFunction{T}(f)
coefftype(::Type{QuadraticFunction{T}}) where {T} = T
Base.zero(::Type{QuadraticFunction{T}}) where {T} = QuadraticFunction(QuadraticTerm{T}[], zero(AffineFunction{T}))
zero!(f::QuadraticFunction) = (empty!(f.quadratic); zero!(f.affine); f)

# TODO: conversions

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
add!(f::AffineFunction, x::LinearTerm) = (push!(f.linear, x); f)
add!(f::AffineFunction{T}, x::Variable) where {T} = add!(f, LinearTerm{T}(x))
add!(f::AffineFunction, x::AffineFunction) = (append!(f.linear, x.linear); f.constant[] += x.constant[]; f)

add!(f::QuadraticFunction, x::Union{<:Number, <:LinearTerm, Variable, <:AffineFunction}) = (add!(f.affine, x); f)
add!(f::QuadraticFunction, x::QuadraticTerm) = (push!(f.quadratic, x); f)
add!(f::QuadraticFunction, x::QuadraticFunction) = (append!(f.quadratic, x.quadratic); add!(f.affine, x.affine); f)

add!(dest, x, y) = (copyto!(dest, x); add!(dest, y); dest)


# subtract!
subtract!(f::AffineFunction, x::Number) = (f.constant[] -= x; f)
subtract!(f::AffineFunction, x::LinearTerm) = add!(f, -x)
subtract!(f::AffineFunction{T}, x::Variable) where {T} = subtract!(f, LinearTerm{T}(x))
function subtract!(f::AffineFunction, x::AffineFunction)
    offset = length(f)
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
    offset = length(f)
    resize!(f.quadratic, offset + length(x.quadratic))
    @inbounds for i in eachindex(x.quadratic)
        f.quadratic[offset + i] = -x.quadratic[i]
    end
    subtract!(f.affine, x.affine)
    f
end

subtract!(dest, x, y) = (copyto!(dest, x); subtract!(dest, y); dest)


# multiply!
function multiply!(dest::QuadraticFunction, x::AffineFunction, y::Union{Variable, <:LinearTerm})
    zero!(dest)
    resize!(dest.quadratic, length(x.linear))
    @inbounds for i in eachindex(x.linear)
        dest.quadratic[i] = x.linear[i] * y
    end
    add!(dest.affine, x.constant[] * y)
    dest
end

function multiply!(dest::QuadraticFunction, x::AffineFunction, y::AffineFunction)
    zero!(dest)
    xlinear = x.linear
    ylinear = y.linear
    resize!(dest.quadratic, length(xlinear) * length(ylinear))
    k = 1
    @inbounds for i in eachindex(xlinear)
        for j in eachindex(ylinear)
            dest.quadratic[k] = xlinear[i] * ylinear[j]
            k += 1
        end
    end
    destaffine = dest.affine
    xconst = x.constant[]
    yconst = y.constant[]
    resize!(destaffine.linear, length(xlinear) + length(ylinear))
    k = 1
    @inbounds for i in eachindex(xlinear)
        destaffine.linear[k] = xlinear[i] * yconst
        k += 1
    end
    @inbounds for i in eachindex(ylinear)
        destaffine.linear[k] = ylinear[i] * xconst
        k += 1
    end
    destaffine.constant[] = xconst * yconst
    dest
end

# Operators
for (op, fun!) in [(:+, add!), (:-, subtract!)]
    @eval begin
        Base.$op(x::LinearTerm{T}, y::LinearTerm{T}) where {T} = AffineFunction([x, $op(y)], zero(T))
        Base.$op(x::LinearTerm, y::LinearTerm) = +(promote(x, y)...)
        Base.$op(x::Variable, y::Variable) = $op(LinearTerm{Int}(x), LinearTerm{Int}(y))
        Base.$op(x::AffineFunction{T}, y::AffineFunction{S}) where {T, S} =
            $fun!(AffineFunction{promote_type(T, S)}(x), y)
        Base.$op(x::LinearTerm, y::Number) = AffineFunction([x], $op(y))
        Base.$op(x::Number, y::LinearTerm) = AffineFunction([$op(y)], x)
        Base.$op(x::Variable, y::T) where {T<:Number} = $op(LinearTerm{T}(x), y)
        Base.$op(x::T, y::Variable) where {T<:Number} = $op(x, LinearTerm{T}(y))
        Base.$op(x::AffineFunction, y::Number) = $fun!(AffineFunction(x), y)
        Base.$op(x::Number, y::AffineFunction) = $op(y, x)

        Base.$op(x::QuadraticTerm{T}, y::QuadraticTerm{T}) where {T} = QuadraticFunction([x, $op(y)], zero(AffineFunction{T}))
        Base.$op(x::QuadraticTerm, y::QuadraticTerm) = +(promote(x, y)...)
        Base.$op(x::QuadraticFunction{T}, y::QuadraticFunction{S}) where {T, S} =
            $fun!(QuadraticFunction{promote_type(T, S)}(x), y)
        Base.$op(x::QuadraticFunction{T}, y::Union{S, LinearTerm{S}, QuadraticTerm{S}, AffineFunction{S}}) where {S <: Number, T} =
            $fun!(QuadraticFunction{promote_type(S, T)}(x), y)
        Base.$op(x::Union{S, LinearTerm{S}, QuadraticTerm{S}, AffineFunction{S}}, y::QuadraticFunction{T}) where {S <: Number, T} =
            $fun!(QuadraticFunction{promote_type(S, T)}(y), x)
    end
end

Base.:*(x::AffineFunction{T}, y::Variable) where {T} = multiply!(zero(QuadraticFunction{T}), x, y)
Base.:*(y::Variable, x::AffineFunction{T}) where {T} = multiply!(zero(QuadraticFunction{T}), x, y)
Base.:*(x::AffineFunction{T}, y::Union{AffineFunction{S}, LinearTerm{S}}) where {T, S} = multiply!(zero(QuadraticFunction{promote_type(T, S)}), x, y)
Base.:*(y::Union{AffineFunction{S}, LinearTerm{S}}, x::AffineFunction{T}) where {T, S} = multiply!(zero(QuadraticFunction{promote_type(T, S)}), x, y)
Base.:*(x::AffineFunction{T}, y::AffineFunction{S}) where {T, S} = multiply!(zero(QuadraticFunction{promote_type(T, S)}), x, y)


# Number-like interface
const FunctionsTypes = Union{Variable, <:LinearTerm, <:QuadraticTerm, <:AffineFunction, <:QuadraticFunction}
Base.transpose(x::FunctionsTypes) = x
Base.dot(x::FunctionsTypes, y::FunctionsTypes) = x * y
Base.zero(::T) where {T<:FunctionsTypes} = zero(T)
Base.one(::T) where {T<:FunctionsTypes} = one(T)
if VERSION >= v"0.7-"
    Base.adjoint(x::FunctionsTypes) = x
end


# Array operations
# TODO: reduce code duplication
function vecdot!(dest, x::DenseVector, y::DenseVector)
    # fallback
    vecdot(x, y)
end

function vecdot!(dest::AffineFunction, x::DenseVector, y::DenseVector)
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
        x::DenseVector{<:Union{<:Number, Variable, <:LinearTerm}},
        y::DenseVector{<:Union{<:Number, Variable, <:LinearTerm}})
    zero!(dest)
    @boundscheck length(x) == length(y) || throw(DimensionMismatch())
    quadratic = dest.quadratic
    resize!(quadratic, length(x))
    @inbounds for i in eachindex(x)
        quadratic[i] = x[i] * y[i]
    end
    dest
end

function vecdot!(dest::QuadraticFunction,
        xs::DenseVector{<:AffineFunction},
        ys::DenseVector{<:Union{Variable, <:LinearTerm}})
    zero!(dest)
    @boundscheck length(xs) == length(ys) || throw(DimensionMismatch())
    for i in eachindex(xs)
        x = xs[i]
        y = ys[i]
        xlinear = x.linear
        @inbounds for j in eachindex(xlinear)
            push!(dest.quadratic, xlinear[j] * y)
        end
        push!(dest.affine.linear, x.constant[] * y)
    end
    dest
end
function vecdot!(dest::QuadraticFunction,
    xs::DenseVector{<:Union{Variable, <:LinearTerm}},
    ys::DenseVector{<:AffineFunction})
    vecdot!(dest, ys, xs)
end

function vecdot!(dest::QuadraticFunction, xs::DenseVector{<:AffineFunction}, ys::DenseVector{<:AffineFunction})
    zero!(dest)
    @boundscheck length(xs) == length(ys) || throw(DimensionMismatch())
    for i in eachindex(xs)
        x = xs[i]
        y = ys[i]
        xlinear = x.linear
        ylinear = y.linear
        @inbounds for j in eachindex(xlinear)
            for k in eachindex(ylinear)
                push!(dest.quadratic, xlinear[j] * ylinear[k])
            end
        end
        xconst = x.constant[]
        yconst = y.constant[]
        @inbounds for j in eachindex(xlinear)
            push!(dest.affine.linear, xlinear[j] * yconst)
        end
        @inbounds for j in eachindex(ylinear)
            push!(dest.affine.linear, ylinear[j] * xconst)
        end
        dest.affine.constant[] += xconst * yconst
    end
    dest
end

function vecadd!(dest::DenseVector{AffineFunction{T}}, x::DenseVector, y::DenseVector) where T
    n = length(x)
    @boundscheck n == length(y) || throw(DimensionMismatch())
    resize!(dest, n)
    @inbounds for i in eachindex(dest)
        if isassigned(dest, i)
            zero!(dest[i])
        else
            dest[i] = zero(AffineFunction{T})
        end
        add!(dest[i], x[i], y[i])
    end
    dest
end

function vecsubtract!(dest::DenseVector{AffineFunction{T}}, x::DenseVector, y::DenseVector) where T
    n = length(x)
    @boundscheck n == length(y) || throw(DimensionMismatch())
    resize!(dest, n)
    @inbounds for i in eachindex(dest)
        if isassigned(dest, i)
            zero!(dest[i])
        else
            dest[i] = zero(AffineFunction{T})
        end
        subtract!(dest[i], x[i], y[i])
    end
    dest
end

function affine_matvecmul!(
        y::DenseVector{AffineFunction{T}},
        A::DenseMatrix{T},
        x::DenseVector{Variable}) where {T<:LinearAlgebra.BlasFloat}
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

function quadratic_bilinear_mul!(
        dest::QuadraticFunction,
        Q::DenseMatrix,
        x::TransposeVector{Variable, <:DenseVector{Variable}},
        y::DenseVector{Variable})
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
    affine_matvecmul!(similar(x, AffineFunction{T}, size(A, 1)) , A, x)
end


if VERSION >= v"0.7-"
    # TODO
else
    LinearAlgebra.At_mul_B(A::StridedMatrix{T}, x::StridedVector{Variable}) where {T <: LinearAlgebra.BlasFloat} =
        At_mul_B!(Vector{AffineFunction{T}}(undef, size(A, 2)), transpose(A), x)
    LinearAlgebra.Ac_mul_B(A::StridedMatrix{T}, x::StridedVector{Variable}) where {T <: LinearAlgebra.BlasFloat} =
        Ac_mul_B!(Vector{AffineFunction{T}}(undef, size(A, 2)), adjoint(A), x)

    LinearAlgebra.A_mul_B!(y::StridedVector{AffineFunction{T}}, A::StridedMatrix{T}, x::StridedVector{Variable}) where {T <: LinearAlgebra.BlasFloat} =
        affine_matvecmul!(y, A, x)
    LinearAlgebra.At_mul_B!(y::StridedVector{AffineFunction{T}}, A::StridedMatrix{T}, x::StridedVector{Variable}) where {T <: LinearAlgebra.BlasFloat} =
        affine_matvecmul!(y, transpose(A), x)
    LinearAlgebra.Ac_mul_B!(y::StridedVector{AffineFunction{T}}, A::StridedMatrix{T}, x::StridedVector{Variable}) where {T <: LinearAlgebra.BlasFloat} =
        affine_matvecmul!(y, adjoint(A), x)
end
# TODO: A_mul_B!, etc.

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

# TODO: mul! in 0.7

end
