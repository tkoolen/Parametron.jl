module Functions

export
    Variable,
    Constant,
    LinearTerm,
    QuadraticTerm,
    Scaled,
    Sum,
    AffineFunction,
    QuadraticFunction

import ..SimpleQP: quad

const SparseSymmetric64 = Symmetric{Float64,SparseMatrixCSC{Float64,Int}}

abstract type Fun end

# Variable
struct Variable
    index::Int
end


# Constant
struct Constant <: Fun
    v::Vector{Float64}
end
outputdim(f::Constant) = length(f.v)
@inline isquadratic(::Type{Constant}) = false
(f::Constant)() = f.v


# LinearTerm
struct LinearTerm <: Fun
    A::Matrix{Float64}
    x::Vector{Variable}

    function LinearTerm(A::Matrix{Float64}, x::Vector{Variable})
        size(A, 2) === length(x) || error()
        new(A, x)
    end
end
outputdim(f::LinearTerm) = size(f.A, 1)
@inline isquadratic(::Type{LinearTerm}) = false
(f::LinearTerm)(vals::Associative{Variable, Float64}) = f.A * getindex.(vals, f.x)


# QuadraticTerm
struct QuadraticTerm <: Fun
    Q::SparseSymmetric64
    x::Vector{Variable}
end
quad(Q::SparseSymmetric64, x::Vector{Variable}) = QuadraticTerm(Q, x)
QuadraticTerm() = QuadraticTerm(Symmetric(spzeros(0, 0)), Variable[])
outputdim(f::QuadraticTerm) = 1
@inline isquadratic(::Type{QuadraticTerm}) = true
(f::QuadraticTerm)(vals::Associative{Variable, Float64}) = quad(f.Q, getindex.(vals, f.x))


# Scaled
struct Scaled{T<:Fun} <: Fun
    scalar::Float64
    val::T
    Scaled{T}(scalar::Float64, val::T) where {T<:Fun} = new{T}(scalar, val)
end
Scaled(scalar::Float64, val::T) where {T<:Fun} = Scaled{T}(scalar, val)
outputdim(f::Scaled) = outputdim(f.val)
@inline isquadratic(::Type{Scaled{T}}) where {T<:Fun} = isquadratic(T)
(f::Scaled)(x...) = f.scalar * f.val(x...)


# Sum
struct Sum{T<:Fun} <: Fun
    terms::Vector{T}

    function Sum{T}(terms::Vector{T}) where {T}
        isempty(terms) && throw(ArgumentError())
        m = outputdim(terms[1])
        for i = 2 : length(terms)
            outputdim(terms[i]) == m || throw(DimensionError())
        end
        new{T}(terms)
    end
end
Sum(terms::Vector{T}) where {T<:Fun} = Sum{T}(terms)
outputdim(f::Sum) = isempty(f.terms) ? 0 : outputdim(f.terms[1])
@inline isquadratic(::Type{Sum{T}}) where {T<:Fun} = isquadratic(T)
(f::Sum)(x...) = sum(term -> term(x...), f.terms)


# AffineFunction
struct AffineFunction <: Fun
    linear::Sum{Scaled{LinearTerm}}
    constant::Sum{Scaled{Constant}}

    function AffineFunction(linear::Sum{Scaled{LinearTerm}}, constant::Sum{Scaled{Constant}})
        outputdim(linear) === outputdim(constant) || throw(ArgumentError())
        new(linear, constant)
    end
end
outputdim(f::AffineFunction) = outputdim(f.constant)
@inline isquadratic(::Type{AffineFunction}) = false
(f::AffineFunction)(vals::Associative{Variable, Float64}) = f.linear(vals) .+ f.constant()


# QuadraticFunction
struct QuadraticFunction <: Fun
    quadratic::Sum{Scaled{QuadraticTerm}}
    affine::AffineFunction

    function QuadraticFunction(quadratic::Sum{Scaled{QuadraticTerm}}, affine::AffineFunction)
        outputdim(quadratic) === outputdim(affine) || throw(ArgumentError())
        new(quadratic, affine)
    end
end
outputdim(f::QuadraticFunction) = 1
@inline isquadratic(::Type{QuadraticFunction}) = true
(f::QuadraticFunction)(vals::Associative{Variable, Float64}) = f.quadratic(vals) .+ f.affine(vals)


# Promotion
Base.promote_rule(::Type{T}, ::Type{T}) where {T<:Fun} = T
Base.promote_rule(::Type{T}, ::Type{Scaled{T}}) where {T<:Fun} = Scaled{T}
Base.promote_rule(::Type{T}, ::Type{Sum{T}}) where {T<:Fun} = Sum{T}
Base.promote_rule(::Type{T}, ::Type{<:Fun}) where {T<:Fun} = isquadratic(T) ? QuadraticFunction : AffineFunction


# Conversion
Base.convert(::Type{Scaled{T}}, x::Fun) where {T<:Fun} = Scaled{T}(1.0, convert(T, x))
Base.convert(::Type{Scaled{T}}, x::Scaled{T}) where {T<:Fun} = x
Base.convert(::Type{Sum{T}}, x::Fun) where {T<:Fun} = Sum{T}(T[convert(T, x)])
Base.convert(::Type{Sum{T}}, x::Sum{T}) where {T<:Fun} = x
Base.convert(::Type{AffineFunction}, x::Constant) = convert(AffineFunction, convert(Scaled{Constant}, x))
Base.convert(::Type{AffineFunction}, x::LinearTerm) = convert(AffineFunction, convert(Scaled{LinearTerm}, x))
Base.convert(::Type{AffineFunction}, x::Scaled{Constant}) = convert(AffineFunction, convert(Sum{Scaled{Constant}}, x))
Base.convert(::Type{AffineFunction}, x::Scaled{LinearTerm}) = convert(AffineFunction, convert(Sum{Scaled{LinearTerm}}, x))
Base.convert(::Type{AffineFunction}, x::Sum{Scaled{LinearTerm}}) =
    AffineFunction(x, convert(Sum{Scaled{Constant}}, Constant(zeros(outputdim(x)))))
Base.convert(::Type{AffineFunction}, x::Sum{Scaled{Constant}}) =
    AffineFunction(convert(Sum{Scaled{LinearTerm}}, LinearTerm(zeros(outputdim(x), 0), Vector{Variable}())), x)
Base.convert(::Type{QuadraticFunction}, x::QuadraticFunction) = x
Base.convert(::Type{QuadraticFunction}, x::QuadraticTerm) = convert(QuadraticFunction, convert(Scaled{QuadraticTerm}, x))
Base.convert(::Type{QuadraticFunction}, x::Scaled{QuadraticTerm}) =
    convert(QuadraticFunction, convert(Sum{Scaled{QuadraticTerm}}, x))
Base.convert(::Type{QuadraticFunction}, x::Sum{Scaled{QuadraticTerm}}) =
    QuadraticFunction(x, convert(AffineFunction, Constant(zeros(outputdim(x)))))
Base.convert(::Type{QuadraticFunction}, x::Fun) =
    QuadraticFunction(convert(Sum{Scaled{QuadraticTerm}}, QuadraticTerm()), convert(AffineFunction, x))


# Operations
Base.:*(x::Real, f::T) where {T<:Fun} = simplify(Scaled{T}(Float64(x), simplify(f)))
Base.:*(f::Fun, x::Real) = x * f
Base.:-(f::Fun) = -1.0 * f
Base.:+(f1::Fun, f2::Fun) = +(promote(f1, f2)...)
Base.:+(f1::T, f2::T) where {T<:Fun} = simplify(Sum{T}(T[simplify(f1), simplify(f2)]))
Base.:-(f1::Fun, f2::Fun) = f1 + -f2


# Simplification
simplify(f::Fun) = f
simplify(f::Scaled{<:Scaled}) = simplify(Scaled(f.scalar * f.val.scalar, f.val.val))
simplify(f::Scaled{<:Sum}) = simplify(Sum([Scaled(f.scalar, term) for term in f.val.terms]))

function simplify(f::Scaled{AffineFunction})
    linear = simplify(Scaled(f.scalar, f.val.linear))
    constant = simplify(Scaled(f.scalar, f.val.constant))
    AffineFunction(linear, constant)
end

function simplify(f::Scaled{QuadraticFunction})
    quadratic = simplify(Scaled(f.scalar, f.val.quadratic))
    affine = simplify(Scaled(f.scalar, f.val.affine))
    QuadraticFunction(quadratic, affine)
end

function simplify(f::Sum{Sum{T}}) where {T}
    terms = Vector{T}()
    for sum in f.terms
        for term in sum.terms
            push!(terms, term)
        end
    end
    simplify(Sum(terms))
end

function simplify(f::Sum{AffineFunction})
    linear = simplify(Sum([term.linear for term in f.terms]))
    constant = simplify(Sum([term.constant for term in f.terms]))
    AffineFunction(linear, constant)
end

function simplify(f::Sum{QuadraticFunction})
    quadratic = simplify(Sum([term.quadratic for term in f.terms]))
    affine = simplify(Sum([term.affine for term in f.terms]))
    QuadraticFunction(quadratic, affine)
end

end
