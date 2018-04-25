module Functions

export
    Variable,
    Constant,
    LinearTerm,
    QuadraticTerm,
    Scaled,
    Sum,
    AffineFunction

using MathOptInterface

const MOI = MathOptInterface
const SparseSymmetric64 = Symmetric{Float64,SparseMatrixCSC{Float64,Int}}

struct Variable
    index::MOI.VariableIndex
end
Variable(index::Int) = Variable(MOI.VariableIndex(index))

abstract type Fun end
# Base.show(io::IO, f::Fun) = print(io, typeof(f)) # FIXME

struct Constant <: Fun
    v::Vector{Float64}
end
outputdim(f::Constant) = length(f.v)
(f::Constant)() = f.v


struct LinearTerm <: Fun
    A::Matrix{Float64}
    x::Vector{Variable}

    function LinearTerm(A::Matrix{Float64}, x::Vector{Variable})
        size(A, 2) === length(x) || error()
        new(A, x)
    end
end
outputdim(f::LinearTerm) = size(f.A, 1)
(f::LinearTerm)(vals::Associative{Variable, Float64}) = f.A * getindex.(vals, f.x)


struct QuadraticTerm <: Fun
    Q::SparseSymmetric64
    x::Vector{Variable}
end
(f::QuadraticTerm)(vals::Associative{Variable, Float64}) = quad(f.Q, getindex.(vals, f.x))


struct Scaled{T<:Fun} <: Fun
    scalar::Float64
    val::T
    Scaled{T}(scalar::Float64, val::T) where {T<:Fun} = new{T}(scalar, val)
end
Scaled(scalar::Float64, val::T) where {T<:Fun} = Scaled{T}(scalar, val)
outputdim(f::Scaled) = outputdim(f.val)
(f::Scaled)(x...) = f.scalar * f.val(x...)


struct Sum{T<:Fun} <: Fun
    terms::Vector{T}

    function Sum{T}(terms::Vector{T}) where {T}
        isempty(terms) && throw(ArgumentError())
        m = outputdim(terms[1])
        for i = 2 : length(terms)
            outputdim(terms[i]) == m || throw(ArgumentError())
        end
        new{T}(terms)
    end
end
Sum(terms::Vector{T}) where {T<:Fun} = Sum{T}(terms)
outputdim(f::Sum) = isempty(f.terms) ? 0 : outputdim(f.terms[1])
(f::Sum)(x...) = sum(term -> term(x...), f.terms)


struct AffineFunction <: Fun
    linear::Sum{Scaled{LinearTerm}}
    constant::Sum{Scaled{Constant}}

    function AffineFunction(linear::Sum{Scaled{LinearTerm}}, constant::Sum{Scaled{Constant}})
        outputdim(linear) === outputdim(constant) || throw(ArgumentError())
        new(linear, constant)
    end
end
outputdim(f::AffineFunction) = outputdim(f.constant)
(f::AffineFunction)(vals::Associative{Variable, Float64}) = f.linear(vals) .+ f.constant()

Base.:*(A::Matrix{Float64}, x::Vector{Variable}) = LinearTerm(A, x)
Base.:*(A::AbstractMatrix, x::Vector{Variable}) = error("Only Matrix{Float64} is supported.")
quad(Q::SparseSymmetric64, x::Vector{Variable}) = QuadraticTerm(Q, x)

# Promotion
Base.promote_rule(::Type{T}, ::Type{T}) where {T<:Fun} = T
Base.promote_rule(::Type{T}, ::Type{Scaled{T}}) where {T<:Fun} = Scaled{T}
Base.promote_rule(::Type{T}, ::Type{Sum{T}}) where {T<:Fun} = Sum{T}
Base.promote_rule(::Type{<:Fun}, ::Type{<:Fun}) = AffineFunction

# Conversion
Base.convert(::Type{Scaled{T}}, x) where {T<:Fun} = Scaled{T}(1.0, convert(T, x))
Base.convert(::Type{Scaled{T}}, x::Scaled{T}) where {T<:Fun} = x
Base.convert(::Type{Sum{T}}, x) where {T<:Fun} = Sum{T}(T[convert(T, x)])
Base.convert(::Type{Sum{T}}, x::Sum{T}) where {T<:Fun} = x
Base.convert(::Type{AffineFunction}, x::Constant) = convert(AffineFunction, convert(Scaled{Constant}, x))
Base.convert(::Type{AffineFunction}, x::LinearTerm) = convert(AffineFunction, convert(Scaled{LinearTerm}, x))
Base.convert(::Type{AffineFunction}, x::Scaled{Constant}) = convert(AffineFunction, convert(Sum{Scaled{Constant}}, x))
Base.convert(::Type{AffineFunction}, x::Scaled{LinearTerm}) = convert(AffineFunction, convert(Sum{Scaled{LinearTerm}}, x))
Base.convert(::Type{AffineFunction}, x::Sum{Scaled{LinearTerm}}) = AffineFunction(x, convert(Sum{Scaled{Constant}}, Constant(zeros(outputdim(x)))))
Base.convert(::Type{AffineFunction}, x::Sum{Scaled{Constant}}) = AffineFunction(convert(Sum{Scaled{LinearTerm}}, LinearTerm(zeros(outputdim(x), 0), Vector{Variable}())), x)

# Operations
Base.:*(f::Fun, x::Real) = x * f
Base.:*(x::Real, f::T) where {T<:Fun} = simplify(Scaled{T}(Float64(x), f))
Base.:-(f::Fun) = -1.0 * f
Base.:+(f1::Fun, f2::Fun) = +(promote(f1, f2)...)
Base.:+(f1::T, f2::T) where {T<:Fun} = simplify(Sum{T}(T[f1, f2]))
Base.:-(f1::Fun, f2::Fun) = f1 + -f2

# Simplification
simplify(f::Fun) = f
simplify(f::Scaled{<:Scaled}) = simplify(Scaled(f.scalar * f.val.scalar, f.val.val))
simplify(f::Scaled{<:Sum}) = simplify(Sum([Scaled(f.scalar, term) for term in f.val.terms]))

function simplify(f::Scaled{AffineFunction})
    linear = simplify(Scaled(f.scalar, f.linear))
    constant = simplify(Scaled(f.scalar, f.constant))
    AffineFunction(linear, constant)
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

end
