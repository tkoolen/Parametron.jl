module SimpleQP

export
    Variable,
    LinearFunction,
    AffineFunction

import MathOptInterface

const MOI = MathOptInterface

struct Variable
    index::MOI.VariableIndex
end
Variable(index::Int) = Variable(MOI.VariableIndex(index))

struct LinearFunction
    scales::Vector{Vector{Base.RefValue{Float64}}}
    As::Vector{Matrix{Float64}}
    xs::Vector{Vector{Variable}}

    function LinearFunction(
            scales::Vector{Vector{Base.RefValue{Float64}}},
            As::Vector{Matrix{Float64}},
            xs::Vector{Vector{Variable}})
        nterms = length(scales)
        length(As) === nterms || throw(ArgumentError())
        length(xs) === nterms || throw(ArgumentError())
        @boundscheck begin
            m = isempty(As) ? 0 : size(As[1], 1)
            for i = 1 : nterms
                size(As[i], 1) === m || throw(DimensionMismatch())
                size(As[i], 2) === length(xs[i]) || throw(DimensionMismatch())
            end
        end
        new(scales, As, xs)
    end
end

LinearFunction(scale::Base.RefValue{Float64}, A::Matrix{Float64}, x::Vector{Variable}) = LinearFunction([[scale]], [A], [x])
LinearFunction(scale::Number, A::Matrix{Float64}, x::Vector{Variable}) = LinearFunction(Ref(Float64(scale)), A, x)
LinearFunction(A::Matrix{Float64}, x::Vector{Variable}) = LinearFunction(1.0, A, x)

outputdim(f::LinearFunction) = isempty(f.As) ? 0 : size(f.As[1], 1)
numterms(f::LinearFunction) = length(f.xs)

Base.:*(scale::Float64, A::Matrix{Float64}, x::Vector{Variable}) = LinearFunction(scale, A, x)
Base.:*(A::Matrix{Float64}, x::Vector{Variable}) = LinearFunction(A, x)

function Base.:*(scale::Base.RefValue{Float64}, f::LinearFunction)
    scales = [push!(copy(f.scales[i]), scale) for i = 1 : numterms(f)]
    LinearFunction(scales, f.As, f.xs)
end

Base.:*(scale::Number, f::LinearFunction) = Ref(Float64(scale)) * f
Base.:*(f::LinearFunction, scale::Union{Base.RefValue{Float64}, Number}) = scale * f

function Base.:+(f1::LinearFunction, f2::LinearFunction)
    scales = append!(copy(f1.scales), f2.scales)
    As = append!(copy(f1.As), f2.As)
    xs = append!(copy(f1.xs), f2.xs)
    LinearFunction(scales, As, xs)
end

function Base.:-(f::LinearFunction)
    n = numterms(f)
    scales = Vector{Vector{Base.RefValue{Float64}}}(n)
    for i = 1 : n
        scales[i] = push!(copy(f.scales[i]), Ref(-1.0))
    end
    LinearFunction(scales, f.As, f.xs)
end

Base.:-(f1::LinearFunction, f2::LinearFunction) = f1 + -f2

function (f::LinearFunction)(vals::Associative{Vector{Variable}, Vector{Float64}})
    ret = zeros(outputdim(f))
    for i = 1 : numterms(f)
        scale = prod(s -> s[], f.scales[i])
        ret += scale * f.As[i] * vals[f.xs[i]]
    end
    ret
end

struct AffineFunction
    linear::LinearFunction
    constant::Vector{Float64}

    function AffineFunction(linear::LinearFunction, constant::Vector{Float64})
        @boundscheck outputdim(linear) === length(constant) || throw(DimensionMismatch())
        new(linear, constant)
    end
end

outputdim(f::AffineFunction) = length(f.constant)

Base.:+(f1::AffineFunction, f2::LinearFunction) = AffineFunction(push!(copy(f1.linearterms), f2), f1.constant)
Base.:+(f1::LinearFunction, f2::AffineFunction) = f2 + f1
Base.:+(f1::AffineFunction, f2::AffineFunction) =
    AffineFunction(append!(copy(f1.linearterms), f2.linearterms), f1.constant + f2.constant)

end # module
