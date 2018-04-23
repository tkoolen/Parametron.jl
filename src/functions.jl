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

Base.:*(scale::Union{Float64, Base.RefValue{Float64}}, A::Matrix{Float64}, x::Vector{Variable}) = LinearFunction(scale, A, x)
Base.:*(A::Matrix{Float64}, x::Vector{Variable}) = LinearFunction(A, x)

function (f::LinearFunction)(vals::Associative{Vector{Variable}, Vector{Float64}})
    ret = zeros(outputdim(f))
    for i = 1 : numterms(f)
        scale = prod(s -> s[], f.scales[i])
        ret += scale * f.As[i] * vals[f.xs[i]]
    end
    ret
end


struct QuadraticForm
    scales::Vector{Vector{Base.RefValue{Float64}}}
    As::Vector{SparseSymmetric64}
    xs::Vector{Vector{Variable}}

    function QuadraticForm(
            scales::Vector{Vector{Base.RefValue{Float64}}},
            As::Vector{SparseSymmetric64},
            xs::Vector{Vector{Variable}})
        nterms = length(scales)
        length(As) === nterms || throw(ArgumentError())
        length(xs) === nterms || throw(ArgumentError())
        for i = 1 : nterms
            As[i].uplo == 'U' || throw(ArgumentError())
            @boundscheck begin
                ni = length(xs[i])
                size(As[i]) === (ni, ni) || throw(DimensionMismatch())
            end
        end
        new(scales, As, xs)
    end
end

QuadraticForm(scale::Base.RefValue{Float64}, Q::SparseSymmetric64, x::Vector{Variable}) = QuadraticForm([[scale]], [Q], [x])
QuadraticForm(scale::Number, Q::SparseSymmetric64, x::Vector{Variable}) = QuadraticForm(Ref(Float64(scale)), Q, x)
QuadraticForm(Q::SparseSymmetric64, x::Vector{Variable}) = QuadraticForm(1.0, Q, x)
QuadraticForm() = QuadraticForm(Vector{Vector{Base.RefValue{Float64}}}(), Vector{SparseSymmetric64}(), Vector{Vector{Variable}}())

numterms(f::QuadraticForm) = length(f.xs)
quad(Q::SparseSymmetric64, x::Vector{Variable}) = QuadraticForm(Q, x)

function (f::QuadraticForm)(vals::Associative{Vector{Variable}, Vector{Float64}})
    ret = zero(Float64)
    for i = 1 : numterms(f)
        scale = prod(s -> s[], f.scales[i])
        ret += scale * quad(f.As[i], vals[f.xs[i]])
    end
    ret
end


for Function in [LinearFunction, QuadraticForm]
    @eval begin
        function Base.:*(scale::Base.RefValue{Float64}, f::$Function)
            scales = [push!(copy(f.scales[i]), scale) for i = 1 : numterms(f)]
            $Function(scales, f.As, f.xs)
        end

        Base.:*(scale::Number, f::$Function) = Ref(Float64(scale)) * f
        Base.:*(f::$Function, scale::Union{Base.RefValue{Float64}, Number}) = scale * f

        function Base.:+(f1::$Function, f2::$Function)
            scales = append!(copy(f1.scales), f2.scales)
            As = append!(copy(f1.As), f2.As)
            xs = append!(copy(f1.xs), f2.xs)
            $Function(scales, As, xs)
        end

        function Base.:-(f::$Function)
            n = numterms(f)
            scales = Vector{Vector{Base.RefValue{Float64}}}(n)
            for i = 1 : n
                scales[i] = push!(copy(f.scales[i]), Ref(-1.0))
            end
            $Function(scales, f.As, f.xs)
        end

        Base.:-(f1::$Function, f2::$Function) = f1 + -f2

        function Compat.copyto!(dest::$Function, src::$Function)
            nterms = numterms(src)
            resize!(dest.scales, nterms)
            resize!(dest.As, nterms)
            resize!(dest.xs, nterms)
            for i = 1 : nterms
                dest.scales[i] = copy(src.scales[i])
                dest.As[i] = src.As[i]
                dest.xs[i] = copy(src.xs[i])
            end
            dest
        end
    end
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

Base.:+(f::LinearFunction, c::Vector{Float64}) = AffineFunction(f, c)
Base.:+(f1::AffineFunction, f2::LinearFunction) = AffineFunction(f1.linear + f2, f1.constant)
Base.:+(f1::LinearFunction, f2::AffineFunction) = f2 + f1
Base.:+(f1::AffineFunction, f2::AffineFunction) = AffineFunction(f1.linear + f2.linear, f1.constant + f2.constant)

(f::AffineFunction)(vals::Associative{Vector{Variable}, Vector{Float64}}) = f.linear(vals) + f.constant
