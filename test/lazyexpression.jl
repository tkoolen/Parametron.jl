module LazyExpressionTest

using Test
using Random
using LinearAlgebra
using Parametron
using StaticArrays

import Random
import Parametron: setdirty!, mock_model
import Parametron.Functions: canonicalize

@testset "basics" begin
    model = mock_model()
    a = 2
    b = Parameter{Int}(() -> 3, model)
    c = 4
    expr = @expression(a + b * c)
    @test expr() == a + b() * c
end

@testset "bad expression" begin
    @test_throws(ArgumentError, @expression x ? y : z)
end

@testset "parameter" begin
    model = mock_model()
    a = 3
    b = 4.0
    cval = Ref(0)
    c = Parameter{Int}(() -> cval[], model)
    expr = @expression(a + b * c)

    @test expr() == a + b * cval[]
    cval[] = 4
    setdirty!(c)
    @test expr() == a + b * cval[]
end

@testset "nested" begin
    model = mock_model()
    a = 3
    b = 4.0
    cval = Ref(5)
    c = Parameter{Int}(() -> cval[], model)
    expr1 = @expression(a + b * c)
    expr2 = @expression(4 * expr1)
    @test expr2() == 4 * expr1()
    show(devnull, expr1)
end

module M
export SpatialMat, angular, linear
struct SpatialMat
    angular::Matrix{Float64}
    linear::Matrix{Float64}
end
angular(mat::SpatialMat) = mat.angular
linear(mat::SpatialMat) = mat.linear
end

using .M

@testset "user functions" begin
    mat = SpatialMat(rand(3, 4), rand(3, 4))
    scalar = Ref(1.0)
    updatemat! = let scalar = scalar # https://github.com/JuliaLang/julia/issues/15276
        mat -> (mat.angular .= scalar[]; mat.linear .= scalar[]; mat)
    end
    model = mock_model()
    pmat = Parameter(updatemat!, mat, model)
    pmat_angular = @expression angular(pmat)
    result = pmat_angular()
    @test result === angular(mat)
    @test all(result .== scalar[])

    setdirty!(model)
    allocs = @allocated begin
        setdirty!(model)
        pmat_angular()
    end
    @test allocs == 0
end

@testset "matvecmul!" begin
    m = mock_model()
    A = Parameter(rand!, zeros(3, 4), m)
    x = Variable.(1 : 4)
    expr = @expression A * x
    @test expr() == A() * x
    setdirty!(m)
    allocs = @allocated begin
        setdirty!(m)
        expr()
    end
    @test allocs == 0

    wrapped = Parametron.wrap(expr)
    @test wrapped isa Parametron.WrappedExpression{Vector{AffineFunction{Float64}}}
    setdirty!(m)
    @test wrapped() == expr()
    allocs = @allocated begin
        setdirty!(m)
        wrapped()
    end
    @test allocs == 0
    wrapped² = Parametron.wrap(wrapped)
    @test wrapped² === wrapped
end

@testset "StaticArrays" begin
    m = mock_model()
    A = Parameter{SMatrix{3, 3, Int, 9}}(m) do
        @SMatrix ones(Int, 3, 3)
    end
    x = Variable.(1 : 3)

    expr1 = @expression A * x
    @test canonicalize.(expr1()) == canonicalize.(A() * x)
    setdirty!(A)
    allocs = @allocated expr1()
    @test_broken allocs == 0
    @test expr1() isa SVector{3, AffineFunction{Int}}

    y = Parameter(m, val=SVector{3}(x))
    expr2 = @expression y + y
    @test expr2() == y() + y()
    allocs = @allocated expr2()
    @test_broken allocs == 0
    @test expr2() isa SVector{3, AffineFunction{Int}}

    expr3 = @expression y - y
    @test expr3() == y() - y()
    allocs = @allocated expr3()
    @test_broken allocs == 0
    @test expr3() isa SVector{3, AffineFunction{Int}}
end

@testset "mul! optimization" begin
    m = mock_model()
    weight = Parameter(() -> 3, m)
    x = Variable.(1 : 3)
    expr = @expression weight * (x ⋅ x)
    vals = Dict(zip(x, [1, 2, 3]))
    xvals = map(var -> vals[var], x)
    @test expr()(vals) == 3 * xvals ⋅ xvals
    allocs = @allocated expr()
    @test_broken allocs == 0
end

@testset "scale! optimization" begin
    model = mock_model()
    x = [Variable(model) for _ in 1 : 3]
    dt = Parameter{Float64}(model) do
        2.0
    end
    @testset "Number" begin
        expr = @expression dt * [1, 2, 3, 4, 5]
        @test expr() == dt() * [1, 2, 3, 4, 5]
        @test @allocated(expr()) == 0

        expr = @expression [1, 2, 3, 4, 5] * dt
        @test expr() == [1, 2, 3, 4, 5] * dt()
        @test @allocated(expr()) == 0
    end
    @testset "Variable" begin
        expr = @expression dt * x
        @test expr() == dt() * x
        @test @allocated(expr()) == 0

        expr = @expression x * dt
        @test expr() == x * dt()
        @test @allocated(expr()) == 0
    end

    @testset "AffineFunction" begin
        Ax = ones(3, 3) * x
        expr = @expression dt * Ax
        @test expr() == dt() * Ax
        @test @allocated(expr()) == 0

        expr = @expression Ax * dt
        @test expr() == Ax * dt()
        @test @allocated(expr()) == 0
    end
end


@testset "vcat optimization" begin
    Random.seed!(42)
    m = mock_model()
    A = Parameter(rand!, zeros(3, 4), m)
    B = Parameter(rand!, zeros(3, 3), m)
    x = Variable.(1 : 4)
    y = Variable.(5 : 7)
    f1 = @expression A * x
    f2 = @expression B * y

    v1 = @expression vcat(f1, f2)
    @test v1() == vcat(f1(), f2())
    setdirty!(m)
    @test (@allocated begin
        setdirty!(m)
        v1()
    end) == 0

    # Make sure we expand vcat expressions
    v2 = @expression [f1; f2]
    @test v2() == vcat(f1(), f2())
    @test (@allocated begin
        setdirty!(m)
        v2()
    end) == 0

    # Make sure static arrays still work
    C = Parameter{SMatrix{2, 2, Int, 4}}(m) do
        @SMatrix [1 2; 3 4]
    end
    z = SVector(Variable(8), Variable(9))
    f3 = @expression C * z
    v3 = @expression [f3; f3]
    @test v3() == vcat(f3(), f3())
    @test v3() isa SVector{4, AffineFunction{Int}}
    @test_broken (@allocated begin
        setdirty!(m)
        v3()
    end) == 0

    # Other numbers of arguments
    v4 = @expression vcat(f3)
    @test v4() == f3()
    @test_broken (@allocated begin
        setdirty!(m)
        v4()
    end) == 0

    v5 = @expression vcat(f1, f2, f3)
    @test v5() == vcat(f1(), f2(), f3())
    @test_broken (@allocated begin
        setdirty!(m)
        v5()
    end) == 0

    # Generic fallbacks that allocate memory but should still give the right answer
    @test (@expression vcat(f1, 3))() == vcat(f1(), 3)
end

@testset "vect optimization" begin
    m = mock_model()
    x = Variable.(1:2)
    p = Parameter(m) do
        SVector(1., 2.)
    end
    expr = @expression [dot(p, x)]
    @test expr() == [dot(p(), x)]
    @test_broken @allocated(expr()) == 0

    p2 = Parameter(m) do
        2.0
    end
    ex2 = @expression [p2]
    @test ex2() == [2.0]
    @test @allocated(ex2()) == 0
end

@testset "convert optimization" begin
    m = mock_model()
    A = Parameter{SMatrix{3, 3, Int, 9}}(m) do
        @SMatrix ones(Int, 3, 3)
    end
    x = Variable.(1 : 3)
    expr = @expression convert(Vector, A * x)

    @test canonicalize.(expr()) == canonicalize.(A() * x)
    @test expr() isa Vector{AffineFunction{Int}}
    allocs = @allocated expr()
    @test_broken allocs == 0
end

@testset "vecdot optimization" begin
    m = mock_model()
    x = Variable.(1:3)
    p = Parameter(m) do
        SVector(1., 2., 3.)
    end
    ex1 = @expression dot(x, p)
    ex2 = @expression dot(p, x)
    @test ex1() == ex2() == dot(p(), x)
    @test ex1() isa AffineFunction
    @test_broken @allocated(ex1()) == 0
    @test_broken @allocated(ex2()) == 0

    ex3 = @expression dot(x + p, p)
    ex4 = @expression dot(p, x + p)
    @test ex3() == ex4() == dot(p(), x + p())
    @test ex3() isa AffineFunction
    @test_broken @allocated(ex3()) == 0
    @test_broken @allocated(ex4()) == 0
end

@testset "adjoint optimization" begin
    m = mock_model()
    x = Variable.(1:2)
    p = Parameter(m) do
        @SMatrix [1. 2.; 3. 4.]
    end
    ex = @expression adjoint(p) * x
    @test canonicalize.(ex()) == canonicalize.(adjoint(p()) * x)
    @test_broken @allocated(ex()) == 0

    p2 = Parameter(rand!, zeros(2, 2), m)
    ex2 = @expression adjoint(p2) * x
    @test ex2() == adjoint(p2()) * x
    @test @allocated(ex2()) == 0
end

@testset "issue 26" begin
    model = mock_model()
    n = 2
    x = Variable.(1 : n)
    vals = Dict(zip(x, rand(n)))
    q = Parameter(zeros(n), model) do q
        q[1] = 1
        q[2] = 2
    end
    expr1 = @expression transpose(x) * Matrix(I, n, n) * x + q' * x
    expr2 = @expression q' * x + transpose(x) * Matrix(I, 2, 2) * x
    @test expr1()(vals) ≈ expr2()(vals) atol = 1e-14
end

@testset "issue 32" begin
    model = mock_model()
    v = [Variable(model) for _ in 1:2]
    v0 = zeros(2)
    Δt = 0.01
    u = [Variable(model) for _ in 1:2]
    H = Parameter(zeros(2, 2), model) do H
        H[1, 1] = rand()
        H[2, 2] = rand()
    end
    c = Parameter(zeros(2), model) do c
        c[1] = rand()
        c[2] = rand()
    end
    expr = @expression(H * (v - v0) - Δt * (u - c))
    @test expr() == H() * (v - v0) - Δt * (u - c())
    @test @allocated(expr()) == 0
end

@testset "issue 23, fully qualified names" begin
    model = mock_model()
    p = Parameter{Int}(() -> 2, model)
    expr = @expression Base.:+(2, p)
    @test expr() == 4
end

struct MyWrapper{T}
    x::T
end

function allocations_in_local_scope(x)
    @allocated x()
end

@testset "getfield optimization" begin
    model = mock_model()
    p = Parameter(model) do
        MyWrapper(rand())
    end
    ex1 = @expression p.x
    @test ex1() == p().x
    setdirty!(p)
    @test ex1() == p().x
    setdirty!(p)
    @test allocations_in_local_scope(ex1) == 0

    ex2 = @expression p.x + 1
    @test ex2() == p().x + 1
    setdirty!(p)
    @test ex2() == p().x + 1
    setdirty!(p)
    @test allocations_in_local_scope(ex2) == 0
end

@testset "hcat" begin
    model = mock_model()
    p = Parameter(model) do
        rand(3, 3)
    end
    hcat_expr = @expression [p p]
    @test hcat_expr() == [p() p()]
    setdirty!(p)
    @test hcat_expr() == [p() p()]
end

@testset "getindex" begin
    model = mock_model()
    p = Parameter(model) do
        rand(3, 3)
    end
    getindex_expr = @expression p[:,2]
    @test getindex_expr() == p()[:,2]
    setdirty!(p)
    @test getindex_expr() == p()[:,2]
end

@testset "A' * x (https://github.com/tkoolen/Parametron.jl/pull/79#issuecomment-423732633)" begin
    model = mock_model()
    n, m = 5, 15
    Adata = randn(n, m)
    A = Parameter(model, val=Adata)
    x = [Variable(model) for _ = 1:n]
    resid = @expression A' * x
    @test resid() == Adata' * x
end

@testset "tuple (https://github.com/tkoolen/Parametron.jl/issues/88)" begin
    A = [1, 2, 3, 4]
    @test(@expression(reshape(A, (2, 2))) == reshape(A, (2, 2)))
end

@testset "scale SVector" begin
    model = mock_model()
    x = Parameter(model, val=SVector(Variable(1), Variable(2)))
    expr = @expression x * 3
    @test expr() == SVector(3 * x.val[][1], 3 * x.val[][2])
    expr = @expression 3 * x
    @test expr() == SVector(x.val[][1] * 3, x.val[][2] * 3)
end

end # module
