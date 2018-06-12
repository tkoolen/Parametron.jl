module FunctionsTest

using Compat
using Compat.Test
using StaticArrays: @SVector
using SimpleQP
using SimpleQP.Functions

@testset "LinearTerm" begin
    x = Variable(1)
    @test LinearTerm(4.5, x) === 4.5 * x
    @test LinearTerm(4.5, x) === x * 4.5
    @test convert(LinearTerm{Float64}, x) === 1.0 * x
    @test promote(LinearTerm(4.5, x), x) === (4.5 * x, 1.0 * x)
    @test +(2 * x) === 2 * x
    @test -(2 * x) === -2 * x
    @test 2 * (3 * x) === x * 6
    @test (2 * x) * 3 === 6 * x
    buf = IOBuffer()
    show(buf, LinearTerm(4.5, x))
    @test String(take!(buf)) == "4.5 * x1"
end

@testset "QuadraticTerm" begin
    x = Variable(1)
    y = Variable(2)
    @test x * y === QuadraticTerm(1, x, y)
    @test (2 * x) * y === QuadraticTerm(2, x, y)
    @test x * (2 * y) === QuadraticTerm(2, x, y)

    @test promote(QuadraticTerm(4.5, x, y), QuadraticTerm(3, x, y)) === (QuadraticTerm(4.5, x, y), QuadraticTerm(3.0, x, y))
    @test +(2 * x * y) === 2 * x * y
    @test -(2 * x * y) === -2 * x * y
    @test 2 * (3 * x * y) === x * 6 * y
    @test (y * 2 * x) * 3 === 6 * y * x
    buf = IOBuffer()
    show(buf, QuadraticTerm(1, x, y))
    @test String(take!(buf)) == "1 * x1 * x2"
end

@testset "AffineFunction" begin
    x = Variable(1)
    y = Variable(2)

    vals1 = Dict(x => 1.0, y => 2.0)
    f1 = 2 * x + 3 * y + 5
    @test f1(vals1) == 2.0 + 3.0 * 2.0 + 5
    buf = IOBuffer()
    show(buf, f1)
    @test String(take!(buf)) == "2 * x1 + 3 * x2 + 5"

    vals2 = Dict(x => 2, y => -1)
    f2 = 2.5 * x + 4 * y + 1
    @test f2(vals2) == 2.5 * 2 - 4 + 1

    f3 = zero(typeof(f2))
    @test typeof(f3) == typeof(f2)
    @test f3(vals2) === 0.0

    f4 = f1 + f2
    @test f4(vals1) == f1(vals1) + f2(vals1)
    @test f4(vals2) == f1(vals2) + f2(vals2)

    f5 = f1 + 4.0
    @test f5(vals1) == f1(vals1) + 4.0

    f6 = 1 + f5
    @test f6(vals1) == 1 + f5(vals1)
end

@testset "QuadraticFunction" begin
    x = Variable.(1 : 2)
    vals = Dict(zip(x, [1.0, 2.0]))
    xvals = getindex.(vals, x)
    y = zero(QuadraticFunction{Int})

    @test (x ⋅ x)(vals) == xvals ⋅ xvals
    @test Functions.vecdot!(y, x, x)(vals) == xvals ⋅ xvals

    a = [4.0, 5.0]
    @test (x ⋅ (a .* x))(vals) == xvals ⋅ (a .* xvals)
    @test Functions.vecdot!(y, x , (a .* x))(vals) == xvals ⋅ (a .* xvals)
    @test ((a .* x) ⋅ x)(vals) == (a .* xvals) ⋅ xvals
    @test Functions.vecdot!(y, (a .* x), x)(vals) == (a .* xvals) ⋅ xvals
    @test ((a .* x) ⋅ (a .* x))(vals) == (a .* xvals) ⋅ (a .* xvals)
    @test Functions.vecdot!(y, (a .* x), (a .* x))(vals) == (a .* xvals) ⋅ (a .* xvals)

    b = [6.0, 7.0]
    @test (x ⋅ (a .* x .+ b))(vals) == xvals ⋅ (a .* xvals .+ b)
    @test Functions.vecdot!(y, x, (a .* x .+ b))(vals) == xvals ⋅ (a .* xvals .+ b)
    @test ((a .* x .+ b) ⋅ x)(vals) == xvals ⋅ (a .* xvals .+ b)
    @test Functions.vecdot!(y, (a .* x .+ b), x)(vals) == xvals ⋅ (a .* xvals .+ b)
    @test ((b .* x .+ a) ⋅ (a .* x .+ b))(vals) == (b .* xvals .+ a) ⋅ (a .* xvals .+ b)
    @test Functions.vecdot!(y, b .* x .+ a, a .* x .+ b)(vals) == (b .* xvals .+ a) ⋅ (a .* xvals .+ b)
end

@testset "matvecmul!" begin
    A = ones(Int, 3, 4)
    x = Variable.(1 : size(A, 2))
    y = Vector{AffineFunction{Int}}(undef, size(A, 1))
    Functions.matvecmul!(y, A, x)
    @test y == fill(sum(x), size(A, 1))
    allocs = @allocated Functions.matvecmul!(y, A, x)
    @test allocs == 0
end

@testset "mul!" begin
    x = map(Variable, 1 : 3)
    aff = [1, 2, 3]' * x + 4
    quad = x[1]^2 + 2 * x[1] * x[3] + 3 * x[2] + 4

    @testset "affine" begin
        dest = zero(aff)
        for i = 1 : 2
            Functions.mul!(dest, aff, 2)
            @test dest == 2 * x[1] + 4 * x[2] + 6 * x[3] + 8

            Functions.mul!(dest, 3, aff)
            @test dest == 3 * x[1] + 6 * x[2] + 9 * x[3] + 12
        end
    end
    @testset "quadratic" begin
        dest = zero(quad)
        for i = 1 : 2
            Functions.mul!(dest, quad, 2)
            @test dest == 2 * x[1]^2 + 4 * x[1] * x[3] + 6 * x[2] + 8

            Functions.mul!(dest, 2, quad)
            @test dest == 2 * x[1]^2 + 4 * x[1] * x[3] + 6 * x[2] + 8

            Functions.mul!(dest, aff, x[1])
            @test dest == [1, 2, 3]' * (x .* x[1]) + 4 * x[1]

            Functions.mul!(dest, x[1], aff)
            @test dest == [1, 2, 3]' * (x .* x[1]) + 4 * x[1]

            Functions.mul!(dest, aff, aff)
            @test dest == ([1, 2, 3]' * x + 4)^2
        end
    end
end

@testset "Matrix operations" begin
    x = Variable.(1 : 2)
    A1 = [1.0 2.0; 3.0 4.0]
    fs = A1 * x
    @test fs isa Vector{AffineFunction{Float64}}
    vals = Dict(zip(x, [2.0, 5.0]))
    fvals = [f(vals) for f in fs]
    xvals = getindex.(vals, x)
    @test fvals == A1 * xvals

    a = [1, 2]
    @test a ⋅ x == x ⋅ a == a[1] * x[1] + a[2] * x[2]

    gs = fs .+ a
    gvals = [g(vals) for g in gs]

    @test (gs ⋅ gs)(vals) == gvals ⋅ gvals
    @test (gs ⋅ x)(vals) == gvals ⋅ xvals
    @test (x ⋅ gs)(vals) == xvals ⋅ gvals

    h = x ⋅ x
    @test h(vals) == xvals ⋅ xvals
    buf = IOBuffer()
    show(buf, h)
    @test String(take!(buf)) == "1 * x1 * x1 + 1 * x2 * x2 + 0"
end

@testset "concatenation" begin
    x = Variable(1)
    y = Variable(2)
    f1 = 3x + 10
    f2 = 0.1 * y - 0.5
    f3 = x + y
    v1 = @SVector [f1, f2]
    v2 = @SVector [f2, f3]
    v3 = [f3, f2, f1]
    @test Functions.vcat!(deepcopy(v1), v1) == v1
    @test Functions.vcat!(deepcopy(vcat(v1, v2)), v1, v2) == vcat(v1, v2)
    @test Functions.vcat!(deepcopy(vcat(v1, v2, v3)), v1, v2, v3) == vcat(v1, v2, v3)
    @test_throws DimensionMismatch Functions.vcat!(zeros(AffineFunction{Float64}, 1), v1)
    @test_throws DimensionMismatch Functions.vcat!(zeros(AffineFunction{Float64}, 3), v1)
    @test_throws DimensionMismatch Functions.vcat!(zeros(AffineFunction{Float64}, 3), v1, v2)
    @test_throws DimensionMismatch Functions.vcat!(zeros(AffineFunction{Float64}, 5), v1, v2)
    @test_throws DimensionMismatch Functions.vcat!(zeros(AffineFunction{Float64}, 6), v1, v2, v3)
    @test_throws DimensionMismatch Functions.vcat!(zeros(AffineFunction{Float64}, 8), v1, v2, v3)
end

end # module
