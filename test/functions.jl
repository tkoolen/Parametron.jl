module FunctionsTest

using Compat
using Compat.Test
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

# TODO:
# @testset "QuadraticFunction" begin
#     x = Variable.(1 : 2)

# end

@testset "matvecmul!" begin
    A = ones(Int, 3, 4)
    x = Variable.(1 : size(A, 2))
    y = Vector{AffineFunction{Int}}(undef, size(A, 1))
    Functions.matvecmul!(y, A, x)
    @test y == fill(sum(x), size(A, 1))
    allocs = @allocated Functions.matvecmul!(y, A, x)
    @test allocs == 0
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

end # module
