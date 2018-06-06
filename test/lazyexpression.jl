module LazyExpressionTest

using Compat
using Compat.Test
using SimpleQP

import ..MockModel
import SimpleQP: setdirty!

@testset "basics" begin
    a = 2
    b = 3
    c = 4
    expr = @expression(a + b * c)
    @test expr() == a + b * c
end

@testset "bad expression" begin
    @test_throws(ArgumentError, @expression x ? y : z)
end

@testset "parameter" begin
    model = MockModel()
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
    model = MockModel()
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
    updatemat! = let scalar = scalar
        mat -> (mat.angular .= scalar[]; mat.linear .= scalar[]; mat)
    end
    model = MockModel()
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
    m = MockModel()
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

    wrapped = SimpleQP.WrappedExpression{Vector{AffineFunction{Float64}}}(expr)
    setdirty!(m)
    @test wrapped() == expr()
    allocs = @allocated begin
        setdirty!(m)
        wrapped()
    end
    @test allocs == 0
end

end
