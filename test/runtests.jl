module SimpleQPTest

using Compat
using Compat.Test
using SimpleQP
using SimpleQP.Functions

@testset "Linear function basics" begin
    x = [Variable(1), Variable(2)]
    A1 = [1.0 2.0; 3.0 4.0]
    f1 = LinearTerm(A1, x)
    f2 = 1.0 * LinearTerm(A1, x)
    vals = Dict(zip(x, [2.0, 5.0]))
    @test f1(vals) == f2(vals)
    @test f1(vals) == A1 * getindex.(vals, x)

    f3 = 2 * f2
    @test f3(vals) == 2 * f2(vals)

    f4 = f2 * 2.0
    @test f4(vals) == f3(vals)

    A2 = [4.0 5.0; 6.0 7.0]
    f5 = 1.5 * LinearTerm(A1, x) + 2.5 * LinearTerm(A2, x)
    @test f5(vals) == 1.5 * A1 * getindex.(vals, x) + 2.5 * A2 * getindex.(vals, x)

    f6 = -f5
    @test f6(vals) == -f5(vals)

    f7 = LinearTerm(A1, x) - 0.5 * LinearTerm(A2, x)
    @test f7(vals) == A1 * getindex.(vals, x) - 0.5 * A2 * getindex.(vals, x)

    y = [Variable(3), Variable(4)]
    A2 = [4.0 5.0; 6.0 7.0]
    vals = merge(Dict(zip(x, [2.0, 5.0])), Dict(zip(y, [-1.0, 2.0])))

    f1 = 0.3 * LinearTerm(A1, x)
    f2 = 2.0 * LinearTerm(A2, y)
    f3 = f1 + f2
    @test f3(vals) == f1(vals) + f2(vals)
end

@testset "Linear function modification" begin
    x = [Variable(1), Variable(2)]
    A1 = [1.0 2.0; 3.0 4.0]
    vals = Dict(zip(x, [2.0, 5.0]))
    f1 = 0.5 * LinearTerm(A1, x)
    f1val = f1(vals)
    A1 .*= 2
    @test f1(vals) == 2.0 * f1val
end

@testset "AffineFunction" begin
    x = Variable.(1 : 3)
    y = Variable.(4 : 6)
    A1 = [1.0 2.0 3.0; 3.0 4.0 5.0]
    A2 = [4.0 5.0 5.0; 6.0 7.0 8.0]
    vals = merge(Dict(zip(x, [2.0, 5.0, 0.5])), Dict(zip(y, [-1.0, 2.0, 4.5])))
    c = [0.1, 0.8]

    f = LinearTerm(A1, x) - 3.0 * LinearTerm(A2, y) + Constant(c)
    @test f(vals) == A1 * getindex.(vals, x) - 3.0 * A2 * getindex.(vals, y) + c
end

@testset "Quadratic" begin
    n = 4
    A1 = Symmetric(sparse(triu(reshape(1.0 : n^2, n, n))))
    x = [Variable(i) for i = 1 : n]

    f1 = quad(A1, x)
    xval = collect(1.0 : n)
    vals = Dict(zip(x, xval))
    @test f1(vals) ≈ dot(xval, A1 * xval) atol = 1e-15

    f2 = f1 + 2 * f1
    @test f2(vals) ≈ 3 * f1(vals) atol = 1e-15

    m = 2
    y = [Variable(i) for i = n + 1 : n + m]
    yval = collect(7.0 : 6.0 + m)
    merge!(vals, Dict(zip(y, yval)))
    A2 = Symmetric(sparse(triu(reshape(5.0 : 4.0 + m^2, m, m))))
    f3 = 0.5 * quad(A1, x) + 2 * quad(A2, y)
    @test f3(vals) ≈ 0.5 * dot(xval, A1 * xval) + 2 * dot(yval, A2 * yval) atol = 1e-15

    f4 = f3 + [3.0]
    @test f4(vals) ≈ f3(vals) + 3.0 atol = 1e-15

    @test_throws DimensionMismatch quad(A2, x)

    @test QuadraticTerm()(vals) == 0.0
end

using Compat
using Compat.Test
using SimpleQP
using OSQP.MathOptInterfaceOSQP
import MathOptInterface
const MOI = MathOptInterface

function defaultoptimizer()
    optimizer = OSQPOptimizer()
    MOI.set!(optimizer, OSQPSettings.Verbose(), false)
    MOI.set!(optimizer, OSQPSettings.EpsAbs(), 1e-8)
    MOI.set!(optimizer, OSQPSettings.EpsRel(), 1e-16)
    MOI.set!(optimizer, OSQPSettings.MaxIter(), 10000)
    MOI.set!(optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) # required for deterministic behavior
    optimizer
end

@testset "Model basics" begin
    optimizer = defaultoptimizer()
    model = Model(optimizer)
    n = 2
    x = [Variable(model) for _ = 1 : n]
    x1, x2 = x
    @test x1.index == 1
    @test x2.index == 2

    Q = Symmetric(speye(n))
    objective = quad(Q, x)
    setobjective!(model, Senses.Min, objective)

    # cf = eye(n) * x + -ones(n)
    # add_nonnegative_constraint!(model, cf)


    SimpleQP.initialize!(model)
    @test model.initialized


    solve!(model)

    SimpleQP.update!(model)
end

end
