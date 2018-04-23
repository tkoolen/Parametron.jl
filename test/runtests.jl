module SimpleQPTest

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

@testset "LinearFunction basics" begin
    x = [Variable(1), Variable(2)]
    A1 = [1.0 2.0; 3.0 4.0]
    f1 = A1 * x
    f2 = 1.0 * A1 * x
    vals = Dict(x => [2.0, 5.0])
    @test f1(vals) == f2(vals)
    @test f1(vals) == A1 * vals[x]

    f3 = 2 * f2
    @test f3(vals) == 2 * f2(vals)

    f4 = f2 * 2.0
    @test f4(vals) == f3(vals)

    A2 = [4.0 5.0; 6.0 7.0]
    f5 = 1.5 * A1 * x + 2.5 * A2 * x
    @test f5(vals) == 1.5 * A1 * vals[x] + 2.5 * A2 * vals[x]

    f6 = -f5
    @test f6(vals) == -f5(vals)

    f7 = A1 * x - 0.5 * A2 * x
    @test f7(vals) == A1 * vals[x] - 0.5 * A2 * vals[x]
end

@testset "LinearFunction modification" begin
    x = [Variable(1), Variable(2)]
    A1 = [1.0 2.0; 3.0 4.0]
    vals = Dict(x => [2.0, 5.0])
    p = Ref(0.5)
    f1 = p * A1 * x
    f1val = f1(vals)
    p[] = 2.0
    @test f1(vals) == 4.0 * f1val

    f1val = f1(vals)
    A1 .*= 2
    @test f1(vals) == 2 * f1val
end

@testset "LinearFunction sum" begin
    x = [Variable(1), Variable(2)]
    y = [Variable(3), Variable(4)]
    A1 = [1.0 2.0; 3.0 4.0]
    A2 = [4.0 5.0; 6.0 7.0]
    vals = Dict(x => [2.0, 5.0], y => [-1.0, 2.0])

    f1 = 0.3 * A1 * x
    f2 = 2.0 * A2 * y
    f3 = f1 + f2
    @test f3(vals) == f1(vals) + f2(vals)

    p = Ref(2.0)
    f4 = 0.3 * A1 * x + p * A2 * y
    @test f4(vals) == f3(vals)

    copyto!(f1, f4)
    @test f1(vals) == f4(vals)
end

@testset "AffineFunction" begin
    x = [Variable(1), Variable(2), Variable(3)]
    y = [Variable(3), Variable(4), Variable(5)]
    A1 = [1.0 2.0 3.0; 3.0 4.0 5.0]
    A2 = [4.0 5.0 5.0; 6.0 7.0 8.0]
    vals = Dict(x => [2.0, 5.0, 0.5], y => [-1.0, 2.0, 4.5])
    c = [0.1, 0.8]

    f = A1 * x - 3.0 * A2 * y + c
    @test f(vals) == A1 * vals[x] - 3.0 * A2 * vals[y] + c
end

@testset "QuadraticForm" begin
    n = 4
    A1 = Symmetric(sparse(triu(reshape(1.0 : n^2, n, n))))
    x = [Variable(i) for i = 1 : n]

    f1 = quad(A1, x)
    vals = Dict(x => collect(1.0 : n))
    @test f1(vals) ≈ dot(vals[x], A1 * vals[x]) atol = 1e-15

    f2 = f1 + 2 * f1
    @test f2(vals) ≈ 3 * f1(vals) atol = 1e-15

    m = 2
    y = [Variable(i) for i = 1 : m]
    vals[y] = collect(7.0 : 6.0 + m)
    A2 = Symmetric(sparse(triu(reshape(5.0 : 4.0 + m^2, m, m))))
    f3 = 0.5 * quad(A1, x) + 2 * quad(A2, y)
    @test f3(vals) ≈ 0.5 * dot(vals[x], A1 * vals[x]) + 2 * dot(vals[y], A2 * vals[y]) atol = 1e-15

    @test_throws DimensionMismatch quad(A2, x)

    @test QuadraticForm()(vals) == 0.0

    copyto!(f1, f3)
    @test f1(vals) == f3(vals)
end

@testset "Model basics" begin
    optimizer = defaultoptimizer()
    model = Model(optimizer)
    n = 2
    x = [Variable(model) for _ = 1 : n]
    x1, x2 = x
    @test x1.index == MOI.VariableIndex(1)
    @test x2.index == MOI.VariableIndex(2)

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
