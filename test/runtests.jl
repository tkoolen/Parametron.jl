module SimpleQPTest

include("functions.jl")

using Compat
using Compat.Test
using SimpleQP
using SimpleQP.Functions
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

function test_unconstrained(model, x, Q, r, s; atol=1e-15)
    solve!(model)
    xval = value.(model, x)
    expected = -2 * Q \ r'
    @test xval ≈ expected atol = atol
    @test objectivevalue(model) ≈ dot(xval, Q * xval) + dot(vec(r), xval) + s[1] atol = atol
end

@testset "Model: unconstrained" begin
    optimizer = defaultoptimizer()
    model = Model(optimizer)
    n = 2
    x = [Variable(model) for _ = 1 : n]
    x1, x2 = x
    @test x1.index == 1
    @test x2.index == 2

    Q = Symmetric(speye(n))
    r = zeros(1, n)
    s = [1.0]
    objective = QuadraticTerm(Q, x) + LinearTerm(r, x) + s

    setobjective!(model, Senses.Min, objective)
    SimpleQP.initialize!(model)
    @test model.initialized

    test_unconstrained(model, x, Q, r, s)

    srand(1)
    for i = 1 : 10
        Q[1, 1] = rand()
        Q[2, 2] = rand()
        r[1] = rand() - 0.5
        r[2] = rand() - 0.5
        test_unconstrained(model, x, Q, r, s)
    end

    s[1] = 2.0
    @test_throws ArgumentError solve!(model) # can't modify constant offset
    s[1] = 1.0
end

end
