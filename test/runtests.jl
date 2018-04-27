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

    rand_data = function (rng)
        Q[1, 1] = rand(rng)
        Q[2, 2] = rand(rng)
        r[1] = rand(rng) - 0.5
        r[2] = rand(rng) - 0.5
    end
    rng = MersenneTwister(1)
    for i = 1 : 10
        rand_data(rng)
        test_unconstrained(model, x, Q, r, s)
    end

    rand_data(rng)
    let model = model
        @allocated solve!(model)
    end

    s[1] = 2.0
    @test_throws ArgumentError solve!(model) # can't modify constant offset
    s[1] = 1.0
end

end # module
