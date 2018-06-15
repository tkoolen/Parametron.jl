module ModelTests

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

function test_unconstrained(model, x, Q, r, s; atol=1e-15)
    solve!(model)
    xval = value.(model, x)
    expected = -2 * Q() \ r()
    @test xval ≈ expected atol = atol
    @test objectivevalue(model) ≈ dot(xval, Q() * xval) + dot(r(), xval) + s() atol = atol
end

@testset "Model: unconstrained" begin
    optimizer = defaultoptimizer()
    model = Model(optimizer)
    n = 2
    x = [Variable(model) for _ = 1 : n]
    x1, x2 = x
    @test x1.index == 1
    @test x2.index == 2

    rng = MersenneTwister(1)

    sval = Ref(1.0)
    Q = let rng = rng # https://github.com/JuliaLang/julia/issues/15276
        Parameter(eye(n), model) do Q
            Q[1, 1] = rand(rng)
            Q[2, 2] = rand(rng)
        end
    end
    r = let rng = rng # https://github.com/JuliaLang/julia/issues/15276
        Parameter(zeros(n), model) do r
            rand!(rng, r)
        end
    end
    s = let sval = sval # https://github.com/JuliaLang/julia/issues/15276
        Parameter{Float64}(model) do
            sval[]
        end
    end

    @objective(model, Minimize, transpose(x) * Q * x + r ⋅ x + s)
    SimpleQP.initialize!(model)
    @test model.initialized

    test_unconstrained(model, x, Q, r, s)

    allocs = @allocated solve!(model)
    @test_broken allocs == 0
    @test allocs <= 16

    sval[] = 2.0
    @test_throws ArgumentError solve!(model) # can't modify constant offset
end

@testset "Model: equality constrained" begin
    # Minimize ||A x - b||^2 = x' A' A x - (2 * A' * b)' x + b' * b
    # subject to C x = d

    test_equality_constrained = function(A, b, C, d, x; rtol = 1e-4)
        C⁺ = pinv(C)
        Q = I - C⁺ * C
        expected = Q * (pinv(A * Q) * (b - A * C⁺ * d)) + C⁺ * d # note: can be quite badly conditioned
        @test x ≈ expected rtol = rtol
    end

    n = 8
    m = 2

    optimizer = defaultoptimizer()
    model = Model(optimizer)
    x = [Variable(model) for _ = 1 : n]

    rng = MersenneTwister(1234)
    randrng! = let rng = rng
        x -> rand!(rng, x)
    end
    A = Parameter(randrng!, zeros(n, n), model)
    b = Parameter(randrng!, zeros(n), model)
    C = Parameter(randrng!, zeros(m, n), model)
    d = Parameter(randrng!, zeros(m), model)

    residual = @expression A * x - b

    # TODO: currently need to subtract b ⋅ b to make it so that constant is zero (due to limitation in MOI 0.3)
    @objective(model, Minimize, residual ⋅ residual - b ⋅ b)
    @constraint(model, C * x == d)
    @test_throws ArgumentError @constraint(model, C * x ≈ d)

    for i = 1 : 100
        allocs = @allocated solve!(model)
        if i > 1
            @test allocs == 0
        end
        test_equality_constrained(A(), b(), C(), d(), value.(model, x))
    end
end

@testset "Model: box constrained" begin
    # Minimize || x - p ||^2 = x' I x - 2 * p' * x + p ' * p
    # subject to l <= x <= p
    # where p is outside {x : l <= x <= p}

    n = 10

    optimizer = defaultoptimizer()
    model = Model(optimizer)
    x = [Variable(model) for _ = 1 : n]

    rng = MersenneTwister(1234)
    l = let rng = rng # https://github.com/JuliaLang/julia/issues/15276
        Parameter(zeros(n), model) do l
            @. l = -rand(rng)
            l
        end
    end
    u = let rng = rng # https://github.com/JuliaLang/julia/issues/15276
        Parameter(u -> rand!(rng, u), zeros(n), model)
    end
    p = let l = l, u = u, rng = rng # https://github.com/JuliaLang/julia/issues/15276
        Parameter(zeros(n), model) do p
            if rand(rng, Bool)
                p .= 2 .* l()
            else
                p .= 2 .* u()
            end
            p
        end
    end

    residual = @expression x - p
    # TODO: currently need to subtract p ⋅ p to make it so that constant is zero (due to limitation in MOI 0.3)
    @objective(model, Minimize, residual ⋅ residual - p ⋅ p)
    @constraint(model, x >= l)
    @constraint(model, x <= u)

    for testnum = 1 : 100
        allocs = @allocated solve!(model)
        expected = p() ./ 2
        @test value.(model, x) ≈ expected rtol = 1e-4
        testnum > 1 && @test allocs == 0
    end
end

end # module
