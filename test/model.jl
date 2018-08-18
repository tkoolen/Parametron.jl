module ModelTests

using Compat
using Compat.Test
using Compat.Random
using Compat.LinearAlgebra
using Parametron
using OSQP.MathOptInterfaceOSQP
using GLPK
using StaticArrays: SVector

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

function test_unconstrained(model, x, Q, r, s; atol=1e-8)
    solve!(model)
    @test terminationstatus(model) == MOI.Success
    @test primalstatus(model) == MOI.FeasiblePoint
    @test dualstatus(model) == MOI.FeasiblePoint
    xval = value.(Ref(model), x)
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
        Parameter(Matrix(1.0I, n, n), model) do Q
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

    Parametron.initialize!(model)
    @test model.initialized

    test_unconstrained(model, x, Q, r, s)
    solve!(model)
    allocs = @allocated solve!(model)
    @test allocs == 0

    # constant modification
    sval[] = 2.0
    test_unconstrained(model, x, Q, r, s)
    allocs = @allocated solve!(model)
    @test allocs == 0
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

    @objective(model, Minimize, residual ⋅ residual)
    @constraint(model, C * x == d)
    @test_throws ArgumentError @constraint(model, C * x ≈ d)

    for i = 1 : 100
        allocs = @allocated solve!(model)
        @test terminationstatus(model) == MOI.Success
        @test primalstatus(model) == MOI.FeasiblePoint
        if i > 1
            @test allocs == 0
        end
        test_equality_constrained(A(), b(), C(), d(), value.(Ref(model), x))
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
            l .= .-rand.(Ref(rng))
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
        @test value.(Ref(model), x) ≈ expected rtol = 1e-4
        testnum > 1 && @test allocs == 0
    end
end

@testset "Issue 30 #1" begin
    optimizer = defaultoptimizer()
    model = Model(optimizer)
    x = Variable(model)
    @constraint(model, [x] <= [-3.])
    @objective model Minimize x^2
    solve!(model)
    @test value(model, x) ≈ -3. atol = 1e-8
    allocs = @allocated solve!(model)
    @test allocs == 0
end

@testset "Issue 30 #2" begin
    optimizer = defaultoptimizer()
    model = Model(optimizer)
    constraintexpr(x, upperbound) = @expression [x] - [upperbound]
    x = Variable(model)
    @constraint model constraintexpr(x, -3.0) <= [0.0]
    @objective model Minimize x^2
    solve!(model)
    @test value(model, x) ≈ -3. atol = 1e-8
    allocs = @allocated solve!(model)
    @test allocs == 0
end

@testset "Issue 29" begin
    model = Model(defaultoptimizer())
    x = Variable(model)
    @constraint model [x] >= [0]
    @objective model Minimize x^2 + 1
    solve!(model)
    @test value(model, x) ≈ 0 atol=1e-8
    @test objectivevalue(model) ≈ 1 atol=1e-8
end

@testset "MOI Issue 426" begin
    model = Model(defaultoptimizer())
    x = Variable(model)
    for vector_set_type in (MOI.Nonnegatives, MOI.Nonpositives, MOI.Zeros)
        @test(length(@inferred(MOI.get(model.backend, MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64}, vector_set_type}()))) == 0)
    end
    @constraint model [x] >= [0]
    @constraint model [x] <= [1]
    @constraint model [x] == [0.5]
    for vector_set_type in (MOI.Nonnegatives, MOI.Nonpositives, MOI.Zeros)
        @test(length(@inferred(MOI.get(model.backend, MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64}, vector_set_type}()))) == 1)
    end
end

@testset "boolean basics" begin
    optimizer = GLPK.Optimizer()
    model = Model(optimizer)
    x = Variable(model)
    @constraint model x ∈ {0, 1}
    @objective model Maximize x

    solve!(model)

    @test terminationstatus(model) == MOI.Success
    @test primalstatus(model) == MOI.FeasiblePoint
    @test value(model, x) ≈ 1.0 atol=1e-8
end

@testset "integer basics" begin
    optimizer = GLPK.Optimizer()
    model = Model(optimizer)
    x = Variable(model)
    @constraint model x ∈ ℤ
    @constraint model [x] >= [0.5]
    @objective model Minimize x

    solve!(model)

    @test terminationstatus(model) == MOI.Success
    @test primalstatus(model) == MOI.FeasiblePoint
    @test value(model, x) ≈ 1.0 atol=1e-8
end

@testset "scalar constraints 1" begin
    optimizer = defaultoptimizer()
    model = Model(optimizer)
    x = Variable(model)
    @constraint(model, x <= -3)
    @objective model Minimize x^2
    solve!(model)
    @test value(model, x) ≈ -3. atol = 1e-8
    allocs = @allocated solve!(model)
    @test allocs == 0
end

@testset "scalar constraints 2" begin
    # from https://github.com/JuliaOpt/MathOptInterface.jl/blob/575647cfa344ed7656e0eee16537017184a011a8/src/Test/contquadratic.jl#L3
    optimizer = defaultoptimizer()
    model = Model(optimizer)
    x = Variable(model)
    y = Variable(model)
    z = Variable(model)
    @objective model Minimize x^2 + x * y + y^2 + y * z + z^2
    @constraint model x + 2y + 3z >= 4
    @constraint model x +  y      >= 1
    solve!(model)

    @test terminationstatus(model) == MOI.Success
    @test primalstatus(model) == MOI.FeasiblePoint
    @test objectivevalue(model) ≈ 13/7 atol=1e-6
    @test value.(Ref(model), [x, y, z]) ≈ [4/7, 3/7, 6/7] atol=1e-6

    allocs = @allocated solve!(model)
    @test allocs == 0
end

@testset "scalar constraints 3" begin
    # from https://github.com/JuliaOpt/MathOptInterface.jl/blob/575647cfa344ed7656e0eee16537017184a011a8/src/Test/contquadratic.jl#L160
    optimizer = defaultoptimizer()
    model = Model(optimizer)
    x = Variable(model)
    y = Variable(model)
    @objective model Minimize 2 * x^2 + y^2 + x * y + x + y + 1
    @constraint model x >= 0
    @constraint model -y <= 0
    @constraint model x + y == 1
    solve!(model)

    @test terminationstatus(model) == MOI.Success
    @test primalstatus(model) == MOI.FeasiblePoint
    @test objectivevalue(model) ≈ 2.875 atol=1e-6
    @test value.(Ref(model), [x, y]) ≈ [0.25, 0.75] atol=1e-6

    allocs = @allocated solve!(model)
    @test allocs == 0
end

@testset "AbstractVector constraints" begin
    model = Model(defaultoptimizer())
    x = SVector(Variable(model), Variable(model))
    @constraint model x == SVector(1.0, 2.0)
end

@testset "default objective (issue #62)" begin
    model = Model(defaultoptimizer())
    x = Variable(model)
    @constraint model x == 1
    solve!(model)
    @test terminationstatus(model) == MOI.Success
    @test primalstatus(model) == MOI.FeasiblePoint
    @test objectivevalue(model) == 0.0
    @test value(model, x) ≈ 1.0 atol=1e-6
end

end # module
