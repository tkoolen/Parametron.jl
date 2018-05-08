module SimpleQPTest

include("functions.jl")

# using Compat
# using Compat.Test
# using SimpleQP
# using OSQP.MathOptInterfaceOSQP
# import MathOptInterface
# const MOI = MathOptInterface

# function defaultoptimizer()
#     optimizer = OSQPOptimizer()
#     MOI.set!(optimizer, OSQPSettings.Verbose(), false)
#     MOI.set!(optimizer, OSQPSettings.EpsAbs(), 1e-8)
#     MOI.set!(optimizer, OSQPSettings.EpsRel(), 1e-16)
#     MOI.set!(optimizer, OSQPSettings.MaxIter(), 10000)
#     MOI.set!(optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) # required for deterministic behavior
#     optimizer
# end

# function test_unconstrained(model, x, Q, r, s; atol=1e-15)
#     solve!(model)
#     xval = value.(model, x)
#     expected = -2 * Q \ r'
#     @test xval ≈ expected atol = atol
#     @test objectivevalue(model) ≈ dot(xval, Q * xval) + dot(vec(r), xval) + s[1] atol = atol
# end

# @testset "Model: unconstrained" begin
#     optimizer = defaultoptimizer()
#     model = Model(optimizer)
#     n = 2
#     x = [Variable(model) for _ = 1 : n]
#     x1, x2 = x
#     @test x1.index == 1
#     @test x2.index == 2

#     Q = Symmetric(speye(n))
#     r = zeros(1, n)
#     s = [1.0]
#     objective = QuadraticTerm(Q, x) + LinearTerm(r, x) + Constant(s)

#     setobjective!(model, Senses.Min, objective)
#     SimpleQP.initialize!(model)
#     @test model.initialized

#     test_unconstrained(model, x, Q, r, s)

#     rand_data = function (rng)
#         Q[1, 1] = rand(rng)
#         Q[2, 2] = rand(rng)
#         r[1] = rand(rng) - 0.5
#         r[2] = rand(rng) - 0.5
#     end
#     rng = MersenneTwister(1)
#     for i = 1 : 10
#         rand_data(rng)
#         test_unconstrained(model, x, Q, r, s)
#     end

#     rand_data(rng)
#     allocs = @allocated solve!(model)
#     @test allocs == 0

#     s[1] = 2.0
#     @test_throws ArgumentError solve!(model) # can't modify constant offset
#     s[1] = 1.0
# end

# @testset "Model: equality constrained" begin
#     # Minimize ||A x - b||^2 = x' A' A x - (2 * A' * b)' x + b' * b
#     # subject to C x = d

#     rand_data! = function (A, b, C, d, rng)
#         rand!(rng, A)
#         rand!(rng, b)
#         rand!(rng, C)
#         rand!(rng, d)
#     end

#     update_objective_matrices! = function (P, qt, r, A, b)
#         P.data[:] = triu(A' * A)
#         qt[:] = -2 * b' * A
#         r[:] = b' * b
#     end

#     test_equality_constrained = function(A, b, C, D, x; rtol = 1e-4)
#         C⁺ = pinv(C)
#         Q = I - C⁺ * C
#         expected = Q * (pinv(A * Q) * (b - A * C⁺ * d)) + C⁺ * d # note: can be quite badly conditioned
#         @test x ≈ expected rtol = rtol
#     end

#     n = 8
#     m = 2

#     A = zeros(n, n)
#     b = zeros(n)
#     C = zeros(m, n)
#     d = zeros(m)

#     P = Symmetric(sparse(ones(n, n)))
#     qt = zeros(1, n)
#     r = [0.0]

#     optimizer = defaultoptimizer()
#     model = Model(optimizer)
#     x = [Variable(model) for _ = 1 : n]
#     objective = QuadraticTerm(P, x) + LinearTerm(qt, x) #+ r
#     setobjective!(model, Senses.Min, objective)
#     @constraint(model, LinearTerm(C, x) == Constant(d))

#     rng = MersenneTwister(1234)
#     for i = 1 : 100
#         rand_data!(A, b, C, d, rng)
#         update_objective_matrices!(P, qt, r, A, b)
#         allocs = @allocated solve!(model)
#         if i > 1
#             @test allocs == 0
#         end
#         test_equality_constrained(A, b, C, d, value.(model, x))
#     end
# end

# @testset "Model: box constrained" begin
#     # Minimize || x - p ||^2 = x' I x - 2 * p' * x + p ' * p
#     # subject to l <= x <= p
#     # where p is outside {x : l <= x <= p}

#     n = 10

#     l = [[0.0] for _ = 1 : n]
#     u = [[0.0] for _ = 1 : n]

#     optimizer = defaultoptimizer()
#     model = Model(optimizer)
#     x = [Variable(model) for _ = 1 : n]

#     P = Symmetric(speye(n))
#     qt = zeros(1, n)
#     r = [0.0]

#     objective = QuadraticTerm(P, x) + LinearTerm(qt, x) #+ r
#     setobjective!(model, Senses.Min, objective)
#     for i = 1 : n
#         @constraint(model, LinearTerm(fill(1.0, 1, 1), [x[i]]) >= Constant(l[i]))
#         @constraint(model, LinearTerm(fill(1.0, 1, 1), [x[i]]) <= Constant(u[i]))
#     end

#     rng = MersenneTwister(1234)
#     for testnum = 1 : 100
#         for i = 1 : n
#             l[i] .= -rand(rng)
#             u[i] .= rand(rng)
#         end
#         lvec = [l[i][1] for i = 1 : n]
#         uvec = [u[i][1] for i = 1 : n]

#         for vertex in [lvec, uvec]
#             p = 2 .* vertex
#             qt .= -2 .* p'
#             allocs = @allocated solve!(model)
#             @test value.(model, x) ≈ vertex rtol = 1e-4
#             testnum > 1 && @test allocs == 0
#         end
#     end
# end

end # module
