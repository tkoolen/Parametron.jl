# Parametron

[![Build Status](https://travis-ci.org/tkoolen/Parametron.jl.svg?branch=master)](https://travis-ci.org/tkoolen/Parametron.jl)
[![codecov.io](http://codecov.io/github/tkoolen/Parametron.jl/coverage.svg?branch=master)](http://codecov.io/github/tkoolen/Parametron.jl?branch=master)

Parametron makes it easy to set up and efficiently (ideally, with *zero* allocation) solve instances of a **parameterized family** of optimization problems.

As an example, we'll use the [OSQP](https://github.com/oxfordcontrol/OSQP.jl) solver to solve the following quadratic program:

```
Minimize ||A x - b||^2
subject to C x = d
```

with decision variable vector `x`, and where `A`, `b`, `C`, and `d` are parameters with random values, to be re-sampled each time the problem is re-solved.

Here's the basic problem setup:

```julia
# create a MathOptInterface optimizer instance
using OSQP.MathOptInterfaceOSQP
optimizer = OSQPOptimizer()

# create a Parametron.Model, which holds problem information
using Parametron
model = Model(optimizer)

# create decision variables and parameters
n = 8; m = 2
x = [Variable(model) for _ = 1 : n]
A = Parameter(rand!, zeros(n, n), model)
b = Parameter(rand!, zeros(n), model)
C = Parameter(rand!, zeros(m, n), model)
d = Parameter(zeros(m), model) do d
    # do syntax makes it easy to create custom Parameters
    rand!(d)
    d .*= 2
end

# the @expression macro can be used to create 'lazy' expressions,
# which can be used in constraints or the objective function, and
# can be evaluated at a later time, automatically updating the
# Parameters in the process (if needed).
residual = @expression A * x - b

# set the objective function
@objective(model, Minimize, residual ⋅ residual)

# add the constraints. You could have multiple @constraint calls
# as well. ==, <=, and >= are supported.
@constraint(model, C * x == d)
```

Now that the problem is set up, we can solve and obtain the solution as follows:

```julia
julia> solve!(model)
-----------------------------------------------------------------
           OSQP v0.3.0  -  Operator Splitting QP Solver
              (c) Bartolomeo Stellato,  Goran Banjac
        University of Oxford  -  Stanford University 2017
-----------------------------------------------------------------
problem:  variables n = 8, constraints m = 2
          nnz(P) + nnz(A) = 88
settings: linear system solver = suitesparse ldl,
          eps_abs = 1.0e-03, eps_rel = 1.0e-03,
          eps_prim_inf = 1.0e-04, eps_dual_inf = 1.0e-04,
          rho = 1.00e-01 (adaptive),
          sigma = 1.00e-06, alpha = 1.60, max_iter = 4000
          check_termination: on (interval 25),
          scaling: on, scaled_termination: off
          warm start: on, polish: off

iter   objective    pri res    dua res    rho        time
   1  -7.8949e-01   9.57e-01   1.02e+03   1.00e-01   1.34e-04s
  25  -2.0032e+00   2.87e-04   4.82e-03   1.00e-01   1.76e-04s

status:               solved
number of iterations: 25
optimal objective:    -2.0032
run time:             1.81e-04s
optimal rho estimate: 5.16e-02

julia> value.(model, x)
8-element Array{Float64,1}:
 -0.365181
 -0.119036
 -0.267222
  1.41655
  0.69472
  0.993475
 -0.631194
 -1.02733
```

Note that the next time `solve!` is called, the update functions of our parameters (`A`, `b`, `C`, and `d`) will be called again (*once* for each parameter), resulting in a different optimum:

```julia
julia> solve!(model)
iter   objective    pri res    dua res    rho        time
   1  -1.4419e+00   2.57e-01   5.79e+02   1.00e-01   1.53e-05s
  25  -3.2498e+00   1.34e-04   2.74e-03   1.00e-01   3.10e-05s

status:               solved
number of iterations: 25
optimal objective:    -3.2498
run time:             3.63e-05s
optimal rho estimate: 7.79e-02

julia> value.(model, x)
8-element Array{Float64,1}:
  0.220836
 -0.462071
  0.509146
  0.667908
 -0.850638
  0.7877
  1.01581
 -0.992135
```

Note that the solver is warm-started. Also note that updating the parameters and solving a new QP instance is quite fast:

```julia
julia> MathOptInterface.set!(optimizer, OSQPSettings.Verbose(), false) # silence the optimizer

julia> using BenchmarkTools

julia> @btime solve!($model)
  85.077 μs (0 allocations: 0 bytes)
```

The performance and lack of allocations stems from the fact that the 'lazy expressions' used for the constraints and objective function are automatically optimized to calls to in-place functions.
