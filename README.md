# SimpleQP

[![Build Status](https://travis-ci.org/tkoolen/SimpleQP.jl.svg?branch=master)](https://travis-ci.org/tkoolen/SimpleQP.jl)
[![codecov.io](http://codecov.io/github/tkoolen/SimpleQP.jl/coverage.svg?branch=master)](http://codecov.io/github/tkoolen/SimpleQP.jl?branch=master)

SimpleQP makes it easy to set up and efficiently (ideally, with *zero* allocations) solve instances of a **parameterized family** of quadratic programs.

As an example, we'll use the [OSQP](https://github.com/oxfordcontrol/OSQP.jl) solver is used to solve the following problem:

```
Minimize ||A x - b||^2
subject to C x = d
```

with decision variable vector `x`, and where `A`, `b`, `C`, and `d` are parameters with random values, to be re-sampled each time the problem is re-solved.

Here's the basic problem setup:
```julia
using OSQP.MathOptInterfaceOSQP
optimizer = OSQPOptimizer()

using SimpleQP
model = Model(optimizer)
n = 8; m = 2

x = [Variable(model) for _ = 1 : n]
A = Parameter(rand!, zeros(n, n), model)
b = Parameter(rand!, zeros(n), model)
C = Parameter(rand!, zeros(m, n), model)
d = Parameter(rand!, zeros(m), model)

residual = @expression A * x - b

# currently need to subtract b ⋅ b to make it so that constant is zero (due to limitation in MOI 0.3)
@objective(model, Minimize, residual ⋅ residual - b ⋅ b) 
@constraint(model, C * x == d)
```

after which we can solve the problem and obtain the solution as follows:

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

Now note that the next time `solve!` is called, the `rand!` function will be called again to update the parameters `A`, `b`, `C`, and `d`, resulting in a different optimum:

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

Note the solver warm start. Also note that setting up and solving the problem is very efficient:

```julia
julia> MathOptInterface.set!(optimizer, OSQPSettings.Verbose(), false) # silence the optimizer

julia> using BenchmarkTools

julia> @btime solve!($model)
  85.077 μs (0 allocations: 0 bytes)
```


