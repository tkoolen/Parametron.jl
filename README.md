# Parametron

[![Build Status](https://travis-ci.org/tkoolen/Parametron.jl.svg?branch=master)](https://travis-ci.org/tkoolen/Parametron.jl)
[![codecov.io](http://codecov.io/github/tkoolen/Parametron.jl/coverage.svg?branch=master)](http://codecov.io/github/tkoolen/Parametron.jl?branch=master)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://tkoolen.github.io/Parametron.jl/latest)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://tkoolen.github.io/Parametron.jl/stable)

Parametron makes it easy to set up and efficiently (ideally, with *zero* allocation) solve instances of a **parameterized family** of optimization problems.

## Example 1

As an example, we'll use the [OSQP](https://github.com/oxfordcontrol/OSQP.jl) solver to solve the following quadratic program:

```
Minimize ||A x - b||^2
subject to C x = d
```

with decision variable vector `x`, and where `A`, `b`, `C`, and `d` are parameters with random values, to be re-sampled each time the problem is re-solved.

Here's the basic problem setup:

```julia
# create a MathOptInterface optimizer instance
using OSQP
optimizer = OSQP.Optimizer()

# create a Parametron.Model, which holds problem information
using Parametron
using Random, LinearAlgebra
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
julia> using MathOptInterface; using OSQP.MathOptInterfaceOSQP: OSQPSettings; MathOptInterface.set(optimizer, OSQPSettings.Verbose(), false) # silence the optimizer

julia> using BenchmarkTools

julia> @btime solve!($model)
  51.863 μs (0 allocations: 0 bytes)
```

The performance and lack of allocations stems from the fact that the 'lazy expressions' used for the constraints and objective function are automatically optimized to calls to in-place functions.

## Example 2

Of course, in many real-world problems you are unlikely to update your parameters with random values.
Here's an illustration showing how you might control these values more directly, fitting a vector
`g` in a model

```julia
g' * X[:,i] ≈ p[i]
```

for a set of vectors in columns of `X`.

This example also demonstrates a different style of updating parameters. Whereas in the previous example we
supplied an 'update function' (e.g., `rand!`) that is automatically called when `solve!` is
called, in this example we use the syntax

```julia
Parameter(model, val=some_manually_updated_mutable_object)
```

to create a `Parameter` whose value may be updated manually between calls to the `solve!` function.

```julia
using Parametron, OSQP.MathOptInterfaceOSQP
using Random

n, m = 5, 15
Xdata = randn(n, m)
pdata = Vector{Float64}(undef, m);
model = Model(OSQP.Optimizer())
X = Parameter(model, val=Xdata)
p = Parameter(model, val=pdata)
g = [Variable(model) for _ = 1:n]
resid = @expression X'*g - p
@objective(model, Minimize, resid'*resid)

# Try with a specific ground-truth `g`
ggt = randn(n)
pdata .= Xdata'*ggt  # set the values in-place using `.=`
solve!(model)

julia> value.(model, g)
5-element Array{Float64,1}:
  0.6710700783457044
  1.3999896266657308
  0.5666247642146109
 -1.018123491138979
 -0.7464853233374451

julia> ggt
5-element Array{Float64,1}:
  0.671068170731507
  1.399985646860983
  0.5666231534233734
 -1.0181205969900424
 -0.7464832010803155
```

You can re-fit the model after updating `pdata` and/or `Xdata` in-place.
