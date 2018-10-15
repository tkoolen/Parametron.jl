# Parametron.jl

Parametron makes it easy to set up and efficiently (ideally, with zero allocation) solve instances of a parameterized family of optimization problems.

Parametron interfaces with various solvers through MathOptInterface.jl, similar to JuMP. However, unlike JuMP, the focus of Parametron is on efficient and user-friendly problem modification.

## Installation

### Installing Julia

Download links and more detailed instructions are available on the [Julia website](http://julialang.org/). The latest version of Parametron.jl requires Julia 0.7, but we recommend downloading 1.0 (the latest stable Julia release at the time of writing).

!!! warning

    Do **not** use `apt-get` or `brew` to install Julia, as the versions provided by these package managers tend to be out of date.

### Installing Parametron

To install the latest tagged release of Parametron, start Julia and enter `Pkg` mode by pressing `]`. Then simply run

```julia
add Parametron
```

To use the latest master version and work on the bleeding edge (generally, not recommended), instead run

```julia
add Parametron#master
```

A third option is to clone the repository (to the directory printed by `julia -e 'import Pkg; println(Pkg.devdir())'`):

```julia
dev Parametron
```

## Usage

See the [README](https://github.com/tkoolen/Parametron.jl/blob/master/README.md) for usage examples.

## API

See the API section for detailed documentation of exported functions.
