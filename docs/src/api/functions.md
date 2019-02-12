# Functions

```@meta
CurrentModule = Parametron.Functions
```

```@docs
Functions
```

## Types

```@docs
Variable
LinearTerm
QuadraticTerm
AffineFunction
QuadraticFunction
```

## Exported functions

```@docs
canonicalize
canonicalize(::QuadraticTerm)
canonicalize(::AffineFunction)
canonicalize(::QuadraticFunction)
canonicalize!
```

```@docs
prune_zero
prune_zero!
```

## In-place math functions (unexported)

```@docs
add!
subtract!
muladd!
vecdot!
vecadd!
vecsubtract!
matvecmul!
bilinearmul!
scale!
vcat!
```
