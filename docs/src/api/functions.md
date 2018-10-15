# Functions

```@meta
CurrentModule = Parametron.Functions
```

```@docs
Functions
```

## Types

```@docs
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
