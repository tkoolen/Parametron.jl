# Model

## Creating a `Model`

```@docs
Model
```

## Creating Decision variables

```@docs
Variable(::Model)
```

## Adding constraints and an objective function

```@docs
@constraint
@objective
setobjective!
```

## Solving

```@docs
solve!
Parametron.setdirty!
Parametron.initialize!
Parametron.update!
```

## Accessing solver results

```@docs
value
objectivevalue
terminationstatus
primalstatus
dualstatus
```
