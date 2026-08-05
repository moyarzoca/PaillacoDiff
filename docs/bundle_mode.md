# Bundle Mode (Association-based)

In bundle mode you encapsulate everything in an Association. Tensors are stored as sparse Associations of independent components and computed on demand.

## Expected Bundle Keys

### Metric bundle (coordinate basis)

| Key | Required | Type | Purpose |
|---|---|---|---|
| `"ds2"` | Yes | Expression | Line element string (e.g. `-f[r] d[t]^2 + ...`) |
| `"coord"` | Yes | List | Coordinate variables |
| `"constants"` | No | List | Symbols with zero exterior derivative |
| `"assum"` | No | Expression | Assumptions for simplification (e.g. `r > 0`) |

### Vielbein bundle (orthonormal frame)

| Key | Required | Type | Purpose |
|---|---|---|---|
| `"basis"` | Yes | List | Vielbein 1-form symbols (e.g. `{e[1], e[2], e[3], e[4]}`) |
| `"eU"` | Yes | List | Explicit vielbein 1-forms as differential expressions |
| `"coordbasis"` | Yes | List | Coordinate 1-forms `d[coord]` |
| `"signature"` | Yes | List | Metric signature, e.g. `{-1, 1, 1, 1}` |
| `"constants"` | No | List | Constant symbols |
| `"assum"` | No | Expression | Simplification assumptions |

## Initialization

```mathematica
bundle = <|
  "ds2" -> ds2,
  "coord" -> {t, r, θ, φ},
  "constants" -> {M, Q}
|>;
InitializeBundle[bundle]
```

`InitializeBundle` converts the list to an Association and sets `d[c] = 0` for every constant.

## Tensor Levels

`PaiComputeBundleTensors[bundle, level]` computes tensors up to the requested level, skipping anything already present:

```
level: "metric" → "ChrisUdd" → "Rdddd" → "Rdd" → "RicciScalar"
```

| Level | Computes | Stored in `bundle["Tensors"]` |
|---|---|---|
| `"metric"` | `gdd`, `gUU` | `<|"gdd" → Assoc, "gUU" → Assoc|>` |
| `"ChrisUdd"` | Christoffel symbols | `<|"ChrisUdd" → Assoc|>` |
| `"Rdddd"` | Riemann tensor | `<|"Rdddd" → Assoc|>` |
| `"Rdd"` | Ricci tensor | `<|"Rdd" → Assoc|>` |
| `"RicciScalar"` | Ricci scalar | `<|"RicciScalar" → scalar|>` |

Tensors are stored as sparse Associations of independent components (exploiting symmetries), not full arrays. Use `GetTensorArray[bundle, name]` to retrieve a dense array.

```mathematica
PaiComputeBundleTensors[bundle, "RicciScalar"]  (* computes all dependencies *)
GetTensorArray[bundle, "gdd"]                   (* returns Dim × Dim dense array *)
```

## Building Operators

### Covariant derivative

```mathematica
nablaT = PaiCovD[bundle, Tdd, "dd"]
```

`PaiCovD` computes the coordinate-basis covariant derivative of a tensor array. The `indices` string gives the tensor index variance with `"U"` for upper and `"d"` for lower indices. The returned array has the derivative index first.

### Hodge star

```mathematica
hodge = BuildHodge[bundle, Simplify];
hodge[F]  (* returns *F *)
```

Dispatches to `BuildHodgeMetric` or `BuildHodgeVielbein` based on whether `bundle["UseVielbein"]` is `True`.

### Form squares

```mathematica
sq = BuildSquaresTools[bundle, Simplify];
sq["FormSquare"][F]     (* F_{μ1…μp} F^{μ1…μp} *)
sq["FormSquaredd"][F]   (* F_{μ λ2…λp} F_ν^{ λ2…λp} *)
```

## Vielbein Bundle Workflow

```mathematica
eIN = {eT[t, r] d[t], eR[t, r] d[r], r d[θ], r Sin[θ] d[φ]};
vielbeinBundle = <|
  "basis" -> {e[1], e[2], e[3], e[4]},
  "eU" -> eIN,
  "coordbasis" -> {d[t], d[r], d[θ], d[φ]},
  "signature" -> {-1, 1, 1, 1},
  "constants" -> {M, Q}
|>;
InitVielbein[vielbeinBundle];
PaiComputeBundleTensors[vielbeinBundle, "RicciScalar", Simplify]
```

The vielbein director (`PaiComputeBundleTensorsVielbein`) follows this dependency chain:

```
spinConnection → curvatureForm → RddddFlat → RddFlat → RicciScalarFlat
```

| Function | Computes | Stored in bundle |
|---|---|---|
| `PaiComputeSpinConnection[bundle]` | Spin connection 1-form `ω_{ab}` | `bundle["spinConnection", "dd"]` |
| `PaiComputeCurvatureForm[bundle]` | Curvature 2-form `Ω^a_b` | `bundle["curvatureForm", "Ud"]` |
| `PaiComputeRddddFlat[bundle]` | Riemann `R_{abcd}` | `bundle["FlatTensors", "Rdddd"]` |
| `PaiComputeRddFlat[bundle]` | Ricci `R_{ab}` | `bundle["FlatTensors", "Rdd"]` |
| `PaiComputeRicciScalarFlat[bundle]` | Ricci scalar | `bundle["FlatTensors", "RicciScalar"]` |

## Utility Functions

| Function | Purpose |
|---|---|
| `TensorProductContract[t1, t2, {{i1,j1}, ...}]` | Contracted tensor product |
| `RaiseIndices[sparse, bundle, {pos1, pos2, ...}]` | Raise specified indices using bundle metric |
| `PaiCovD[bundle, tensor, indices]` | Compute the coordinate-basis covariant derivative; derivative index is first |
| `GetTensorArray[bundle, name]` | Return dense array for a tensor (computes if needed) |
| `InitializeBundle[bundle]` | Turn list into Association, set d[constants] = 0 |
