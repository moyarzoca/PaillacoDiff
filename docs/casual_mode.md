# Casual Mode (Global Variables)

In casual mode you set global variables and call functions directly. Best for quick interactive exploration.

## Global Variables You Set

| Variable | Type | Purpose | Used By |
|---|---|---|---|
| `coord` | List | Coordinate variables `{t, r, θ, φ}` | `DiffToMatrix`, `ComputeChrisUdd`, `ComputeRdd`, `ComputeRicciScalar`, `FormSquare`, `FormSquaredd`, `Hstar`, `DNAofForm`, `FormToSparse`, `coordcontraction` |
| `Dim` | Integer | Spacetime dimension | `MyHStarE`, `Contractione`, `GenerateGamma` |
| `gdd` | Matrix | Metric tensor `g_{μν}` | `ComputeChrisUdd`, `ComputeRdd`, `ComputeRicciScalar` |
| `gUU` | Matrix | Inverse metric `g^{μν}` | `FormSquare`, `FormSquaredd`, `Hstar` |
| `sqrtdetg` | Scalar | `Sqrt[-Det[gdd]]` | `Hstar` |
| `ChrisUdd` | Array | Christoffel symbols (set by `ComputeChrisUdd`) | `ComputeRdd`, `ComputeRicciScalar` |
| `Rdd` | Array | Ricci tensor (set by `ComputeRdd`) | `ComputeRicciScalar` |
| `RicciScalar` | Scalar | Ricci scalar (set by `ComputeRicciScalar`) | — |

### Vielbein-specific globals

| Variable | Type | Purpose | Set By |
|---|---|---|---|
| `\[Eta]dd` | Matrix | Flat (Minkowski) metric | `SetVielbein` |
| `\[Eta]UU` | Matrix | Inverse flat metric | `SetVielbein` |
| `e` | Symbol | Vielbein 1-forms `e^a` | `SetVielbein` |
| `eBasis` | List | `{e[1], ..., e[Dim]}` | `SetVielbein` |
| `eTodx` | Rule | `e^a → e^a_μ dx^μ` | `SetVielbein` |
| `dxToe` | Rule | `dx^μ → dx^μ(e) e^a` | `SetVielbein` |
| `eamuUd` | Matrix | Vielbein `e^a_μ` | `SetVielbein` |
| `eamudU` | Matrix | Inverse vielbein `e_a^μ` | `SetVielbein` |
| `\[Omega]Ud` | Array | Spin connection `ω^a_b` | `ComputeSpinConnection` |
| `\[Omega]dd` | Array | Spin connection `ω_{ab}` | `ComputeSpinConnection` |

### Clifford/Pauli globals

| Variable | Type | Purpose |
|---|---|---|
| `idGam` | Symbol | Identity in Clifford algebra |
| `id2` | Symbol | 2×2 identity in Pauli sector |
| `GU` | List | Explicit Gamma matrices (set by `GenerateGamma`) |

## Typical Workflow

```mathematica
(* 1. Set coordinates and dimension *)
coord = {t, r, θ, φ};
Dim = 4;

(* 2. Define the line element *)
ds2 = -f[r] d[t]^2 + d[r]^2/f[r] + r^2 (d[θ]^2 + Sin[θ]^2 d[φ]^2);
f[r_] := 1 - 2 M/r + Q^2/r^2/2;
d[Q] = d[M] = 0;

(* 3. Extract metric and its inverse *)
gdd = DiffToMatrix[ds2];
gUU = Inverse[gdd];
sqrtdetg = Sqrt[-Det[gdd]];

(* 4. Compute curvature *)
ComputeRicciScalar[]

(* 5. Work with forms *)
A = Q/r d[t];
F = d[A];
Hstar[F]
FormSquare[F]
FormSquaredd[F]

(* 6. Verify Einstein equations *)
EEdd = Rdd - 1/2 gdd RicciScalar - (FormSquaredd[F] - 1/4 gdd FormSquare[F]);
Simplify[EEdd]
```

## Functions that Default to Globals

| Function | Default globals |
|---|---|
| `DiffToMatrix[ds2]` | `coord` |
| `FormSquare[X]` | `gUU`, `coord` |
| `FormSquaredd[X]` | `gUU`, `coord` |
| `Hstar[X]` | `gUU`, `sqrtdetg`, `coord` |
| `DNAofForm[X]` | `coord` |
| `FormToSparse[X]` | `coord` |
| `FormsToMatrix[X]` | `coord` |
| `coordcontraction[X]` | `coord` |
| `ComputeChrisUdd[]` | `gdd`, `coord` |
| `ComputeRdd[]` | `gdd`, `coord` |
| `ComputeRicciScalar[]` | `gdd`, `coord` |

All these also accept explicit arguments to override the globals:

```mathematica
ComputeChrisUdd[Simplify, {myGdd, myCoord}]
FormSquare[X, mygUU, Simplify, myBasis]
```
