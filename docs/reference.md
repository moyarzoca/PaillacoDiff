# PaillacoDiff Reference

## 1. Form Manipulation (lines 51–112)

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `FormDegree` | **Public** | Returns the degree of a differential form | 0 for scalars, 1 for one-forms, additive for Wedge/Times. Handles `e[n]`, `Wedge`, `Times`, `Plus`, `List`. |
| `Wedge` | **Public** | Exterior (wedge) product | Flat, OneIdentity. Handles associativity, zero factors, scalars, antisymmetry (identical odd-degree → 0), reordering with `Signature`. Routes Gamma factors through `CenterDot`. |
| `WedgeDot` | Private | Inner product of form lists | Used internally for `Wedge[A_List, B_List]`. |
| `d` | **Public** | Exterior derivative | Listable. Graded Leibniz rule on Wedge/Times. Acts on scalars via partial derivatives. Zero on numeric/d constants. |
| `IsThisGood` | **Dead** | Debugging helper | Checks if a polyform has uniform degree. Candidate for removal. |

## 2. Gamma and Sigma Matrices (lines 64–287)

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `Gam` | **Public** | Abstract Dirac gamma matrix `Gam[i]` | Indexed symbol used in abstract Clifford algebra computations. |
| `GamWeight` | Private | Returns weight of an expression (0 = non-Gamma, 1 = Gamma, 2 = Sigma) | Internal routing for CenterDot vs Dot selection. |
| `GamQ` | **Public** | Predicate: is expression Gamma-weighted? | Used to decide routing in `Wedge` and `CenterDot`. |
| `idGam` | **Global** | Identity element in Clifford algebra | Used as neutral element in Gamma products. |
| `CenterDot` | **Public** | Abstract multiplication for Gamma matrices | Flat, OneIdentity. Distributes over Plus, handles scalars. |
| `auxCenterDot` | Private | Helper for CenterDot on TensorProduct | Routes Gamma→CenterDot, Pauli→Dot inside TensorProduct factors. |
| `auxproduct` | Private | Selects CenterDot or Dot per factor | Used by `auxCenterDot`. |
| `TensorProduct` | **Public** | Overloaded Kronecker product | Handles scalars, Plus distribution. Also used for abstract tensor products of Gamma/Pauli sectors. |
| `Dot` | **Public** | Overloaded `.` operator | Extended to handle Gamma-weighted scalars. |
| `simpGamma` | Private | Simplifies Gamma matrix products | Pairs identical Gamma indices via Minkowski metric `\[Eta]dd`. Core of Gamma simplification. |
| `simpGam2` | **Dead** | Alternative Gamma simplifier | Uses `simpGamList`. Has debug `Print` statements. Unused in the rest of the code. |
| `simpGamList` | **Dead** | Helper for simpGam2 | Splits sorted index list into pairs for contraction. |
| `\[Sigma]rules` | Private | Pauli matrix multiplication rules | `σ[i].σ[i] → id2`, cyclic products `σ[1].σ[2] → I σ[3]`, etc. |
| `prodGamSig` | Private | Safe product of Gamma × Sigma expressions | Decomposes TensorProduct via `RuleFromTensorProd`, simplifies each pair. |
| `RuleFromTensorProd` | Private | Decomposes a TensorProduct into basis→coefficient rules | Used by `prodGamSig`. |
| `getCoefTensorProd` | Private | Extracts a single rule from a TensorProduct | Helper for `RuleFromTensorProd`. |
| `GenerateGamma` | **Public** | Builds explicit Gamma matrices | Constructs Gamma matrices as Kronecker products of Pauli matrices. Supports `"Symbolic"` (abstract TensorProduct) and numeric modes. Sets global `GU`. |
| `CliffordMap` | **Public** | Maps a form in the vielbein basis to Clifford algebra | `X /. dxToe /. Wedge → Inactive[CenterDot] /. e[i] :> Gam[i] // Activate`. |

### Global Variables (Section 2)

| Variable | Purpose |
|---|---|
| `idGam` | Identity in Clifford algebra |
| `id2` | 2×2 identity (Pauli sector) |
| `GU` | Explicit Gamma matrix representation (set by `GenerateGamma`) |
| `\[Eta]dd` | Flat Minkowski metric (used by `simpGamma`) |

## 3. Tensor Algebra and Contractions (lines 292–376)

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `Extractor` | **Public** | Contraction operator: given `F = F1 + F2 ^ A`, returns `F2` | Extracts a 1-form `A` from a polyform. Supports `"right"` and `"left"` sides. |
| `Extractorleft` | **Public** | Left-side contraction operator | Like `Extractor` with `side → "left"` by default. |
| `coordcontraction` | **Public** | Contracts a form with all coordinate 1-forms | Returns a list of (p-1)-forms. Defaults to global `coord`. |
| `GenerateGamma` | **Public** | (see Section 2) | |

## 4. Hodge Star and Related Operations (lines 389–803)

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `PolyFormQ` | **Public** | Tests if an expression is a polyform (mixed degrees) | Returns `True` if terms have different `FormDegree`. |
| `coeffBaseElement` | Private | Extracts coefficient and index list from a single form term | Internal building block for `DNAofForm`. |
| `DNAofForm` | **Public** | Decomposes a form into DNA: list of `{coefficient, sorted-index-list}` | The core component decomposition. Uses global `coord` as fallback. |
| `SparseFromDNA` | **Public** | Converts DNA representation to a `SparseArray` | Antisymmetrizes by permuting indices with signatures. |
| `DNAFromSparse` | **Public** | Converts a `SparseArray` back to DNA | Inverse of `SparseFromDNA`. |
| `nonzeroCompsFirst` | Private | Finds unique sorted index sets with a fixed first index | Helper for `FormSquaredd`. |
| `nonzeroComps` | Private | Finds all unique sorted index sets | Helper for `FormSquare` and `Hstar`. |
| `getCompRule` | Private | Looks up a component from a rules list | Helper. |
| `RaiseAllSparse` | Private | Raises all indices of a sparse form tensor using `gUU` | Used by `FormSquare` and `Hstar`. |
| `BuildSquaresTools` | **Public** | Constructs `FormSquare` and `FormSquaredd` closures for a bundle | Returns `<|"FormSquare" → fn, "FormSquaredd" → fn|>`. Handles metric and vielbein bundles. |
| `FormSquare` | **Public** | Computes the scalar squared norm of a p-form | `F_{μ1…μp} F^{μ1…μp}`. Uses DNA/SparseArray machinery. |
| `FormSquaredd` | **Public** | Computes the symmetric 2-tensor `F_{μ λ2…λp} F_ν^{ λ2…λp}` | Returns a `Dim × Dim` matrix. |
| `Hstar` | **Public** | Hodge star operator in coordinate basis | `(*F)_{ν1…ν_{D-p}} = (√\|g\| / p!(D-p)!) F^{μ1…μp} ε_{μ1…μp ν1…ν_{D-p}} dx^{ν1} ^ …`. |
| `FormToSparse` | **Public** | Converts a form to a `SparseArray` | Uses `DNAofForm` + `SparseFromDNA`. |
| `FormsToMatrix` | **Public** | Converts a form to a dense matrix (Normal array) | Wrapper around `FormToSparse`. |
| `MyHStarE` | Private | Hodge star in the vielbein basis | Superseded by `BuildHodgeVielbein`. Uses `\[Eta]UU` and `e`. |

### Global Variables (Section 4)

| Variable | Purpose |
|---|---|
| `gdd` | Metric tensor (used by `Hstar` as fallback) |
| `gUU` | Inverse metric (used by `FormSquare`, `FormSquaredd`, `Hstar`) |
| `coord` | Coordinates (used by `DNAofForm`, `Hstar`, etc.) |
| `sqrtdetg` | √\|det(g)\| (used by `Hstar` as fallback) |
| `\[Eta]UU` | Inverse flat metric (used by `MyHStarE`) |
| `e` | Vielbein basis 1-forms (used by `MyHStarE`) |
| `Dim` | Dimension (used by `MyHStarE`) |

## 5. Riemannian Geometry Tools (lines 806–991)

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `ClearGeometric` | **Public** | Clears global `ChrisUdd`, `Rdd`, `RicciScalar` | Convenience utility. |
| `DiffToMatrix` | **Public** | Extracts metric tensor `g_{μν}` from a `ds^2` line element | Parses the symbolic squared line element. |
| `Computegdd` | **Public** | Computes `gdd` and `sqrtdetg` for a bundle | Mutates a bundle Association. |
| `ComputeChrisUdd` | **Public** | Computes Christoffel symbols `Γ^μ_{νρ}` | Defaults to globals `{gdd, coord}`. Stores result in global `ChrisUdd`. |
| `ComputeRdd` | **Public** | Computes the Ricci tensor | Defaults to globals. Uses `ComputeChrisUdd` internally. Stores result in global `Rdd`. |
| `ComputeRicciScalar` | **Public** | Computes the Ricci scalar `R = g^{μν} R_{μν}` | Defaults to globals. Stores result in global `RicciScalar`. |
| `SetVielbein` | **Public** | Sets up the vielbein basis | Defines global `\[Eta]dd`, `\[Eta]UU`, `eTodx`, `dxToe`, `eBasis`, `gdd`, `gUU`, `eamuUd`, `eamudU`. |
| `ComputeSpinConnection` | **Public** | Computes spin connection 1-form `ω^a_b` | Calls `SetVielbein` then `ComputeChrisUdd`. Stores global `\[Omega]Ud`, `\[Omega]dd`. |
| `PrintIndex` | Private | Pretty-prints a tensor with subscripts/superscripts | Internal formatting. |
| `PrintIndices` | Private | Pretty-prints a tensor with a list of dn/up indicators | Internal formatting. |
| `$chrisdef` | Private | Display string for the Christoffel symbol definition | Internal formatting. |
| `$defRicciTensor` | Private | Display string for the Ricci tensor definition | Internal formatting. |
| `$defRicciScalar` | Private | Display string for the Ricci scalar definition | Internal formatting. |

### Global Variables (Section 5)

| Variable | Purpose |
|---|---|
| `coord` | Coordinate list (default for tensor computations) |
| `Dim` | Spacetime dimension |
| `gdd` | Metric tensor `g_{μν}` |
| `gUU` | Inverse metric `g^{μν}` |
| `ChrisUdd` | Christoffel symbols `Γ^μ_{νρ}` |
| `Rdd` | Ricci tensor `R_{μν}` |
| `RicciScalar` | Ricci scalar `R` |
| `\[Eta]dd` | Flat metric (Lorentzian) |
| `\[Eta]UU` | Inverse flat metric |
| `e` | Vielbein 1-forms `e^a` |
| `eTodx` | Replacement rule `e^a → e^a_μ dx^μ` |
| `dxToe` | Replacement rule `dx^μ → dx^μ(e) e^a` |
| `eamuUd` | Vielbein matrix `e^a_μ` |
| `eamudU` | Inverse vielbein `e_a^μ` |
| `eBasis` | List of `e[1], …, e[Dim]` |
| `\[Omega]Ud` | Spin connection 1-form `ω^a_b` |
| `\[Omega]dd` | Spin connection with both indices down `ω_{ab}` |

## 6. Bundle Framework (lines 1020–1612)

### Core Bundle Functions

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `InitializeBundle` | **Public** | Initializes a bundle Association | Sets `d[constant] = 0` for all listed constants. |
| `TensorProductContract` | **Public** | Contracted tensor product | Wrapper around `TensorContract[Inactive[TensorProduct][…], …]` with `Activate`. |
| `RaiseIndices` | **Public** | Raises specified indices of a sparse tensor using bundle metric | Uses `GetTensorArray[bundle, "gUU"]`. |
| `PaiCovD` | **Public** | Computes the coordinate-basis covariant derivative of a tensor array | Uses `GetTensorArray[bundle, "ChrisUdd"]`. The `indices` string uses `U`/`d`; the derivative index is first. |
| `GetTensorArray` | **Public** | Retrieves a tensor array from a bundle, computing it on demand | Orchestrates via `PaiComputeBundleTensors`. Returns dense array. |
| `CleanZeros` | Private | Removes zero entries from an Association | `KeySelect[X, X[#] =!= 0 &]`. |
| `PaiComponent` | Private | Dispatcher for component lookup by tensor name | Routes to symmetric/Riemann lookups. |
| `tensorRank` | Private | Association mapping tensor names → rank | `{"gdd" → 2, "gUU" → 2, "Rdd" → 2, "ChrisUdd" → 3, "Rdddd" → 4}`. |
| `PaiComponent2sym` | Private | 2-index symmetric component lookup | Exploits symmetry: `{i,j}` → `{Min[i,j], Max[i,j]}`. |
| `PaiComponent3symLast` | Private | 3-index component lookup symmetric in last two indices | Used for Christoffel. |
| `PaiComponent4Riem` | Private | 4-index Riemann component lookup | Uses pair-exchange symmetries. |
| `failRequirements` | Private | Checks if a bundle has all required keys | `Not[And @@ Map[KeyExistsQ[bundle, #] &, req]]`. |

### Tensor Computation Functions

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `PaiComputeMetric` | **Public** | Computes `gdd` and `gUU` Associations from a bundle's `"ds2"` | Stored as sparse Associations of independent components. |
| `PaiComputeChrisUdd` | **Public** | Computes Christoffel symbols from bundle metric | Uses Accessor functions `gdd[i,j]`, `gUU[i,j]` with `PaiComponent`. |
| `PaiComputeRdddd` | **Public** | Computes the Riemann tensor `R_{μνρσ}` | Assumes `ChrisUdd` is already in bundle. |
| `PaiComputeRdd` | **Public** | Computes the Ricci tensor from the Riemann tensor | Contracts `g^{μα} R_{αβνρ}`. |
| `PaiComputeRicciScalar` | **Public** | Computes the Ricci scalar from the Ricci tensor | `Sum[g^{ij} R_{ij}]`. |

### Bundle Orchestrator

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `PaiComputeBundleTensors` | **Public** | Main director: computes tensors up to a requested level | Dispatches to metric or vielbein path. Levels: `metric`, `ChrisUdd`, `Rdddd`, `Rdd`, `RicciScalar`. |
| `PaiComputeBundleTensorsMetric` | **Public** | Metric-branch director | Lazy step-by-step computation: metric → ChrisUdd → Rdddd → Rdd → RicciScalar. |
| `PaiComputeBundleTensorsVielbein` | **Public** | Vielbein-branch director | Steps: spinConnection → curvatureForm → Rdddd → Rdd → RicciScalar. |

### Array Extraction

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `GetArrayRdd` | **Public** | Returns dense Ricci tensor array from bundle | Convenience bridge between bundle and array form. |
| `NewComputeRdd` | **Public** | Alias/HoldFirst wrapper | (Used in tuple with `GetArrayRdd`) |

### Hodge Builders

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `BuildHodge` | **Public** | Dispatches Hodge star construction | Routes to `BuildHodgeMetric` or `BuildHodgeVielbein`. |
| `BuildHodgeMetric` | **Public** | Builds a Hodge star function for a coordinate-basis bundle | Returns `Function[{X}, Hstar[X, gUU, sqrtdetg, d[coord], simp]]`. |
| `BuildHodgeVielbein` | **Public** | Builds a Hodge star function for a vielbein bundle | Returns `Function[{X}, Hstar[X, etainv, sqrtdeteta, basis, simp]]`. |

### Vielbein Bundle Functions

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `ConstructContraction` | Private | Builds a contraction operator for a vielbein bundle | Returns `Function[{X}, Table[Extractor[X, eIter, "left"], {eIter, symbs}]]`. |
| `InitVielbein` | **Public** | Initializes a vielbein bundle | Sets up `eTodx`, `dxToe`, `contraction`, assigns `d[eIter]`. |
| `PaiComputeSpinConnection` | **Public** | Computes spin connection 1-form `ω_{ab}` in bundle | Uses contraction on `d[e]`. Stores in `bundle["spinConnection"]`. |
| `PaiComputeCurvatureForm` | **Public** | Computes curvature 2-form `Ω^a_b = dω^a_b + ω^a_c ∧ ω^c_b` | Stores in `bundle["curvatureForm"]`. |
| `PaiComputeRddddFlat` | **Public** | Computes Riemann in flat (vielbein) basis `R_{abcd}` | From curvature 2-form. |
| `PaiComputeRddFlat` | **Public** | Computes Ricci in flat basis `R_{ab}` | Contracted from `R_{abcd}`. |
| `PaiComputeRicciScalarFlat` | **Public** | Computes Ricci scalar in flat basis | `Tr[η^{ac} R_{cb}]`. |

### Additional Utilities

| Function | Visibility | Purpose | Notes |
|---|---|---|---|
| `inP` | Private | Interior product with a vielbein 1-form `e^a` | Superseded by `ConstructContraction`. |
| `Contractione` | Private | Contracts a form with all vielbein 1-forms | Superseded by `ConstructContraction`. |

## 7. Functions Removed

| Function | Reason |
|---|---|
| `IsThisGood` | Debugging leftover. Deleted. |
| `simpGam2` / `simpGamList` | Alternative Gamma simplifier, unused. Superseded by `simpGamma`. Deleted. |
| `GetArrayRdd` / `NewComputeRdd` | Unused bundle bridge. Deleted. |

## 8. Final Public API (after BeginPackage refactor)

### Public Functions (50)

`FormDegree`, `Wedge`, `d`, `PolyFormQ`, `Gam`, `GamQ`, `CenterDot`, `TensorProduct`, `Dot`, `simpGamma`, `GenerateGamma`, `CliffordMap`, `Extractor`, `Extractorleft`, `coordcontraction`, `DNAofForm`, `SparseFromDNA`, `DNAFromSparse`, `BuildSquaresTools`, `FormSquare`, `FormSquaredd`, `Hstar`, `FormToSparse`, `FormsToMatrix`, `ClearGeometric`, `DiffToMatrix`, `Computegdd`, `ComputeChrisUdd`, `ComputeRdd`, `ComputeRicciScalar`, `SetVielbein`, `ComputeSpinConnection`, `InitializeBundle`, `TensorProductContract`, `RaiseIndices`, `PaiCovD`, `GetTensorArray`, `PaiComputeMetric`, `PaiComputeChrisUdd`, `PaiComputeRdddd`, `PaiComputeRdd`, `PaiComputeRicciScalar`, `PaiComputeBundleTensors`, `BuildHodge`, `BuildHodgeMetric`, `BuildHodgeVielbein`, `InitVielbein`, `PaiComputeSpinConnection`, `PaiComputeCurvatureForm`, `PaiComputeRddddFlat`, `PaiComputeRddFlat`, `PaiComputeRicciScalarFlat`

### Public Globals (21)

`coord`, `Dim`, `gdd`, `gUU`, `ChrisUdd`, `Rdd`, `RicciScalar`, `sqrtdetg`, `\[Eta]dd`, `\[Eta]UU`, `e`, `eTodx`, `dxToe`, `eamuUd`, `eamudU`, `eBasis`, `\[Omega]Ud`, `\[Omega]dd`, `idGam`, `id2`, `GU`

### Moved to Private

`MyHStarE`, `inP`, `Contractione`, `PrintIndex`, `PrintIndices`, `$chrisdef`, `$defRicciTensor`, `$defRicciScalar`

### Deleted

`IsThisGood`, `simpGam2`, `simpGamList`, `GetArrayRdd`, `NewComputeRdd`
