# PaillacoDiff Package Documentation

## Author
Marcelo Oyarzo

---

## Package Structure

The code is organized into the following sections:

### 1. Form Manipulation
- `FormDegree`
- `Wedge`
- `d`
- `PolyFormQ`, etc.

### 2. Gamma and Sigma Matrices
- `GamWeight`
- `GamQ`
- `CenterDot`
- `simpGamma`
- `prodGamSig`, etc.

### 3. Tensor Algebra and Contractions
- `TensorProduct`
- `Extractor`
- `coordcontraction`
- `RaiseIndices`

### 4. Hodge Star and Related Operations
- `Hstar`
- `MyHStar`
- `MyHStarE`
- `DNAofForm`
- `SparseFromDNA`, etc.

### 5. Riemannian Geometry Tools
- `DiffToMatrix`
- `Computegdd`
- `ComputeChrisUdd`
- `ComputeRdd`

### 6. Utilities
- `FormsToMatrix`
- `DNAtoMatrix`
- `DNAtoForms`

---
Tree of dependencies
---
#### `Wedge`
- Dependencies: `FormDegree`, `GamQ`

#### `d`
- Dependencies: `Wedge`, `FormDegree`

#### `PolyFormQ`
- Dependencies: `FormDegree`

#### `DNAofForm`
- Dependencies: `FormDegree`, `PolyFormQ`, `coeffBaseElement`
- Inputs: `X` (p-form), `coord` (coordinate list)
- Output: DNA representation of the p-form

#### `SparseFromDNA`
- Inputs: DNA of p-form, form degree `p`, `Dim` (dimension)
- Output: Sparse array representation of the form


#### `DNAFromSparse`
- Inputs: Sparse Array of a p-form
- Output: DNA of the p-form


#### `RaiseAllSparse`
- Inputs: Sparse array of a p-form, `gUU`, `p`
- Outputs: Sparse array of the p-form with upper indices

#### nonzeroComps
- Input: Sparse Array of a p-form
- Output: non-zero positions of the sparse array

#### getCompRule
- Input: tensor as Association or ArrayRules of sparse Array, basis label
- Output: component along the basis label

#### `FormSquare`
- Dependencies: `FormDegree`, `DNAofForm`, `SparseFromDNA`, `RaiseAllSparse`, `nonzeroComps`, `getcompRule`
- Input: p-form, `gUU`, simplification rule, `coord`
- Output: Square of the p-form (scalar)

#### `nonzeroCompsFirst`
- Input: Sparse array of a p-form, first index `mu` (integer)
- Output: List of the non-zero components of the tensor with first index `mu` and `mu` removed

#### `FormSquaredd`
- Dependencies: `FormDegree`, `DNAofForm`, `SparseFromDNA`, `nonzeroCompsFirst`, `getcompRule`
- Input: p-form, `gUU`, simplification rule, `coord`
- Output: symmetric rank 2 tensor square of the p-form
