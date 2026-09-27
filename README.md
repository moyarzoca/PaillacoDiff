# PaillacoDiff

PaillacoDiff is a Wolfram Language package for differential geometry and exterior algebra. It provides a unified set of tools for working with tensors and differential forms. It includes curvature tensors, Hodge duals, contractions, and related geometric operations, with conventions specified in [Conventions](https://moyarzoca.github.io/paillacodiff/#conventions). It also provides a compact interface for defining tensors using an explicit GRTensor-like notation.


## Installation

Clone the repository:

```bash
git clone https://github.com/moyarzoca/PaillacoDiff.git
```

and load `PaillacoDiff.wl` directly from Mathematica:

```mathematica
Get["/path/to/PaillacoDiff/PaillacoDiff.wl"];
```

Alternatively, copy `PaillacoDiff.wl` into the same directory as your notebook and load it with:

```mathematica
Get[FileNameJoin[{NotebookDirectory[], "PaillacoDiff.wl"}]];
```

## Two Usage Modes

PaillacoDiff supports two complementary ways of working with a geometry.

- **Global Mode** — geometric data such as `ds2` and `coord` are defined globally, and functions are called directly:

  ```mathematica
  Hstar[X]
  FormSquare[X]
  Paillaco["R(dn,dn)"]
  ```

  This mode is convenient for quick and interactive calculations.

- **Bundle Mode** — the geometric data are stored in a bundle, and the same operations receive the bundle explicitly:

  ```mathematica
  Hstar[bundle][X]
  FormSquare[bundle][X]
  Paillaco[bundle]["R(dn,dn)"]
  ```

  This mode is convenient when working with several geometries or when a calculation should be self-contained and reproducible.

Whenever an operation does not depend on geometric data, its syntax is identical in both modes, for example:

```mathematica
d[X]
FormDegree[X]
```

## Example: Reissner–Nordström

Here we combine the differential-form utilities and tensor machinery in a single example to verify the Einstein–Maxwell equations for the Reissner–Nordström solution. We first define the geometry and the electromagnetic field:

```mathematica
bundle = <|
    "ds2" -> -f[r]*d[t]^2 + d[r]^2/f[r] 
             + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2),
    "coord" -> {t, r, theta, phi}
|>;

d[Q] = d[M] = 0;

A = 2*Q/r*d[t];
F = d[A];

PaiDef[bundle]["F(dn,dn)", F];

PaiDef[bundle][
    "E{a b} := R{a b} - 1/2*g{a b}*Ricciscalar 
     + 1/2*(F{a c}*F{b ^c} - 1/4*g{a b}*F{c d}*F{^c ^d})"
];

PaiCalc[bundle]["E(dn,dn)"];
```
Then, we verify the equations in the solution

```mathematica
(* Maxwell equation *)
Simplify[d[Hstar[bundle][F]]]

(* Einstein equation *)
Simplify[
    PaiComponents[bundle]["E(dn,dn)"] /. 
    f -> Function[{r}, 1 - 2*M/r + Q^2/r^2]
]
```

## Tests

Run all tests using `wolframscript` from the repository root with

```bash
wolframscript -file tests/run_tests.wls
```

Use

```bash
wolframscript -file tests/run_tests.wls --verbose
```

to show the full test output.

## Reference

- [Differential forms](https://moyarzoca.github.io/paillacodiff/#differential-forms)
- [Tensor manipulation](https://moyarzoca.github.io/paillacodiff/#tensor-manipulation)
