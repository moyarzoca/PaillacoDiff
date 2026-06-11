# PaillacoDiff



A Wolfram Mathematica package for differential geometry and exterior algebra.

## Installation

Put `PaillacoDiff.wl` in the same folder as your `.nb` file and load it with:

```mathematica
Get[FileNameJoin[{NotebookDirectory[], "PaillacoDiff.wl"}]];
```

## Two Usage Modes

PaillacoDiff supports two complementary workflows:

- **Casual Mode** — set global variables and call functions directly. Best for quick interactive work.
  → [docs/casual_mode.md](docs/casual_mode.md)

- **Bundle Mode** — Association-based with on-demand computation and caching. Best for structured, reproducible computations.
  → [docs/bundle_mode.md](docs/bundle_mode.md)

## Full Example: Reissner–Nordström

```mathematica
coord = {t, r, \[Theta], \[Phi]};
Dim = 4;
f[r_] := 1 - 2 M/r + Q^2/r^2/2;
ds2 = -f[r] d[t]^2 + d[r]^2/f[r] + r^2 (d[\[Theta]]^2 + Sin[\[Theta]]^2 d[\[Phi]]^2);

d[Q] = d[M] = 0;
A = Q/r d[t];
F = d[A];

gdd = DiffToMatrix[ds2];
gUU = Inverse[gdd];
sqrtdetg = Sqrt[-Det[gdd]];

ComputeRicciScalar[]

Print["Maxwell equation d(*F) = 0:"];
Simplify[d[Hstar[F]]]

Print["Einstein equation R_{mn} - 1/2 g_{mn} R = T_{mn}:"];
EEdd = Rdd - 1/2 gdd RicciScalar - (FormSquaredd[F] - 1/4 gdd FormSquare[F]);
Simplify[EEdd]
```

## Reference

See [docs/reference.md](docs/reference.md) for a complete function catalog.
