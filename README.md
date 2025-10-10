# PaillacoDiff
A Wolfram Mathematica package for differential geometry and exterior algebra.

## Definitions

Let $F$ be a $p-{\rm form}$

$$F = \frac{1}{p!} F_{\mu_1 \dots \mu_p} dx^{\mu_1} \wedge \dots \wedge dx^{\mu_p}$$

The package supports

- Wedge product. ([esc] ^ [esc])
- Exterior derivatives. (d)
- Hodge star on coordinates basis. (MyHStar)
  
$$\star F = \frac{\sqrt{-\det g}}{p!(D-p)!}F^{\mu_1 \dots \mu_p} \epsilon_{\mu_1 \dots \mu_p \nu_1 \dots \nu_{D-p}} dx^{\nu_1} \wedge \dots \wedge dx^{\nu_{D-p}}$$

  $\qquad$ The normalization for the symbol in the coordiantes basis is $\epsilon_{1 \dots D} = 1$.

- Christoffel symbols.
- Ricci tensor.
- Ricci scalar. 
- Square of differential forms  (FormSquare)

  
$$F_{\lambda_1 \dots \lambda_p} F^{\lambda_1 \dots \lambda_p}$$


- Square of differential forms with two symmetric indices (FormSquaredd)

$$F_{\mu \lambda_2 \dots \lambda_p} F_{\nu}{}^{\lambda_2 \dots \lambda_p}$$
  
- Contraction operator  
- Supports vielbein basis
- Abstract Gamma matrices to compute Killings spinors
- ...

## Load PaillacoDiff
As a .wl file, you can load it using Get function following the directory of the file.

As a first attempt you can put PaillacoDiff.wl in the same folder as you .nb file and load it with

```mathematica
Get[FileNameJoin[{NotebookDirectory[], "PaillacoDiff.wl"}]];
```

## Example of usage: Reissner–Nordström

```mathematica
coord = {t, r, theta, phi};
Dim = 4;
ds2 = -f[r]*d[t]^2 + d[r]^2/f[r] + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2);
f[r_] := 1 - 2*M/r + Q^2/r^2/2
d[Q] = d[M] = 0;
A = Q/r*d[t];
F = d[A];

gdd = DiffToMatrix[ds2];
gUU = Inverse[gdd];
sqrtdetg = Sqrt[-Det[gdd]];

ComputeRicciScalar[]

ME = Simplify[d[MyHStar[F]]];
EEdd = Rdd - 1/2*gdd*RicciScalar - (FormSquaredd[F] - 1/4*gdd*FormSquare[F]);

Print[Simplify[ME]];
Print[Simplify[EEdd]];
```










