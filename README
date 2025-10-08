# PaillacoDiff
A Wolfram Mathematica package for differential geometry and exterior algebra.

## Definitions
- Wedge product  
- Exterior derivatives  
- Hodge star  
- Christoffel symbols  
- Ricci tensor  
- Ricci scalar  
- Square of differential forms  
- Square of differential forms with two symmetric indices  
- Contraction operator  
- Supports vielbein basis  

## Example of usage: Reissner–Nordström

```mathematica
Get["PaillacoDiff.wl"]

coord = {t, r, theta, phi};
ds2 = -f[r]*d[t]^2 + d[r]^2/f[r] + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2);
A = Q/r*d[t];
F = d[A];

gdd = DiffToMatrix[ds2];
sqrtdetg = Sqrt[-Det[gdd]];
ComputeRicciScalar[]

ME = d[MyHstar[F]] == 0
EEdd = Rdd - 1/2*gdd*RicciScalar - FormSquaredd[F] - 1/4*gdd*FormSquare[F]
```










