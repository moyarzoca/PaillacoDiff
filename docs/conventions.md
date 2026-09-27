# Conventions

This document summarizes the geometric conventions used by PaillacoDiff.

Greek indices $\mu,\nu,\rho,\ldots$ denote coordinate indices, while Latin indices $a,b,c,\ldots$ denote vielbein indices.

In component notation, `dn` and `up` denote lower and upper coordinate indices, while `vdn` and `vup` denote lower and upper vielbein indices.

PaillacoDiff uses the following names for the standard geometric tensors:

- `g` -- metric tensor, e.g. `g{a b}` and `g(dn,dn)`.
- `R` -- Riemann tensor, e.g. `R{a b c d}` and `R(dn,dn,dn,dn)`.
- `R` -- Ricci tensor, e.g. `R{a b}` and `R(dn,dn)`.
- `Ricciscalar` -- Ricci scalar.
- `omega` -- spin connection 1-form, e.g. `omega(vdn,vdn)`
- `Rform` -- curvature 2-form, e.g. `Rform(vdn,vdn)`

### Metric and Christoffel symbols

The metric and its inverse satisfy

$$
g^{\mu\rho}g_{\rho\nu}=\delta^\mu{}_\nu.
$$

PaillacoDiff uses the Levi-Civita connection,

$$
\Gamma^\rho{}_{\mu\nu}
=
\frac{1}{2}g^{\rho\sigma}
\left(
\partial_\mu g_{\sigma\nu}
+
\partial_\nu g_{\sigma\mu}
-
\partial_\sigma g_{\mu\nu}
\right).
$$

### Riemann tensor

The Riemann tensor convention is

$$
R^\rho{}_{\sigma\mu\nu}
=
\partial_\mu \Gamma^\rho{}_{\sigma\nu}
-
\partial_\nu \Gamma^\rho{}_{\sigma\mu}
+
\Gamma^\rho{}_{\mu\lambda}
\Gamma^\lambda{}_{\sigma\nu}
-
\Gamma^\rho{}_{\nu\lambda}
\Gamma^\lambda{}_{\sigma\mu}.
$$

### Ricci tensor

The Ricci tensor is obtained by contracting the first and third indices of the Riemann tensor,

$$
R_{\mu\nu}
=
R^\rho{}_{\mu\rho\nu}
=
g^{\rho\sigma}R_{\rho\mu\sigma\nu}.
$$

### Ricci scalar

The Ricci scalar is

$$
R
=
g^{\mu\nu}R_{\mu\nu}.
$$

### Vielbein

The vielbein one-forms are defined by

$$
e^a
=
e^a{}_\mu\,dx^\mu.
$$

The spacetime metric is related to the flat metric $\eta_{ab}$ by

$$
g_{\mu\nu}
=
\eta_{ab}\,
e^a{}_\mu e^b{}_\nu,
$$

The inverse vielbein satisfies

$$
e_a{}^\mu e^a{}_\nu
=
\delta^\mu{}_\nu,
\qquad
e_a{}^\mu e^b{}_\mu
=
\delta_a{}^b.
$$

### Spin connection

PaillacoDiff uses the torsion-free spin connection. It is defined by the first Cartan structure equation,

$$
de^a
+
\omega^a{}_b\wedge e^b
=
0.
$$

With both flat indices lowered,

$$
\omega_{ab}
=
\eta_{ac}\omega^c{}_b,
$$

and metric compatibility implies

$$
\omega_{ab}
=
-\omega_{ba}.
$$

### Curvature two-form

The curvature two-form is defined by the second Cartan structure equation,

$$
\mathcal{R}^a{}_b
=
d\omega^a{}_b
+
\omega^a{}_c\wedge\omega^c{}_b.
$$

With both flat indices lowered,

$$
\mathcal{R}_{ab}
=
d\omega_{ab}
+
\omega_{ac}\wedge\omega^c{}_b.
$$

Its relation to the Riemann tensor in the vielbein basis is

$$
\mathcal{R}_{ab}
=
\frac{1}{2}
R_{abcd}\,
e^c\wedge e^d.
$$

Thus, the Riemann, Ricci tensor, and Ricci scalar computed in the coordinate and vielbein formulations use the same curvature convention.
