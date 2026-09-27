# Conventions

This document summarizes the geometric conventions used by PaillacoDiff.

Greek indices \(\mu,\nu,\rho,\ldots\) denote coordinate indices, while Latin indices \(a,b,c,\ldots\) denote vielbein indices.

## Metric and Christoffel symbols

The metric and its inverse satisfy

\[
g^{\mu\rho}g_{\rho\nu}=\delta^\mu{}_\nu.
\]

PaillacoDiff uses the Levi-Civita connection,

\[
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
\]

The connection is torsion-free,

\[
\Gamma^\rho{}_{\mu\nu}
=
\Gamma^\rho{}_{\nu\mu}.
\]

## Riemann tensor

The Riemann tensor convention is

\[
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
\]

The fully covariant Riemann tensor is

\[
R_{\rho\sigma\mu\nu}
=
g_{\rho\lambda}
R^\lambda{}_{\sigma\mu\nu}.
\]

With this convention,

\[
R_{\rho\sigma\mu\nu}
=
-R_{\sigma\rho\mu\nu}
=
-R_{\rho\sigma\nu\mu}
=
R_{\mu\nu\rho\sigma}.
\]

## Ricci tensor

The Ricci tensor is obtained by contracting the first and third indices of the Riemann tensor,

\[
R_{\mu\nu}
=
R^\rho{}_{\mu\rho\nu}
=
g^{\rho\sigma}R_{\rho\mu\sigma\nu}.
\]

## Ricci scalar

The Ricci scalar is

\[
R
=
g^{\mu\nu}R_{\mu\nu}.
\]

## Vielbein

The vielbein one-forms are defined by

\[
e^a
=
e^a{}_\mu\,dx^\mu.
\]

The spacetime metric is related to the flat metric \(\eta_{ab}\) by

\[
g_{\mu\nu}
=
\eta_{ab}\,
e^a{}_\mu e^b{}_\nu,
\]

or equivalently,

\[
ds^2
=
\eta_{ab}\,e^a e^b.
\]

The inverse vielbein satisfies

\[
e_a{}^\mu e^a{}_\nu
=
\delta^\mu{}_\nu,
\qquad
e_a{}^\mu e^b{}_\mu
=
\delta_a{}^b.
\]

## Spin connection

PaillacoDiff uses the torsion-free spin connection. It is defined by the first Cartan structure equation,

\[
de^a
+
\omega^a{}_b\wedge e^b
=
0.
\]

With both flat indices lowered,

\[
\omega_{ab}
=
\eta_{ac}\omega^c{}_b,
\]

and metric compatibility implies

\[
\omega_{ab}
=
-\omega_{ba}.
\]

Equivalently, the vielbein postulate is

\[
\partial_\mu e^a{}_\nu
+
\omega^a{}_{b\mu}e^b{}_\nu
-
\Gamma^\rho{}_{\mu\nu}e^a{}_\rho
=
0.
\]

## Curvature two-form

The curvature two-form is defined by the second Cartan structure equation,

\[
\mathcal{R}^a{}_b
=
d\omega^a{}_b
+
\omega^a{}_c\wedge\omega^c{}_b.
\]

With both flat indices lowered,

\[
\mathcal{R}_{ab}
=
d\omega_{ab}
+
\omega_{ac}\wedge\omega^c{}_b.
\]

Its relation to the Riemann tensor in the vielbein basis is

\[
\mathcal{R}_{ab}
=
\frac{1}{2}
R_{abcd}\,
e^c\wedge e^d.
\]

Thus, the Riemann, Ricci tensor, and Ricci scalar computed in the coordinate and vielbein formulations use the same curvature convention.
