

# Differential forms

We normalize a $p$ -form as

$$
F = \frac{1}{p!} F_{\mu_1 \dots \mu_p} dx^{\mu_1} \wedge \dots \wedge dx^{\mu_{p}}
$$

---

### `FormDegree`

```mathematica
FormDegree[expr]
```

Returns the degree of a differential form. Scalars have degree `0`.

---

### `Wedge`

```mathematica
Wedge[x, y, ...]
```

Computes the exterior product of differential forms.

---

### `Hstar`

```mathematica
Hstar[F]
Hstar[bundle][F]
```

Computes the Hodge dual of `F`.

`Hstar[F]` uses the current global geometry, while `Hstar[bundle][F]` uses the geometry stored in `bundle`.

$$
\star F = \frac{\sqrt{-g}}{p!(D-p)!}F^{\mu_1 \dots \mu_p} \epsilon_{\mu_1 \dots \mu_p \nu_1 \dots \nu_{D-p}} dx^{\nu_1} \wedge \dots \wedge dx^{\nu_{D-p}}
$$

---

### `FormSquare`

```mathematica
FormSquare[F]
FormSquare[bundle][F]
```
Computes the square of `F`

$$
F_{\mu_1 \dots \mu_p} F^{\mu_1 \dots \mu_p}
$$

---

### `FormSquaredd`

```mathematica
FormSquaredd[F]
FormSquaredd[bundle][F]
```

Computes the contraction `F` with itself, leaving two free indices

$$
F_{\mu \rho_1 \dots \rho_{p-1}} F_{\nu}{}^{\rho_1 \dots \rho_{p-1}}
$$

---

# Tensor manipulation

PaillacoDiff uses two index notations. In symbolic tensor expressions, indices are written inside braces: lower indices are written as `a`, while upper indices are written as `^a`. Covariant derivatives are specified after `;`, for example `xi{a ;b}`.

When requesting tensor components, indices are written inside parentheses: `dn` and `up` denote lower and upper coordinate indices, while `vdn` and `vup` denote lower and upper vielbein indices.

---

### `PaiDef`

Tensors can be defined in two ways.

From a symbolic tensor expression:

```mathematica
PaiDef["T{a b} := xi{a ;b}"]
PaiDef[bundle]["T{a b} := xi{a ;b}"]
```

The first definition is global, while the second is local to `bundle`.

Alternatively, a tensor can be defined from an explicit differential form, bilinear expression, or array:

```mathematica
PaiDef[bundle]["F(dn,dn)", d[r] \[Wedge] d[t]]
PaiDef[bundle]["S(dn,dn)", d[t]^2 + d[r]^2]
PaiDef[bundle]["V(up)", {1, 0, 0, 0}]
```

Here `F` is defined from a 2-form, `S` from a bilinear expression in the differentials, and `V` directly from its components.

The golden rule is to always specify the position of the indices, since PaillacoDiff converts the input expression into an array representation with a definite index structure.

---

### `PaiCalc`

Computes the requested tensor components.

```mathematica
PaiCalc["T(dn,dn)"]
PaiCalc[bundle]["T(dn,dn)"]
```
The first uses the global geometry, while the second uses the geometry stored in bundle.

---

### `PaiComponents`

Returns the components of a previously computed tensor as an array.

```mathematica
PaiComponents["T(dn,dn)"]
PaiComponents[bundle]["T(dn,dn)"]
```
The first uses global mode, while the second uses bundle mode.


---


