

# Differential forms

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
Hstar[X]
Hstar[bundle][X]
```

Computes the Hodge dual of `X`.

`Hstar[X]` uses the current global geometry, while `Hstar[bundle][X]`
uses the geometry stored in `bundle`.

---

# Tensor manipulation

### `PaiDef`

```mathematica
PaiDef["T{indices} := expression"]
PaiDef[bundle]["T{indices} := expression"]
```

Defines a tensor using PaillacoDiff's explicit index notation.

Example:

```mathematica
PaiDef[bundle][
    "E{a b} := R{a b} - 1/2*g{a b}*Ricciscalar"
]
```

---

### `PaiCalc`

```mathematica
PaiCalc["T(dn,dn)"]
PaiCalc[bundle]["T(dn,dn)"]
```

Computes the requested tensor components.

---

### `PaiComponents`

```mathematica
PaiComponents["T(dn,dn)"]
PaiComponents[bundle]["T(dn,dn)"]
```

Returns previously computed tensor components.




