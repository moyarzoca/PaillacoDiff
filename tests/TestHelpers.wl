ClearAll[MakeSparseTestBundle];

MakeSparseTestBundle[] := <|
    "ds2" -> -f[r] d[t]^2 + d[r]^2/f[r] + r^2 (d[theta]^2 + Sin[theta]^2 d[phi]^2),
    "coord" -> {t, r, theta, phi},
    "A" -> p[r]*d[t]
|>;

MakeSparseTestBundleVielbein5[] := <|
    "eU"->{Q[r]*d[t], d[r]/P[r], M[r]*d[theta], M[r]*Sin[theta]*(d[phi1]), M[r]*Cos[theta]*d[phi2]},
    "coord" -> {t,r,theta,phi1, phi2},
    "basis"->Array[V, 5]
|>;

MakeVielbeinBundleAdS5[] := <|
    "eU" -> {Sqrt[r^2+1]*d[t], d[r]/Sqrt[r^2+1], r*d[theta], r*Sin[theta]*d[phi1], r*Cos[theta]*d[phi2]},
    "coord" -> {t,r,theta,phi1, phi2},
    "basis"->Array[Z, 5]
|>;
