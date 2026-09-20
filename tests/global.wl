
globalTests = {
    Hold[VerificationTest[
        Module[{},
            ds2=-f[r]*d[t]^2 + d[r]^2/f[r] + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2);
            coord = {t,r,theta,phi};

            PaiCompute["R(dn,dn)", Simplify];
            Simplify[Simplify[PaiComponents["R(dn,dn)"]] /. f->Function[{r}, 1-2*M/r]]
        ],
        ConstantArray[0, {4, 4}],
        TestID -> "<global> Schwarzschild"
    ]]
    ,
    Hold[VerificationTest[
        Module[{},
            ds2=-f[r]*d[t]^2 + d[r]^2/f[r] + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2);
            coord = {t,r,theta,phi};
            PaiDef["G{a b} := R{a b} -1/2*g{a b}*Ricciscalar"];
            PaiDef["div{a} := G{a b ;^b}"];
            PaiCompute["div(dn)"];
            PaiComponents["div(dn)"]//Simplify
        ],
        ConstantArray[0, {4}],
        TestID -> "<global> Divergence of Einstein tensor"
    ]]
    ,
    Hold[VerificationTest[
        Module[{},
            ds2=-f[r]*d[t]^2 + d[r]^2/f[r] + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2);
            coord = {t,r,theta,phi};
            PaiCompute["R(dn,dn)"];

            ds2 = -d[t]^2 + d[r]^2 + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2);
            coord = {t,r,theta,phi};

            Paillaco["R(dn, dn)", Simplify]
        ],
        ConstantArray[0, {4, 4}],
        TestID -> "<global> consistent change from (ds2, coord) to (ds2', coord')"
    ]]
    ,
    Hold[VerificationTest[
        Module[{},
            ds2=-f[r]*d[t]^2 + d[r]^2/f[r] + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2);
            coord = {t,r,theta,phi};

            Simplify[Simplify[Paillaco["R(up,dn)"]] /. f->Function[{r}, 1-2*M/r]]
        ],
        ConstantArray[0, {4, 4}],
        TestID -> "<global> R(up,dn) on demand"
    ]]
    ,
    Hold[VerificationTest[
        Module[{A, F},
            ds2=-f[r]*d[t]^2 + d[r]^2/f[r] + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2);
            coord = {t,r,theta,phi};

            A = r^2*d[theta]/2;
            F = d[A];

            Simplify[Hstar[F] - r*Sin[theta]*f[r]*d[t] \[Wedge] d[phi], {r>0, 0<theta<Pi}]

        ],
        0,
        TestID -> "<global> Hstar computation on demand"
    ]]
    ,
    Hold[VerificationTest[
        Module[{A, F},
            ds2=-f[r]*d[t]^2 + d[r]^2/f[r] + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2);
            coord = {t,r,theta,phi};

            A = r^2*d[theta]/2;
            F = d[A];

            Simplify[FormSquare[F], {r>0, 0<theta<Pi}]

        ],
        2*f[r],
        TestID -> "<global> FormSquare computation on demand"
    ]]
    ,
    Hold[VerificationTest[
        Module[{A, F},
            ds2=-f[r]*d[t]^2 + d[r]^2/f[r] + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2);
            coord = {t,r,theta,phi};

            A = r^2*d[theta]/2;
            F = d[A];

            Simplify[FormSquaredd[F], {r>0, 0<theta<Pi}]

        ],
        DiagonalMatrix[{0, 1, r^2*f[r], 0}],
        TestID -> "<global> FormSquaredd computation on demand"
    ]]
    ,
    Hold[VerificationTest[
        Module[{A, F},
            ds2=-f[r]*d[t]^2 + d[r]^2/f[r] + r^2*(d[theta]^2 + Sin[theta]^2*d[phi]^2);
            coord = {t,r,theta,phi};

            A = r^2*d[theta]/2;
            F = d[A];

            Contraction[F]

        ],
        {0, r*d[theta], -r*d[r], 0},
        TestID -> "<global> Contraction computation on demand"
    ]]

}
