testDir = DirectoryName[$InputFileName];
repoDir = DirectoryName[testDir];

Get[FileNameJoin[{repoDir, "PaillacoDiff.wl"}]];
Get[FileNameJoin[{testDir, "TestHelpers.wl"}]];


VerificationTest[
    Module[{bund, A, F},
        bund = MakeSparseTestBundle[];
        A = bund["A"];
        F = d[A];

        PaiDef[bund]["F2tensor(dn,dn)", FormSquaredd[bund][F]];
        PaiDef[bund]["F2", FormSquare[bund][F]];
        PaiDef[bund]["TfromFromSquare{a b} := F2tensor{a b} - 1/4*g{a b}*F2"];
        PaiDef[bund]["F(dn,dn)", F];

        PaiDef[bund]["TfromF{a b} := F{a c}*F{b ^c} - 1/4*g{a b}*F{c d}*F{^c ^d}"];
        PaiDef[bund]["MustBeZero{a b} := TfromF{a b} - TfromFromSquare{a b}"];
        PaiCompute[bund]["MustBeZero(dn,dn)"];
        Simplify[PaiComponents[bund]["MustBeZero(dn,dn)"]]
    ],
    ConstantArray[0, {4, 4}],
    TestID -> "energy-momentum tensor"
];

VerificationTest[
    Module[{bund, A, F},
        bund = MakeSparseTestBundle[];
        A = bund["A"];
        F = d[A];

        PaiDef[bund]["F(dn,dn)", F];
        PaiDef[bund]["T{a b} := F{a c}*F{b ^c} - 1/4*g{a b}*F{c d}*F{^c ^d}"];
        PaiDef[bund]["E{a b} := R{a b}-1/2*g{a b}*Ricciscalar -1/2*T{a b}"];
        PaiCompute[bund]["E(dn,dn)"];
        Simplify[PaiComponents[bund]["E(dn, dn)"]
            /.f->Function[{r}, 1 - 2*m/r + Q^2/r^2/4]
            /. p->Function[{r}, Q/r]
        ]
    ],
    ConstantArray[0, {4, 4}],
    TestID -> "Reissner-Nordstrom field equations"
];

VerificationTest[
    Module[{bund, A, F, xi1, xi2, xi3, xi4, xi},
        bund = MakeSparseTestBundle[];

        xi1 = {0,0, Cos[phi], -Cot[theta]*Sin[phi]};
        xi2 = {0,0, -Sin[phi], -Cot[theta]*Cos[phi]};
        xi3 = {0, 0, 0, 1};
        xi4 = {1, 0, 0, 0};
        xi = c1*xi1 + c2*xi2 + c3*xi3 + c4*xi4;

        PaiDef[bund]["xi(up)", xi];
        PaiDef[bund]["KE{a b} := xi{a ;b} + xi{b ;a}"];
        PaiCompute[bund]["KE(dn, dn)"];

        Simplify[PaiComponents[bund]["KE(dn, dn)"]]
    ],
    ConstantArray[0, {4, 4}],
    TestID -> "Times translation + rotation Killing equation"
];

VerificationTest[
    Module[{bundV, omegaDD, omegaUD, gUU},
        bundV = MakeSparseTestBundleVielbein5[];

        PaiCompute[bundV]["omega(dn,dn)"];
        PaiCompute[bundV]["omega(up,dn)"];
        PaiCompute[bundV]["g(up,up)"];

        omegaDD = PaiComponents[bundV]["omega(dn,dn)"];
        omegaUD = PaiComponents[bundV]["omega(up,dn)"];
        gUU = PaiComponents[bundV]["g(up,up)"];

        Simplify[Normal[omegaUD] - gUU . Normal[omegaDD]]

    ],
    ConstantArray[0, {5, 5}],
    TestID -> "Vielbein - mixed operations flat/curve indices"
];

VerificationTest[
    Module[{bundV, xivdn},
        bundV = MakeSparseTestBundleVielbein5[];
        PaiDef[bundV]["xi(up)", {1,0,0,0,0}];
        PaiCompute[bundV]["xi(vdn)"];
        xivdn = PaiComponents[bundV]["xi(vdn)"];
        PaiDef[bundV]["zi(vdn)", xivdn];
        Paillaco[bundV]["zi(up)"]
    ],
    {1,0,0,0,0},
    TestID -> "Vielbein - flat/coord index transformations"
];


VerificationTest[
    Module[{bundV2},
        bundV2 = MakeVielbeinBundleAdS5[];
        PaiCompute[bundV2]["Rform(vup, vup)"];
        Simplify[ PaiComponents[bundV2]["Rform(vup,vup)"] + Table[Z[a]\[Wedge]Z[b], {a, 5}, {b, 5}]]
    ],
    ConstantArray[0, {5, 5}],
    TestID -> "Vielbein - curvature 2-form check"
];

