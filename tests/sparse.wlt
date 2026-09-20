testDir = DirectoryName[$InputFileName];
repoDir = DirectoryName[testDir];

Get[FileNameJoin[{repoDir, "PaillacoDiff.wl"}]];
Get[FileNameJoin[{testDir, "TestHelpers.wl"}]];

bund = MakeSparseTestBundle[];



dim = Length[bund["coord"]];

Xdd = SparseArray[
    {
        {1, 2} -> 2,
        {2, 1} -> -1,
        {2, 3} -> 3,
        {4, 4} -> 5
    },
    {dim, dim}
]//Normal;

Ydd = SparseArray[
    {
        {1, 1} -> 4,
        {2, 3} -> -2,
        {3, 2} -> 1,
        {4, 1} -> 3
    },
    {dim, dim}
]//Normal;

Hdd = SparseArray[
    {
        {1, 2} -> 1,
        {2, 1} -> -1,
        {3, 4} -> 2,
        {4, 3} -> -2
    },
    {dim, dim}
]//Normal;

PaiDef[bund]["X(dn,dn)", Xdd];
PaiDef[bund]["Y(dn,dn)", Ydd];
PaiDef[bund]["H(dn,dn)", Hdd];

definitions = {
    "A1{a b} := X{a b}",
    "A2{a b} := -X{a b}",
    "A3{a b} := X{a b} + Y{a b}",
    "A4{a b} := X{a b} - Y{a b}",
    "A5{a b} := 3*X{a b}",
    "A6{a b} := 2*X{a b} - 3*Y{a b}",
    "A7{a b} := 2*(X{a b} - Y{a b})",
    "A8{a b} := (X{a b} - Y{a b}) + (Y{a b} + X{a b})",

    "S1 := g{a b}*g{^a ^b}",
    "S2 := X{a ^a}",
    "S3 := X{a b}*X{^a ^b}",

    "A9{a b} := X{a c}*Y{b ^c}",
    "A10{a b} := X{b a}",
    "A11{a b} := X{a b} + Y{b a}",

    "A12{a b} := (X{a c} - Y{a c})*X{b ^c}",

    "S4 := (X{a b} - Y{a b})*(X{^a ^b} + Y{^a ^b})",

    "A13{a b} := g{a b}*X{c ^c}",

    "A14{a b} := H{a c}*H{b ^c} - 1/4*g{a b}*H{c d}*H{^c ^d}",

    "A15{a b} := 2*(X{a b} - (Y{a b} - X{a b}))",

    "A16{a b} := (X{a c} - 2*Y{a c})*(X{b ^c} + Y{b ^c})- 1/4*g{a b}*X{c d}*X{^c ^d}"
};

        PartialContract[Xdd - 2*Ydd, Xdd + Ydd] - 1/4*gddTest*Contract2[Xdd, Xdd]
Map[PaiDef[#]&, definitions];
PaiCompute[bund]["g(dn,dn)"];
PaiCompute[bund]["g(up,up)"];
gddTest = PaiComponents[bund]["g(dn,dn)"];

gUUTest = PaiComponents[bund]["g(up,up)"];

Raise2[M_] := gUUTest . M . gUUTest;

Mixed[M_] := M . gUUTest;

Contract2[A_, B_] :=
    Tr[A.Transpose[Raise2[B]]];

PartialContract[A_, B_] :=
    A . Transpose[Mixed[B]];

expected = <|
    "A1(dn,dn)"  -> Xdd,
    "A2(dn,dn)"  -> -Xdd,
    "A3(dn,dn)"  -> Xdd + Ydd,
    "A4(dn,dn)"  -> Xdd - Ydd,
    "A5(dn,dn)"  -> 3*Xdd,
    "A6(dn,dn)"  -> 2*Xdd - 3*Ydd,
    "A7(dn,dn)"  -> 2*(Xdd - Ydd),
    "A8(dn,dn)"  -> 2*Xdd,

    "S1" -> dim,
    "S2" -> Tr[Mixed[Xdd]],
    "S3" -> Contract2[Xdd, Xdd],

    "A9(dn,dn)" ->
        PartialContract[Xdd, Ydd],

    "A10(dn,dn)" ->
        Transpose[Xdd],

    "A11(dn,dn)" ->
        Xdd + Transpose[Ydd],

    "A12(dn,dn)" ->
        PartialContract[Xdd - Ydd, Xdd],

    "S4" ->
        Contract2[Xdd - Ydd, Xdd + Ydd],

    "A13(dn,dn)" ->
        gddTest*Tr[Mixed[Xdd]],

    "A14(dn,dn)" ->
        PartialContract[Hdd, Hdd]
        - 1/4*gddTest*Contract2[Hdd, Hdd],

    "A15(dn,dn)" ->
        4*Xdd - 2*Ydd,

    "A16(dn,dn)" ->
        PartialContract[Xdd - 2*Ydd, Xdd + Ydd] - 1/4*gddTest*Contract2[Xdd, Xdd]
|>;

tests = KeyValueMap[
    Function[{spec, result},
        VerificationTest[
            Module[{},

                PaiCompute[bund][spec];
                Expand[Factor[Normal[PaiComponents[bund][spec]]]]

            ],
            Expand[Factor[Normal[result]]],
            TestID -> spec
        ]
    ],
    expected
];

TestReport[tests];
