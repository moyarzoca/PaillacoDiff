#!/usr/bin/env wolframscript

globalFile = $ScriptCommandLine[[2]];
testNumber = ToExpression[$ScriptCommandLine[[3]]];

testDir = DirectoryName[globalFile];
repoDir = DirectoryName[testDir];

Get[FileNameJoin[{repoDir, "PaillacoDiff.wl"}]];
Get[globalFile];

report = Block[
    {$Output = {}},
    TestReport[{
        ReleaseHold[globalTests[[testNumber]]]
    }]
];

results = report["Results"];

successQ[test_] := test["Outcome"] === "Success";

Scan[
    Function[test,
        If[
            successQ[test],

            Print["[PASS] ", test["TestID"]],

            Print["[FAIL] ", test["TestID"]];
            Print["       File:     ", FileNameTake[test["TestFileName"]]];
            Print["       Type:     ", test["FailureType"]];
            Print["       Expected: ", test["ExpectedOutput"]];
            Print["       Actual:   ", test["ActualOutput"]];

            If[
                test["ActualMessages"] =!= {},
                Print["       Messages: ", test["ActualMessages"]]
            ];
        ]
    ],
    results
];
Exit[
    If[TrueQ[report["ReportSucceeded"]], 0, 1]
];
