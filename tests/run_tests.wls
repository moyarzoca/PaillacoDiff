#!/usr/bin/env wolframscript

testDir = DirectoryName[$InputFileName];
testFiles = Sort[FileNames["*.wlt", testDir]];
globalFile = FileNameJoin[{testDir, "global.wl"}];


If[
    testFiles === {},
    Print["No .wlt test files found in: ", testDir];
    Exit[2];
];

verboseQ = MemberQ[$ScriptCommandLine, "--verbose"];

report = If[verboseQ,
            TestReport[testFiles],
                Block[{$Output = {}},TestReport[testFiles]]
];

results = report["Results"];

successQ[test_] := test["Outcome"] === "Success";

cleanOutput[expr_] := expr /. HoldForm[x_] :> x;

passed = Select[results, successQ];
failed = Select[results, Function[test, Not[successQ[test]]]];

Print[""];
Print["================ PaillacoDiff Tests ================"];

Scan[
    Function[test,
        If[
            successQ[test],

            Print["[PASS] ", test["TestID"]],

            Print["[FAIL] ", test["TestID"]];
            Print["       File:     ", FileNameTake[test["TestFileName"]]];
            Print["       Type:     ", test["FailureType"]];
            Print["       Expected: ", cleanOutput[test["ExpectedOutput"]]];
            Print["       Actual:   ", cleanOutput[test["ActualOutput"]]];

            If[
                test["ActualMessages"] =!= {},
                Print["       Messages: ", test["ActualMessages"]]
            ];
        ]
    ],
    results
];


Print["===================================================="];
Print[
    Length[passed], " passed, ",
    Length[failed], " failed",
    "   (",
    Round[
    	QuantityMagnitude[UnitConvert[report["AbsoluteTimeUsed"], "Seconds"]]
    , 0.0001],
    " s)"
];

Print["===================================================="];

Get[globalFile];
nGlobalTests = Length[globalTests];

clockGlobal = AbsoluteTiming[
    globalResults = Map[
        Function[i,
            RunProcess[{
                "wolframscript",
                "-file",
                FileNameJoin[{testDir, "run_fresh.wls"}],
                globalFile,
                ToString[i]
            }]
        ],
        Range[nGlobalTests]
    ];
]

Scan[
    Function[result,
        Print[StringTrim[result["StandardOutput"]]]
    ],
    globalResults
];

globalSucceeded = AllTrue[
    globalResults,
    #["ExitCode"] === 0 &
];


globalPassed = Count[globalResults, x_ /; x["ExitCode"] === 0];
globalFailed = Length[globalResults] - globalPassed;

Print["===================================================="];

Print[
    globalPassed, " global passed, ",
    globalFailed, " global failed" ,
    "(", clockGlobal[[1]], " s)"
];

Print["===================================================="];

If[
    TrueQ[report["ReportSucceeded"]] && globalSucceeded,
    Exit[0],
    Exit[1]
];

