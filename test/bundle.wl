(* ::Package:: *)

Get["/home/marcelo/src/PaillacoDiff/PaillacoDiff.wl"];
Print["[ Testing ]
    PaiComputeBundleTensors -> PaiComputegddAgUU
                            -> PaiComputeChrisUdd
                            -> PaiComputeRdddd
                            -> PaiComputeRdd
                            -> PaiComputeRicciScalar
    no simplification in computing
"]
Print["Global: 
	bundle
	sqrtdetg
"]

fr = 1-2*M/r;
bundle = 
	<|
		"ds2" -> -fr*d[t]^2 + d[r]^2/fr + r^2*(d[th]^2 + Sin[th]^2*d[ph]^2),
		"coord" -> {t,r,th,ph}
	|>;
	
bundle = PaiComputeBundleTensors[bundle, "ChrisUdd"];

Print["[ OK ]   ChrisUdd Computed"];

bundle = PaiComputeBundleTensors[bundle, "Rdddd"];

Print["[ OK ]   Rdddd Computed"];

bundle = PaiComputeBundleTensors[bundle, "Rdd"];

Print["[ OK ]   Rdd Computed"];

bundle = PaiComputeBundleTensors[bundle, "RicciScalar"];

Print["[ OK ]   RicciScalar Computed"];

Print["\nSchwarzschild?\n"]

bool = Simplify[bundle["Tensors","RicciScalar"]] === 0;
Print["[ "<>ToString[bool]<>" ]   RicciScalar"];

bool = CleanZeros[Association[Simplify[Normal[bundle["Tensors", "Rdd"]]]]]===<||>;
Print["[ "<>ToString[bool]<>" ]   Rdd"];

hstar = BuildHodge[bundle];

Print["[ OK ] BuildHodge"];

bool = Simplify[d[hstar[d[r]\[Wedge]d[t]/r^2]]] === 0;
Print["[ "<>ToString[bool]<>" ]   Maxwell equations"];

sqrtdetg = Sqrt[-Simplify[Det[DiffToMatrix[bundle["ds2"], bundle["coord"]]]]];
bool = Simplify[hstar[Y] - Y*sqrtdetg*Apply[Wedge,d[bundle["coord"]]]] === 0;
Print["[ "<>ToString[bool]<>" ]   Hodge dual of scalar"]



