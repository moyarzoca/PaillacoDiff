(* ::Package:: *)

Get["/home/marcelo/src/PaillacoDiff/PaillacoDiff.wl"];
Print["** Testing 
	ComputeRicciScalar -> ComputeRdd -> ComputeChrisUdd
		no simplification
	DiffToMatrix
"]
Print["Global: coord, ds2, gdd"]
coord = {t,r,th,ph};
fr = 1-2*M/r;
ds2 = -fr*d[t]^2 + d[r]^2/fr + r^2*(d[th]^2 + Sin[th]^2*d[ph]^2);
gdd = DiffToMatrix[ds2];
Block[{$Output = {}}, ComputeRicciScalar[]];
bool1 = Simplify[Rdd]===Rdd*0;
bool2 = Simplify[RicciScalar]===0;
Print[And[bool1,bool2]];

Print["\n-------------------------------------------\n"]



ClearAll["Global`*"];
Get["/home/marcelo/src/PaillacoDiff/PaillacoDiff.wl"];
Print["** Testing 
	ComputeSpinConnection -> SetVielbein -> FormsToMatrix
		no simplification
"];

Print["Global: 
	coord, Dim,
	eU, eta,
	d[M]=0"];
	
coord = {t,r,th,ph};
Dim = Length[coord];
d[M]=0;
eU = {-F[r]*d[t],d[r]/G[r],r*d[th],r*Sinh[th]*d[ph]};
eta = IdentityMatrix[Length[coord]];
eta[[1,1]] = -1;
Block[{$Output = {}},
	ComputeSpinConnection[eU, eta]
];

Print["
** Checking anti-symmetry of \[Omega]dd\n 
", Simplify[\[Omega]dd == - Transpose[\[Omega]dd]]]



