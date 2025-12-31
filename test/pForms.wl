(* ::Package:: *)

Get["/home/marcelo/src/PaillacoDiff/PaillacoDiff.wl"];
Print["[ Testing Reisner-Nordstrom D=4 ]
    MyHStar(s) -> DNAtoForms(s)
               -> DNAofHStarU(s) -> DNAofMatrix(s)
		
    Hstar -> RaiseAllSparse
          -> nonzeroComps
          -> SparseFromDNA
          -> getCompRule
          -> DNAofForm
	
    ComputeRicciScalar
    DiffToMatrix
    
    no simplification in computing
"]
Print["Global: 
	coord, gdd, gUU, sqrtdetg
	ds2, A, F
	d[Q]=0
"]
coord = {t,r,th,ph};
fr = 1-2*M/r + Q^2/r^2/4;
ds2 = -fr*d[t]^2 + d[r]^2/fr + r^2*(d[th]^2 + Sin[th]^2*d[ph]^2);
d[Q] = 0;
d[M] = 0;
A = Q/r*d[t];
F = d[A];
gdd = DiffToMatrix[ds2];
gUU = Inverse[gdd];
sqrtdetg = Simplify[Sqrt[-Det[gdd]]];


bool1 = ((Simplify[d[Simplify[MyHStar[F]]]] === 0)
	&&(Simplify[MyHStar[MyHStar[F]]+F]===0)
	&&((MyHStar[Y]-Y*sqrtdetg*Apply[Wedge, d[coord]])===0)
	);
Print[ "[ ", ToString[bool1], " ]", "   MyHStar - Maxwell equations"];


bool2 = ((Simplify[d[Simplify[Hstar[F]]]] === 0)
	&&((Simplify[MyHStar[Hstar[F]]+F])===0)
	&&((Hstar[Y]-Y*sqrtdetg*Apply[Wedge, d[coord]])===0)); 
Print[ "[ ", ToString[bool2], " ]", "   Hstar - Maxwell equations"];


bool3 = (
	((MyHStar[A]-Hstar[A])===0)
	&&(Simplify[MyHStar[F]-Hstar[F]]===0)
	);
Print["[ ", bool3, " ]","   They both agree"];



FF = FormSquare[F];
Print["[ OK ]   FormSquare 2-form"];
FFdd = FormSquaredd[F];
Print["[ OK ]   FormSquaredd 2-form"];

bool = Dimensions[FFdd] === {4,4};
Print["[ "<>ToString[bool]<>" ]   Array dimensions"];

AA = FormSquare[A];
Print["[ OK ]   FormSquare 1-form"];
AAdd = FormSquaredd[A];
Print["[ OK ]   FormSquaredd 1-form"];

bool = Dimensions[AAdd] === {4,4};
Print["[ "<>ToString[bool]<>" ]   Array dimensions"];

Block[{$Output = {}},ComputeRicciScalar[]];
bool = DeleteDuplicates[Flatten[Simplify[Rdd - 1/2*RicciScalar*gdd - 1/2*(FFdd -1/4*gdd*FF)]]] === {0};
Print["[ "<>ToString[bool]<>" ]   Einstein equation"];




