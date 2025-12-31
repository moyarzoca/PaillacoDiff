(* ::Package:: *)

Get["/home/marcelo/src/PaillacoDiff/PaillacoDiff.wl"];
Print["[ Testing ]
    MyHStar(s) -> DNAtoForms(s)
               -> DNAofHStarU(s) -> DNAofMatrix(s)
		
    Hstar -> RaiseAllSparse
          -> nonzeroComps
          -> SparseFromDNA
          -> getCompRule
          -> DNAofForm
    DiffToMatrix
    * no simplification in computing
"]
Print["Global: 
	coord, gdd, gUU, sqrtdetg
	ds2, A, F
	d[Q]=0
"]
coord = {t,r,th,ph};
fr = 1-2*M/r;
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
Print[ "[ ", ToString[bool1], " ]", "   MyHStar"];
bool2 = ((Simplify[d[Simplify[Hstar[F]]]] === 0)
	&&((Simplify[MyHStar[Hstar[F]]+F])===0)
	&&((Hstar[Y]-Y*sqrtdetg*Apply[Wedge, d[coord]])===0)); 
Print[ "[ ", ToString[bool2], " ]", "   Hstar"];
bool3 = (
	((MyHStar[A]-Hstar[A])===0)
	&&(Simplify[MyHStar[F]-Hstar[F]]===0)
	);
Print["[ ", bool3, " ]","   They both agree"];




