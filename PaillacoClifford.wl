
(* ::Package:: *)

BeginPackage["PaillacoClifford`", {"PaillacoDiff`"}]

(* ---------- Public functions ---------- *)

Gam::usage = "Gam[i] is an abstract Dirac gamma matrix."
GamQ::usage = "GamQ[expr] tests whether expr contains gamma matrices."
CenterDot::usage = "CenterDot[x, y] is the abstract multiplication of gamma matrices."
simpGamma::usage = "simpGamma[expr] simplifies products of gamma matrices."
GenerateGamma::usage = "GenerateGamma[dim, type] constructs explicit gamma matrix representations."
CliffordMap::usage = "CliffordMap[expr] maps a form in the vielbein basis to Clifford algebra."

(* ---------- Public globals ---------- *)

idGam::usage = "Identity element in the Clifford algebra."
id2::usage = "2x2 identity matrix in the Pauli sector."
GU::usage = "Explicit gamma matrix representation."

Begin["`Private`"]

(*============== Gamma matrices initialization ==============*)

ClearAll[GamWeight]
ClearAll[Gam]
GamWeight[X_]:=If[MatchQ[X,Gam[___]],1,0];
GamWeight[idGam]:=1;
GamWeight[x_CenterDot]:=Times@@Map[GamWeight,List@@x];
GamWeight[x_Times]:=Plus@@Map[GamWeight,List@@x];
GamWeight[x_Plus]:=If[MemberQ[Map[GamWeight,List@@x],1],1,0];
GamWeight[x_List]:=If[MemberQ[GamWeight/@x,1],1,0];
GamWeight[x_TensorProduct]:=Plus@@Map[GamWeight,List@@x]/Length[x];
GamQ[X_]:=If[GamWeight[X]>=1,True,False];
GamWeight[\[Sigma][i_?IntegerQ]]:=2;
GamWeight[id2]:=2;
d[Gam[x_?IntegerQ]]:=0;

PaiRegisterNonCommutativeScalarQ[GamQ];

Wedge[Times[gam_,x_],y__]:=CenterDot[gam,Wedge[x,y]]/;GamQ[gam]

(*============== Gamma matrix abstract product ============== *)

ClearAll[CenterDot];
Default[CenterDot] := 1;
SetAttributes[CenterDot,{Flat,OneIdentity}];
CenterDot[0,y__] := 0;
CenterDot[x__,0] := 0;
CenterDot[x__,y_Plus] := Plus@@(CenterDot[x,#]&/@(List@@y));
CenterDot[x_Plus,y__] := Plus@@(CenterDot[#,y]&/@(List@@x));
CenterDot[x__,Times[sca_,y_]] := Times[sca,CenterDot[x,y]]/;(NumericQ[sca]||!GamQ[sca]);
CenterDot[Times[sca_,x_],y__] := Times[sca,CenterDot[x,y]]/;(NumericQ[sca]||!GamQ[sca]);
CenterDot[y_,x_^n_] := x^n*y;
CenterDot[x_^n_,y_] := x^n*y;
CenterDot[Times[gam_,x_],y__] := CenterDot[gam,Wedge[x,y]]/;GamQ[gam];

ClearAll[auxCenterDot];
CenterDot[X_TensorProduct,Y__TensorProduct] := Fold[auxCenterDot,X,{Y}];
auxproduct[X_,Y_] := If[GamWeight[X]===1,CenterDot[X,Y],Dot[X,Y]];
auxCenterDot[X_TensorProduct,Y_TensorProduct] := TensorProduct[Sequence@@Table[auxproduct[(List@@X)[[i]] ,(List@@Y)[[i]]],{i,Length@X}]];

CliffordMap[X_] := X/.dxToe/.Wedge->Inactive[CenterDot]/.e[i_]:>Gam[i]//Activate;
(*============ Abstract tensor product between Gamma and Pauli Matrices ============*)

Unprotect[TensorProduct]
ClearAll[TensorProduct]
SetAttributes[TensorProduct,{Flat,OneIdentity}];
TensorProduct[Times[const_,x_],y__]:=const*TensorProduct[x,y]/;NumberQ[const]||!GamQ[const]
TensorProduct[x__,Times[const_,y_]]:=const*TensorProduct[x,y]/;NumberQ[const]||!GamQ[const]
TensorProduct[x__,y_Plus]:=Plus@@(TensorProduct[x,#]&/@(List@@y));
TensorProduct[x_Plus,y__]:=Plus@@(TensorProduct[#,y]&/@(List@@x));
Unprotect[Dot];
Dot[Times[const_,x_],y__]:=const*Dot[x,y]/;NumberQ[const]||!GamQ[const];
Dot[x__,Times[const_,y_]]:=const*Dot[x,y]/;NumberQ[const]||!GamQ[const];

(*======= Simplification rule for Gamma matrices  and sigma rules======*)

simpGamma[exp_]:=Module[{auxgamsimp,outputlist,collectedexpression},
		auxgamsimp[cd_]:=cd/.CenterDot[y__]:>outputlist[y];
		
	outputlist[y__] := Module[{ourlist,jumpsign,oursign,refinelist,RepPos,initialRepPos,newlist,ourfactorjumps,whowilldie,thegammafactor},
						ourlist=If[DeleteDuplicates[((#/.Gam[x_]:>x )&/@{y})]==={idGam},
									Return[idGam],DeleteCases[((#/.Gam[x_]:>x )&/@{y}),
									idGam]];
		
		RepPos[aList_]:=Module[{auxPos,gathered,PosEven,Rep},
						gathered = GatherBy[Range[Length@aList],aList[[#]]&];
						auxPos = Select[gathered,Length[#]>1&];
						PosEven = If[OddQ[Length[#]],Most[#],#]&/@auxPos;
						Rep = Table[aList[[First@index]],{index,auxPos}];Return[{Rep,PosEven}];
						];
						
		jumpsign := (-1)^(#1-#2-1)&;
		initialRepPos = RepPos[ourlist];
		newlist = ourlist;
		ourfactorjumps = 1;
		
		If[Last[initialRepPos]==={},
			Return[Signature[newlist]*(CenterDot[Sequence@@(Gam[#]&/@Sort[newlist])]/.{CenterDot[z_]:>z,
			CenterDot[]:>1})]];
		whowilldie = First@Partition[Last[RepPos[newlist]][[1]],2];
		
		Do[whowilldie=First@Partition[Last[RepPos[newlist]][[1]],2];
			ourfactorjumps = (ourfactorjumps*jumpsign[Sequence@@(whowilldie)]*
				Power[\[Eta]dd[[First[RepPos[newlist]][[1]],First[RepPos[newlist]][[1]]]],((Length[#]/2)&@whowilldie)]);
			newlist=Delete[newlist,Partition[whowilldie,1]];
		,{ii,Length[Partition[Flatten@Last@initialRepPos,2]]}];
		
		thegammafactor=(CenterDot[Sequence@@(Gam[#]&/@Sort[newlist])]/.{CenterDot[z_]:>z,CenterDot[]:>1});
		
		Return[ourfactorjumps*Signature[newlist]*If[NumberQ[thegammafactor],thegammafactor*idGam,thegammafactor]]];
	
	collectedexpression=Collect[Expand[exp],{idGam,CenterDot[Y__],Gam[JJ_]}];
	Return[If[Head[collectedexpression]===Plus,auxgamsimp[#]&/@collectedexpression,auxgamsimp[#]&@collectedexpression]];
]

d[\[Sigma][i_?IntegerQ]]:=0;d[id2]=0;

\[Sigma]rules = 
	{
		\[Sigma][j_?IntegerQ] . \[Sigma][j_?IntegerQ]:>id2,
		Dot[id2,X_]:>X,Dot[X_,id2]:>X,
		id2 . id2:>id2,
		\[Sigma][1] . \[Sigma][2]:>I \[Sigma][3],
		\[Sigma][2] . \[Sigma][3]:>I \[Sigma][1],
		\[Sigma][1] . \[Sigma][3]:>-I \[Sigma][2],
		\[Sigma][2] . \[Sigma][1]:>-I \[Sigma][3],
		\[Sigma][3] . \[Sigma][2]:>-I \[Sigma][1],
		\[Sigma][3] . \[Sigma][1]:>I \[Sigma][2]
		};

(*=== Safe Production between Gam Sigma matrices ===*)

RuleFromTensorProd[X_] := 
Module[{collected,listTerms},
	collected = Collect[X,_TensorProduct];
	listTerms = 
		If[
		Head[collected]===Plus,
			Apply[List,collected],
				{collected}
		];
	Map[getCoefTensorProd,listTerms]
];

Clear[getCoefTensorProd]
getCoefTensorProd[X_] := 
	Module[{basis,coef},
		basis = Cases[X, _TensorProduct,{0,\[Infinity]}];
		If[Length[basis]===1,
			basis = basis[[1]];
			coef  = Coefficient[X,basis]
		];
		Return[basis->coef]
	];
	
prodGamSig[0,X2_]:=0;
prodGamSig[X2_,0]:=0;
prodGamSig[X1_,X2_]:=
	Module[{rules1, rules2,alltuples,thesum},
		rules1 = Association[RuleFromTensorProd[X1]];
		rules2 = Association[RuleFromTensorProd[X2]];
		alltuples = Tuples[{Keys[rules1],Keys[rules2]}];
			
		thesum = Sum[
			rules1[tup[[1]]]*rules2[tup[[2]]]*(simpGamma[CenterDot[tup[[1]],tup[[2]]]]//.\[Sigma]rules)
		,
		{tup,alltuples}];
		
		Return[thesum];
	];


GenerateGamma[GU_,type_:"Symbolic"] := 
Module[{iterationsigma,pair,auxlastelement,beforegam,gu,readytoKronProd,sigmarule},
	Clear[GamTosigma];
	GamTosigma:={Gam[i_]:>GU[[i]],idGam:>TensorProduct@@Table[id2,{j,IntegerPart[Dim/2]}]};
	Print[Style["GamToSigma   ",Bold],"rule to map abstract \[CapitalGamma] matrices Gam[i] to explicit \[CapitalGamma] matrices",
		" in terms of tensor product"," of abstract Pauli matrices \[Sigma][i]"];
	Clear[\[Sigma]];
	d[\[Sigma][i_?IntegerQ]]:=0;d[id2]=0;
	sigmarule={\[Sigma][j_?IntegerQ] . \[Sigma][j_?IntegerQ]:>id2,
		Dot[id2,X_]:>X,Dot[X_,id2]:>X,
		id2 . id2:>id2,\[Sigma][1] . \[Sigma][2]:>I \[Sigma][3],\[Sigma][2] . \[Sigma][3]:>I \[Sigma][1],\[Sigma][1] . \[Sigma][3]:>-I \[Sigma][2],
		\[Sigma][2] . \[Sigma][1]:>-I \[Sigma][3],\[Sigma][3] . \[Sigma][2]:>-I \[Sigma][1],\[Sigma][3] . \[Sigma][1]:>I \[Sigma][2]};
	Print[Style["simpsigma  ",Bold],"Function that simplifies Dot[\[Sigma][i],\[Sigma][j]] Pauli matrices product"];
	Explicitsigma=(#/.TensorProduct->KroneckerProduct)/.{id2->IdentityMatrix[2],\[Sigma][i_?IntegerQ]:>PauliMatrix[i]}&;
	Print[Style["Explicitsigma  ",Bold],"Function that converts TensorProduct into KroneckerProduct and abstract Pauli matrices into explicit Pauli matrices"];
	simpsigma[x_]:=x//.sigmarule;
	
	Table[pair[i]=Table[\[Sigma][3],i-1],{i,IntegerPart[Dim/2]}];
	Do[beforegam[2i-1]=Append[pair[i],\[Sigma][1]];beforegam[2i]=Append[pair[i],\[Sigma][2]];,{i,IntegerPart[Dim/2]}];
	Do[beforegam[i]=Flatten[Append[beforegam[i],Table[id2,IntegerPart[Dim/2]-Length[beforegam[i]]]]],{i,2*IntegerPart[Dim/2]}];
	If[OddQ[Dim],auxlastelement=beforegam[Dim-1];auxlastelement[[-1]]=\[Sigma][3];beforegam[Dim]=auxlastelement];
	If[type==="Symbolic",
	Return[GU=Table[If[i===1,I,1]TensorProduct[Sequence@@(Array[beforegam,{Dim}][[i]])],{i,Dim}];]
	]
	;
	readytoKronProd=Array[beforegam,{Dim}]/.{\[Sigma][1]->PauliMatrix[1],\[Sigma][2]->PauliMatrix[2],\[Sigma][3]->PauliMatrix[3],id2->IdentityMatrix[2]};
	GU=Table[If[i===1,I,1]KroneckerProduct[Sequence@@(readytoKronProd[[i]])],{i,Dim}];
];

End[]

EndPackage[]
