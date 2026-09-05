
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

Clear[GamQ];

GamQ[expr_] := Not[FreeQ[expr, _Gam | idGam]];d[Gam[x_?IntegerQ]]:=0;

PaiRegisterNonCommutativeScalarQ[GamQ];

Wedge[Times[gam_,x_],y__]:=CenterDot[gam,Wedge[x,y]]/;GamQ[gam]

(*============== Gamma matrix abstract product ============== *)

Clear[CenterDot];
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

BuildCliffordMap[bundle_] := Module[
    {dxToeRules, basisToGamma},

    dxToeRules = bundle["dxToe"];

    basisToGamma = Normal[AssociationThread[bundle["basis"], Map[Gam, Range[Length[bundle["basis"]]]]]];

    Function[{X},
        Activate[X /. dxToeRules /. Wedge -> Inactive[CenterDot] /. basisToGamma]]
];

CliffordMap[X_, bundle_] := BuildCliffordMap[bundle][X];

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

End[]

EndPackage[]
