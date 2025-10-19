(* ::Package:: *)

"
Author : Marcelo Oyarzo
Acknowledgement: I learned some of properties of FormDegree, Wedge and d[] from the RGTC source code.
So I want to thank the author of RGTC package for keeping it open source. 
Also I thank Ruggero Noris and Stefano Maurelli the code and 
for poiting my out issues and helped me to improve the code.
"

ClearAll[FormDegree]
FormDegree[d[x_]]:=1+FormDegree[x];
FormDegree[x_]:=0;
FormDegree[e[n_Integer]]:=1;
FormDegree[x_Wedge]:=Plus@@Map[FormDegree,List@@x];
FormDegree[x_Times]:=Plus@@Map[FormDegree,List@@x];
FormDegree[x_Plus]:=FormDegree[First@x];
IsThisGood[x_Plus]:=If[Length@DeleteDuplicates[FormDegree/@List@@x]==1,Print["Yes"],Print["No"]];
FormDegree[x_List]:=FormDegree/@x;
IsThisGood[x_List]:=If[Length@DeleteDuplicates[FormDegree/@Flatten@x]==1,Print["Yes"],Print["No"]];

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

(*============== Definition Wedge ==============*)

ClearAll[Wedge]
Default[Wedge]:=1;
Wedge/:Wedge[]:=1;
Wedge/:Wedge[arg_/;!(Head[arg]===Pattern)]:=arg;
SetAttributes[Wedge,{Flat,OneIdentity}];
Wedge[0,y__]:=0;
Wedge[x__,0]:=0;
Wedge[x__,y_Plus]:=Plus@@(Wedge[x,#]&/@(List@@y));
Wedge[x_Plus,y__]:=Plus@@(Wedge[#,y]&/@(List@@x));
Wedge[x__,y_List]:=(Wedge[x,#]&/@y);
Wedge[x_List,y__]:=(Wedge[#,y]&/@x);
Wedge[y_,x_^n_.]:=x^n*Wedge[y]/;FormDegree[x]===0;
Wedge[x_^n_.,y_]:=x^n*Wedge[y]/;FormDegree[x]===0;
Wedge[x__,Times[sca_,y_]]:=Times[sca,Wedge[x,y]]/;NumericQ[sca]||(FormDegree[sca]===0&&(!GamQ[sca]||!GamQ[{x}]));
Wedge[Times[sca_,x_],y__]:=Times[sca,Wedge[x,y]]/;NumericQ[sca]||(FormDegree[sca]===0&&!GamQ[sca]);
Wedge[x_,y___,x_]:=0/;OddQ[FormDegree[x]]&&!GamQ[x]&&!GamQ[{y}];
Wedge[y__]:=Signature[{y}]*Wedge@@Sort[{y}]/;Sort[{y}]=!={y}&&Union[FormDegree[{y}]]==={1}&&!GamQ[{y}];
Wedge[Times[gam_,x_],y__]:=CenterDot[gam,Wedge[x,y]]/;GamQ[gam]

Wedge[A_List,B_List]:=Module[{matrixdimension},
If[Dimensions[A]!=Dimensions[B],Return[Print["Incompatible matrix dimensions for multiplication"]]];
matrixdimension=Dimensions[A][[1]];
Table[Sum[Wedge[A[[ii,kk]],B[[kk,jj]]],{kk,matrixdimension}],{ii,matrixdimension},{jj,matrixdimension}]]

(*===== Exterior derivative =====*)

SetAttributes[d,{Listable}];
d[x_Wedge/;Length@x===2] := Wedge[d[First[x]],Last[x]]+(-1)^FormDegree[First[x]]*Wedge[First[x],d[Last[x]]];
d[x_Times|x_Wedge] := Wedge[d[First[x]],Rest[x]]+(-1)^FormDegree[First[x]]*Wedge[First[x],d[Rest[x]]];
d[x_Plus]:=d/@x;
d[x_?NumericQ|x_d] := 0;
nodHeads={Pattern,Blank,Condition,RuleDelayed,SeriesData};
d[h_[y__]/;FreeQ[nodHeads,h]] := 
	Sum[
		Derivative[Sequence@@ReplacePart[Table[0,{Length[{y}]}],i->1]][h][y]*d[{y}[[i]]]
	,{i,Length[{y}]}]/;(FormDegree[h[y]]===0 && FreeQ[{Integer,Blank,Pattern,Condition},Head[First[{y}]]]);

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

\[Sigma]rules={\[Sigma][j_?IntegerQ] . \[Sigma][j_?IntegerQ]:>id2,
	Dot[id2,X_]:>X,Dot[X_,id2]:>X,
	id2 . id2:>id2(*,
	\[Sigma][i_?IntegerQ].\[Sigma][j_?IntegerQ]:>\[Sigma][j].\[Sigma][i]/;i>j *),\[Sigma][1] . \[Sigma][2]:>I \[Sigma][3],\[Sigma][2] . \[Sigma][3]:>I \[Sigma][1],\[Sigma][1] . \[Sigma][3]:>-I \[Sigma][2],
	\[Sigma][2] . \[Sigma][1]:>-I \[Sigma][3],\[Sigma][3] . \[Sigma][2]:>-I \[Sigma][1],\[Sigma][3] . \[Sigma][1]:>I \[Sigma][2]}
	

(*====== Towards contraction operator =====*)

Clear[Extractor];
Extractor[F3_List,y_]:=Extractor[#,y]&/@F3;
Extractor[x__,oneform_]:=
Module[{listx,listxref,killz},
	killz[y__,KK_]:=
	Module[{pos},pos=Flatten[Position[{y},KK]];If[pos==={},Return[0]];
		Return[(Wedge@@DeleteCases[{y},KK])*(-1)^(Length[{y}]-pos)]
		];
	listx = If[Head[x]===Plus, List@@x,{x}];
	listxref = Select[listx,(!FreeQ[#,oneform])&];
	Return[Plus@@Flatten[#/.Wedge[y__]:>killz[y,oneform]/.oneform->1&/@listxref]];
	];
Extractor::usage="Given a polyform F = F1+ F2\[Wedge]A,  with A a 1-form,  Extractor[F,A] returns F2.";

Clear[Extractorleft];
Clear[coordcontraction];
Extractorleft[F3_List,y_]:=Extractor[#,y]&/@F3;
Extractorleft[x__,oneform_]:=
Module[{listx,listxref,killzleft},
	killzleft[y__,KK_]:=
	Module[{pos},pos=Flatten[Position[{y},KK]];If[pos==={},Return[0]];
		Return[(Wedge@@DeleteCases[{y},KK])*(-1)^(pos-1)]
		];
	listx=If[Head[x]===Plus, List@@x,{x}];
	listxref=Select[listx,(!FreeQ[#,oneform])&];
	Return[Plus@@Flatten[#/.Wedge[y__]:>killzleft[y,oneform]/.oneform->1&/@listxref]];
];
coordcontraction[X_List,coord_:coord]:=Map[coordcontraction[#,coord]&,X];
coordcontraction[X_,coord_:coord]:=Extractorleft[X,#]&/@d[coord];
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

(*
			---- DNAofForm ----
*)

PolyFormQ[expr_] := Module[{terms, degs,exprExpand,degsDiff},
	exprExpand = Expand[expr];
	terms = 
		If[
		Head[exprExpand] === Plus,
			Apply[List,exprExpand],
				{exprExpand}
		];
	degs = Map[FormDegree, terms];
	degsDiff = DeleteDuplicates[degs];
	If[
	Length[degsDiff]>1,
		Return[True],
			Return[False]
	];
];

Clear[coeffBaseElement];

coeffBaseElement[pform_, coord_]:=
	Module[{baseElement, coordBasisQ,deg, coeff, 
		baseAlong, Dim, mapcoord,baseNumb, sign},
		
		deg = FormDegree[pform];
		Dim = Length[coord];
		mapcoord = AssociationThread[coord -> Range[Dim]];
		
		coordBasisQ = FreeQ[pform, e[_]];
		
		Which[
		(deg===1)&&coordBasisQ,
			baseElement = Cases[pform, _d, {0, Infinity}];
			baseNumb = baseElement/.d[y_]:>{mapcoord[y]}
			,
		(deg===1)&&Not[coordBasisQ],
			baseElement = Cases[pform, _e, {0, Infinity}];
			baseNumb = baseElement/.e[y_]:>{y}
			,
		(deg > 1)&&coordBasisQ,
			baseElement = Cases[pform, _Wedge, {0, Infinity}];
			baseNumb = baseElement/.Wedge[YY__]:>{YY}/.d[y_]:>{mapcoord[y]};
			,
		(deg > 1)&&Not[coordBasisQ],
			baseElement = Cases[pform, _Wedge, {0, Infinity}];
			baseNumb = baseElement/.Wedge[YY__]:>{YY}/.e[y_]:>{y};
		];
		
		baseNumb = Flatten[baseNumb];
		sign = Signature[baseNumb];
		baseNumb = Sort[baseNumb];
		
		If[
		Length[baseElement]===1,
			baseAlong = baseElement[[1]];
			coeff = Coefficient[pform, baseAlong],
				Message[DNAofForm::noBaseFound, baseElement];
				Return[$Failed]
		];
	Return[{sign*coeff, baseNumb}];
	];

DNAofForm::nocord = 
	"No coordinates provided to computed DNAofForm";
DNAofForm::noBaseFound = 
	"No base element found";
	
DNAofForm /: DNAofForm[X_,coordIN_:Automatic] := 
	Module[{TermsXArray,listtermscoeffs,mapcoord,Dimint,coordint,mappiator,Collected, track},
		coordint = 
			Which[
				coordIN =!= Automatic,
					coordIN,
				ValueQ[coord], 
					coord,
				True,
					Message[DNAofForm::nocord];
					Return[$Failed]
			];
		Dimint = Length[coordint];
		
		Collected = 
			Which[
			FormDegree[X]===1,
				Collect[Expand@X, {d[1forms_], e[num_]}],
			(FormDegree[X]>1)||PolyFormQ[X],
				Collect[Expand@X, Wedge[listforms___]]
			];
		
		TermsXArray = 
			If[
			Head[Collected] === Plus,
				Apply[List,Collected],
					{Collected}
				];
		
		Return[Map[coeffBaseElement[#, coordint]&,TermsXArray]];
	];


Errorcoord[] := Print[Style["Coordinates not defined. ",Red,14], "The code requiers a global variabled called ",
					Style[coord,Bold]," which is an array of the coordiantes."]

FormsToMatrix[X_,formd_:"formdegree"] := 
Module[
{mapcoord, formdegree, listtermsX, listtermsX1, listtermsforms, listtermsnumbers,
listtermsnumbers1, listtermscoeffs, listtermscoeffs1, bigmatrix, indicess, numterms, listforms, thevector},
	If[
	Length@coord>0,
		Do[mapcoord[coord[[ii]]] = ii, {ii, Dim}],
			Return[Errorcoord[]]
	];
	If[
	X == 0,
		If[
		formd!="formdegree",
			Return[ConstantArray[0,ConstantArray[Dim,formd]]],
				Return[0]
		]
	];
	formdegree = FormDegree[X];
	(*The next line is only for 1-forms*)
	If[
	formdegree == 1,
		listtermsX1 = If[Head@Collect[Expand@X, d[YY_]] === Plus, List @@ Collect[Expand@X, d[YY_]], {X}];
		listtermscoeffs1 = Map[# /. d[terms_] :> 1 &, listtermsX1];
		listtermsnumbers1 = Map[# /. d[terms_]*coeff2_. :> {mapcoord[terms]} &, listtermsX1];
		Return[D[X,d[#]]&/@coord];
         ];
      (*for p-forms with p>1 the following code is the correct*)
      listtermsX = If[Head[Collect[Expand@X, Wedge[listforms___]]] === Plus, List @@ Collect[Expand@X, Wedge[listforms__]], {Collect[Expand@X, Wedge[listforms___]]}];(*from the p-form we construct a list with the non-vanishing terms*)
      listtermscoeffs = Map[# /. Wedge[YY__] :> 1 &, listtermsX];
      (*the coefficients*)
      numterms = Length@listtermscoeffs;
      listtermsforms = Map[# /. Wedge[terms__]*coeff2_. :> {terms} &, listtermsX];(*list of forms*)
      listtermsnumbers = Map[Map[mapcoord, #] &, (listtermsforms /. Wedge[XX__] :> XX /. d[YZ_] :> YZ)];(*create a list of numbers related to the non trivial entries. ---The next line is in case of emergence*)
(*Return[{listtermsX,listtermscoeffs,listtermsforms,listtermsnumbers}];*)
      bigmatrix = ConstantArray[0, Table[Dim, formdegree]];(*creats the list*)
      Do[Do[bigmatrix[[Sequence @@ indicess]] = listtermscoeffs[[II]] Signature[indicess] Signature[listtermsnumbers[[II]]] , {indicess, Select[Tuples[listtermsnumbers[[II]], formdegree], Signature[#] != 0 &]}], {II, numterms}];
      Return[bigmatrix];
];

DNAofMatrix[X_] :=
Module[{blockofX,LengthBlocks,nonzeroX,CoeffcientsofX},
	nonzeroX = SparseArray[X]["NonzeroPositions"];
	blockofX = DeleteDuplicates[Sort/@nonzeroX];
	CoeffcientsofX = Table[X[[Sequence@@blockofX[[IIinx]]]],{IIinx,Length@blockofX}];
	Table[{CoeffcientsofX[[IIinx]], blockofX[[IIinx]]}, {IIinx, Length@blockofX}]
];

DNAofHStarU[X_,type_:"ListOfInput",sqrtdetgcoord_:{sqrtdetg,coord}] := 
Module[{DNAofX,DNAofHStarX,coordint,Dimint,sqrtdetgint},
	If[
	sqrtdetgcoord==={gdd,coord},
		coordint=coord;
		sqrtdetgint=sqrtdetg,
			coordint=sqrtdetgcoord[[2]];
			sqrtdetgint=sqrtdetgcoord[[1]]
	];
	Dimint=Length@coordint;
	If[
	ValueQ[sqrtdetgcoord[[1]]]===False,
		Return[Print["You must declare a global variable ",Style["sqrtdetg",Bold]," which is Sqrt[-det[Subscript[g, \[Mu]\[Nu]]]]."]]

		];
	If[
	type==="DNA",
		DNAofX=X,DNAofX=DNAofMatrix[X]];
		DNAofHStarX = 
		Table[{sqrtdetgint*DNAofX[[JJinx,1]]*Signature[{Sequence@@DNAofX[[JJinx,2]],Sequence@@Complement[Range[Dimint],DNAofX[[JJinx,2]]]}]
		,Complement[Range[Dimint],DNAofX[[JJinx,2]]]}
		,{JJinx,Length@DNAofX}];
	Return[DNAofHStarX];
]

"The function eats a matrix with Upper indices and computes the DNA of HStar matrix";
"The idea of the Complement[...] part is that given a non-trivial element of X, 
which is the p-form we want to HStar, it will fill the first entries of the Levi-Civita symbol and 
we fill the rest using Complement[Range[Dim],] function, where the first entry is the full list and 
the second element is a subset of the list, it returns the complement.";

DNAtoMatrix[DNA_]:=
Module[{listtermsnumbers,bigmatrix,formdegree},
	formdegree = Length[DNA[[1,2]]];
	listtermsnumbers = Table[DNA[[IIinx,2]], {IIinx, Length@DNA}];
	bigmatrix = ConstantArray[0,Table[Dim,formdegree]];
	Do[
		Do[
			bigmatrix[[Sequence@@indicess]]=DNA[[IIinx,1]] Signature[indicess] Signature[listtermsnumbers[[IIinx]]] 
		,{indicess,Select[Tuples[listtermsnumbers[[IIinx]],formdegree],Signature[#]!=0&]}]
	,{IIinx, Length@DNA}];
	Return[bigmatrix];
];

DNAtoForms[DNA_,coordIN_:coord] :=
Module[{II,numbertocoord,coordint},
	If[
	coordIN===coord,
		coordint=coord,
			coordint=coordIN
	];
	numbertocoord[ii_]:=coordint[[ii]];
	Return[Sum[DNA[[IIinx,1]]* Wedge[Sequence@@(d[numbertocoord[#]]&/@(DNA[[IIinx,2]]))],{IIinx,Length@DNA}]];
];

RaiseIndices[Xd_] :=
Module[{degreeform},
	degreeform=Length@Dimensions[Xd];
	Activate@TensorContract[Inactive[TensorProduct][Xd, Sequence@@Table[gUU,{IIinx,degreeform}]],Table[{iiinx,degreeform+2*iiinx-1},{iiinx,degreeform}]]
]

(*====== HStar in the coordinates basis ======*)

MyHStar[X_List,simp_:Identity,gddcoord_:{gdd,coord}]:=MyHStar/@X;

MyHStar[X_,simp_:Identity,gddcoordgUUsqrt_:{gdd,coord}]:=
Module[{formdegree,DNA,BADDNAUp,DNAUp,thetuples,auxtiempo,gUUnonZero,gUUint,Dimint,coordint,sqrtdetgint},
	
	If[X==0, Return[0]];
	
	{coordint, Dimint, gUUint, sqrtdetgint} = 
		Which[
			gddcoordgUUsqrt === {gdd, coord},
				{coord, Length[coord], gUU, sqrtdetg},
		
			Length[gddcoordgUUsqrt] == 4,
				{gddcoordgUUsqrt[[2]], Length[gddcoordgUUsqrt[[2]]], gddcoordgUUsqrt[[3]], gddcoordgUUsqrt[[4]]},
		
			True,
				{gddcoordgUUsqrt[[2]], Length[gddcoordgUUsqrt[[2]]], Inverse[gddcoordgUUsqrt[[1]]], simp[Sqrt[-Det[gddcoordgUUsqrt[[1]]]]]}
		];

	gUUnonZero =
		Table[
			#[[2]]&/@Select[SparseArray[gUUint]["NonzeroPositions"],#[[1]]==iiinx&]
		,{iiinx,Dimint}];

	formdegree = FormDegree[X]; 

	If[
	formdegree===0,
		Return[X*sqrtdetgint Wedge@@(d/@coordint)]
	];
	
	DNA = DNAofForm[X,coordint];
	auxtiempo = SessionTime[];
	
	thetuples = 
		Table[
			Select[
				Tuples[
					Table[
						gUUnonZero[[Last[DNA[[IIinx]]][[kkinx]]]]
					,{kkinx,formdegree}]
				]
			,Signature[#]!=0&]
		,{IIinx,Length@DNA}];
	
	BADDNAUp = 
	Flatten[
		Table[
			Table[
				{First[DNA[[IIinx]]]*Product[gUUint[[Last[DNA[[IIinx]]][[llinx]],indexx[[llinx]]]],{llinx,formdegree}],indexx}
			,{indexx,thetuples[[IIinx]]}]
		,{IIinx,Length@DNA}]
	,1];
	
	auxtiempo = SessionTime[];
	DNAUp = DNAofForm[DNAtoForms[BADDNAUp,coordint],coordint];
	Return[DNAtoForms[DNAofHStarU[DNAUp,"DNA",{sqrtdetgint,coordint}],coordint]];
];

HStarT[X_]:=(-1)^(FormDegree[X]*(Dim-FormDegree[X]))*MyHStar[X];

(*====== Squares of differential forms ======*)

Clear[FormSquare];
FormSquare[Xform_,gintUUinput_:gUU,simp_:Identity] :=
Module[{degreeform,Xd,XU,gintUU,listindices},
	If[
	gintUUinput===gUU,
		gintUU=gUU,
			gintUU=gintUUinput
	];
	
	degreeform = FormDegree[Xform];
	Xd = simp[FormsToMatrix[Xform]];
	listindices = Join[Table[{iiinx,2*degreeform+2*iiinx-1},{iiinx,degreeform}], Table[{iiinx+degreeform,2*degreeform+2*iiinx},{iiinx,degreeform}]];
	Return[Activate@TensorContract[Inactive[TensorProduct][Xd,Xd,Sequence@@Table[gintUU,{IIinx,degreeform}]], listindices]]
];

Clear[FormSquaredd];
FormSquaredd[0,simp_:Identity]:=0

FormSquaredd[Xform_,gintUUinput_:gUU,simp_:Identity] :=
Module[{degreeform,Xd,XdU,DNAofXd,gintUU,listindices},
	If[
	gintUUinput===gUU,
		gintUU=gUU,
			gintUU=gintUUinput
	];
	DNAofXd=simp[DNAofForm[Xform]];
	Xd=SparseArray[DNAtoMatrix[DNAofXd]];
	degreeform=FormDegree[Xform];
	listindices = Join[Table[{iiinx+1,2*degreeform+2*iiinx-1},{iiinx,degreeform-1}], Table[{iiinx+1+degreeform,2*degreeform+2*iiinx},{iiinx,degreeform-1}]];
	Return[
		Activate[Simplify[TensorContract[Inactive[TensorProduct][Xd, Xd, Sequence@@Table[gintUU,{IIinx,degreeform-1}]], listindices]]]
	]
];
(*====== Hodge star in the vielbein basis ======*)

ClearAll[MyHStarE];

MyHStarE[X_List,simp_:Identity,gddcoord_:{gdd,coord}]:=MyHStarE/@X;

MyHStarE[X_]:=
Module[{formdegree,DNA,lengthDNA,func\[Eta],mapindices,relevanmatrix\[Eta],signinterm,complementindices,III,term},
	If[
		X==0,
		Return[0]
	];
	formdegree=FormDegree[X]; 
	If[
	formdegree===0,
		Return[Wedge@@(e/@Range[Dim])]
	];
	DNA=DNAofForm[X];
	lengthDNA=Length[DNA];
	func\[Eta]=\[Eta]UU[[#]]&;
	Do[
		relevanmatrix\[Eta] = func\[Eta]/@DNA[[IIIinx,2]];
		mapindices = SparseArray[relevanmatrix\[Eta]]["NonzeroPositions"];
		signinterm = Times@@(relevanmatrix\[Eta][[Sequence@@#]]&/@mapindices);
		complementindices = Complement[Range[Dim],mapindices[[All,2]]];
		term[IIIinx] = DNA[[IIIinx,1]]*signinterm*Signature[{Sequence@@(mapindices[[All,2]]),Sequence@@complementindices}]*(Wedge@@(e/@complementindices))
	,{IIIinx,1,lengthDNA}
	];
	
	Return[Sum[term[IIIinx],{IIIinx,lengthDNA}]];
];
(* ====== Riemann geometry ====== *)

ClearGeometric[]:=Module[{},Clear[ChrisUdd];Clear[Rdd];Clear[RicciScalar];Return[Print["Clear OK - ChrisUdd, Rdd, RicciScalar"]]];
DiffToMatrix[themetric_,coordIn_:coord]:=
	Module[{Dimint},
		Dimint=Length@coordIn;
		Table[
			If[iiinx!=jjinx,
				1/2*Coefficient[Collect[Expand[themetric],d[X_]d[Y_]],d[coordIn[[iiinx]]]d[coordIn[[jjinx]]]],
					Coefficient[Collect[Expand[themetric],d[X_]d[Y_]],d[coordIn[[iiinx]]]d[coordIn[[jjinx]]]]
			]
		,{iiinx,Dimint},{jjinx,Dimint}]
	];


Computegdd[bundle_Association] := 
	Module[{gdd, sqrtdetg, update, copybundle},
		If[
			KeyExistsQ[bundle, "gdd"] ,
				Return[bundle],
					gdd = DiffToMatrix[bundle["ds2"], bundle["coord"]];
					sqrtdetg = Sqrt[-Det[gdd]];
		];
		copybundle = bundle;
		update = AssociateTo[copybundle, <|"gdd" -> gdd, "sqrtdetg" -> sqrtdetg|>];
		Return[update];
		
	]

ComputeChrisUdd[simp_:Identity,gddcoord_:{gdd,coord}] :=
	Module[{dgdd,dGamdd,dChrisdd,initialTime,valuesimp,gddint,coordint,gUUint,Dim},
		If[
		gddcoord=!={gdd,coord},
			gddint=gddcoord[[1]];
			coordint=gddcoord[[2]],
				gddint=gdd;coordint=coord
		];
		gUUint = simp[Inverse[gddint]];
		Dim = Length@coordint;
		initialTime = SessionTime[];
		dgdd[iiinx_,jj_,kk_] := simp[D[gddint[[jj,kk]],coordint[[iiinx]]]];
		dGamdd[kk_,iiinx_,jj_] := simp[1/2 dgdd[iiinx,jj,kk]+1/2 dgdd[jj,iiinx,kk]-1/2 dgdd[kk,iiinx,jj]];
		dChrisdd = Array[dGamdd,{Dim,Dim,Dim}];
		ChrisUdd = simp[gUUint . dChrisdd];
		If[
		simp===Identity,
			valuesimp="without simplification.",
				valuesimp="with simplification."
		];
		Return[
		Print[Style["ChrisUdd",Bold],   "   ",$chrisdef,"  Christoffel symbols computed in ",SessionTime[]-initialTime," sec. ",valuesimp]
		];
	];
Clear[ComputeRdd];
ComputeRdd["conf"] = <|"ComputeChris"->True|>;

ComputeRdd[simp_:Identity,gddcoord_:{gdd,coord}] :=
Module[{useglobalgdd,TrChrisd,term1,term2,term3,term4,auxterm1f,auxterm2f,initialTime,valuesimp,gddint,coordint,Dim},
	initialTime=SessionTime[];
	Clear[Rdd];
	useglobalgdd = True;
	
	If[
	gddcoord=!={gdd,coord},
		Print["** Computing Rdd with metric in arguments"];
		useglobalgdd = False;
		gddint=gddcoord[[1]];
		coordint=gddcoord[[2]];
		Clear[ChrisUdd];
		ComputeChrisUdd[simp,gddcoord];
		,
			gddint=gdd;
			coordint=coord
	];
	
	If[
	(ComputeRdd["conf"]["ComputeChris"])&&useglobalgdd,
		Clear[ChrisUdd,Rdd];
		ComputeChrisUdd[simp,gddcoord];
	];
	
	Dim=Length@coordint;
	
	TrChrisd=simp[TensorContract[ChrisUdd,{{1,2}}]];
	term3=simp[TrChrisd . ChrisUdd];
	term4=simp[Activate@TensorContract[Inactive[TensorProduct][ChrisUdd,ChrisUdd],{{1,5},{3,4}}]];
	auxterm1f[iiinx_,jj_]:=simp[Sum[D[ChrisUdd[[kk,iiinx,jj]],coordint[[kk]]],{kk,Dim}]];
	term1=Array[auxterm1f,{Dim,Dim}];
	auxterm2f[iiinx_,jj_]:=simp[D[TrChrisd[[jj]],coordint[[iiinx]]]];
	term2=Array[auxterm2f,{Dim,Dim}];
	
	Rdd=simp[term1-term2+term3-term4];
	
	If[
	simp===Identity,
		valuesimp="without simplification.",
			valuesimp="with simplification."
	];
	Return[Print[Style["Rdd",Bold],   "   ",$defRicciTensor,"  Ricci tensor computed in ",SessionTime[]-initialTime," sec. ",valuesimp]]
];

ComputeRicciScalar[simp_:Identity,gddcoord_:{gdd,coord}]:=
Module[{InitialTime,valuesimp,gddint,coordint,gUUint},
	If[
	gddcoord=!={gdd,coord},
		gddint=gddcoord[[1]];
		coordint=gddcoord[[2]],
			gddint=gdd;coordint=coord
	];
	Clear[ChrisUdd,Rdd];
	ComputeRdd[simp,gddcoord];
	gUUint=Inverse[gddint];
	InitialTime=SessionTime[];
	RicciScalar=simp[Activate@TensorContract[Inactive[TensorProduct][Rdd,gUUint],{{1,3},{2,4}}]];
	If[
	simp===Identity,
		valuesimp="without simplification.",
			valuesimp="with simplification."
	];
	Return[Print[Style["RicciScalar",Bold],   "   ",$defRicciScalar,"  Ricci scalar computed in ",SessionTime[]-InitialTime," sec. ",valuesimp]];
];

"Here we consider the definition of the contraction operator Contracione that take a p-form in the vielbein basis an 
return a (p-1)-form with a Lorentz index attaced at the beggining."

SetAttributes[inP,Listable](*inP for Inned Product*)
inP[x_,0]=0;inP[x_,y_]:=0/;FormDegree[y]===0;
inP[x_Plus,y_]:=inP[#,y]&/@x
inP[x_,u_*y_]:=u*inP[x,y]/;FormDegree[u]===0;
inP[x_,y_Plus]:=inP[x,#]&/@y
inP[x_.*e[a_],y_.*e[b_]]:=x*y*IdentityMatrix[Dim][[a,b]]
inP[x_.*e[j_],y_.*HoldPattern[Wedge[e[k_],p__]]]:=x*y*(IdentityMatrix[Dim][[j,k]]*Wedge[p]-Wedge[e[k],inP[e[j],Wedge[p]]])
Contractione[X_]:= Table[inP[e[a1111],X],{a1111,Dim}];
ClearAll[SetVielbein]

SetVielbein[eIN_,flatmetric_,simp_:Identity]:=
Module[{},
	ClearAll[\[Eta]dd,\[Eta]UU,eTodx,dxToe,eBasis,gdd,gUU,eamuUd,eamudU];
	\[Eta]dd=flatmetric;
	\[Eta]UU=Inverse[\[Eta]dd];
	eTodx=Table[e[iiinx]->eIN[[iiinx]],{iiinx,Dim}];
	dxToe=Solve[Table[eIN[[iiinx]]==e[iiinx],{iiinx,Dim}],Table[d[coord[[jjinx]]],{jjinx,Dim}]]//Last;
	eBasis=Array[e,{Dim}];
	Do[d[eBasis[[iiinx]]]=simp[(d[eIN]/.dxToe)][[iiinx]],{iiinx,Dim}];
	eamuUd=(FormsToMatrix/@eIN);
	eamudU=Transpose[Inverse[eamuUd]];
	Clear[gdd,gUU];
	gdd=simp[Table[Sum[eamuUd[[a,\[Mu]1]]\[Eta]dd[[a,b]]eamuUd[[b,\[Mu]2]],{a,Dim},{b,Dim}],{\[Mu]1,Dim},{\[Mu]2,Dim}]];
	gUU=simp[Table[Sum[eamudU[[a,\[Mu]1]]\[Eta]UU[[a,b]]eamudU[[b,\[Mu]2]],{a,Dim},{b,Dim}],{\[Mu]1,Dim},{\[Mu]2,Dim}]];
	
	(*-End code-Next only print what is what.*)
	
	Print[
	Style["The following global variables were defined:\n",Purple],Style["eTodx",Bold],
	"  Rule to chage basis from ",PrintIndex["e",{"a"}]," to ",PrintIndex["dx",{"\[Mu]"}],"\n",
	Style["dxToe",Bold],"  Rule to chage basis from ",PrintIndex["dx",{"\[Mu]"}]," to ",PrintIndex["e",{"a"}],"\n",
	Style["eamuUd",Bold],"  Vielbein matrix ",PrintIndex["e",{"a",-"\[Mu]"}],"\n",
	Style["eamudU",Bold],"  Inverse vielbein matrix ",PrintIndex["e",{-"a","\[Mu]"}]," s.t. ",
	PrintIndex["e",{"a",-"\[Mu]"}],PrintIndex["e",{-"a","\[Nu]"}]," = ",PrintIndex["\[Delta]",{"\[Nu]",-"\[Mu]"}]," and ",
	PrintIndex["e",{"a",-"\[Mu]"}],PrintIndex["e",{-"b","\[Mu]"}]," = ",PrintIndex["\[Delta]",{"a",-"b"}]
	,"\n",Style["\[Eta]dd",Bold],"  Flat metric ",PrintIndex["\[Eta]",{-"a",-"b"}],"\n",
	Style["\[Eta]UU",Bold],"  Inverse flat metric ",PrintIndex["\[Eta]",{"a","b"}]
	,"\n",Style["gdd",Bold],"  ",PrintIndex["g",{-"\[Mu]",-"\[Nu]"}]," = ",PrintIndex["\[Eta]",{-"a",-"b"}],
	PrintIndex["e",{"a",-"\[Mu]"}],PrintIndex["e",{"b",-"\[Nu]"}],"\n",Style["gUU",Bold],"  ",PrintIndex["g",{"\[Mu]","\[Nu]"}]," = ",
	PrintIndex["\[Eta]",{"a","b"}],PrintIndex["e",{-"a","\[Mu]"}],PrintIndex["e",{-"b","\[Nu]"}]
	]
];

ComputeSpinConnection[eIN_,flatmetric_,simp_:Identity]:=
Module[{secondterm\[Omega],GUdd},
	SetVielbein[eIN,flatmetric];
	If[
	Dimensions[gdd]!={Dim,Dim},
		Return[Print[Style["Metric ",Red,14],Style["gdd ",Bold,Red,14],Style["not defined\n",Red,14]," The program requieres a global variabled ",
		Style["gdd ",Bold],"being an square array filled by the component of the metric tensor ",PrintIndex["g",{-"\[Mu]",-"\[Nu]"}]]];
	];
	ComputeChrisUdd[simp];
	ClearAll[\[Omega]Ud];ClearAll[\[Omega]dd];
	secondterm\[Omega] = Activate@TensorContract[Inactive[TensorProduct][d[eamuUd],eamudU],{{2,4}}]/.dxToe;
	GUdd = Activate@TensorContract[Inactive[TensorProduct][ChrisUdd,eamuUd,eamudU,eamudU],{{1,5},{2,7},{3,9}}];
	\[Omega]Ud = Activate@TensorContract[Inactive[TensorProduct][GUdd,eBasis],{{2,4}}]-secondterm\[Omega];
	\[Omega]dd = Activate@TensorContract[Inactive[TensorProduct][\[Eta]dd,\[Omega]Ud],{{2,3}}];
	Print[
		Style["\[Omega]Ud",Bold],"  Spin connection 1-form ",PrintIndex["\[Omega]",{"a",-"b"}]," = ",PrintIndex["\[Omega]",{-"c","a",-"b"}],PrintIndex["e",{"c"}],"\n",
		Style["\[Omega]dd",Bold],"  ",PrintIndex["\[Omega]",{-"a",-"b"}]," = ",PrintIndex["\[Eta]",{-"a",-"c"}],PrintIndex["\[Omega]",{"c",-"b"}]
	];
];
(*====== Print of the definition of the objects ======*)

PrintIndex[g_,index_]:=Module[{auxobject},
auxobject=g;
Do[If[Head[index[[II]]]===Times,auxobject=Subscript[auxobject,(-1)*index[[II]]],auxobject=Superscript[auxobject,index[[II]]]],{II,Length@index}];
Return[auxobject]
]

PrintIndices[g_,listdnup_,index_]:=Module[{auxobject},
auxobject=g;
Do[If[listdnup[[II]]===dn,auxobject=Subscript[auxobject,index[[II]]],auxobject=Superscript[auxobject,index[[II]]]],{II,Length@index}];
Return[auxobject]
]

$chrisdef=
Row[{PrintIndices["\[CapitalGamma]",{up,dn,dn},{"\[Mu]","\[Nu]","\[Rho]"}],"=",1/2,PrintIndices["g",{up,up},{"\[Mu]","\[Lambda]"}],PrintIndices["(\[PartialD]",{dn},{"\[Nu]"}],
PrintIndices["g",{dn,dn},{"\[Rho]","\[Lambda]"}],"+",PrintIndices["\[PartialD]",{dn},{"\[Rho]"}],PrintIndices["g",{dn,dn},{"\[Nu]","\[Lambda]"}],"-",PrintIndices["\[PartialD]",{dn},{"\[Lambda]"}],
PrintIndices["g",{dn,dn},{"\[Nu]","\[Rho]"}],")"}];

$defRicciTensor=
Row[{PrintIndices["R",{dn,dn},{"\[Mu]","\[Nu]"}]," = ",PrintIndices["\[PartialD]",{dn},{"\[Rho]"}],PrintIndices["\[CapitalGamma]",{up,dn,dn},{"\[Rho]","\[Mu]","\[Nu]"}],"-",
PrintIndices["\[PartialD]",{dn},{"\[Mu]"}],PrintIndices["\[CapitalGamma]",{up,dn,dn},{"\[Rho]","\[Rho]","\[Nu]"}],"+",PrintIndices["\[CapitalGamma]",{up,dn,dn},{"\[Rho]","\[Rho]","\[Lambda]"}],
PrintIndices["\[CapitalGamma]",{up,dn,dn},{"\[Lambda]","\[Mu]","\[Nu]"}],"-",PrintIndices["\[CapitalGamma]",{up,dn,dn},{"\[Rho]","\[Mu]","\[Lambda]"}],PrintIndices["\[CapitalGamma]",{up,dn,dn},{"\[Lambda]","\[Rho]","\[Nu]"}]}];
$defRicciScalar=
Row[{"R"," = ",PrintIndices["R",{dn,dn},{"\[Mu]","\[Nu]"}],PrintIndices["g",{up,up},{"\[Mu]","\[Nu]"}]}];
(*==========================================================================================================================================*)


(*
|------------------------------------------------------
| Thinking on Association for saving tensor
|------------------------------------------------------

                     --- Utils ---
*)

CleanZeros[X_Association] := KeySelect[X, X[#] =!= 0 &];

ATensorToTensor2sym[Xab_Association, {i_, j_}] := 
	Module[{m,n},
		{m,n} = If[i <= j, {i,j}, {j,i}];
		First[Lookup[Xab, {{m,n}}, 0]]
	];

ATensorToTensor3symLast[Xabc_Association, {p_, i_, j_}] := 
	Module[{q,m,n},
		{q,m,n} = If[i <= j,{p,i,j},{p,j,i}];
		First[Lookup[Xabc, {{q,m,n}}, 0]]
	];

ATensorToTensorRiem[Xabcd_Association, {i_, j_, k_, l_}] := 
	Module[{p,q,m,n,keyMissing,sign=1},
		
		{p,q,m,n} = If[i <= j, {i,j,k,l}, sign = -sign;{j,i,k,l}];
		{p,q,m,n} = If[m <= n, {p,q,m,n}, sign = -sign;{p,q,n,m}];

		keyMissing = Not[KeyExistsQ[Xabcd, {p,q,m,n}]];
		{p,q,m,n} = If[keyMissing, {m,n,p,q}, {p,q,m,n}];
		
		sign*First[Lookup[Xabcd, {{p,q,m,n}}, 0]]
	];

failRequirements[bundle_Association, req_List] :=
	Module[{keyInBundle, keyNegation},
		keyInBundle[key_] := KeyExistsQ[bundle, key];
		keyNegation = Not[Apply[And, Map[keyInBundle, req]]];
		Return[keyNegation];
	];

(*
				---- Metric ----
*)

ComputeAgddAgUU[bundle_Association] := 
	Module[{Dim,gdd,Agdd,gUU,AgUU,buildSymmetric2},
		Dim = Length[bundle["coord"]];
		gdd = DiffToMatrix[bundle["ds2"], bundle["coord"]];
		gUU = Inverse[gdd];
		
		buildSymmetric2[X_] := 
			Association[
				Table[
				{i, j} -> X[[i,j]]
				,{i,Dim}, {j,i, Dim}]
			];
		
		Agdd = buildSymmetric2[gdd];
		AgUU = buildSymmetric2[gUU];

		Agdd = KeySelect[Agdd, Agdd[#] =!= 0 &];
		AgUU = KeySelect[AgUU, AgUU[#] =!= 0 &];
		
		Return[<|"Agdd" -> Agdd, "AgUU" -> AgUU|>];
	];

(*
				---- ChrisUdd ----
*)

ComputeAChrisUdd[bundle_Association] := 
	Module[{coord,Dim,gdd,Agdd,gUU,AgUU,AChrisUdd,dgdd,dGamdd,GamUdd},
		
		coord = bundle["coord"];
		Dim = Length[coord];
		Agdd = bundle["ATensors","Agdd"];
		AgUU = bundle["ATensors","AgUU"];
		gdd[i_,j_] := ATensorToTensor2sym[Agdd, {i,j}];
		gUU[i_,j_] := ATensorToTensor2sym[AgUU, {i,j}];
		
		dgdd[i_,j_,k_] := D[gdd[j,k],coord[[i]]];
		dGamdd[k_,i_,j_] := 1/2*dgdd[i,j,k]+1/2*dgdd[j,i,k]-1/2*dgdd[k,i,j];
		GamUdd[k_,i_,j_] := Sum[gUU[k,l]*dGamdd[l,i,j] ,{l,Dim}];
		
		AChrisUdd = 
			Association[
				Table[
					{i,j,k} -> GamUdd[i,j,k]
				,{i,Dim}, {j, Dim},{k,j,Dim}]
			];
	
		AChrisUdd = CleanZeros[AChrisUdd];
		Return[AChrisUdd];
	];

(*
				---- Riemdddd ----
*)

ComputeARiemdddd[bundle_Association] := 
	Module[{coord,Dim,AChrisUdd, ChrisUdd, dChrisUdd,ChrisChrisUddd,RiemUddd,Agdd,gdd, Riemdddd,seen,list, 
			ARiemdddd,ChrisUddArray,dChrisUddArray,RiemUdddArray,ChrisChrisUdddArray,gddArray,RiemddddArray,ChrisUddArrayToDer},
		
		If[
			failRequirements[bundle["ATensors"], {"AChrisUdd"}],
				Return["Missing Tensors/AChrisUdd"]
		];
		
		coord = bundle["coord"];
		Dim = Length[coord];
		Agdd = bundle["ATensors","Agdd"];
		AChrisUdd = bundle["ATensors","AChrisUdd"];
		gdd[i_,j_] := ATensorToTensor2sym[Agdd, {i,j}];
		gddArray = SparseArray[Array[gdd,{Dim,Dim}]];
		ChrisUdd[i_,j_,k_] := ATensorToTensor3symLast[AChrisUdd, {i,j,k}]; 
		ChrisUddArrayToDer = Array[ChrisUdd, {Dim, Dim, Dim}];
		ChrisUddArray = SparseArray[ChrisUddArrayToDer];
		dChrisUddArray = SparseArray[Transpose[Table[D[ChrisUddArrayToDer,coord[[i]]],{i,Dim}],{3,1,4,2}]];

		(*dChrisUdd[i_,j_,k_,l_] := dChrisUddArray[[i,j,k,l]];*)

		
		(*dChrisUdd[i_,j_,k_,l_] := D[ChrisUdd[j,k,l],coord[[i]]];*)
		ChrisChrisUdddArray = SparseArray[Transpose[ChrisUddArray . ChrisUddArray,{1,3,4,2}]];
		(*Sum[ChrisUdd[i,k,n]*ChrisUdd[n,l,j] ,{n,Dim}];*)
		RiemUdddArray = dChrisUddArray - Transpose[dChrisUddArray, {1,2,4,3}] + ChrisChrisUdddArray - Transpose[ChrisChrisUdddArray, {1,2,4,3}];
		(*RiemUddd[i_,j_,k_,l_] := dChrisUdd[k,i,l,j]-dChrisUdd[l,i,k,j] + ChrisChrisUddd[i,j,k,l]-ChrisChrisUddd[i,j,l,k];
		Riemdddd[i_,j_,k_,l_] := Sum[ gdd[i,n]*RiemUddd[n,j,k,l] ,{n,Dim}];*)
		RiemddddArray = gddArray . RiemUdddArray;
		seen = <||>;
		list = {};
		Print[Dimensions[RiemUdddArray]];
		Do[
			If[KeyExistsQ[seen, {k, l, i, j}], Continue[]];
			AppendTo[list, {i, j, k, l} -> RiemddddArray[[i, j, k, l]]];
			seen[{i, j, k, l}] = True;
			seen[{k, l, i, j}] = True;
		,
		{i, Dim}, {j, i + 1, Dim}, {k, Dim}, {l, k + 1, Dim}
		];
  
		ARiemdddd = Association[list];
		ARiemdddd = CleanZeros[ARiemdddd];
		Return[ARiemdddd];
	];
	
(*
				---- Ricdd ----
*)

	
ComputeARicdd[bundle_Association] := 
	Module[{coord,Dim,AgUU,ARiemdddd,gUU,Rdddd,Ricdd,ARicdd,gUUArray,RddddArray },
		
		If[
			failRequirements[bundle["ATensors"], {"ARiemdddd"}],
				Return["Missing ATensors/ARiemdddd"]
		];
		
		coord = bundle["coord"];
		Dim = Length[coord];
		AgUU = bundle["ATensors","AgUU"];
		ARiemdddd = bundle["ATensors","ARiemdddd"];
		gUU[i_,j_] := ATensorToTensor2sym[AgUU, {i,j}];
		Rdddd[i_,j_,k_,l_] := ATensorToTensorRiem[ARiemdddd,{i,j,k,l}];
		gUUArray = SparseArray[Array[gUU,{Dim,Dim}]];
		RddddArray = SparseArray[Array[Rdddd,{Dim,Dim,Dim,Dim}]];
		(*Ricdd[i_,j_] :=  Sum[gUU[m,n]*Rdddd[m,i,n,j],{m,Dim}, {n,Dim}];*)
		Ricdd = Activate[TensorContract[Inactive[TensorProduct][gUUArray, RddddArray],{{1,3},{2,5}}]];
		ARicdd = 
			Association[
				Table[{i, j} -> Ricdd[[i,j]],{i,Dim}, {j,i, Dim}]
			];
		ARicdd = CleanZeros[ARicdd];
		Return[ARicdd];
	];

(*
				---- RicciScalar ----
*)

ComputeARicciScalar[bundle_Association] := 
	Module[{coord,Dim,AgUU,ARiemdddd,gUU,Rdddd,Ricdd,ARicdd,RicciScalar},
		
		If[
			failRequirements[bundle["ATensors"], {"ARicdd"}],
				Return["Missing ATensors/ARicdd"]
		];
		
		coord = bundle["coord"];
		Dim = Length[coord];
		AgUU = bundle["ATensors","AgUU"];
		ARicdd = bundle["ATensors","ARicdd"];
		gUU[i_,j_] := ATensorToTensor2sym[AgUU, {i,j}];
		Ricdd[i_,j_] := ATensorToTensor2sym[ARicdd,{i,j}];
		RicciScalar = Sum[gUU[i,j]*Ricdd[i,j],{i,Dim},{j,Dim}];
		Return[RicciScalar];
	];
	

(*
		    ------------------------------------
			---       Orcheta Director       ---
			------------------------------------
*)

ComputeBundleTensors[bundleIN_Association, level_: "Rdddd", simp_:Identity] := 
Module[{ATensors, needMetric, needChris, needRiemann, AgddgUU, Agdd, AgUU, 
	AChrisUdd, ARiemdddd, bundle, ARicdd,needRicci, needRicciScalar, RicciScalar},
	
	bundle = bundleIN;
	ATensors = Lookup[bundle, "ATensors", <||>];

	(* --- Check what is missing --- *)
	needMetric   = Not[KeyExistsQ[ATensors, "Agdd"]] || Not[KeyExistsQ[ATensors, "AgUU"]];
	needChris    = Not[KeyExistsQ[ATensors, "AChrisUdd"]] && MemberQ[{"RicciScalar","ChrisUdd", "Rdddd", "Rdd"}, level];
	needRiemann  = Not[KeyExistsQ[ATensors, "ARiemdddd"]] && MemberQ[{"RicciScalar","Rdddd", "Rdd"}, level];
	needRicci    = Not[KeyExistsQ[ATensors, "ARicdd"]] && MemberQ[{"RicciScalar","Rdd"}, level];
	needRicciScalar = Not[KeyExistsQ[ATensors, "ARicciScalar"]] && MemberQ[{"RicciScalar"}, level];
	
	(* --- Metric --- *)
	If[needMetric,
		Print["** Computing metric"];
		AgddgUU = ComputeAgddAgUU[bundle];
		Agdd = Map[simp, AgddgUU["Agdd"]];
		AgUU = Map[simp, AgddgUU["AgUU"]];
		ATensors = Join[ATensors, <|"Agdd" -> Agdd, "AgUU" -> AgUU|>];
		bundle = AssociateTo[bundle, "ATensors" -> ATensors];
	];
	
	If[level === "metric", Return[bundle]];
	
	(* --- Christoffel --- *)
	If[needChris,
		Print["** Computing Christoffel"];
		AChrisUdd = ComputeAChrisUdd[bundle];
		AChrisUdd = Map[simp, AChrisUdd];
		ATensors = AssociateTo[ATensors, "AChrisUdd" -> AChrisUdd];
		bundle = AssociateTo[bundle, "ATensors" -> ATensors];
	];

	If[level === "ChrisUdd", Return[bundle]];

	(* --- Riemann --- *)
	If[needRiemann,
		Print["** Computing Riemann"];
		ARiemdddd = ComputeARiemdddd[bundle];
		Print["    -Simplifying"];
		ARiemdddd = Map[simp, ARiemdddd];
		ATensors = AssociateTo[ATensors, "ARiemdddd" -> ARiemdddd];
		bundle = AssociateTo[bundle, "ATensors" -> ATensors];
	];
	
	If[level === "Rdddd", Return[bundle]];
	
	(* --- Ricci --- *)
	If[needRicci,
		Print["** Computing Ricci tensor"];
		ARicdd = ComputeARicdd[bundle];
		ARicdd = Map[simp, ARicdd];
		ATensors = AssociateTo[ATensors, "ARicdd" -> ARicdd];
		bundle = AssociateTo[bundle, "ATensors" -> ATensors];
	];
	If[level === "Rdd", Return[bundle]];
	
	(* --- RicciScalar --- *)
	If[needRicciScalar,
		Print["** Computing RicciScalar"];
		RicciScalar = simp[ComputeARicciScalar[bundle]];
		ATensors = AssociateTo[ATensors, "ARicciScalar" -> RicciScalar];
		bundle = AssociateTo[bundle, "ATensors" -> ATensors];
	];
	
	If[level === "RicciScalar", Return[bundle]];

];


(*
    ---- "Conneting with previous notation" ----
*)

NewComputeRdd[bundleIN_Association, simp_:Identity] := 
	Module[{bundle,ARicdd, Ricdd, Dim, Rdd},
		bundle = bundleIN;
		Dim = Length[bundle["coord"]];
		bundle = ComputeBundleTensors[bundle, "Rdd", simp];
		ARicdd = bundle["ATensors", "ARicdd"];
		Ricdd[i_,j_] := ATensorToTensor2sym[ARicdd, {i, j}];
		Rdd = Array[Ricdd, {Dim,Dim}];
		Return[{Rdd,bundle}];
	]


Clear[BuildHodge];
BuildHodge[bundleIN_Association, simp_:Simplify] := 
	Module[{gdd,sqrtdetg,needMetric,ATensors,bundle,AgddgUU,Agdd,AgUU,Dim,gddMap,
		gUU,gUUMap,coord},
		bundle = bundleIN;
		coord = bundle["coord"];
		Dim = Length[coord];
		ATensors = Lookup[bundle, "ATensors", <||>];
	
		needMetric   = Not[KeyExistsQ[ATensors, "Agdd"]] || Not[KeyExistsQ[ATensors, "AgUU"]];
		If[needMetric,
				bundle = ComputeBundleTensors[bundle, "metric", simp]
		];
		Agdd = Map[Simplify, bundle["ATensors","Agdd"]];
		AgUU = Map[Simplify, bundle["ATensors","AgUU"]];
		gddMap[i_,j_] := ATensorToTensor2sym[Agdd,{i,j}];
		gUUMap[i_,j_] := ATensorToTensor2sym[AgUU,{i,j}];
		gdd = Array[gddMap,{Dim,Dim}];
		gUU = Array[gUUMap,{Dim,Dim}];
		sqrtdetg = Sqrt[-Det[gdd]];
		If[KeyExistsQ[bundle, "assum"],
				sqrtdetg = Simplify[sqrtdetg,bundle["assum"]],
					sqrtdetg = Simplify[sqrtdetg]
		];
		Return[Function[{X}, 
			MyHStar[X, simp, {gdd, coord, gUU, sqrtdetg}]
			]];
	];

