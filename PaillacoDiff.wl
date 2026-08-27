(* ::Package:: *)

BeginPackage["PaillacoDiff`"]

(* ---------- Public functions ---------- *)

FormDegree::usage = "FormDegree[expr] returns the degree of a differential form (0 for scalars)."
Wedge::usage = "Wedge[x, y, ...] is the exterior (wedge) product of forms."
d::usage = "d[expr] is the exterior derivative."
PolyFormQ::usage = "PolyFormQ[expr] tests whether expr is a sum of forms of different degrees."

Extractor::usage = "Extractor[F, A] extracts the coefficient of 1-form A in polyform F."
Extractorleft::usage = "Extractorleft[F, A] extracts A from the left side of each term."
coordcontraction::usage = "coordcontraction[X, coord] contracts X with all coordinate 1-forms."

DNAofForm::usage = "DNAofForm[X] decomposes form X into {{coeff, indices}, ...}."
SparseFromDNA::usage = "SparseFromDNA[DNA, dim, deg] converts DNA to a SparseArray."
DNAFromSparse::usage = "DNAFromSparse[sparse] converts a SparseArray back to DNA."
BuildSquaresTools::usage = "BuildSquaresTools[bundle] builds FormSquare/FormSquaredd closures."
FormSquare::usage = "FormSquare[X] computes F_{mu1...mup} F^{mu1...mup}."
FormSquaredd::usage = "FormSquaredd[X] computes F_{mu l2...lp} F_nu^{ l2...lp}."
Hstar::usage = "Hstar[X] computes the Hodge dual of form X."
FormToSparse::usage = "FormToSparse[X] converts a form to a SparseArray."
FormsToMatrix::usage = "FormsToMatrix[X, deg, coord] converts a form X to a dense matrix. deg is an Integer and the degree of the form, and coord are the coordinates."

ClearGeometric::usage = "ClearGeometric[] clears global tensors ChrisUdd, Rdd, RicciScalar."
DiffToMatrix::usage = "DiffToMatrix[ds2, coord] extracts the metric tensor from a line element."
Computegdd::usage = "Computegdd[bundle] computes gdd and sqrtdetg for a bundle."
ComputeChrisUdd::usage = "ComputeChrisUdd[] computes Christoffel symbols from global gdd, coord."
ComputeRdd::usage = "ComputeRdd[] computes the Ricci tensor."
ComputeRicciScalar::usage = "ComputeRicciScalar[] computes the Ricci scalar."
SetVielbein::usage = "SetVielbein[eIN, eta] sets up the vielbein basis and defines global variables."
ComputeSpinConnection::usage = "ComputeSpinConnection[eIN, eta] computes the spin connection 1-form."

InitMetricBundle::usage = "InitializeBundle[bundle] initialises a bundle Association."
TensorProductContract::usage = "TensorProductContract[t1, t2, ..., {{i1,j1}, ...}] contracts tensor products."
RaiseIndices::usage = "RaiseIndices[sparse, bundle, positions] raises specified indices."
PaiCovD::usage = "PaiCovD[bundle, tensor, indices] computes the coordinate-basis covariant derivative of tensor. indices is a string of U/d characters describing tensor index variance. For instace for  tensor TUdU indices must be the string UdU. The covariant derivative index is added at the beginning of the tensor"
GetTensorArray::usage = "GetTensorArray[bundle, name] retrieves a tensor array, computing on demand."
PaiComputeMetric::usage = "PaiComputeMetric[bundle] computes metric from bundle's ds2."
PaiComputeChrisUdd::usage = "PaiComputeChrisUdd[bundle] computes Christoffel symbols from bundle."
PaiComputeRdddd::usage = "PaiComputeRdddd[bundle] computes the Riemann tensor."
PaiComputeRdd::usage = "PaiComputeRdd[bundle] computes the Ricci tensor."
PaiComputeRicciScalar::usage = "PaiComputeRicciScalar[bundle] computes the Ricci scalar."
PaiComputeBundleTensors::usage = "PaiComputeBundleTensors[bundle, level] computes tensors up to a requested level."

BuildHodge::usage = "BuildHodge[bundle] builds a Hodge star function for a bundle."
BuildHodgeMetric::usage = "BuildHodgeMetric[bundle] builds a coordinate-basis Hodge star function."
BuildHodgeVielbein::usage = "BuildHodgeVielbein[bundle] builds a vielbein-basis Hodge star function."
PaiComputeSpinConnection::usage = "PaiComputeSpinConnection[bundle] computes spin connection in bundle."
PaiComputeCurvatureForm::usage = "PaiComputeCurvatureForm[bundle] computes curvature 2-form."
PaiComputeRddddFlat::usage = "PaiComputeRddddFlat[bundle] computes Riemann in flat (vielbein) basis."
PaiComputeRddFlat::usage = "PaiComputeRddFlat[bundle] computes Ricci in flat basis."
PaiComputeRicciScalarFlat::usage = "PaiComputeRicciScalarFlat[bundle] computes Ricci scalar in flat basis."

(* ---------- Public globals ---------- *)

PaiNonCommutativeScalarQ::usage =
  "PaiNonCommutativeScalarQ[expr] tests whether expr contains a registered noncommutative scalar coefficient.";
PaiRegisterNonCommutativeScalarQ::usage =
  "PaiRegisterNonCommutativeScalarQ[test] registers a predicate test[expr] used by Wedge to detect noncommutative scalar coefficients.";

coord::usage = "List of coordinate variables."
Dim::usage = "Spacetime dimension."
gdd::usage = "Metric tensor g_{mu nu}."
gUU::usage = "Inverse metric g^{mu nu}."
ChrisUdd::usage = "Christoffel symbols Gamma^mu_{nu rho}."
Rdddd::usage = "Riemann tensor R_{mu nu rho sigma}."
Rdd::usage = "Ricci tensor R_{mu nu}."
RicciScalar::usage = "Ricci scalar R."
sqrtdetg::usage = "Sqrt[-det(g)]."

\[Eta]dd::usage = "Flat (Minkowski) metric."
\[Eta]UU::usage = "Inverse flat metric."
e::usage = "Vielbein basis 1-forms e^a."
eTodx::usage = "Rule mapping e^a to e^a_mu dx^mu."
dxToe::usage = "Rule mapping dx^mu to e^a."
eamuUd::usage = "Vielbein matrix e^a_mu."
eamudU::usage = "Inverse vielbein matrix e_a^mu."
eBasis::usage = "List of vielbein basis symbols {e[1], ..., e[Dim]}."

\[Omega]Ud::usage = "Spin connection 1-form omega^a_b."
\[Omega]dd::usage = "Spin connection omega_{ab} (both indices down)."

Begin["`Private`"]

ClearAll[GlobalRequired];
SetAttributes[GlobalRequired, HoldAll];

GlobalRequired::missing =
	"Global Mode requires `1` to be defined.";

GlobalRequired[x_Symbol] :=
	If[
		!ValueQ[x],
		Message[GlobalRequired::missing, HoldForm[x]];
		Abort[]
	];

GlobalRequired[x_Symbol, xs__Symbol] := (
	GlobalRequired[x];
	GlobalRequired[xs]
);

ClearAll[ResolveGlobal];
SetAttributes[ResolveGlobal, HoldRest];

ResolveGlobal[varIn_, varGlob_Symbol, func_:Identity] := If[
	varIn === "Global",
		GlobalRequired[varGlob];
		func[varGlob],
			varIn
];

(*
	--------- Form Degree --------- 

*)

ClearAll[FormDegree]
FormDegree[d[x_]]:=1+FormDegree[x];
FormDegree[x_]:=0;
FormDegree[e[n_Integer]]:=1;
FormDegree[x_Wedge]:=Plus@@Map[FormDegree,List@@x];
FormDegree[x_Times]:=Plus@@Map[FormDegree,List@@x];
FormDegree[x_Plus]:=FormDegree[First@x];
FormDegree[x_List]:=FormDegree/@x;


(*=============== Wedge ===============*)

$PaiNonCommutativeScalarTests = {};

PaiRegisterNonCommutativeScalarQ[test_] := Module[{},
  If[
    ! MemberQ[$PaiNonCommutativeScalarTests, test],
    AppendTo[$PaiNonCommutativeScalarTests, test]
  ];
  test
];

PaiNonCommutativeScalarQ[expr_] :=  AnyTrue[$PaiNonCommutativeScalarTests, TrueQ[#[expr]] &];

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
Wedge[x__,Times[sca_,y_]]:=Times[sca,Wedge[x,y]]/;NumericQ[sca]||(FormDegree[sca]===0&&(!PaiNonCommutativeScalarQ[sca]||!PaiNonCommutativeScalarQ[{x}]));
Wedge[Times[sca_,x_],y__]:=Times[sca,Wedge[x,y]]/;NumericQ[sca]||(FormDegree[sca]===0&&!PaiNonCommutativeScalarQ[sca]);
Wedge[x_,y___,x_]:=0/;OddQ[FormDegree[x]]&&!PaiNonCommutativeScalarQ[x]&&!PaiNonCommutativeScalarQ[{y}];
Wedge[y__]:=Signature[{y}]*Wedge@@Sort[{y}]/;Sort[{y}]=!={y}&&Union[FormDegree[{y}]]==={1}&&!PaiNonCommutativeScalarQ[{y}];


WedgeDot[a_, b_] := Inner[Wedge, a, b, Plus];
Wedge[A_List,B_List] := WedgeDot[A,B];
(*====== Exterior derivative (d) ======*)

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

(*====== Towards contraction operator =====*)

Clear[Extractor];
(*Extractor[F3_List,y_]:=Extractor[#,y]&/@F3;*)

Extractor[xIn__,oneform_, side_:"right"]:=
Module[{listx,listxref,killz, x},
	
	killz[pform__,form_]:=
		Module[{pos, signJumpRight, signJumpLeft},
			pos = Flatten[Position[{pform},form]];
			If[pos==={},
				Return[0]
			];
			
			signJumpRight = (-1)^(Length[{pform}]-pos);
			signJumpLeft = (-1)^(pos-1);
			
			Which[
			side==="right",
				Return[(Wedge@@DeleteCases[{pform},form])*signJumpRight],
			side==="left",	
				Return[(Wedge@@DeleteCases[{pform},form])*signJumpLeft],
			True,
				Print["pick side for the extraction: left or right. right by default"]
				Return[Apply[Wedge,pform]]
			];
			];
	x = Expand[xIn];	
	listx = 
		If[Head[x]===Plus,
			List@@x,
				{x}];
	
	listxref = Select[listx,(!FreeQ[#,oneform])&];
	
	Return[
		Plus@@Flatten[
			Map[#/.Wedge[y__] :> killz[y,oneform]/.oneform -> 1&, listxref]
		]
	];
	];
Extractor::usage="Given a polyform F = F1+ F2\[Wedge]A,  with A a 1-form,  Extractor[F,A] returns F2.";

Clear[Extractorleft];
Clear[coordcontraction];
Extractorleft[F3_List,y_] := Extractor[#,y]&/@F3;
Extractorleft[x__,oneform_] := Module[
	{listx,listxref,killzleft},
	killzleft[y__,KK_] := Module[
		{pos},
		pos = Flatten[Position[{y},KK]];
		If[pos==={},Return[0]];
		Return[(Wedge@@DeleteCases[{y},KK])*(-1)^(pos-1)]
	];
	listx=If[Head[x]===Plus, List@@x,{x}];
	listxref=Select[listx,(!FreeQ[#,oneform])&];
	Return[Plus@@Flatten[#/.Wedge[y__]:>killzleft[y,oneform]/.oneform->1&/@listxref]];
];
coordcontraction[X_List, coord_:coord]:=Map[coordcontraction[#,coord]&,X];
coordcontraction[X_, coord_:coord]:=Extractorleft[X,#]&/@d[coord];


(*

|------------------------------------------------------------------------------------------------------------
|     p-forms operations (FormSquaredd, FormSquare, Hstar)
|------------------------------------------------------------------------------------------------------------

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

coeffBaseElement[pform_, base_]:=
	Module[{baseElement,deg, coeff, 
		baseAlong, Dim, mapcoord,baseNumb, sign},
		
		deg = FormDegree[pform];
		Dim = Length[base];
		mapcoord = AssociationThread[base -> Range[Dim]];
		
		Which[
		(deg===1),
			baseElement = Cases[pform, Apply[Alternatives, base], {0, Infinity}];
			baseNumb = baseElement /. mapcoord
			,
		(deg > 1),
			baseElement = Cases[pform, _Wedge, {0, Infinity}];
			baseNumb = baseElement/.Wedge[YY__]:>{YY} /. mapcoord;
		];
		
		baseNumb = Flatten[baseNumb];
		sign = Signature[baseNumb];
		baseNumb = Sort[baseNumb];
		
		If[
		Length[baseElement]=!=1,
			Message[DNAofForm::noBaseFound, baseElement];
			Return[$Failed]
		];

		baseAlong = baseElement[[1]];
		coeff = Coefficient[pform, baseAlong];
	Return[{sign*coeff, baseNumb}];
	];

DNAofForm::noBaseFound = 
	"No base element found";
Clear[DNAofForm];

DNAofForm[FormIn_, base_:"Global"] := Module[
	{Collected, baseint, formAsList},

	If[FormDegree[FormIn]===0,
		Return[{FormIn}]
	];
	
	baseint = ResolveGlobal[base, coord, d];

	Collected = 
		Which[
		FormDegree[FormIn]===1,
			Collect[Expand@FormIn, baseint],
		(FormDegree[FormIn]>1)||PolyFormQ[FormIn],
			Collect[Expand@FormIn, _Wedge]
		];

	formAsList = 
		If[
		Head[Collected] === Plus,
			Apply[List,Collected],
				{Collected}
			];
	
	Return[Map[coeffBaseElement[#, baseint]&, formAsList]];
];
	
(*
   ---- Tools DNA and Sparse Array ----
*)

Clear[SparseFromDNA];
SparseFromDNA[DNA_List, Dim_Integer, formdeg_Integer]:=
Module[{rules = <||>, perms, sign, base, comp},
	If[DNA ==={0},
		Return[SparseArray[{}, ConstantArray[Dim, formdeg]]]
	];

	Do[
		comp = CompBase[[1]];
		base = CompBase[[2]];
		perms = Permutations[base];
		Do[
			sign = Signature[perm];
			rules[perm] = sign * comp
		,
		{perm, perms}]
		,
	{CompBase, DNA}];
	SparseArray[Normal[rules], ConstantArray[Dim, formdeg]]
 ];

Clear[DNAFromSparse];
DNAFromSparse[tensor_SparseArray] := 
Module[{tensorRules,indepComps,getDNAcomp},
	tensorRules = ArrayRules[tensor];
	indepComps = 
		Select[
			DeleteDuplicates[
				Map[Sort,Keys[tensorRules]]
			]
		,VectorQ[#,IntegerQ]&];
		
	getDNAcomp[tRules_, comp_] := Apply[List,Reverse[Select[tRules, First[#]===comp&][[1]]]];
	
	Table[
		getDNAcomp[tensorRules, comp]
	,
	{comp, indepComps}]
];

(*====== Squares of differential forms ======*)

(* ---- Tools ---- *)
nonzeroCompsFirst[Tensorsparse_SparseArray, mu_Integer] :=
	Module[{tensorList,keysWithDuplicates},
		tensorList = ArrayRules[Tensorsparse];
		keysWithDuplicates = 
			Map[
			Delete[#,1]&, 
				Select[Keys[tensorList],(#[[1]]===mu)&]
			];

		DeleteDuplicates[Map[Sort,Select[keysWithDuplicates, Signature[#]=!=0&]]]
	];
	
Clear[nonzeroComps];

nonzeroComps[Tensorsparse_SparseArray] :=
	Module[{tensorList,keysWithDuplicates},
		tensorList = ArrayRules[Tensorsparse];
		DeleteDuplicates[Map[Sort,Select[Keys[tensorList],(Signature[#]=!=0)&]]]
	];

Clear[getCompRule];
getCompRule[tensor_List, comp_List] := First[Values[Select[tensor, #[[1]]===comp&]]];

Clear[RaiseAllSparse];
RaiseAllSparse[FformSparse_SparseArray, gUU_SparseArray, formdegree_Integer] :=
	Module[{seqgUU, indicesContract, FtensorSparse},
		seqgUU = Sequence@@Table[gUU, {IIinx,formdegree}];
		indicesContract = Table[{inx, formdegree + 2*inx-1}, {inx, formdegree}];
		Return[Activate[TensorContract[Inactive[TensorProduct][FformSparse,seqgUU], indicesContract]]];
	];
(* Build Form Square and FormSquaredd*)

Clear[BuildSquaresTools];
SetAttributes[BuildSquaresTools, HoldFirst];
BuildSquaresTools[bundle_, simp_:Simplify] := Module[{gUU, eta, etainv, basis, dxToe},
	Which[
	KeyExistsQ[bundle, "ds2"]===True,
		PaiComputeBundleTensors[bundle, "metric", simp];
		gUU = GetTensorArray[bundle, "gUU"];
		basis = d[bundle["coord"]];
		Return[
		   <| "FormSquare" -> Function[{X}, FormSquare[X, gUU, simp,  basis]],
		      "FormSquaredd" -> Function[{X}, FormSquaredd[X, gUU, simp,  basis]]
		    |>
		],
	KeyExistsQ[bundle, "eU"]===True,
		eta = DiagonalMatrix[bundle["signature"]];
		etainv = Inverse[eta];
		basis = bundle["basis"];
		dxToe = bundle["dxToe"];
		Return[
		    <| "FormSquare" -> Function[{X}, FormSquare[X /. dxToe, etainv, simp,  basis]],
		    "FormSquaredd" -> Function[{X}, FormSquaredd[X /. dxToe, etainv, simp,  basis]]
		    |>
		],
	True,
		Print["[ Aborting ] BuildSquaresTools: ds2 nor eU not given"];
		Abort[];
	];
];


Clear[FormSquare];
FormSquare[Xform_, gUUIN_:"Global", simp_:Identity, basisIN_:"Global"] :=
Module[{deg,FformDNA, FformSparse,gintUU,listindices,seqgUU, FtensorSparse,
	indicesContract,FormComps,TensorComps,InterComps,FformRule,FtensorRule,
	FformValues, FtensorValues, basisint,Dim},
	If[
	Xform ===0,
		Return[0]
	];

	gintUU = SparseArray[ResolveGlobal[gUUIN, gUU]];
	basisint = ResolveGlobal[basisIN, coord, d];

	deg = FormDegree[Xform];
	Dim = Length[basisint];
	FformDNA = simp[DNAofForm[Xform, basisint]];
	FformSparse = SparseFromDNA[FformDNA, Dim, deg];
	FtensorSparse = RaiseAllSparse[FformSparse, gintUU, deg];
	
	FormComps   = nonzeroComps[FformSparse];
	TensorComps = nonzeroComps[FtensorSparse];
	InterComps  = Intersection[FormComps , TensorComps];
	
	FformRule   = ArrayRules[FformSparse];
	FtensorRule = ArrayRules[FtensorSparse];
	
	FformValues   = simp[Map[getCompRule[FformRule,  #]&, InterComps]];
	FtensorValues = simp[Map[getCompRule[FtensorRule, #]&, InterComps]];
	
	Return[(FformValues . FtensorValues)*(deg)!]
];


Clear[FormSquaredd];
FormSquaredd[0,__]:=0

FormSquaredd[Xform_, gintUUinput_:"Global", simp_:Identity, basisIN_:"Global"] :=
Module[{deg,FformDNA,FformSparse,gintUU,
	basisint,Dim,seqgUU,indexcontr, nonzeroUp, nonzeroDn, nonzeroInter,
	nonzeroInterUp, nonzeroInterDn,FformRule,FtensorRule,nonzeroXd,nonzeroXdU,Xsqdd,
	Xdmunu, XdUmunu, FtensorSparse},
	
	gintUU = SparseArray[ResolveGlobal[gintUUinput, gUU]];
	basisint = ResolveGlobal[basisIN, coord, d];

	Dim = Length[basisint];
	deg = FormDegree[Xform];
	FformDNA = simp[DNAofForm[Xform, basisint]];
	FformSparse = SparseFromDNA[FformDNA, Dim, deg];
	
	seqgUU = Sequence@@Table[gintUU,{IIinx,deg-1}];
	
	indexcontr = Table[{iiinx+1,deg+2*iiinx-1},{iiinx,deg-1}];
	FtensorSparse = Activate[TensorContract[Inactive[TensorProduct][FformSparse,seqgUU],indexcontr]];
	
	nonzeroDn[mu_] := nonzeroCompsFirst[FformSparse, mu];
	nonzeroUp[nu_] := nonzeroCompsFirst[FtensorSparse, nu];
	nonzeroInter[mu_, nu_] := Intersection[nonzeroDn[mu],nonzeroUp[nu]];
	
	nonzeroInterDn[mu_,nu_] := Map[Prepend[#,mu]&,nonzeroInter[mu,nu]];
	nonzeroInterUp[mu_,nu_] := Map[Prepend[#,nu]&,nonzeroInter[mu,nu]];
	
	FformRule = ArrayRules[FformSparse];
	FtensorRule = ArrayRules[FtensorSparse];
		
	Table[nonzeroXd[mu,nu]  = Map[getCompRule[FformRule,  #]&, nonzeroInterDn[mu,nu]],{mu,Dim},{nu,mu,Dim}];
	Table[nonzeroXdU[mu,nu] = Map[getCompRule[FtensorRule, #]&, nonzeroInterUp[mu,nu]],{mu,Dim},{nu,mu,Dim}];
	Xsqdd = ConstantArray[0,{Dim,Dim}];
	Do[
		Xdmunu = simp[nonzeroXd[mu,nu]];
		XdUmunu = simp[nonzeroXdU[mu,nu]];
		Xsqdd[[mu,nu]] = Xdmunu . XdUmunu;
		Xsqdd[[nu,mu]] = Xsqdd[[mu,nu]];
	,{mu, Dim}, {nu, mu, Dim}];
	
	Return[((deg-1)!)*Xsqdd];

];


Clear[Hstar];
Hstar[Xform_, gintUUIN_:"Global", sqrtdetgIN_:"Global", baseIN_:"Global", simp_:Identity] := 
	Module[{gintUU, coordint, Dim, deg, FformDNA, FformSparse, FtensorSparse,
		TensorComps,FtensorRule,FtensorValues, FtensorDict, compToStar,starF, sqrtdetgint, baseint},

		If[
		Xform ===0,
			Return[0]
		];

		gintUU = SparseArray[ResolveGlobal[gintUUIN, gUU]];
		sqrtdetgint = ResolveGlobal[sqrtdetgIN, sqrtdetg];
		baseint = ResolveGlobal[baseIN, coord, d];
		
		Dim = Length[baseint];
		deg = FormDegree[Xform];
		
		If[
			deg===0,
				Return[Xform*sqrtdetgint Wedge@@(baseint)]
		];
		
		FformDNA = simp[DNAofForm[Xform, baseint]];
		FformSparse = SparseFromDNA[FformDNA, Dim, deg];
		FtensorSparse = RaiseAllSparse[FformSparse, gintUU, deg];
		
		TensorComps = nonzeroComps[FtensorSparse];
		FtensorRule = ArrayRules[FtensorSparse];
		FtensorValues   = simp[Map[getCompRule[FtensorRule,  #]&, TensorComps]];
		
		FtensorDict = AssociationThread[TensorComps, FtensorValues];
		
		compToStar[formcomp_] :=
			Module[{complement, toepsilon},
				complement = Complement[Range[Dim],formcomp];
				toepsilon = Flatten[{formcomp,complement}];
				<|"eps"->toepsilon, "basis" -> Map[baseint[[#]]&, complement]|>
			];
		
		starF = 
			sqrtdetgint*Sum[
				FtensorDict[comp]*Signature[compToStar[comp]["eps"]]*Apply[Wedge, compToStar[comp]["basis"]]
			,
			{comp, TensorComps}];
		Return[starF]
	];



Clear[FormToSparse];
Clear[FormsToMatrix];
FormToSparse[X_, formdegIN_:"deg", coordIN_:"Global"] :=
Module[{coordint, Dimint, formdegint},
	coordint = ResolveGlobal[coordIN, coord];
	Dimint = Length[coordint];
	Which[
		(X===0)&&(formdegIN==="deg"),
			Print["** Not possible to generate array. Null p-from and p is not given."];
			Return[0],
		(X===0),
			Return[SparseArray[{}, ConstantArray[Dimint, formdegIN]]]
	];
	formdegint = FormDegree[X];
	Return[SparseFromDNA[DNAofForm[X, d[coordint]], Dimint,formdegint]];
];

FormsToMatrix[X_, formdegIN_:"deg", coordIN_:"Global"] := Normal[FormToSparse[X, formdegIN, coordIN]];

(* ====== Riemann geometry ====== *)

ClearGeometric[]:=Module[{},Clear[ChrisUdd];Clear[Rdd];Clear[RicciScalar];Return[Print["Clear OK - ChrisUdd, Rdd, RicciScalar"]]];

DiffToMatrix[themetric_, coordIn_:"Global"] := Module[
	{Dimint, coordMod},

	coordMod = ResolveGlobal[coordIn, coord];
	metricCollected = Collect[Expand[themetric], _d];

	Table[
		If[xIter=!=yIter,
			1/2*Coefficient[metricCollected, d[xIter]*d[yIter]],
				Coefficient[metricCollected, d[xIter]*d[yIter]]
		]
	,{xIter, coordMod}, {yIter, coordMod}]
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
		
	];
Clear[GlobalGeometryID];

GlobalGeometryID[gdd_, coord_] := Hash[HoldComplete[{gdd, coord}]];

Clear[BuildGlobalBundle];

BuildGlobalBundle[gdd_, coord_, id_] := Module[
    {ds2, buildSymmetric2, Agdd, AgUU, Dim, gUU},

    ds2 = d[coord].gdd.d[coord];
    Dim = Length[coord];

    buildSymmetric2[X_] := Association[
        Table[{iIter, jIter} -> X[[iIter,jIter]], {iIter, Dim}, {jIter, iIter, Dim}]
    ];

    gUU = Inverse[gdd];

    Agdd = buildSymmetric2[gdd];
    AgUU = buildSymmetric2[gUU];

    Agdd = KeySelect[Agdd, Agdd[#] =!= 0 &];
    AgUU = KeySelect[AgUU, AgUU[#] =!= 0 &];

    globalBundle = <|
    	"GeometryID" -> id,
        "coord" -> coord,
        "ds2" -> ds2,
	"GlobalSync"-> {},
        "Tensors" -> <|
            "gdd" -> Agdd,
            "gUU" -> AgUU
        |>
    |>;
];


Clear[InitGlobalBundle];

InitGlobalBundle[gddIN_:"Global", coordIN_:"Global"] := Module[
    {gddint, coordint, id},

    coordint = ResolveGlobal[coordIN, coord];
    gddint   = ResolveGlobal[gddIN, gdd];

    id = GlobalGeometryID[gddint, coordint];

    If[(!AssociationQ[globalBundle]) || (Lookup[globalBundle, "GeometryID", None] =!= id),
    	Print["Initialize New Global Bundle"];
    	BuildGlobalBundle[gddint, coordint, id]
    ];

    globalBundle
];

Clear[SetGlobalTensor];
SetAttributes[SetGlobalTensor, HoldFirst];

SetGlobalTensor[symbol_, tensorName_] := Module[{},
    If[
        KeyExistsQ[globalBundle["Tensors"], tensorName] && Not[MemberQ[globalBundle["GlobalSync"], tensorName]],
        symbol = GetTensorArray[globalBundle, tensorName];
	AppendTo[globalBundle["GlobalSync"], tensorName];
    ];
];

Clear[SyncGlobalTensors];

SyncGlobalTensors[] := Module[{},
    SetGlobalTensor[ChrisUdd, "ChrisUdd"];
    SetGlobalTensor[Rdd, "Rdd"];
    SetGlobalTensor[RicciScalar, "RicciScalar"];
];

Clear[ComputeChrisUdd];

ComputeChrisUdd[simp_:Identity, gddcoord_:{"Global", "Global"}] := Module[
    {},
    InitGlobalBundle[First[gddcoord], Last[gddcoord]];
    PaiComputeBundleTensors[globalBundle, "ChrisUdd", simp];
    SyncGlobalTensors[];
];

Clear[ComputeRdd];

ComputeRdd[simp_:Identity, gddcoord_:{"Global", "Global"}] :=Module[
	{},
	InitGlobalBundle[First[gddcoord], Last[gddcoord]];
    	PaiComputeBundleTensors[globalBundle, "Rdd", simp];
	SyncGlobalTensors[];
];

ComputeRicciScalar[simp_:Identity, gddcoord_:{"Global", "Global"}] := Module[
	{},
	InitGlobalBundle[First[gddcoord], Last[gddcoord]];
    	PaiComputeBundleTensors[globalBundle, "RicciScalar", simp];
	SyncGlobalTensors[];
	];

"Here we consider the definition of the contraction operator Contracione that take a p-form in the vielbein basis an 
return a (p-1)-form with a Lorentz index attaced at the beggining."

SetAttributes[inP,Listable](*inP for Inned Product*)
inP[x_,0] = 0;
inP[x_,y_] := 0/;FormDegree[y]===0;
inP[x_Plus,y_] := inP[#,y]&/@x;
inP[x_,y_Plus] := inP[x,#]&/@y;
inP[x_,u_*y_] := u*inP[x,y]/;FormDegree[u]===0;
inP[x_.*e[a_],y_.*e[b_]] := x*y*KroneckerDelta[a,b];
inP[x_.*e[j_],y_.*HoldPattern[Wedge[e[k_],p__]]] := x*y*(KroneckerDelta[j,k]*Wedge[p]-Wedge[e[k],inP[e[j],Wedge[p]]])
Contractione[X_, DimIn_:Dim] := Table[inP[e[a1111],X],{a1111,DimIn}];

ClearAll[SetVielbein];
SetVielbein[eIN_,flatmetric_,simp_:Identity] := Module[{},
	ClearAll[\[Eta]dd,\[Eta]UU,eTodx,dxToe,eBasis,gdd,gUU,eamuUd,eamudU];
	\[Eta]dd=flatmetric;
	\[Eta]UU=Inverse[\[Eta]dd];
	eTodx=Table[e[iiinx]->eIN[[iiinx]],{iiinx,Dim}];
	dxToe=Solve[Table[eIN[[iiinx]]==e[iiinx],{iiinx,Dim}],Table[d[coord[[jjinx]]],{jjinx,Dim}]]//Last;
	eBasis=Array[e,{Dim}];
	Do[d[eBasis[[iiinx]]]=simp[(d[eIN]/.dxToe)][[iiinx]],{iiinx,Dim}];
	eamuUd=Map[FormsToMatrix[#, 1, coord]&, eIN];
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
(*====== Build Print definition of the objects ======*)

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
	
(*==========================================================================================================================================*)

(*
|------------------------------------------------------
|     Thinking on Association for saving tensor
|------------------------------------------------------

                     --- Utils ---
*)

TensorProductContract[Tensors__, contractIndices_List] := Activate@TensorContract[Inactive[TensorProduct][Tensors], contractIndices];

RaiseIndices[TensorSparsedown_, bundle_, indicesRaisePosition_] := 
	Module[{gUU, RaisePositions, rank, RaiseRelations, metricSequence, TensorUpPermuted, indexPermutation},
		RaisePositions = Sort[indicesRaisePosition];
		rank = Length[Dimensions[TensorSparsedown]];
		RaiseRelations = Table[{RaisePositions[[n]], rank + 2*n-1}, {n, 1, Length[RaisePositions]}];
		gUU = SparseArray[GetTensorArray[bundle, "gUU"]];
		metricSequence = Sequence @@ ConstantArray[gUU, Length[RaisePositions]];
		TensorUpPermuted = TensorProductContract[TensorSparsedown, metricSequence,RaiseRelations];
		indexPermutation = Join[Complement[Range[rank], RaisePositions], RaisePositions];
		Return[Transpose[TensorUpPermuted, indexPermutation]];
	];

(*
				---- Covariant derivative ----
*)

Clear[PaiCovD];
ClearAll[ComputeTermCovD];

ComputeTermCovD[ChrUdd_, Tensor_, pos_, NLegs_, Uord_] := Module[
	{UnorderedProduct, toTranslate, term},
	Which[
	Uord === "d",
		UnorderedProduct = TensorProductContract[ChrUdd, Tensor, {{1, 3 + pos}}];
		toTranslate = DeleteCases[Range[2, NLegs + 1], 1 + pos];
		term = -Transpose[UnorderedProduct, {1, 1 + pos, Sequence @@ toTranslate}],
	Uord === "U",
		UnorderedProduct = TensorProductContract[ChrUdd, Tensor, {{3, 3 + pos}}];
		toTranslate = DeleteCases[Range[1, NLegs + 1], 1 + pos];
		term = Transpose[UnorderedProduct, {1 + pos, Sequence @@ toTranslate}],
	True,
		Print["[Aborting] character << " <> ToString[Uord] <> " >> is not U or d"];
		Abort[];
	];
	Return[term]
];

SetAttributes[PaiCovD, HoldFirst];
PaiCovD[bundle_, tensor_, indices_String] :=
	Module[{indicesSplit, terms, chrUdd, nLegs, coord},
		indicesSplit = Characters[indices];
		nLegs = Length[indicesSplit];
		chrUdd = GetTensorArray[bundle, "ChrisUdd"];
		terms = Table[
			ComputeTermCovD[chrUdd, tensor, pos, nLegs, indicesSplit[[pos]]]
		,{pos, nLegs}];
		coord = bundle["coord"];
		Return[Table[D[tensor, xIter], {xIter, coord}] + Total[terms]];
	];

SetAttributes[GetTensorArray, HoldFirst];
GetTensorArray[bundle_, tensorName_, simp_:Identity] := 
	Module[{PaiTensor, TensorComponents, Dim, TensorArray, DimensionsTensor},
		Dim = Length[bundle["coord"]];
		PaiComputeBundleTensors[bundle, tensorName, simp];
		
		PaiTensor = bundle["Tensors", tensorName];
		
		If[
		tensorName === "RicciScalar",
			Return[PaiTensor]
		];
		
		TensorComponents[indices__] := PaiComponent[PaiTensor, {indices}, tensorName];
		DimensionsTensor = ConstantArray[Dim, tensorRank[tensorName]];
		TensorArray = Array[TensorComponents, DimensionsTensor];
		Return[TensorArray];
	];

CleanZeros[X_Association] := KeySelect[X, X[#] =!= 0 &];
PaiComponent::unk = "Tensor `1` desconocido."
PaiComponent[PaiTensor_, {indices__}, tensorName_] := 
	Which[
	(tensorName === "Rdd")||(tensorName === "gdd")||(tensorName === "gUU"),
		PaiComponent2sym[PaiTensor, {indices}],
	tensorName === "ChrisUdd",
		PaiComponent3symLast[PaiTensor, {indices}],
	tensorName === "Rdddd",
		PaiComponent4Riem[PaiTensor, {indices}],
	True,
		Message[PaiComponent::unk, tensorName]
	];

tensorRank = <|
	"gdd" -> 2, "gUU" -> 2, "Rdd" -> 2,
	"ChrisUdd" -> 3,
	"Rdddd" -> 4
	|>;

PaiComponent2sym[Xab_Association, {i_, j_}] := 
	Module[{m,n},
		{m,n} = If[i <= j, {i,j}, {j,i}];
		First[Lookup[Xab, {{m,n}}, 0]]
	];

PaiComponent3symLast[Xabc_Association, {p_, i_, j_}] := 
	Module[{q,m,n},
		{q,m,n} = If[i <= j,{p,i,j},{p,j,i}];
		First[Lookup[Xabc, {{q,m,n}}, 0]]
	];

PaiComponent4Riem[Xabcd_Association, {i_, j_, k_, l_}] := 
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

PaiComputeMetric[bundle_Association] := 
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
		
		Return[<|"gdd" -> Agdd, "gUU" -> AgUU|>];
	];

(*
				---- ChrisUdd ----
*)

PaiComputeChrisUdd[bundle_Association] := Module[
	{coord, Dim, gdd, Paigdd, gUU, AgUU,
	AChrisUdd, dgdd, dGamdd, GamUdd, i, j, k, l},

	coord = bundle["coord"];
	Dim = Length[coord];
	Paigdd = bundle["Tensors","gdd"];
	AgUU = bundle["Tensors","gUU"];
	gdd[i_,j_] := PaiComponent[Paigdd, {i,j}, "gdd"];
	gUU[i_,j_] := PaiComponent[AgUU, {i,j}, "gUU"];
	
	dgdd[i_,j_,k_] := dgdd[i,j,k] = D[gdd[j,k], coord[[i]]];
	dGamdd[k_,i_,j_] := dGamdd[k,i,j] = 1/2*dgdd[i,j,k] + 1/2*dgdd[j,i,k] - 1/2*dgdd[k,i,j];
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

Clear[PaiComputeRdddd];

PaiComputeRdddd[bundle_Association] :=
Module[
    {coord, Dim, Agdd, AChrisUdd, gdd, ChrisUdd, dChrisUdd,
    RiemUddd, Riemdddd, pairs, nPairs, rules, p, q, e},

    coord = bundle["coord"];
    Dim = Length[coord];

    Agdd = bundle["Tensors", "gdd"];
    AChrisUdd = bundle["Tensors", "ChrisUdd"];

    gdd[i_, j_] :=
        PaiComponent[Agdd, {i, j}, "gdd"];

    ChrisUdd[a_, b_, c_] :=
        PaiComponent[AChrisUdd, {a, b, c}, "ChrisUdd"];

    dChrisUdd[a_, b_, c_, mu_] := dChrisUdd[a, b, c, mu] = D[ChrisUdd[a, b, c], coord[[mu]]];

    RiemUddd[a_, b_, c_, d_] := RiemUddd[a, b, c, d] = (dChrisUdd[a, b, d, c] - dChrisUdd[a, b, c, d] + 
    	Sum[ChrisUdd[a, c, e] ChrisUdd[e, b, d] - ChrisUdd[a, d, e] ChrisUdd[e, b, c], {e, Dim}]);

    Riemdddd[a_, b_, c_, d_] := Riemdddd[a, b, c, d] = Sum[gdd[a, e] RiemUddd[e, b, c, d], {e, Dim}];

    pairs = Subsets[Range[Dim], {2}];
    nPairs = Length[pairs];

    rules = Table[
                With[
                    {ab = pairs[[p]], cd = pairs[[q]]},
                    Join[ab, cd] -> Apply[Riemdddd, Join[ab, cd]]
                ],
                {p, nPairs},
                {q, p, nPairs}
            ];
    rules = Flatten[rules, 1];

    CleanZeros[Association[rules]]
]

(*
				---- Ricdd ----
*)

	
PaiComputeRdd[bundle_Association] := 
	Module[{coord,Dim,AgUU,ARiemdddd,gUU,Rdddd,Ricdd,ARicdd,gUUArray,RddddArray },
		
		If[
			failRequirements[bundle["Tensors"], {"Rdddd"}],
				Return["Missing Tensors/Rdddd"]
		];
		
		coord = bundle["coord"];
		Dim = Length[coord];
		AgUU = bundle["Tensors","gUU"];
		ARiemdddd = bundle["Tensors","Rdddd"];
		gUU[i_,j_] := PaiComponent[AgUU, {i,j}, "gUU"];
		Rdddd[i_,j_,k_,l_] := PaiComponent[ARiemdddd, {i, j, k, l}, "Rdddd"];
		gUUArray = SparseArray[Array[gUU,{Dim,Dim}]];
		RddddArray = SparseArray[Array[Rdddd,{Dim,Dim,Dim,Dim}]];
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

PaiComputeRicciScalar[bundle_Association] := 
	Module[{coord,Dim,AgUU,ARiemdddd,gUU,Rdddd,Ricdd,ARicdd,RicciScalar},
		
		If[
			failRequirements[bundle["Tensors"], {"Rdd"}],
				Return["Missing Tensors/Rdd"]
		];
		
		coord = bundle["coord"];
		Dim = Length[coord];
		AgUU = bundle["Tensors","gUU"];
		ARicdd = bundle["Tensors","Rdd"];
		gUU[i_,j_] := PaiComponent[AgUU, {i, j}, "gUU"];
		Ricdd[i_,j_] := PaiComponent[ARicdd,{i, j}, "Rdd"];
		RicciScalar = Sum[gUU[i,j]*Ricdd[i,j],{i,Dim},{j,Dim}];
		Return[RicciScalar];
	];
	

(*
		    	------------------------------------
			---       Orchestra Director     ---
			------------------------------------
*)

PaiComputeBundleTensors::usage = 
	"
	PaiComputeBundleTensors[bundle, level, simp_:]
	levels  ->    {metric, ChrisUdd, Rdddd, Rdd, RicciScalar}  
	Default -> Rdddd
	Usage   -> Mutate bundle according to level
	"

ClearAll[PaiComputeBundleTensors];
SetAttributes[PaiComputeBundleTensors, HoldFirst];
PaiComputeBundleTensors[bundleIN_, level_: "RicciScalar", simp_:Identity] := Module[{}, 
	Which[
	KeyExistsQ[bundleIN, "ds2"],
		PaiComputeBundleTensorsMetric[bundleIN, level, simp],
	KeyExistsQ[bundleIN, "eU"],
		PaiComputeBundleTensorsVielbein[bundleIN, level, simp],
	True,
		Print["[ Aborting ] netierh ds2 or eU was provided"];
		Abort[];
	];
];

ClearAll[PaiComputeBundleTensorsVielbein];
SetAttributes[PaiComputeBundleTensorsVielbein, HoldFirst];

PaiComputeBundleTensorsVielbein[bundleIN_, level_: "RicciScalar", simp_:Identity] := 
Module[{bundle, needSpinConnection, needCurvature, needRdddd, needRdd, needRicciScalar},
	bundle = bundleIN;

	InitVielbeinBundle[bundleIN, simp];
	bundle=bundleIN;

	If[level === "basic", Return[]];

	needSpinConnection   = Not[KeyExistsQ[bundle, "spinConnection"]];
	needCurvature  = Not[KeyExistsQ[bundle, "curvatureForm"]];
	FlatTensors = Lookup[bundle,"FlatTensor", <| |>];
	needRdddd  = Not[KeyExistsQ[FlatTensors, "Rdddd"]];
	needRdd = Not[KeyExistsQ[FlatTensors, "Rdd"]];
	needRicciScalar = Not[KeyExistsQ[FlatTensors, "RicciScalar"]];
	
	If[needSpinConnection,
		Print["** Computing spin connection"];
		PaiComputeSpinConnection[bundle];
		bundleIN = bundle;
	];

	If[level === "spinConnection", Return[]];
	
	If[needCurvature,
		Print["** Computing curvature 2-form"];
		PaiComputeCurvatureForm[bundle];
		bundleIN = bundle;
	];

	If[level === "curvatureForm", Return[]];

	If[needRdddd,
		Print["** computing Rdddd flat"];
		PaiComputeRddddFlat[bundle];
		bundleIN = bundle;
	];

	If[level === "Rdddd", Return[]];

	If[needRdd,
		Print["** computing Rdd flat"];
		PaiComputeRddFlat[bundle];
		bundleIN = bundle;
	];

	If[level === "Rdd", Return[]];

	If[needRicciScalar,
		Print["** computing RicciScalar"];
		PaiComputeRicciScalarFlat[bundle];
		bundleIN = bundle;
	];
];


ClearAll[PaiComputeBundleTensorsMetric];
SetAttributes[PaiComputeBundleTensorsMetric, HoldFirst];
PaiComputeBundleTensorsMetric[bundleIN_, level_: "Rdddd", simp_:Identity] := 
Module[{Tensors, needMetric, needChris, needRiemann, AgddgUU, Agdd, AgUU, 
	AChrisUdd, ARiemdddd, bundle, ARicdd,needRicci, needRicciScalar, RicciScalar},
	
	bundle = bundleIN;
	Tensors = Lookup[bundle, "Tensors", <||>];

	needMetric   = Not[KeyExistsQ[Tensors, "gdd"]] || Not[KeyExistsQ[Tensors, "gUU"]];
	needChris    = Not[KeyExistsQ[Tensors, "ChrisUdd"]] && MemberQ[{"RicciScalar","ChrisUdd", "Rdddd", "Rdd"}, level];
	needRiemann  = Not[KeyExistsQ[Tensors, "Rdddd"]] && MemberQ[{"RicciScalar","Rdddd", "Rdd"}, level];
	needRicci    = Not[KeyExistsQ[Tensors, "Rdd"]] && MemberQ[{"RicciScalar","Rdd"}, level];
	needRicciScalar = Not[KeyExistsQ[Tensors, "RicciScalar"]] && MemberQ[{"RicciScalar"}, level];
	
	(* --- Metric --- *)
	If[needMetric,
		Print["** Computing metric"];
		AgddgUU = PaiComputeMetric[bundle];
		Agdd = Map[simp, AgddgUU["gdd"]];
		AgUU = Map[simp, AgddgUU["gUU"]];
		Tensors = Join[Tensors, <|"gdd" -> Agdd, "gUU" -> AgUU|>];
		bundle = AssociateTo[bundle, "Tensors" -> Tensors];
		InitMetricBundle[bundle, simp];
	];
	
	If[level === "metric", 
		bundleIN = bundle;
		Return[]
	];
	
	(* --- Christoffel --- *)
	If[needChris,
		Print["** Computing Christoffel"];
		AChrisUdd = PaiComputeChrisUdd[bundle];
		AChrisUdd = Map[simp, AChrisUdd];
		Tensors = AssociateTo[Tensors, "ChrisUdd" -> AChrisUdd];
		bundle = AssociateTo[bundle, "Tensors" -> Tensors];
	];

	If[level === "ChrisUdd",
		bundleIN = bundle;
		Return[]
	];

	(* --- Riemann --- *)
	If[needRiemann,
		Print["** Computing Riemann"];
		ARiemdddd = PaiComputeRdddd[bundle];
		ARiemdddd = Map[simp, ARiemdddd];
		Tensors = AssociateTo[Tensors, "Rdddd" -> ARiemdddd];
		bundle = AssociateTo[bundle, "Tensors" -> Tensors];
	];
	
	If[level === "Rdddd", 
		bundleIN = bundle;
		Return[]
	];
	
	(* --- Ricci --- *)
	If[needRicci,
		Print["** Computing Ricci tensor"];
		ARicdd = PaiComputeRdd[bundle];
		ARicdd = Map[simp, ARicdd];
		Tensors = AssociateTo[Tensors, "Rdd" -> ARicdd];
		bundle = AssociateTo[bundle, "Tensors" -> Tensors];
	];
	If[level === "Rdd", 
		bundleIN = bundle;
		Return[]
	];
	
	(* --- RicciScalar --- *)
	If[needRicciScalar,
		Print["** Computing RicciScalar"];
		RicciScalar = simp[PaiComputeRicciScalar[bundle]];
		Tensors = AssociateTo[Tensors, "RicciScalar" -> RicciScalar];
		bundle = AssociateTo[bundle, "Tensors" -> Tensors];
	];
	
	If[level === "RicciScalar",
		bundleIN = bundle;
		Return[]
	];

];

(*
    ---- "Conneting with previous notation" ----
*)

(*
    ---- Hodge Builders ----
*)

Clear[BuildHodge];
SetAttributes[BuildHodge, HoldFirst];
BuildHodge[bundle_, simp_:Simplify] := Module[{},
	Which[
	KeyExistsQ[bundle, "ds2"],
		Return[BuildHodgeMetric[bundle, simp]],
	KeyExistsQ[bundle, "eU"],
		Return[BuildHodgeVielbein[bundle, simp]],
	True,
		Print["[ Aborting ] BuildHodge: ds2 or eU are not given."];
	];
];

SetAttributes[BuildHodgeMetric, HoldFirst];
BuildHodgeMetric[bundle_, simp_:Simplify] := 
	Module[{gdd,sqrtdetg,needMetric,Tensors,Dim,gddMap,gUU,gUUMap,coord},
		coord = bundle["coord"];
		Dim = Length[coord];
		Tensors = Lookup[bundle, "Tensors", <||>];
	
		needMetric   = Not[KeyExistsQ[Tensors, "gdd"]] || Not[KeyExistsQ[Tensors, "gUU"]];
		If[needMetric,
				PaiComputeBundleTensors[bundle, "metric", simp]
		];
		gdd = GetTensorArray[bundle, "gdd"];
		gUU = GetTensorArray[bundle, "gUU"];
		sqrtdetg = simp[Sqrt[-Det[gdd]]];
		If[KeyExistsQ[bundle, "assum"],
				sqrtdetg = Simplify[sqrtdetg, bundle["assum"]],
					sqrtdetg = Simplify[sqrtdetg]
		];
		Return[Function[{X}, 
			Hstar[X, gUU, sqrtdetg, d[coord], simp]
			]];
	];

BuildHodgeVielbein[bundle_, simp_:Simplify] := Module[{basis, eta, etainv, sqrtdeteta},
	basis = bundle["basis"];
	eta = DiagonalMatrix[bundle["signature"]]; 
	etainv = Inverse[eta];
	sqrtdeteta = Sqrt[-Det[eta]];
	Return[Function[{X}, 
			Hstar[X /. bundle["dxToe"], etainv, sqrtdeteta, basis, simp]
			]];


]

(*
	"Vielbein Computations"
*)

Clear[ConstructContraction];
ConstructContraction[vielbeinBundle_] := Module[
	{symbs, Contraction, InnerExtractor},
	symbs = vielbeinBundle["basis"];
	SetAttributes[InnerExtractor, Listable];
	InnerExtractor[X_,eIter_]:=Extractor[X, eIter, "left"];
	Contraction[X_] := Table[InnerExtractor[X, eIter], {eIter, symbs}];
	Return[Function[{X}, Contraction[X]]]
	];

SetAttributes[InitVielbeinBundle, HoldFirst];

InitVielbeinBundle[vielbeinBundle_, simp_:Simplify] := Module[
	{eTodx, dxToe, symbs, eU, contraction, coordbasis, hstar},
	vielbeinBundle = Association[vielbeinBundle];
	symbs = vielbeinBundle["basis"];
	eU = vielbeinBundle["eU"];
	coordbasis = d[vielbeinBundle["coord"]];
	Do[FormDegree[eIter] = 1, {eIter, symbs}];
	eTodx = Normal[AssociationThread[symbs -> eU]];
	dxToe = Solve[eU == symbs, coordbasis][[1]];
	contraction = ConstructContraction[vielbeinBundle];

	deU = d[symbs] /. eTodx /. dxToe;
	dictde = AssociationThread[d[symbs], deU];
	Do[d[eIter] = Collect[dictde[d[eIter]], _Wedge, simp],{eIter, symbs}];
	AssociateTo[vielbeinBundle, {"eTodx" -> eTodx, "dxToe" -> dxToe, "contraction"->contraction, "UseVielbein" -> True}];
	hstar = BuildHodge[vielbeinBundle, simp];
	vielbeinBundle["Hstar"] = hstar;
	FormSquareTools = BuildSquaresTools[vielbeinBundle, simp];
	vielbeinBundle["FormSquare"] = FormSquareTools["FormSquare"];
	vielbeinBundle["FormSquaredd"] = FormSquareTools["FormSquaredd"];

	];

SetAttributes[InitMetricBundle, HoldFirst];

InitMetricBundle[bundle_, simp_:Identity] := Module[{ToClear, FormSquareTools},
    bundle = Association[bundle];
	If[KeyExistsQ[bundle, "constants"],
    	Do[d[cIter]=0, {cIter, bundle["constants"]}];
	];
    bundle["Hstar"] = BuildHodge[bundle, simp];
    FormSquareTools = BuildSquaresTools[bundle, simp];
    bundle["FormSquare"] = FormSquareTools["FormSquare"];
    bundle["FormSquaredd"] = FormSquareTools["FormSquaredd"];
];

Clear[PaiComputeSpinConnection];
SetAttributes[PaiComputeSpinConnection, HoldFirst];
PaiComputeSpinConnection[frameBundle_, simp_:Identity] := Module[
	{eta, vielbein, deU, contraction, idedU, idedd, iideddU, iideddd, omegadd, symbs, Dim, evald, omegaUd},
	frameBundle = Association[frameBundle];
	eta = DiagonalMatrix[frameBundle["signature"]];
	Dim = Length[frameBundle["signature"]];
	vielbein = frameBundle;
	symbs = vielbein["basis"];
	deU = Table[d[e], {e, symbs}];
	deU = simp[deU];
	contraction = vielbein["contraction"];
	idedU = simp[contraction[deU]];
	idedd = idedU . eta;
	iideddU = contraction[idedU];
	iideddd = iideddU . eta;
	omegadd = 1/2*iideddd . symbs - 1/2*idedd + 1/2*Transpose[idedd];
	AssociateTo[frameBundle, "spinConnection" -> <|"dd"->omegadd|>];
	];
	
Clear[PaiComputeCurvatureForm];
SetAttributes[PaiComputeCurvatureForm, HoldFirst];
PaiComputeCurvatureForm[frameBundle_, simp_: Identity] := Module[
	{eta, omegaUd, vielbein, RRUd, omegadd},
	eta = DiagonalMatrix[frameBundle["signature"]];
	omegadd = frameBundle["spinConnection", "dd"];
	omegaUd = eta . omegadd;
	vielbein = frameBundle;
	RRUd = (d[omegaUd] + omegaUd\[Wedge]omegaUd)/.vielbein["dxToe"];
	RRUd = simp[RRUd];
	AssociateTo[frameBundle, "curvatureForm" -> <|"Ud" -> RRUd|>];
	];

Clear[PaiComputeRddddFlat];
SetAttributes[PaiComputeRddddFlat, HoldFirst];
PaiComputeRddddFlat[frameBundle_, simp_: Identity] := Module[
	{eta, RUd, contraction, RddUd, RUddd, Rdddd},
	eta = DiagonalMatrix[frameBundle["signature"]];
	RUd = frameBundle["curvatureForm", "Ud"];
	contraction = frameBundle["contraction"];
	RddUd = -contraction[contraction[RUd]];
	RUddd = Transpose[RddUd, {3,4,1,2}];
	Rdddd = eta . RUddd;
	AssociateTo[frameBundle, "FlatTensors" -> <|"Rdddd"-> Rdddd|>]
	];

Clear[PaiComputeRddFlat];
SetAttributes[PaiComputeRddFlat, HoldFirst];
PaiComputeRddFlat[frameBundle_, simp_: Identity] := Module[
	{eta, etainv, Rdddd, RUddd, Rdd, FlatTensors},
	eta = DiagonalMatrix[frameBundle["signature"]];
	etainv = Inverse[eta];
	Rdddd = frameBundle["FlatTensors", "Rdddd"];
	RUddd = etainv . Rdddd;
	Rdd = TensorContract[RUddd, {{1, 3}}];
	AssociateTo[frameBundle["FlatTensors"], <|"Rdd"-> Rdd|>]
	];

Clear[PaiComputeRicciScalarFlat];
SetAttributes[PaiComputeRicciScalarFlat, HoldFirst];
PaiComputeRicciScalarFlat[frameBundle_, simp_: Identity] := Module[
	{eta, etainv, Rdd, RUd, RicciScalar, FlatTensors},
	eta = DiagonalMatrix[frameBundle["signature"]];
	etainv = Inverse[eta];
	Rdd = frameBundle["FlatTensors", "Rdd"];
	RUd = etainv . Rdd;
	RicciScalar = Tr[RUd];
	AssociateTo[frameBundle["FlatTensors"], <|"RicciScalar"-> RicciScalar|>]
	];

End[]

EndPackage[]
