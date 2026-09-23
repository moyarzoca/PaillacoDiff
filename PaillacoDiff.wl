(* ::Package:: *)

BeginPackage["PaillacoDiff`"]

(* ---------- Public functions ---------- *)

FormDegree::usage = "FormDegree[expr] returns the degree of a differential form (0 for scalars)."
Wedge::usage = "Wedge[x, y, ...] is the exterior (wedge) product of forms."
d::usage = "d[expr] is the exterior derivative."
PolyFormQ::usage = "PolyFormQ[expr] tests whether expr is a sum of forms of different degrees."

Extractor::usage = "Extractor[F, A] extracts the coefficient of 1-form A in polyform F."
Extractorleft::usage = "Extractorleft[F, A] extracts A from the left side of each term."

DNAofForm::usage = "DNAofForm[X] decomposes form X into {{coeff, indices}, ...}."
SparseFromDNA::usage = "SparseFromDNA[DNA, dim, deg] converts DNA to a SparseArray."
DNAFromSparse::usage = "DNAFromSparse[sparse] converts a SparseArray back to DNA."
FormSquare::usage = "FormSquare[bundle][X] / FormSquare[X] computes F_{mu1...mup} F^{mu1...mup}."
FormSquaredd::usage = "FormSquaredd[bundle][X] / FormSquaredd[X] computes F_{mu l2...lp} F_nu^{ l2...lp}."
Hstar::usage = "Hstar[bundle][X] / Hstar[X] computes the Hodge dual of form X."
Contraction::usage = "Contraction[bundle][X] / Contraction[X] computes the contraction operation of X in the basis d[coord] for metric vielbein and bundle[\"basis\"] for vielbein mode."
FormToSparse::usage = "FormToSparse[X, deg, coord] converts a form to a SparseArray."
FormToMatrix::usage = "FormToMatrix[X, deg, coord] converts a form X to a dense matrix. deg is an Integer and the degree of the form, and coord are the coordinates."

ClearGeometric::usage = "ClearGeometric[] clears global tensors ChrisUdd, Rdd, RicciScalar."
DiffToMatrix::usage = "DiffToMatrix[ds2, coord] extracts the metric tensor from a line element."

TensorProductContract::usage = "TensorProductContract[t1, t2, ..., {{i1,j1}, ...}] contracts tensor products."
RaiseIndices::usage = "RaiseIndices[sparse, bundle, positions] raises specified indices."
LowerIndices::usage = "LowerIndices[sparse, bundle, positions] lower specified indices."
PaiCovD::usage = "PaiCovD[bundle, tensor, indices] computes the coordinate-basis covariant derivative of tensor. indices is a string of U/d characters describing tensor index variance. For instace for  tensor TUdU indices must be the string UdU. The covariant derivative index is added at the beginning of the tensor"
GetTensorArray::usage = "GetTensorArray[bundle, name] retrieves a tensor array, computing on demand."

PaiDef::usage = "PaiDef[\"T{indices}:=expression\"] defines a tensor using GRTensor-like notation.

Tensor indices are written inside braces, with '^' denoting an upper index.
Covariant derivatives are written after ';'. For example,

    PaiDef[\"H{a b}:=18*R{a ^c}*R{b c}*Ricciscalar\"]

defines H_ab, while

    R{a b ;^c ;c}

denotes a covariantly differentiated Ricci tensor.

PaiDef stores the definition symbolically and does not compute tensor components.
The current implementation supports monomial tensor expressions without
parenthesized sums.";


PaiCalc::usage = "PaiCalc[bundle][\"T(indices)\"] computes the tensor in bundle."

PaiComponents::usage = "PaiComponents[bundle][\"T(indices)\"] returns the components of the tensor previously computed."
PaiCompute::usage =
"PaiCompute[bundle][\"T(indices)\"] computes the components of a tensor
previously defined with PaiDef.
The string \"indices\" specifies the requested index positions using
'dn' for lower indices and 'up' for upper indices. For example,
            PaiCompute[bundle][\"H(dn,dn)\"]"
Paillaco::usage =
"Paillaco[bundle][\"T(indices)\"] computes a tensor if necessary and returns its components.";

(* ---------- Public globals ---------- *)

PaiNonCommutativeScalarQ::usage = "PaiNonCommutativeScalarQ[expr] tests whether expr contains a registered noncommutative scalar coefficient.";
PaiRegisterNonCommutativeScalarQ::usage = "PaiRegisterNonCommutativeScalarQ[test] registers a predicate test[expr] used by Wedge to detect noncommutative scalar coefficients.";

PaiSimplify::usage = "PaiSimplify[expr] applies PaillacoDiff's default lightweight algebraic simplification.";
$UsePaiSimplify::usage = "$UsePaiSimplify controls whether PaiSimplify applies automatic simplification. Default is True.";

coord::usage = "List of coordinate variables."
Dim::usage = "Spacetime dimension."
ds2::usage = "Metric ds2 expresed in the coordinate basis d[xmu]*d[xnu]"
gdd::usage = "Metric tensor g_{mu nu}."
gUU::usage = "Inverse metric g^{mu nu}."
ChrisUdd::usage = "Christoffel symbols Gamma^mu_{nu rho}."
Rdddd::usage = "Riemann tensor R_{mu nu rho sigma}."
Rdd::usage = "Ricci tensor R_{mu nu}."
RicciScalar::usage = "Ricci scalar R."
sqrtdetg::usage = "Sqrt[-det(g)]."

eTodx::usage = "Rule mapping e^a to e^a_mu dx^mu."
dxToe::usage = "Rule mapping dx^mu to e^a."
eamuUd::usage = "Vielbein matrix e^a_mu."
eamudU::usage = "Inverse vielbein matrix e_a^mu."

Begin["`Private`"]

globalBundle = <| |>;

ClearAll[GlobalRequired];
SetAttributes[GlobalRequired, HoldAll];

GlobalRequired::missing = "Global Mode requires `1` to be defined.";

GlobalRequired[x_Symbol] :=
	If[!ValueQ[x],
		Message[GlobalRequired::missing, HoldForm[x]];
		Abort[],
            True
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

Clear[MetricDifferentialDegreeQ];

MetricQuadraticInDiffQ[ds2_] := Module[
    {lambda, scaled},
    scaled = Expand[ds2 /. d[_] :> lambda];
    (Exponent[scaled, lambda, Min] === 2) && (Exponent[scaled, lambda, Max] === 2)
];
Clear[ValidateMetricBundle];

ValidateMetricBundle[bundle_] := Module[
    {coord, metricCoord},

    If[!KeyExistsQ[bundle, "coord"],
        Print["[ Aborting ] Bundle requires key \"coord\""];
        Abort[]
    ];

    If[KeyExistsQ[bundle, "ds2"],
        coord = bundle["coord"];
        metricCoord = DeleteDuplicates[Cases[bundle["ds2"], d[x_] :> x, Infinity]];

        If[Sort[coord] =!= Sort[metricCoord],
            Print[
                "[ Aborting ] Coordinate mismatch",
                "\nIn coord but not in metric: ", Complement[coord, metricCoord],
                "\nIn metric but not in coord: ", Complement[metricCoord, coord]
            ];
            Abort[]
        ];
    ];

    If[Not[MetricQuadraticInDiffQ[bundle["ds2"]]],
        Print["[ Aborting ] Metric must be quadratic in your coordinates differentials ", d[coord]];
        Abort[]
    ];

    True
];

Clear[ValidateForm];

ValidateForm[bundle_][X_] := Module[
    {allowedDifferentials, differentialForms, invalid},

    allowedDifferentials = d[bundle["coord"]];

    differentialForms =DeleteDuplicates[Cases[X, d[some_] :> d[some], {0, Infinity}]];

    invalid = Complement[differentialForms, allowedDifferentials];

    If[
        invalid =!= {},
        Print["[ Aborting ] Forms outside the bundle basis: ", invalid,
        "\nDid you forget to declare a constant?"];
        Abort[]
    ];

    True
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

DNAofForm::noBaseFound = "No base element found";
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
BuildSquaresTools[bundle_, simp_:PaiSimplify] := Module[{gUU, eta, etainv, basis, dxToe},
	Which[

	VielbeinBundleQ[bundle],
		eta = GetFlatMetric[bundle];
		etainv = Inverse[eta];
		basis = bundle["basis"];
		dxToe = bundle["dxToe"];
		Return[
		    <| "FormSquare" -> Function[{X}, FormSquareCore[X /. dxToe, etainv, simp,  basis]],
		    "FormSquaredd" -> Function[{X}, FormSquareddCore[X /. dxToe, etainv, simp,  basis]]
		    |>
		],

	KeyExistsQ[bundle, "ds2"],
		PaiComputeBundleTensors[bundle, "metric", simp];
		gUU = GetTensorArray[bundle, "gUU"];
		basis = d[bundle["coord"]];
		Return[
		   <| "FormSquare" -> Function[{X}, FormSquareCore[X, gUU, simp,  basis]],
		      "FormSquaredd" -> Function[{X}, FormSquareddCore[X, gUU, simp,  basis]]
		    |>
		],

	True,
		Print["[ Aborting ] BuildSquaresTools: ds2 nor eU not given"];
		Abort[];
	];
];


Clear[FormSquareCore];
Clear[FormSquare];
FormSquareCore[Xform_, gUU_, simp_:PaiSimplify, basis_] :=
Module[{deg,FformDNA, FformSparse,gintUU,listindices,seqgUU, FtensorSparse,
    indicesContract,FormComps,TensorComps,InterComps,FformRule,FtensorRule,
    FformValues, FtensorValues,Dim},

    If[Xform ===0, Return[0]];

    gintUU = SparseArray[gUU];

    deg = FormDegree[Xform];
    Dim = Length[basis];
	FformDNA = simp[DNAofForm[Xform, basis]];
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

SetAttributes[FormSquare, HoldFirst];

FormSquare[bundle_][X_] /; AssociationQ[bundle] := Module[{},
    If[!KeyExistsQ[bundle, "FormSquare"],
        PaiComputeBundleTensors[bundle, "basicTools"]
    ];
    ValidateForm[bundle][X];
    bundle["FormSquare"][X]
];

FormSquare[X_] /; !AssociationQ[X] := Module[{},
    If[!KeyExistsQ[globalBundle, "FormSquare"],
        InitGlobalBundle[];
        PaiComputeBundleTensors[globalBundle, "basicTools"]
    ];
    ValidateForm[globalBundle][X];
    globalBundle["FormSquare"][X]
];


Clear[FormSquareddCore];
FormSquareddCore[0,__]:=0

FormSquareddCore[Xform_, gUU_, simp_:PaiSimplify, basis_] :=
Module[{deg,FformDNA,FformSparse,gintUU,
	Dim,seqgUU,indexcontr, nonzeroUp, nonzeroDn, nonzeroInter,
	nonzeroInterUp, nonzeroInterDn,FformRule,FtensorRule,nonzeroXd,nonzeroXdU,Xsqdd,
	Xdmunu, XdUmunu, FtensorSparse},
	
	gintUU = SparseArray[gUU];

	Dim = Length[basis];
	deg = FormDegree[Xform];
	FformDNA = simp[DNAofForm[Xform, basis]];
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

Clear[FormSquaredd];
SetAttributes[FormSquaredd, HoldFirst];

FormSquaredd[bundle_][X_] /; AssociationQ[bundle] := Module[{},
    If[!KeyExistsQ[bundle, "FormSquaredd"],
        PaiComputeBundleTensors[bundle, "basicTools"]
    ];
    ValidateForm[bundle][X];
    bundle["FormSquaredd"][X]
];

FormSquaredd[X_] /; !AssociationQ[X] := Module[{},
    If[!KeyExistsQ[globalBundle, "FormSquaredd"],
        InitGlobalBundle[];
        PaiComputeBundleTensors[globalBundle, "basicTools"]
    ];
    ValidateForm[globalBundle][X];
    globalBundle["FormSquaredd"][X]
];

Clear[Hstar, HstarCore];
HstarCore[Xform_, gintUUIN_, sqrtdetg_, base_, simp_:PaiSimplify] := 
	Module[{gintUU, coordint, Dim, deg, FformDNA, FformSparse, FtensorSparse,
		TensorComps,FtensorRule,FtensorValues, FtensorDict, compToStar,starF},

		If[
		Xform ===0,
			Return[0]
		];

		gintUU = SparseArray[gintUUIN];
		
		Dim = Length[base];
		deg = FormDegree[Xform];
		
		If[
			deg===0,
				Return[Xform*sqrtdetg Wedge@@(base)]
		];
		
		FformDNA = simp[DNAofForm[Xform, base]];
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
				<|"eps"->toepsilon, "basis" -> Map[base[[#]]&, complement]|>
			];
		
		starF = 
			sqrtdetg*Sum[
				FtensorDict[comp]*Signature[compToStar[comp]["eps"]]*Apply[Wedge, compToStar[comp]["basis"]]
			,
			{comp, TensorComps}];
		Return[starF]
    ];

SetAttributes[Hstar, HoldFirst];

Hstar[bundle_][X_] /; AssociationQ[bundle] := Module[{},
    If[!KeyExistsQ[bundle, "Hstar"],
        PaiComputeBundleTensors[bundle, "basicTools"]
    ];
    ValidateForm[bundle][X];
    bundle["Hstar"][X]
];

Hstar[X_] /; !AssociationQ[X] := Module[{},
    If[!KeyExistsQ[globalBundle, "Hstar"],
        InitGlobalBundle[];
        PaiComputeBundleTensors[globalBundle, "basicTools"]
    ];
    ValidateForm[globalBundle][X];
    globalBundle["Hstar"][X]
];

Clear[Contraction];

SetAttributes[Contraction, HoldFirst];

Contraction[bundle_][X_] /; AssociationQ[bundle] := Module[{},
    If[!KeyExistsQ[bundle, "contraction"],
        PaiComputeBundleTensors[bundle, "basicTools"]
    ];
    bundle["contraction"][X]
];

Contraction[X_] /; !AssociationQ[X] := Module[{},
    If[!KeyExistsQ[globalBundle, "contraction"],
        InitGlobalBundle[];
        PaiComputeBundleTensors[globalBundle, "basicTools"]
    ];
    globalBundle["contraction"][X]
];


Clear[NotAssociationQ];
NotAssociationQ[x_] := !AssociationQ[x];



Clear[FormToSparse];
Clear[FormToMatrix];
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

FormToMatrix[X_, formdegIN_:"deg", coordIN_:"Global"] := Normal[FormToSparse[X, formdegIN, coordIN]];

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

$UsePaiSimplify = True;

trigSimp = {
    x_.*Cos[h_]^2 + x_.*Sin[h_]^2 :> x,
    x_.*Cosh[h_]^2 - x_.*Sinh[h_]^2 :> x
};

PaiSimplify[expr_] := If[
	TrueQ[$UsePaiSimplify],
		Factor[expr /. trigSimp] /. trigSimp,
			expr
    ];


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
    	"GlobalID" -> id,
        "coord" -> coord,
        "ds2" -> ds2,
	"GlobalSync"-> {},
        "Tensors" -> <|
            "gdd" -> Agdd,
            "gUU" -> AgUU
        |>
    |>;
];

Clear[globalHash];
globalHash[a_, b_] := Hash[HoldComplete[{a, b}]];

Clear[InitGlobalBundle];

InitGlobalBundle[] := Module[
    {gddint, coordint, id, initFrom},

    GlobalRequired[coord];
    coordint = coord;

    Which[
    GlobalRequired[ds2],
        initFrom = "ds2";
        gddint = DiffToMatrix[ds2, coordint],
    GlobalRequired[gdd],
        initFrom = "gdd";
        gddint = ResolveGlobal[gddIN, gdd];
        Print["** Initializing Global bundle from gdd"],
    True,
        Print["[ Aborting ] neither ds2 or gdd provided as global variables"]
    ];

    id = globalHash[gddint, coordint];

    If[(!AssociationQ[globalBundle]) || (Lookup[globalBundle, "GlobalID", None] =!= id),
        Print["** Initializing  new global bundle from "<>initFrom];
    	BuildGlobalBundle[gddint, coordint, id];
        CleanComputedTensors[globalBundle];
    ];

    globalBundle
];

Clear[CleanComputedTensors];

CleanComputedTensors[bundle_] := Module[
    {id},

    If[!KeyExistsQ[bundle, "id"],
        Return[]
    ];
    id = bundle["id"];
    $ComputedTensors[id] = <||>;
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

	
(*==========================================================================================================================================*)

(*
|------------------------------------------------------
|     Thinking on Association for saving tensor
|------------------------------------------------------

                     --- Utils ---
*)

TensorProductContract[Tensors__, contractIndices_List] := Activate@TensorContract[Inactive[TensorProduct][Tensors], contractIndices];

TensorProductContract[tensor_, contractIndices_List] := TensorContract[tensor, contractIndices];

SetAttributes[ApplyIndexChange, HoldRest];
SetAttributes[MatrixForIndexChange, HoldFirst];

MoveIndicesWithMatrix[tensor_, metric_, positions_] := Module[
    {sortedPositions, rank, relations, metricSequence,
     contracted, permutation},

    sortedPositions = Sort[positions];
    rank = Length[Dimensions[tensor]];

    relations = Table[
                    {sortedPositions[[n]], rank + 2 n - 1}
                , {n, Length[sortedPositions]}
                ];

    metricSequence = Sequence @@ ConstantArray[SparseArray[metric], Length[sortedPositions]];

    contracted = TensorProductContract[tensor, metricSequence, relations];

    permutation = Join[
        Complement[Range[rank], sortedPositions],
        sortedPositions
    ];

    Transpose[contracted, permutation]
];

MatrixForIndexChange[bundle_, {"dn", "up"}] := GetTensorArray[bundle, "gUU"];
MatrixForIndexChange[bundle_, {"up", "dn"}] := GetTensorArray[bundle, "gdd"];
MatrixForIndexChange[bundle_, {"vdn", "vup"}] := Inverse[GetFlatMetric[bundle]];
MatrixForIndexChange[bundle_, {"vup", "vdn"}] := GetFlatMetric[bundle];

MatrixForIndexChange[bundle_, {"vdn", "dn"}] := GetTensorArray[bundle, "eamuUd"];
MatrixForIndexChange[bundle_, {"dn", "vdn"}] := Transpose[GetTensorArray[bundle, "eamudU"]];
MatrixForIndexChange[bundle_, {"up", "vup"}] := Transpose[GetTensorArray[bundle,"eamuUd"]];
MatrixForIndexChange[bundle_, {"vup", "up"}] := GetTensorArray[bundle, "eamudU"];

ApplyIndexChange[tensor_, bundle_, change_, positions_] := MoveIndicesWithMatrix[tensor, MatrixForIndexChange[bundle, change], positions];

SetAttributes[RaiseIndices, HoldRest];
SetAttributes[LowerIndices, HoldRest];

RaiseIndices[tensor_, bundle_, positions_] := ApplyIndexChange[tensor, bundle, {"dn", "up"}, positions];

LowerIndices[tensor_, bundle_, positions_] := ApplyIndexChange[tensor, bundle, {"up", "dn"}, positions];

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

GetTensorArray[bundle_, tensorName_, simp_:Automatic] := Module[
	{PaiTensor, TensorComponents, Dim, TensorArray, DimensionsTensor,
	sector, name, level},

    If[MemberQ[{"eamuUd", "eamudU"}, tensorName],
        If[!KeyExistsQ[bundle, tensorName],
            PaiComputeBundleTensors[bundle, "basicTools", simp]
        ];
        Return[bundle[tensorName]]
    ];

	Dim = Length[bundle["coord"]];

	{sector, name, level} = Switch[tensorName,
		"gdd" | "gUU",
			{"Tensors", tensorName, "metric"},

		"Rflatdd",
			{"FlatTensors", "Rdd", "Rdd"},

		"Rflatdddd",
			{"FlatTensors", "Rdddd", "Rdddd"},
		"omegadd",
			{"Forms", "omegadd", "spinConnection"},
		"Rformdd",
			{"Forms", "Rdd", "curvatureForm"},

		"RicciScalar",
			Which[
			KeyExistsQ[Lookup[bundle, "Tensors", <||>], "RicciScalar"],
				{"Tensors", "RicciScalar", "RicciScalar"},

			KeyExistsQ[Lookup[bundle, "FlatTensors", <||>], "RicciScalar"],
				{"FlatTensors", "RicciScalar", "RicciScalar"},

			VielbeinBundleQ[bundle],
				{"FlatTensors", "RicciScalar", "RicciScalar"},

            KeyExistsQ[bundle, "ds2"],
                {"Tensors", "RicciScalar", "RicciScalar"}
			],
		_,
			{"Tensors", tensorName, tensorName}
	];

	If[
		Not[KeyExistsQ[Lookup[bundle, sector, <||>], name]],
		PaiComputeBundleTensors[bundle, level, simp];
	];

	PaiTensor = bundle[sector, name];

	If[name === "RicciScalar",
		Return[PaiTensor]
	];

	TensorComponents[indices__] := PaiComponent[PaiTensor, {indices},
		If[tensorName === "Rformdd", "RFormdd", name]
	];
	DimensionsTensor =
		ConstantArray[Dim, tensorRank[name]];

	TensorArray =
		Array[TensorComponents, DimensionsTensor];

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

	tensorName === "omegadd",
		PaiComponent2anti[PaiTensor, {indices}],

	tensorName === "RFormdd",
		PaiComponent2anti[PaiTensor, {indices}],

	True,
		Message[PaiComponent::unk, tensorName]
	];
tensorRank = <|
	"gdd" -> 2, "gUU" -> 2, "Rdd" -> 2,
	"ChrisUdd" -> 3,
	"Rdddd" -> 4,
	"omegadd"->2
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

SetAttributes[VielbeinBundleQ, HoldFirst];

Clear[VielbeinBundleQ];

VielbeinBundleQ[bundle_] := KeyExistsQ[bundle, "eU"];

Clear[PaiComputeBundleTensors];

SetAttributes[PaiComputeBundleTensors, HoldFirst];

PaiComputeBundleTensors[bundleIN_, level_: "RicciScalar", simp_:Automatic] := Module[
	{}, 
    Which[
        VielbeinBundleQ[bundleIN],
            PaiComputeBundleTensorsVielbein[bundleIN, level, simp],

        KeyExistsQ[bundleIN, "ds2"],
            PaiComputeBundleTensorsMetric[bundleIN, level, simp],

        True,
            Print["[ Aborting ] neither ds2 nor eU was provided"];
            Abort[]
    ]
];

ClearAll[PaiComputeBundleTensorsVielbein];
SetAttributes[PaiComputeBundleTensorsVielbein, HoldFirst];

PaiComputeBundleTensorsVielbein[bundleIN_, level_: "RicciScalar", simp_:PaiSimplify] := Module[
	{bundle, needSpinConnection, needCurvature, needRdddd, needRdd, needRicciScalar, FlatTensors, Forms,
	simpVielbein, simpSpin, simpCurv, simpRiem, simpRicci, simpR},
	bundle = bundleIN;

	If[simp === Automatic,
		simpVielbein = PaiSimplify;
		simpSpin  = PaiSimplify;
		simpCurv  = Identity;
		simpRiem  = Identity;
		simpRicci = Identity;
		simpR     = Identity,
			simpVielbein = simp;
			simpSpin  = simp;
			simpCurv  = simp;
			simpRiem  = simp;
			simpRicci = simp;
			simpR     = simp
	];
	Print["** Constructing bundle Tools: Hstar, FormSquare, FormSquaredd, Contraction"];
	Print[AbsoluteTiming[InitVielbeinBundle[bundleIN, simpVielbein];]];
	bundle=bundleIN;

	If[level === "basicTools", Return[]];

	Forms = Lookup[bundle,"Forms", <| |>];
	needSpinConnection   = Not[KeyExistsQ[Forms, "omegadd"]];
	needCurvature  = Not[KeyExistsQ[Forms, "Rdd"]];
	FlatTensors = Lookup[bundle,"FlatTensors", <| |>];
	needRdddd  = Not[KeyExistsQ[FlatTensors, "Rdddd"]];
	needRdd = Not[KeyExistsQ[FlatTensors, "Rdd"]];
	needRicciScalar = Not[KeyExistsQ[FlatTensors, "RicciScalar"]];
	
	If[needSpinConnection,
		Print["** Computing spin connection"];
		Print[AbsoluteTiming[PaiComputeSpinConnection[bundle, simpSpin];]];
		bundleIN = bundle;
	];

	If[level === "spinConnection", Return[]];
	
	If[needCurvature,
		Print["** Computing curvature 2-form"];
		Print[AbsoluteTiming[PaiComputeCurvatureForm[bundle, simpCurv];]];
		bundleIN = bundle;
	];

	If[level === "curvatureForm", Return[]];

	If[needRdddd,
		Print["** computing Rdddd flat"];
		Print[AbsoluteTiming[PaiComputeRddddFlat[bundle, simpRiem];]];
		bundleIN = bundle;
	];

	If[level === "Rdddd", Return[]];

	If[needRdd,
		Print["** computing Rdd flat"];
		Print[AbsoluteTiming[PaiComputeRddFlat[bundle, simpRicci];]];
		bundleIN = bundle;
	];

	If[level === "Rdd", Return[]];

	If[needRicciScalar,
		Print["** computing RicciScalar"];
		Print[AbsoluteTiming[PaiComputeRicciScalarFlat[bundle, simpR];]];
		bundleIN = bundle;
	];
];


ClearAll[PaiComputeBundleTensorsMetric];
SetAttributes[PaiComputeBundleTensorsMetric, HoldFirst];
PaiComputeBundleTensorsMetric[bundleIN_, level_: "Rdddd", simp_:Automatic] := Module[
	{Tensors, needMetric, needMetricTools, needChris, needRiemann, AgddgUU, Agdd, AgUU, 
	AChrisUdd, ARiemdddd, bundle, ARicdd,needRicci, needRicciScalar, RicciScalar,
	simpMetric, simpChris, simpRiem, simpRicci, simpR},

    ValidateMetricBundle[bundleIN];

	If[simp === Automatic,
		simpMetric = PaiSimplify;
		simpChris  = PaiSimplify;
		simpRiem   = Identity;
		simpRicci  = Identity;
		simpR      = Identity,
			simpMetric = simp;
			simpChris  = simp;
			simpRiem   = simp;
			simpRicci  = simp;
			simpR      = simp
	];
	
	bundle = bundleIN;
	Tensors = Lookup[bundle, "Tensors", <||>];

	needMetric   = Not[KeyExistsQ[Tensors, "gdd"]] || Not[KeyExistsQ[Tensors, "gUU"]];
	needMetricTools = needMetric || Not[KeyExistsQ[bundle, "Hstar"]] || Not[KeyExistsQ[bundle, "FormSquare"]] || Not[KeyExistsQ[bundle, "FormSquaredd"]];
	needChris    = Not[KeyExistsQ[Tensors, "ChrisUdd"]] && MemberQ[{"RicciScalar","ChrisUdd", "Rdddd", "Rdd"}, level];
	needRiemann  = Not[KeyExistsQ[Tensors, "Rdddd"]] && MemberQ[{"RicciScalar","Rdddd", "Rdd"}, level];
	needRicci    = Not[KeyExistsQ[Tensors, "Rdd"]] && MemberQ[{"RicciScalar","Rdd"}, level];
	needRicciScalar = Not[KeyExistsQ[Tensors, "RicciScalar"]] && MemberQ[{"RicciScalar"}, level];
	
	(* --- Metric --- *)
	If[needMetric,
		Print["** Computing metric"];
		AgddgUU = PaiComputeMetric[bundle];
		Agdd = Map[simpMetric, AgddgUU["gdd"]];
		AgUU = Map[simpMetric, AgddgUU["gUU"]];
		Tensors = Join[Tensors, <|"gdd" -> Agdd, "gUU" -> AgUU|>];
		bundle = AssociateTo[bundle, "Tensors" -> Tensors];
	];
	
	If[level === "metric", 
		bundleIN = bundle;
		Return[]
	];

	(* --- Metric bundle tools --- *)
	If[needMetricTools,
		Print["** Computing basicTools : Hstar, FormSquare, FormSquaredd, Contraction"];
		InitMetricTools[bundle, simpMetric];
	];

	If[level === "basicTools", 
		bundleIN = bundle;
		Return[]
	];

	(* --- Christoffel --- *)
	If[needChris,
		Print["** Computing Christoffel"];
		AChrisUdd = PaiComputeChrisUdd[bundle];
		AChrisUdd = Map[simpChris, AChrisUdd];
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
		ARiemdddd = Map[simpRiem, ARiemdddd];
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
		ARicdd = Map[simpRicci, ARicdd];
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
		RicciScalar = simpR[PaiComputeRicciScalar[bundle]];
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
BuildHodge[bundle_, simp_:PaiSimplify] := Module[{},
	Which[

    VielbeinBundleQ[bundle],
        Return[BuildHodgeVielbein[bundle, simp]],

    KeyExistsQ[bundle, "ds2"],
        Return[BuildHodgeMetric[bundle, simp]],

	True,
		Print["[ Aborting ] BuildHodge: ds2 or eU are not given."];
	];
];

SetAttributes[BuildHodgeMetric, HoldFirst];
BuildHodgeMetric[bundle_, simp_:PaiSimplify] := 
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
			HstarCore[X, gUU, sqrtdetg, d[coord], simp]
			]];
	];

BuildHodgeVielbein[bundle_, simp_:PaiSimplify] := Module[{basis, eta, etainv, sqrtdeteta},
	basis = bundle["basis"];
	eta = GetFlatMetric[bundle]; 
	etainv = Inverse[eta];
	sqrtdeteta = Sqrt[-Det[eta]];
	Return[Function[{X}, 
			HstarCore[X /. bundle["dxToe"], etainv, sqrtdeteta, basis, simp]
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

Clear[ValidateVielbeinBundle];

SetAttributes[ValidateVielbeinBundle, HoldFirst];

ValidateVielbeinBundle[bundle_] := Module[
    {coord, basis, eU, allowedDiffs, diffs, invalid, lambda, scaled, homogen},

    If[
        !And[
            KeyExistsQ[bundle, "coord"],
            KeyExistsQ[bundle, "basis"],
            KeyExistsQ[bundle, "eU"]
        ],
        Print[
            "[ Aborting ] Vielbein bundle requires keys ",
            "\"coord\", \"basis\" and \"eU\""
        ];
        Abort[]
    ];

    coord = bundle["coord"];
    basis = bundle["basis"];
    eU = bundle["eU"];

    If[Length[eU] =!= Length[coord] || Length[basis] =!= Length[coord],
        Print[
            "[ Aborting ] Vielbein dimension mismatch",
            "\nLength[coord]: ", Length[coord],
            "\nLength[basis]: ", Length[basis],
            "\nLength[eU]: ", Length[eU]
        ];
        Abort[]
    ];

    allowedDiffs = d[coord];
    diffs = DeleteDuplicates[Cases[eU, _d, Infinity]];
    invalid = Complement[diffs, allowedDiffs];

    If[invalid =!= {},
        Print[
            "[ Aborting ] Unexpected coordinate differentials in eU: ",
            invalid
        ];
        Abort[]
    ];

    scaled = eU /. d[_] :> lambda;

    homogen = !AllTrue[scaled, (Exponent[#, lambda, Min] === 1 && Exponent[#, lambda, Max] === 1)&];

    If[homogen,
        Print[
            "[ Aborting ] Each vielbein must be linear in d[coord]"
        ];
        Abort[]
    ];

    If[!KeyExistsQ[bundle, "eta"] && !KeyExistsQ[bundle, "signature"],
            Print[
            "** [ Observation ] No \"eta\" metric or flat signature provided. ",
            "Using mostly-plus signature (-,+,...,+)"
            ];  
    ];

    True
];

SetAttributes[InitVielbeinBundle, HoldFirst];

InitVielbeinBundle[bundle_, simp_:PaiSimplify] := Module[
	{eTodx, dxToe, symbs, eU, contraction, coordbasis, hstar,
	deU, dictde, FormSquareTools, aIter, muIter, eta,
    etaUU, gdd, gUU, Agdd, AgUU, eamuUd, eamudU, buildSymmetric2, Tensors},

	bundle = Association[bundle];

    ValidateVielbeinBundle[bundle];

	symbs = bundle["basis"];
	eU = bundle["eU"];
	coordbasis = d[bundle["coord"]];
	Do[FormDegree[eIter] = 1, {eIter, symbs}];

	eTodx = Normal[AssociationThread[symbs -> eU]];

    dxToe = Quiet[
        Check[Solve[eU == symbs, coordbasis][[1]],
            Print["[ Aborting ] Could not construct dxToe"];
            Abort[]
        ]
    ];

	contraction = ConstructContraction[bundle];

    eamuUd = Table[
        Coefficient[eU[[aIter]], coordbasis[[muIter]]],
        {aIter, Length[symbs]},
        {muIter, Length[coordbasis]}
    ];

    eamudU = Transpose[Inverse[eamuUd]];

    eta = GetFlatMetric[bundle];
    etaUU = Inverse[eta];

    gdd = Transpose[eamuUd] . eta . eamuUd;
    gUU = Transpose[eamudU] . etaUU . eamudU;

    buildSymmetric2[X_] := Association[Table[{iIter, jIter} -> X[[iIter, jIter]], {iIter, Length[coordbasis]}, {jIter, iIter, Length[coordbasis]}]];

    Agdd = CleanZeros @ Map[simp, buildSymmetric2[gdd]];
    AgUU = CleanZeros @ Map[simp, buildSymmetric2[gUU]];

    Tensors = Lookup[bundle, "Tensors", <||>];

    Tensors = Join[Tensors, <|"gdd" -> Agdd, "gUU" -> AgUU|>];

	deU = d[eU] /. dxToe;
    dictde = AssociationThread[d[symbs], deU];
	Do[d[eIter] = Collect[dictde[d[eIter]], _Wedge, simp],{eIter, symbs}];
	AssociateTo[bundle, {
        "eTodx" -> eTodx,
        "dxToe" -> dxToe,
        "eamuUd" -> eamuUd,
        "eamudU" -> eamudU,
        "contraction"->contraction,
        "UseVielbein" -> True,
        "Tensors" -> Tensors
    }];
	hstar = BuildHodge[bundle, simp];
	bundle["Hstar"] = hstar;
	FormSquareTools = BuildSquaresTools[bundle, simp];
	bundle["FormSquare"] = FormSquareTools["FormSquare"];
	bundle["FormSquaredd"] = FormSquareTools["FormSquaredd"];

	];

SetAttributes[InitMetricTools, HoldFirst];

InitMetricTools[bundle_, simp_:PaiSimplify] := Module[{ToClear, FormSquareTools},
    bundle = Association[bundle];
	If[KeyExistsQ[bundle, "constants"],
    	Do[d[cIter]=0, {cIter, bundle["constants"]}];
	];
    bundle["basis"] = d[bundle["coord"]];
    bundle["contraction"] = ConstructContraction[bundle];
    bundle["Hstar"] = BuildHodge[bundle, simp];
    FormSquareTools = BuildSquaresTools[bundle, simp];
    bundle["FormSquare"] = FormSquareTools["FormSquare"];
    bundle["FormSquaredd"] = FormSquareTools["FormSquaredd"];
];

PaiComponent2anti[X_Association, {i_, j_}] := Which[
	i < j,
		First[Lookup[X, {{i,j}}, 0]],
	i > j,
		-First[Lookup[X, {{j,i}}, 0]],
	True,
		0
];

GetFlatMetric[bundle_] := Which[
    KeyExistsQ[bundle, "eta"],
        bundle["eta"],

    KeyExistsQ[bundle, "signature"],
        DiagonalMatrix[bundle["signature"]],

    True,
        DiagonalMatrix[
            Join[{-1}, ConstantArray[1, Length[bundle["basis"]] - 1]]
        ]
];

Clear[PaiComputeSpinConnection];
SetAttributes[PaiComputeSpinConnection, HoldFirst];
PaiComputeSpinConnection[frameBundle_, simp_:PaiSimplify] := Module[
	{eta, deU, symbs, Dim, a, b, c, k,
	Aomegadd, Forms, omega, deDNA, CUdd, Cddd},

 	frameBundle = Association[frameBundle];
	eta = GetFlatMetric[frameBundle];
	Dim = Length[frameBundle["basis"]];
	symbs = frameBundle["basis"];
	
	deU = Table[d[e], {e, symbs}];
	deU = simp[deU];

	deDNA[a_] := deDNA[a] = If[deU[[a]]===0, 
		<| |>,
		Association[
		Map[#[[2]] -> #[[1]] &, DNAofForm[deU[[a]], symbs]]]
        ];

	CUdd[a_,b_,c_] := Which[
	    b < c,
		First[Lookup[deDNA[a], {{b,c}}, 0]],

	    b > c,
		-First[Lookup[deDNA[a], {{c,b}}, 0]],

	    True,
		0
	];

	Cddd[a_,b_,c_] :=
	    Cddd[a,b,c] =
		Sum[
		    eta[[a,k]]*CUdd[k,b,c],
		    {k,Dim}
		];

	omega[a_,b_] := (1/2 Sum[
        ( Cddd[a,b,c] + Cddd[b,c,a] - Cddd[c,a,b])*symbs[[c]],{c,Dim}]);

	Aomegadd =
		Association[
			Flatten[
				Table[
					{a,b} -> omega[a,b],
					{a,Dim}, {b,a+1,Dim}
				],
				1
			]
		];
	Forms = Lookup[frameBundle, "Forms", <||>];
	Forms = AssociateTo[Forms, "omegadd" -> Aomegadd];
	frameBundle = AssociateTo[frameBundle, "Forms" -> Forms];
];
	
Clear[PaiComputeCurvatureForm];
SetAttributes[PaiComputeCurvatureForm, HoldFirst];
PaiComputeCurvatureForm[frameBundle_, simp_: PaiSimplify] := Module[
	{eta, etaUU, Dim, Aomegadd, vielbein, ARdd, Rformdd, omegadd,
	aIter, bIter, cIter, dIter, Forms},
	eta = GetFlatMetric[frameBundle]; 
	Aomegadd = frameBundle["Forms", "omegadd"];
	Dim = Length[frameBundle["basis"]];

	omegadd[a_,b_] := PaiComponent2anti[Aomegadd, {a,b}];
	etaUU = Inverse[eta];

	Rformdd[a_,b_] := (d[omegadd[a,b]]
		+ Sum[
			etaUU[[cIter,dIter]]*
				Wedge[omegadd[a,cIter], omegadd[dIter,b]],
			{cIter,Dim}, {dIter,Dim}
		]) /. frameBundle["dxToe"];

	ARdd =
		Association[
			Flatten[
				Table[
					{aIter,bIter} -> Rformdd[aIter,bIter],
					{aIter,Dim}, {bIter,aIter+1,Dim}
				],
				1
			]
		];
	ARdd = Map[simp, ARdd];
	ARdd = CleanZeros[ARdd];
	Forms = Lookup[frameBundle, "Forms", <||>];
	Forms = AssociateTo[Forms, "Rdd" -> ARdd];
	frameBundle = AssociateTo[frameBundle, "Forms" -> Forms];
	];

Clear[PaiComputeRddddFlat];
SetAttributes[PaiComputeRddddFlat, HoldFirst];

PaiComputeRddddFlat[frameBundle_, simp_:PaiSimplify] := Module[
	{Dim, ARdd, RFormdd, contraction, Rdddd,
	pairs, nPairs, rules, p, q, ARdddd, FlatTensors,
	symbs, RFormDNA, a, b, c, d},

	Dim = Length[frameBundle["basis"]];
	ARdd = frameBundle["Forms", "Rdd"];
	symbs = frameBundle["basis"];

	RFormdd[a_,b_] := PaiComponent2anti[ARdd, {a,b}];

	RFormDNA[a_,b_] := RFormDNA[a,b] = If[RFormdd[a,b]===0,
		<| |>, Association[Map[
			#[[2]] -> #[[1]] &,
			DNAofForm[RFormdd[a,b], symbs]]]
		];

	Rdddd[a_,b_,c_,d_] := First[Lookup[RFormDNA[a,b], {{c,d}}, 0]];

	pairs = Subsets[Range[Dim], {2}];
	nPairs = Length[pairs];

	rules = Table[
		With[
			{ab = pairs[[p]], cd = pairs[[q]]},
			Join[ab,cd] -> Apply[Rdddd, Join[ab,cd]]
		],
		{p,nPairs},
		{q,p,nPairs}
	];

	rules = Flatten[rules, 1];

	ARdddd = CleanZeros[Association[rules]];
	ARdddd = Map[simp, ARdddd];

	FlatTensors = Lookup[frameBundle, "FlatTensors", <||>];
	FlatTensors = AssociateTo[FlatTensors, "Rdddd" -> ARdddd];
	frameBundle = AssociateTo[frameBundle, "FlatTensors" -> FlatTensors];
];

Clear[PaiComputeRddFlat];
SetAttributes[PaiComputeRddFlat, HoldFirst];

PaiComputeRddFlat[frameBundle_, simp_:PaiSimplify] := Module[
	{etaUU, Dim, ARdddd, Rdddd, Rdd, ARdd,
	FlatTensors, aIter, bIter, cIter, dIter},

	etaUU = Inverse[GetFlatMetric[frameBundle]];
	Dim = Length[frameBundle["basis"]];
	ARdddd = frameBundle["FlatTensors", "Rdddd"];

	Rdddd[a_,b_,c_,d_] :=
		PaiComponent4Riem[ARdddd, {a,b,c,d}];

	Rdd[b_,d_] :=
		Sum[
			etaUU[[aIter,cIter]]*Rdddd[aIter,b,cIter,d],
			{aIter,Dim}, {cIter,Dim}
		];

	ARdd = Association[
		Table[
			{bIter,dIter} -> Rdd[bIter,dIter],
			{bIter,Dim}, {dIter,bIter,Dim}
		]
	];

	ARdd = Map[simp, ARdd];
	ARdd = CleanZeros[ARdd];

	FlatTensors = Lookup[frameBundle, "FlatTensors", <||>];
	FlatTensors = AssociateTo[FlatTensors, "Rdd" -> ARdd];
	frameBundle = AssociateTo[frameBundle, "FlatTensors" -> FlatTensors];
];

Clear[PaiComputeRicciScalarFlat];
SetAttributes[PaiComputeRicciScalarFlat, HoldFirst];

PaiComputeRicciScalarFlat[frameBundle_, simp_:PaiSimplify] := Module[
	{etaUU, Dim, ARdd, Rdd, RicciScalar, FlatTensors, a, b},

	etaUU = Inverse[GetFlatMetric[frameBundle]];
	Dim = Length[frameBundle["basis"]];
	ARdd = frameBundle["FlatTensors", "Rdd"];

	Rdd[a_,b_] :=
		PaiComponent2sym[ARdd, {a,b}];

	RicciScalar =
		Sum[
			etaUU[[a,b]]*Rdd[a,b],
			{a,Dim}, {b,Dim}
		];

	RicciScalar = simp[RicciScalar];

	FlatTensors = Lookup[frameBundle, "FlatTensors", <||>];
	FlatTensors = AssociateTo[
		FlatTensors,
		"RicciScalar" -> RicciScalar
	];
	frameBundle = AssociateTo[
		frameBundle,
		"FlatTensors" -> FlatTensors
	];
];

(* ========================================================== *)

(*                  PaiTensor / GRTensor layer                *)

(* ========================================================== *)

Clear[TensorSignIndices, TensorSignDerivatives, TensorSignRank,
    TensorSignDerivativeOrder, TensorSignAllIndices];

TensorSignIndices[sign_] := sign[[2]];

TensorSignDerivatives[sign_] := sign[[3]];

TensorSignRank[sign_] := Length[TensorSignIndices[sign]];

TensorSignDerivativeOrder[sign_] := Length[TensorSignDerivatives[sign]];

TensorSignAllIndices[sign_] := Join[TensorSignDerivatives[sign],
    TensorSignIndices[sign]];

TensorSignHead[sign_] := sign[[1]];

MakeTensorSign[head_, indices_, derivatives_] := {head, indices, derivatives};

Clear[InitComputedTensors];

$ComputedTensors = <||>;

SetAttributes[InitComputedTensors, HoldFirst];

InitComputedTensors[bundle_] := Module[
	{id},
	

	If[KeyExistsQ[bundle, "id"],
        id = bundle["id"],
            id = CreateUUID["PaiBundle-"];
            AssociateTo[bundle, "id" -> id]
	];

    If[Not[KeyExistsQ[$ComputedTensors, id]],
        AssociateTo[$ComputedTensors, id -> <| |>]
    ];

    If[Not[KeyExistsQ[$DefTensors, id]],
        AssociateTo[$DefTensors, id -> <| |>]
    ];

    id
];

NativeTensorSources[bundle_] := If[
    VielbeinBundleQ[bundle],

    <|
        MakeTensorSign["g", {"dn", "dn"}, {}] -> "gdd",
        MakeTensorSign["g", {"up", "up"}, {}] -> "gUU",
        MakeTensorSign["R", {"vdn", "vdn"}, {}] -> "Rflatdd",
        MakeTensorSign["R", {"vdn", "vdn", "vdn", "vdn"}, {}] -> "Rflatdddd",
        MakeTensorSign["Ricciscalar", {}, {}] -> "RicciScalar",
        MakeTensorSign["Rform", {"vdn", "vdn"}, {}] -> "Rformdd",
        MakeTensorSign["omega", {"vdn", "vdn"}, {}] -> "omegadd"
    |>,

    <|
        MakeTensorSign["g", {"dn", "dn"}, {}] -> "gdd",
        MakeTensorSign["g", {"up", "up"}, {}] -> "gUU",
        MakeTensorSign["Chris", {"up", "dn", "dn"}, {}] -> "ChrisUdd",
        MakeTensorSign["R", {"dn", "dn"}, {}] -> "Rdd",
        MakeTensorSign["R", {"dn", "dn", "dn", "dn"}, {}] -> "Rdddd",
        MakeTensorSign["Ricciscalar", {}, {}] -> "RicciScalar"
    |>
];
Clear[ComputeNativeTensorSeed];
SetAttributes[ComputeNativeTensorSeed, HoldRest]

ComputeNativeTensorSeed[tensorSign_, bundle_, simp_:PaiSimplify] := Module[
    {candidates, best, nativeSign, tensorName},

    candidates = KeySelect[
        NativeTensorSources[bundle],
        AcceptableSeedTensorQ[tensorSign]
    ];

    If[candidates === <||>,
        Return[False]
    ];

    best = FindMostSimilarTensor[candidates, tensorSign];

    nativeSign = First[Keys[best]];
    tensorName = First[Values[best]];

    StoreComputedTensor[
        bundle,
        nativeSign,
        GetTensorArray[bundle, tensorName, simp],
        simp
    ];

    True
];

Clear[StoreComputedTensor];

StoreComputedTensor[bundle_, tensorSign_, tensor_, simp_:PaiSimplify] := Module[{id},
	id = bundle["id"];
	$ComputedTensors[id] = Append[
		$ComputedTensors[id],
		tensorSign -> simp[tensor]
	];
]

ClearAll[ParseIndex, TensorIndices];

ParseIndex[s_String] :=
    If[
        StringStartsQ[s, "^"],
        {StringDrop[s, 1], "up"},
        {s, "dn"}
    ];

TensorIndices[tensor_String] := Module[
    {inside, pieces, tensorIndices, derivativeIndices},

    If[Not[StringContainsQ[tensor, "{"]],
        Return[{}]
    ];

    inside = First@StringCases[tensor, "{" ~~ x___ ~~ "}" :> x];
    pieces = StringSplit[inside, ";"];
    tensorIndices = Map[ParseIndex, StringSplit[First[pieces]]];

    derivativeIndices =
        Map[
            ParseIndex,
            Flatten@Map[StringSplit, Rest[pieces]]
        ];

    Join[Reverse[derivativeIndices], tensorIndices]
];

ClearAll[IndexStructure, ReadTensorSignature];

IndexStructure[s_String] :=
    Map[
        If[StringStartsQ[#, "^"], "up", "dn"] &,
        StringSplit[StringTrim[s]]
    ];

ReadTensorSignature[tensor_String] := Module[
    {head, inside, posCD, tensorPart, derivativePart},

    If[Not[StringContainsQ[tensor, "{"]],
        Return[MakeTensorSign[StringTrim[tensor], {}, {}]]
    ];

    head = StringTrim[First[StringSplit[tensor, "{"]]];

    inside = First[StringCases[tensor, "{" ~~ x___ ~~ "}" :> x]];

    posCD = StringPosition[inside, ";"];

    If[posCD === {},
        tensorPart = inside;
        derivativePart = "",
            tensorPart = StringTake[inside, posCD[[1, 1]] - 1];
            derivativePart = StringDrop[inside, posCD[[1, 1]]]
    ];

    MakeTensorSign[head,
    IndexStructure[tensorPart],
    Reverse[IndexStructure[StringReplace[derivativePart, ";" -> " "]]]
    ]
];


Clear[PaiDef, $DefTensors];
$DefTensors=<| "shared" -> <| |> |>;


SetAttributes[PaiDef, HoldFirst];

DefineStringTensor[storage_, tensorDef_String] := Module[
	{splitDef, tensor, def, TensorSign, previous},
	splitDef = StringSplit[tensorDef, ":="];
	tensor = StringTrim[splitDef[[1]]];
	def = StringTrim[splitDef[[2]]];
    TensorSign = ReadTensorSignature[tensor];

    previous = Select[
        Keys[$DefTensors[storage]], 
            TensorSignHead[#] === TensorSignHead[TensorSign] &&
            TensorSignRank[#] === TensorSignRank[TensorSign] &
    ];

    If[
        previous =!= {},
        Print[
            "[ Aborting ] Tensor < ",
            TensorSignHead[TensorSign],
            " > with rank ",
            TensorSignRank[TensorSign],
            " is already defined"
        ];
        Abort[]
    ];

	AssociateTo[$DefTensors[storage], TensorSign -><|tensor->def|>];
    Print["** Definition created ", TensorSignToString[TensorSign], " in "<>
        If[storage === "shared",
            "shared definitions",
                "local bundle"
        ]
    ];
];



PaiDef[tensorDef_String] := DefineStringTensor["shared", tensorDef];

PaiDef[bundle_][tensorDef_String] := Module[
    {id},
    id = InitComputedTensors[bundle];
    DefineStringTensor[id, tensorDef]
];


PaiObjectToArray[bundle_, tensorSign_, object_] := Module[
    {rank, indices},

    rank = TensorSignRank[tensorSign];
    indices = TensorSignIndices[tensorSign];

    Which[
        ListQ[object] && ArrayDepth[object] === rank,
            Print["** Recognized object as Array"];
            object,

        rank === 0 && Not[ListQ[object]],
            Print["** Recognized object as scalar"];
            object,

        rank === 2 && DeleteDuplicates[indices]==={"dn"} && MetricQuadraticInDiffQ[object] && FormDegree[object]===0,
            Print["** Recognized object as quadratic form"];
            DiffToMatrix[object, bundle["coord"]],

        FormDegree[object] === rank && DeleteDuplicates[indices]==={"dn"},
            Print["** Recognized object as ",rank, "-form"];
            FormToMatrix[object, rank, bundle["coord"]],

        True,
            Print[
                "[ Aborting ] Could not interpret object as tensor ",
                TensorSignToString[tensorSign]
            ];
            Abort[]
    ]
];

PaiDef[bundle_][tensor_String, object_] := Module[
    {tensorSign, rank, Dim, expectedDimensions, tensorArray},

    InitComputedTensors[bundle];

    tensorSign = TensorStringToSign[tensor];

    Dim = Length[bundle["coord"]];
    rank = TensorSignRank[tensorSign];
    expectedDimensions = ConstantArray[Dim, rank];

    tensorArray = PaiObjectToArray[bundle, tensorSign, object];

    If[rank>0 && ListQ[tensorArray] && (Dimensions[tensorArray] =!= expectedDimensions),
        Print[
            "[ Aborting ] Tensor ", tensor,
            " has dimensions ", Dimensions[tensorArray],
            ", expected ", expectedDimensions
        ];
        Abort[]
    ];

    StoreComputedTensor[bundle, tensorSign, tensorArray];

    Print["** Tensor registered ", TensorSignToString[tensorSign]];
];

PaiDef[tensor_String, tensorArray_] := Module[{},
    InitGlobalBundle[];
    PaiDef[globalBundle][tensor, tensorArray]
    ];

$IndexDistanceGraph = Graph[
    {
        Property["vdn" -> "vup", EdgeWeight -> 1],
        Property["vup" -> "vdn", EdgeWeight -> 1],

        Property["dn" -> "up", EdgeWeight -> 2],
        Property["up" -> "dn", EdgeWeight -> 2],

        Property["dn" -> "vdn", EdgeWeight -> 3],
        Property["vdn" -> "dn", EdgeWeight -> 3],

        Property["up" -> "vup", EdgeWeight -> 3],
        Property["vup" -> "up", EdgeWeight -> 3]
    }
];

IndexDistance[from_, to_] := GraphDistance[$IndexDistanceGraph, from, to];

FindIndexPath[from_, to_] := FindShortestPath[$IndexDistanceGraph, from, to];

TruncateTargetIndices[candidate_, target_] := Module[
    {nDerOrder},
    nDerOrder = TensorSignDerivativeOrder[candidate];
    Join[
        Take[TensorSignDerivatives[target], -nDerOrder],
        TensorSignIndices[target]
    ]
];

FindTensorSignPath[candidate_, target_] := Module[{},
    MapThread[
        FindIndexPath,
        {TensorSignAllIndices[candidate], TruncateTargetIndices[candidate, target]}
    ]
];

TensorSignDistance[candidate_, target_] := Module[{},
    Total[
        MapThread[
            IndexDistance,
            {TensorSignAllIndices[candidate], TruncateTargetIndices[candidate, target]}
        ]
    ]
];

Clear[FindMostSimilarTensor];
FindMostSimilarTensor[closests_Association, tensorSign_List] := Module[
    {signs, bestSign},

    signs = Keys[closests];
    bestSign = First[MinimalBy[signs, TensorSignDistance[#, tensorSign]& ]];

    KeyTake[closests, {bestSign}]
];

CurrentIndexChangeBatch[paths_] := Module[
    {indexedChanges},

    indexedChanges = MapIndexed[
        Function[{path, position},
            If[
                Length[path] > 1,
                {Take[path, 2], First[position]},
                Nothing
            ]
        ],
        paths
    ];

    GroupBy[indexedChanges, First -> Last]
];

NextIndexPaths[paths_] := Map[If[Length[#] > 1, Rest[#], #] &, paths];

SetAttributes[FollowTensorSignPaths, HoldRest];

FollowTensorSignPaths[tensor_, bundle_, paths_] := Module[
    {result, remaining, batch},

    result = tensor;
    remaining = paths;

    While[AnyTrue[remaining, Length[#] > 1 &],

        batch = CurrentIndexChangeBatch[remaining];

        KeyValueMap[
            (result = ApplyIndexChange[result, bundle, #1, #2]) &,
            batch
        ];

        remaining = NextIndexPaths[remaining];
    ];

    result
];

Clear[ComputeCovDTensor];
SetAttributes[ComputeCovDTensor, HoldRest];
ComputeCovDTensor[best_, bundle_, simp_:PaiSimplify] := Module[
	{bestSign, bestSparse, bestIndCD, sparseCD, newSign},

    bestSign = Keys[best][[1]];
    bestSparse = Values[best][[1]];
    bestIndCD = StringJoin[TensorSignAllIndices[bestSign] /. {"up"->"U", "dn"-> "d"}];

    newSign = MakeTensorSign[
        TensorSignHead[bestSign],
        TensorSignIndices[bestSign],
        Prepend[TensorSignDerivatives[bestSign], "dn"]
    ];

	Print["** Computing ", TensorSignToString[bestSign]];

	sparseCD = PaiCovD[bundle, bestSparse, bestIndCD];
	
	StoreComputedTensor[bundle, newSign, sparseCD, simp]
];

AcceptableSeedTensorQ[tensorSign_] := And[
    TensorSignHead[#]===TensorSignHead[tensorSign],
	TensorSignRank[#]===TensorSignRank[tensorSign],
    TensorSignDerivativeOrder[#]<=TensorSignDerivativeOrder[tensorSign]
]&;

SetAttributes[ComputeSingleRequiredTensors, HoldRest]

ComputeSingleRequiredTensors[tensorSign_, bundle_, simp_:PaiSimplify] := Module[
    {usefullComputed, closestDerivatives, best, CompTensors,
    defCandidates, defSign, storage},

    CompTensors = $ComputedTensors[bundle["id"]];
	usefullComputed = KeySelect[CompTensors, AcceptableSeedTensorQ[tensorSign]];

	If[usefullComputed === <||>,

        If[ComputeNativeTensorSeed[tensorSign, bundle, simp],
            Return[ComputeSingleRequiredTensors[tensorSign, bundle, simp]]
        ];

        defCandidates = KeySelect[$DefTensors[bundle["id"]], AcceptableSeedTensorQ[tensorSign]];
        storage = bundle["id"];

        If[defCandidates ===<| |>,
            defCandidates = KeySelect[$DefTensors["shared"], AcceptableSeedTensorQ[tensorSign]];
            storage = "shared";
        ];

        If[defCandidates === <||>,
            Print[
                "[ Aborting ] Tensor ", TensorSignToString[tensorSign], " is neither computed nor defined"
            ];
            Abort[]
        ];

		If[Length[defCandidates] > 1,
			Print[
				"[ Aborting ] Multiple definitions can seed tensor ", tensorSign, ": ", Keys[defCandidates]
			];
			Abort[]
		];

		defSign = First[Keys[defCandidates]];

		ComputeFreshTensor[defSign, storage, bundle, simp];

		Return[ComputeSingleRequiredTensors[tensorSign, bundle, simp]];

    ];

	closestDerivatives = KeyTake[usefullComputed, MaximalBy[Keys[usefullComputed], TensorSignDerivativeOrder]
    ];

	best = FindMostSimilarTensor[closestDerivatives, tensorSign];
    bestSign = First[Keys[best]];

    If[
        TensorSignDerivativeOrder[bestSign] ===
            TensorSignDerivativeOrder[tensorSign],

        If[bestSign === tensorSign,
            Return[]
        ];

        sparse = First[Values[best]];

        Print["** Computing ", TensorSignToString[tensorSign]];

        sparse = FollowTensorSignPaths[
            sparse,
            bundle,
            FindTensorSignPath[bestSign, tensorSign]
        ];

        StoreComputedTensor[
            bundle,
            tensorSign,
            sparse,
            simp
        ];

        Return[]
    ];

    ComputeCovDTensor[best, bundle, simp];

    ComputeSingleRequiredTensors[tensorSign, bundle, simp]
];

Clear[ComputeRequiredTensors];

SetAttributes[ComputeRequiredTensors, HoldRest];
ComputeRequiredTensors[requiredTensors_, bundle_, simp_:PaiSimplify]:= Module[{},
	If[Length[$ComputedTensors[bundle["id"]]] === 0,
		InitComputedTensors[bundle]];
	Do[
	ComputeSingleRequiredTensors[tensor, bundle, simp]
	, {tensor, requiredTensors}]
];

FindRequiredScalars[scalars_, bundle_] := Module[
    {CompTensors, scalarSigns, namesInScalars},

    CompTensors = $ComputedTensors[bundle["id"]];

    scalarSigns = DeleteDuplicates @ Join[
        Keys @ KeySelect[CompTensors, ScalarTensorQ],
        Keys @ KeySelect[$DefTensors[bundle["id"]], ScalarTensorQ],
        Keys @ KeySelect[$DefTensors["shared"], ScalarTensorQ],
        Keys @ KeySelect[NativeTensorSources[bundle], ScalarTensorQ]
    ];

    namesInScalars =
        DeleteDuplicates @ Flatten @
            StringCases[
                scalars,
                RegularExpression["[A-Za-z$][A-Za-z0-9$]*"]
            ];

    Select[
        scalarSigns,
        MemberQ[namesInScalars, TensorSignHead[#]] &
    ]
];

ScalarTensorQ[sign_] := TensorSignRank[sign] === 0 && TensorSignDerivativeOrder[sign] === 0;

Clear[EvalScalarQuantities];
EvalScalarQuantities[ComputedTensors_] := Normal[KeyMap[ToExpression[TensorSignHead[#]]&, KeySelect[ComputedTensors, ScalarTensorQ]]]

Clear[PaiCompute];
Clear[PaiCalc];

SetAttributes[PaiCompute, HoldFirst];
SetAttributes[PaiCalc, HoldFirst];
PaiCalc[y___]:=PaiCompute[y];
PaiCompute[bundle_][spec_String, simp_:PaiSimplify] := Module[
    {tensorSign},
    tensorSign = TensorStringToSign[spec];
    InitComputedTensors[bundle];

    If[
        KeyExistsQ[$ComputedTensors[bundle["id"]], tensorSign],
        Print["** Tensor ", spec, " already computed"];
        Return[]
    ];

    ComputeSingleRequiredTensors[tensorSign, bundle, simp];
];

PaiCompute[spec_, simp_:PaiSimplify] /; StringQ[spec] := Module[
    {},
	InitGlobalBundle[];
    PaiCompute[globalBundle][spec, simp];
];

Clear[Paillaco];
SetAttributes[Paillaco, HoldFirst];

Paillaco[bundle_][spec_String, simp_:PaiSimplify] := Module[{},
    PaiCompute[bundle][spec, simp];
    PaiComponents[bundle][spec]
];

Paillaco[spec_String, simp_:PaiSimplify] := Module[{},
    PaiCompute[spec, simp];
    PaiComponents[spec]
];

Clear[PaiComponents];
SetAttributes[PaiComponents, HoldFirst];

PaiComponents[bundle_][spec_String] := Module[
    {tensorSign},

    tensorSign = TensorStringToSign[spec];

    If[
        !KeyExistsQ[$ComputedTensors[bundle["id"]], tensorSign],
        Print["[ Aborting ] Tensor ", spec, " has not been computed"];
        Abort[]
    ];

    $ComputedTensors[bundle["id"]][tensorSign]
];


PaiComponents[spec_] /; StringQ[spec] := PaiComponents[globalBundle][spec];

(*
====================================================
        Towards tree-like Compute tensors
        generalization of Monomial computation
====================================================

    *)

Clear[DecomposeDefinition];

DecomposeDefinition[expr_String] := Module[
    {str, terms, factors},

    str = StripOuterParentheses[StringTrim[expr]];

    If[
        TensorLeafQ[str],
        Return[str]
    ];

    terms = SplitTensorSum[str];

    If[
        Length[terms] > 1,
        Return[
            <|"plus" -> Map[DecomposeDefinition, terms]|>
        ]
    ];

    factors = SplitTensorTimes[str];

    If[
        Length[factors] > 1,
        Return[
            <|"times" -> Map[DecomposeDefinition, factors]|>
        ]
    ];

    str
(*
    Print["[ Aborting ] Could not decompose expression: ", str];
    Abort[]*)
]

Clear[TopLevelOperatorPositions];

TopLevelOperatorPositions[expr_String, ops_List] := Module[
    {chars, par = 0, cur = 0, bra = 0, positions = {}, ch},

    chars = Characters[expr];

    Do[
        ch = chars[[i]];

        If[
            par === 0 && cur === 0 && bra === 0 && MemberQ[ops, ch],
                AppendTo[positions, i]
        ];

        Switch[ch,
            "(", par++,
            ")", par--,
            "{", cur++,
            "}", cur--,
            "[", bra++,
            "]", bra--
        ],
        {i, Length[chars]}
    ];

    positions
];

Clear[BinarySignQ];

BinarySignQ[chars_, i_] := Module[
    {prev, prevChars},
	prevChars = Reverse[Take[chars, i - 1]];
    prev = SelectFirst[prevChars, StringTrim[#] =!= "" &, None];
    prev =!= None && Not@MemberQ[{"+", "-", "*", "/", "^", "(", "[", "{", ",", ";"}, prev]
];

Clear[SplitTensorSum];

Clear[RemoveLeadingPlus];

RemoveLeadingPlus[str_String] := If[
    StringStartsQ[str, "+"],
    StringTrim[StringDrop[str, 1]],
    str
];

SplitTensorSum[expr_String] := Module[
    {str, chars, signPositions, starts, ends, terms, ranges},

    str = StringTrim[expr];
    chars = Characters[str];

    signPositions = Select[
        TopLevelOperatorPositions[str, {"+", "-"}],
        BinarySignQ[chars, #] &
    ];

    If[signPositions === {},
        Return[{str}]
    ];

    starts = Join[{1}, signPositions];
    ends = Join[signPositions - 1, {StringLength[str]}];
    ranges = Transpose[{starts, ends}];

	terms = Map[StringTake[str, #] &, ranges];
    terms = Map[StringTrim, terms];

    Map[RemoveLeadingPlus, terms]
];

Clear[SplitTensorTimes];

SplitTensorTimes[expr_String] := Module[
    {str, TimesPositions, starts, ends, ranges, terms},

    str = StringTrim[expr];

    TimesPositions = TopLevelOperatorPositions[str, {"*"}];

    If[TimesPositions === {},
        If[StringStartsQ[str, "-"] && str =!= "-1",
            Return[
                {"-1", StringTrim[StringDrop[str, 1]]}
            ]
        ];
            Return[{str}]
    ];

    starts = Join[{1}, TimesPositions + 1];
    ends = Join[TimesPositions - 1, {StringLength[str]}];
    ranges = Transpose[{starts, ends}];
	terms = Map[StringTake[str, #] &, ranges];
    Map[StringTrim, terms];

    terms
];

Clear[OuterParenthesizedQ];

OuterParenthesizedQ[expr_String] := Module[
    {str, chars, depth = 0, closesEarly = False},
    str = StringTrim[expr];

    If[StringLength[str] < 2 || StringTake[str, 1] =!= "(" || StringTake[str, -1] =!= ")",
        Return[False]
    ];

    chars = Characters[str];

    Do[
        Switch[
            chars[[i]],
            "(", depth++,
            ")", depth--
        ];

        If[depth === 0 && i < Length[chars],
            closesEarly = True;
            Break[]
        ],
    {i, Length[chars]}
    ];

    depth === 0 && !closesEarly
];

Clear[StripOuterParentheses];

StripOuterParentheses[expr_String] := Module[
    {str},

    str = StringTrim[expr];

    While[
        OuterParenthesizedQ[str],
        str = StringTrim[StringTake[str, {2, -2}]]
    ];

    str
];

Clear[TensorLeafQ];

TensorLeafQ[s_String] := Module[
    {str, open, close, head, openPositions, closePositions, altern},

    str = StringTrim[s];

    openPositions  = StringPosition[str, "{"];
    closePositions = StringPosition[str, "}"];

    If[
        Length[openPositions] =!= 1 || Length[closePositions] =!= 1,
        Return[False]
    ];
    open = openPositions[[1, 1]];
    close = closePositions[[1, 1]];

    If[
        open >= close || close =!= StringLength[str],
        Return[False]
    ];

    head = StringTrim[StringTake[str, open - 1]];
    altern = Alternatives["+", "-", "*", "/", "^","(", ")", "[", "]"];

    head =!= "" && StringFreeQ[head, altern]
]


(*

Evaluation tools

    *)
SetAttributes[EvaluateLeaf, HoldRest];
SetAttributes[EvaluateDefinitionTree, HoldRest];
SetAttributes[EvaluateScalarLeaf, HoldRest];
SetAttributes[EvaluateTensorLeaf, HoldRest];

Clear[EvaluateLeaf];

EvaluateLeaf[leaf_String, bundle_, simp_:PaiSimplify] := If[TensorLeafQ[leaf],
                                          EvaluateTensorLeaf[leaf, bundle, simp],
                                              EvaluateScalarLeaf[leaf, bundle, simp]
   ];

Clear[EvaluateDefinitionTree];

EvaluateDefinitionTree[node_String, bundle_, simp_:PaiSimplify] := EvaluateLeaf[node, bundle, simp];

EvaluateDefinitionTree[<|"plus" -> children_|>, bundle_, simp_:PaiSimplify] := EvaluatePlus[
           Map[EvaluateDefinitionTree[#, bundle, simp] &, children]
       ];

EvaluateDefinitionTree[<|"times" -> children_|>, bundle_, simp_:PaiSimplify] := EvaluateTimes[
           Map[EvaluateDefinitionTree[#, bundle, simp] &, children], simp
       ];

Clear[EvaluateScalarLeaf];

EvaluateScalarLeaf[leaf_String, bundle_, simp_:PaiSimplify] := Module[
       {requiredScalars, value},
       requiredScalars = FindRequiredScalars[{leaf}, bundle];
       ComputeRequiredTensors[requiredScalars, bundle, simp];
       value = ToExpression[
            StringReplace[
                leaf,
                "dimInter" -> ToString[Length[bundle["coord"]]]
            ]
        ] /. EvalScalarQuantities[$ComputedTensors[bundle["id"]]];
       <|"value" -> value, "indices" -> {}|>
];


Clear[EvaluateTensorLeaf];

EvaluateTensorLeaf[leaf_String, bundle_, simp_:PaiSimplify] := Module[
    {tensorSign, value, indices, contractions,
    freePositions},

    tensorSign = ReadTensorSignature[leaf];
    indices = TensorIndices[leaf];

    ComputeRequiredTensors[{tensorSign}, bundle, simp];

    value = $ComputedTensors[bundle["id"]][tensorSign];

    contractions = GetContractionsFromIndices[indices];

    If[contractions =!= {},
        value = TensorProductContract[value, contractions]
    ];

    freePositions = Complement[Range[Length[indices]], Flatten[contractions]];

    <|"value" -> value, "indices" -> indices[[freePositions]]|>
];

Clear[GetContractionsFromIndices];

GetContractionsFromIndices[allIndices_List] := Module[
    {groups, repeated, badMultiplicity, badUpDownPair},

    groups = GatherBy[Range[Length[allIndices]], allIndices[[#, 1]]&];

    badMultiplicity = Select[groups, Length[#] > 2 &];

    If[badMultiplicity =!= {},
        Print[
            "[ Aborting ] Index appears more than twice: ",
            Map[allIndices[[First[#], 1]] &, badMultiplicity]
        ];
        Abort[]
    ];

    repeated = Select[groups, Length[#] == 2 &];

    badUpDownPair = Select[repeated, Length[DeleteDuplicates[allIndices[[#, 2]]]] =!= 2 &];

    If[badUpDownPair =!= {},
        Print[
            "[ Aborting ] Contracted indices must appear once up and once down: ",
            Map[allIndices[[First[#], 1]] &, badUpDownPair]
        ];
        Abort[]
    ];

    repeated
];

Clear[EvaluateTimes];

EvaluateTimes[children_List, simp_:PaiSimplify] := Module[
    {scalars, tensors, scalarFactor, tensorValues,
    allIndices, contractions, freePositions, value},

    scalars = Select[children, #["indices"] === {} &];

    tensors = Select[children, #["indices"] =!= {} &];

    scalarFactor = If[scalars === {},
                       1, 
                           Apply[Times, Lookup[scalars, "value"]]
                   ];

    If[tensors === {},
        Return[
            <|"value" -> scalarFactor, "indices" -> {}|>
        ]
    ];

    tensorValues = simp[Lookup[tensors, "value"]];

    allIndices = Flatten[Lookup[tensors, "indices"], 1];

    contractions = GetContractionsFromIndices[allIndices];

    value = Apply[TensorProductContract, Append[tensorValues, contractions]] /. TensorProduct[aaI_, bbI_]:>aaI*bbI;

    freePositions = Complement[Range[Length[allIndices]], Flatten[contractions]];

    <|
        "value" -> scalarFactor value,
        "indices" -> allIndices[[freePositions]]
    |>
];

Clear[AlignEvaluatedIndices];

AlignEvaluatedIndices[term_Association, targetIndices_List] := Module[
    {indices, value, permutation},

    indices = term["indices"];
    value = term["value"];

    If[Sort[indices] =!= Sort[targetIndices],
        Print[
            "[ Aborting ] Incompatible free indices in sum: ",
            indices,
            " and ",
            targetIndices
        ];
        Abort[]
    ];

    If[targetIndices === {},
        Return[value]
    ];

    permutation = Map[
        First[FirstPosition[targetIndices, #]] &,
        indices
    ];

    If[permutation === Range[Length[permutation]],
        value,
            Transpose[value, permutation]
    ]
];

Clear[EvaluatePlus];

EvaluatePlus[children_List] := Module[
    {targetIndices, values},

    targetIndices = First[children]["indices"];

    values = Map[
        AlignEvaluatedIndices[#, targetIndices] &,
        children
    ];

    <|
        "value" -> Apply[Plus, values],
        "indices" -> targetIndices
    |>
];


Clear[ComputeFreshTensor];
SetAttributes[ComputeFreshTensor, HoldRest];

ComputeFreshTensor[defSign_, storage_, bundle_, simp_:PaiSimplify] := Module[
    {allDef, tensor, def, tree, result,
    targetIndices, tensorSparse},

    allDef = $DefTensors[storage][defSign];

    tensor = First[Keys[allDef]];
    def = First[Values[allDef]];

    tree = DecomposeDefinition[def];

    result = EvaluateDefinitionTree[
        tree,
        bundle,
        simp
    ];

    targetIndices = TensorIndices[tensor];

    tensorSparse = simp[AlignEvaluatedIndices[
        result,
        targetIndices
    ]];

    StoreComputedTensor[bundle, defSign, tensorSparse, simp];

    tensorSparse
];


Clear[TensorSignToString];

TensorSignToString[{head_, indices_, derivatives_}] := Module[
    {ind, der},

    If[indices==={} && derivatives==={},
        Return[head]
    ];

    der = Reverse[derivatives] /. {"up" -> "Dup", "dn" -> "Ddn"};

    head <> "(" <> StringRiffle[Join[der, indices], ","] <> ")"
];


Clear[ParseComputeSpec];

TensorStringToSign[spec_String] := Module[
    {head, inside, indices, str},


    str = StringTrim[spec];

    If[
        StringFreeQ[str, {"(", ")"}],
        Return[MakeTensorSign[str, {}, {}]]
    ];

    head = StringTrim[First[StringSplit[str, "("]]];
    inside = First[StringCases[str, "(" ~~ x___ ~~ ")" :> x]];
    indices = StringTrim /@ StringSplit[inside, ","];

    If[
        Not[AllTrue[indices, MemberQ[{"up", "dn", "vup", "vdn"}, #] &]],
        Print["[ Aborting ] Invalid tensor indices in ", spec];
        Abort[]
    ];


    {head, indices, {}}
];

DefineStringTensor[
    "shared",
    "Weyl{a b c d} := R{a b c d} - 1/(dimInter-2)*(g{a c}*R{d b} - g{b c}*R{d a}-g{a d}*R{c b} + g{b d}*R{c a}) + 1/(dimInter-1)/(dimInter-2)*Ricciscalar*(g{a c}*g{d b} - g{a d}*g{c b})"
];


End[]

EndPackage[]
