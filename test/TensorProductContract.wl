(* ::Package:: *)

Quit


Get["/home/marcelo/src/PaillacoDiff/PaillacoDiff.wl"];

bundle = 
	<|
		"ds2" -> -f[r]*d[t]^2 + d[r]^2/f[r] + r^2/4*(d[th]^2 + Sin[th]^2*d[ph]^2 + (d[psi] + Cos[th]*d[ph])^2),
		"coord" -> {t,r,th,ph,psi}
	|>;
PaiComputeBundleTensors[bundle,"RicciScalar", Simplify]


gdd = SparseArray[GetTensorArray[bundle,"gdd"]];
Rdddd = SparseArray[GetTensorArray[bundle,"Rdddd"]];
Rdd = SparseArray[GetTensorArray[bundle,"Rdd"]];
RUd = RaiseIndices[Rdd, bundle, {1}];
RUddd = RaiseIndices[Rdddd, bundle, {1}];
RdUUU = RaiseIndices[Rdddd, bundle, {2, 3, 4}];
RUUdd = RaiseIndices[Rdddd, bundle, {1,2}];
RicciScalar = GetTensorArray[bundle,"RicciScalar"];

LGB = Simplify[TensorProductContract[RUUdd, RUUdd, {{3,5},{4,6}, {1,7},{2, 8}}] - 4*TensorProductContract[RUd, RUd, {{1,4},{2,3}}] + RicciScalar^2];
RicRiemdd = TensorProductContract[RUd, RUddd, {{2,2 + 1}, {1, 2 + 3}}];
RiemRiemdd = TensorProductContract[Rdddd, RdUUU, {{2,4 + 2}, {3, 4 + 3}, {4, 4 + 4}}];

Hdd = 2*RicciScalar*Rdd - 4*Rdd . RUd - 4*RicRiemdd + 2*RiemRiemdd -1/2*LGB*gdd;
Gdd = Rdd -1/2*gdd*RicciScalar;
Edd = Simplify[Normal[Gdd + Lamb*gdd + alph*Hdd]];

Print[Simplify[Edd[[1,1]]/.{f -> Function[{r}, 1 + r^2/4/alph*(1 - Sqrt[1 + 4*alph*Lamb/3 + mu*alph/r^4])]}] === 0]
