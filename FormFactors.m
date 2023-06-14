(* ::Package:: *)

BeginPackage["FormFactors`"];


f0D::usage = "f0 form factor for B -> D transitions."
fpD::usage = "f+ form factor for B -> D transitions."
fTD::usage = "fT form factor for B -> D transitions."
A0Ds::usage = "A0 form factor for B -> D* transitions."
A1Ds::usage = "A1 form factor for B -> D* transitions."
A2Ds::usage = "A2 form factor for B -> D* transitions."
A12Ds::usage = "A12 form factor for B -> D* transitions."
A3Ds::usage = "A3 form factor for B -> D* transitions."
VDs::usage = "V form factor for B -> D* transitions."
T1Ds::usage = "T1 form factor for B -> D* transitions."
T2Ds::usage = "T2 form factor for B -> D* transitions."
T23Ds::usage = "T23 form factor for B -> D* transitions."
T3Ds::usage = "T3 form factor for B -> D* transitions."


Begin["`Private`"]
	tp[MB_, MM_] := (MB + MM)^2
	tm[MB_, MM_] := (MB-MM)^2
	t0[MB_, MM_] := tp[MB, MM] (1-Sqrt[tm[MB, MM]/tp[MB, MM]])
	z[t_, MB_, MM_] := (Sqrt[tp[MB, MM]-t]-Sqrt[tp[MB, MM]-t0[MB, MM]])/(Sqrt[tp[MB, MM]-t]+Sqrt[tp[MB, MM]-t0[MB, MM]])
	F[q2_, MRF_, \[Alpha]0_, \[Alpha]1_, \[Alpha]2_, MB_, MM_] := 1/(1-q2/MRF^2) (\[Alpha]0 + \[Alpha]1 (z[q2, MB, MM]-z[0, MB, MM])+\[Alpha]2 (z[q2, MB, MM]-z[0, MB, MM])^2)
	MRFf0D = 5.540;
	MRFfpD = 5.325;
	MRFA0Ds = 5.279;
	MRFA1Ds = 5.724;
	MB0 = 5.27966;
	MD = 1.86966;
	MDs = 2.01026;
	\[Alpha]1f0D = 1.845829677892925;
	\[Alpha]2f0D = 0.04744873027689452;
	\[Alpha]0fpD = 0.6494022492599675;
	\[Alpha]1fpD = -1.348121551696309;
	\[Alpha]2fpD = 1.460191255411654;
	\[Alpha]0fTD = 0.5656604515314846;
	\[Alpha]1fTD = -2.492010124810672;
	\[Alpha]2fTD = 10.54601528980545;
	\[Alpha]0A0Ds = 0.670050553248689;
	\[Alpha]1A0Ds = -2.199987732264871;
	\[Alpha]2A0Ds = 4.142108872726158;
	\[Alpha]0A1Ds = 0.599986664237791;
	\[Alpha]1A1Ds = 0.9206797471679162;
	\[Alpha]2A1Ds = 0.804639222042595;
	\[Alpha]0A12Ds = (MB0^2-MDs^2)/(8 MB0 MDs) \[Alpha]0A0Ds;
	\[Alpha]1A12Ds = 0.2290328090903999;
	\[Alpha]2A12Ds = -0.2245847149849823;
	\[Alpha]0VDs = 0.6915010412736222;
	\[Alpha]1VDs = -1.199753864310677;
	\[Alpha]2VDs = -3.773903692351789;
	\[Alpha]0T1Ds = 0.6345457207774636;
	\[Alpha]1T1Ds = -1.415987126407053;
	\[Alpha]2T1Ds = -0.7044948165322704;
	\[Alpha]1T2Ds = 1.585024376626186;
	\[Alpha]2T2Ds = -0.4904037070379432;
	\[Alpha]0T23Ds = 0.8092462257968681;
	\[Alpha]1T23Ds = 0.6096479114969915;
	\[Alpha]2T23Ds = 4.87422542072693;
	f0D[q2_] := F[q2, MRFf0D, \[Alpha]0fpD, \[Alpha]1f0D, \[Alpha]2f0D, MB0, MD]
	fpD[q2_] := F[q2, MRFfpD, \[Alpha]0fpD, \[Alpha]1fpD, \[Alpha]2fpD, MB0, MD]
	fTD[q2_] := F[q2, MRFfpD, \[Alpha]0fTD, \[Alpha]1fTD, \[Alpha]2fTD, MB0, MD]
	A0Ds[q2_] := F[q2, MRFA0Ds, \[Alpha]0A0Ds, \[Alpha]1A0Ds, \[Alpha]2A0Ds, MB0, MDs]
	A1Ds[q2_] := F[q2, MRFA1Ds, \[Alpha]0A1Ds, \[Alpha]1A1Ds, \[Alpha]2A1Ds, MB0, MDs]
	A12Ds[q2_] := F[q2, MRFA1Ds, \[Alpha]0A12Ds, \[Alpha]1A12Ds, \[Alpha]2A12Ds, MB0, MDs]
	\[Lambda][q2_, MB_, MM_] := ((MB+MM)^2-q2)((MB-MM)^2-q2)
	A2Ds[q2_] := ((MB0+MDs)^2(MB0^2-MDs^2-q2) A1Ds[q2] - 16 MB0 MDs^2(MB0+MDs) A12Ds[q2])/\[Lambda][q2, MB0, MDs]
	A3Ds[q2_] := (MB0+MDs)/(2MDs) A1Ds[q2] - (MB0-MDs)/(2MDs) A2Ds[q2]
	VDs[q2_] := F[q2, MRFfpD, \[Alpha]0VDs, \[Alpha]1VDs, \[Alpha]2VDs, MB0, MDs]
	T1Ds[q2_] := F[q2, MRFfpD, \[Alpha]0T1Ds, \[Alpha]1T1Ds, \[Alpha]2T1Ds, MB0, MDs]
	T2Ds[q2_] := F[q2, MRFA1Ds, \[Alpha]0T1Ds, \[Alpha]1T2Ds, \[Alpha]2T2Ds, MB0, MDs]
	T23Ds[q2_] := F[q2, MRFA1Ds, \[Alpha]0T23Ds, \[Alpha]1T23Ds, \[Alpha]2T23Ds, MB0, MDs]
	T3Ds[q2_] := ((MB0^2-MDs^2)(MB0^2+2MDs^2-q2)T2Ds[q2] -8 MB0 MDs^2 (MB0-MDs) T23Ds[q2]) /\[Lambda][q2, MB0, MDs]


End[]
EndPackage[]



