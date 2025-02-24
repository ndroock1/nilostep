

(*
    RHpack: a Mathematica package to go with the book "Equivalents of the Riemann Hypothesis, Vol One Arithmetic Equivalents."
 by Kevin Broughan.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                    !
    ! IF THIS TEXT APPEARS IN YOUR BROWSER WINDOW, USE   !
    ! THE FILE MENU BUTTON AT THE TOP OF YOUR BROWSER TO !
    ! SAVE THE DATA TO A FILE, WITH NAME RHpack.m, IN ANY   !
    ! SUBDIRECTORY ON YOUR DIRECTORY TREE IN WHICH YOU   !
    ! ARE ABLE TO SAVE FILES.                            !
    !                                                    !
    ! THEN READ THE FIRST FEW PARAGRAPHS OF THE RHpack
    ! DOCUMENTATION, OR THE BOOK APPENDIX B, TO SEE HOW TO !
    ! LOAD THE PACKAGE RHpack.m INTO MATHEMATICA.           !
    !                                                    !
    ! TO FIND OUT ABOUT DIRECTORY TREES CONSULT ONE OF   !
    ! THE TUTORIALS ON THE RHpack HOME PAGE OR TALK   !
    ! TO YOUR COMPUTER SUPPORT PERSON.                   !
    ! IF THERE ARE PROBLEMS EMAIL kab@waikato.ac.nz .    !
    !                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should obtain a copy of the GNU General Public License
    along with this program from its website; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


    Contact details: Kevin Broughan, 
                     Department of Mathematics,
                     University of Waikato, 
                     Hamilton Private Bag 3105,
                     New Zealand.
                     Email: kab@waikato.ac.nz

Package: Name RHpack.
Author: Kevin Broughan
Copyright (c) Kevin Broughan, 2017.
Version: 16 October 2017.
*)

Print[StyleForm["RHpack version December 2016 loading ....
Copyright(c) 2016 Kevin A. Broughan.  
RHpack comes with ABSOLUTELY NO WARRANTY. This is free software, 
and you are welcome to redistribute it under the terms of the GNU Public Licence. 
For assistance consult kab@waikato.ac.nz .....",FontWeight->"Bold"]];
coq=CompositeQ;

BeginPackage["RHpack`"]



(* Functions: 16 October 2017
Parameters:*)
$a1::usage = "";
 $a2::usage = ""; 
$a3::usage = "";
 $H::usage = "";
 $R::usage = "";
$RosserSchoenfeldR::usage = "";
$FordR::usage = "";

(*Chapter 2: zeta *)
Psi0::usage = "[x,nz,d]";
PlotVonMangoldtPsi::usage = "[a,b]";pvm=PlotVonMangoldtPsi;
BacklundB::usage = "[T,d]";blb=BacklundB;
ZetaZeroHeight::usage = "[n,d]";zzh=ZetaZeroHeight;
LargestZetaZero::usage = "[H,d]";lzz=LargestZetaZero;
NextZetaZero::usage = "[h,d]";nzz=NextZetaZero;

(* Chapter 3: Estimates *)
ComputeTheta::usage = "[x,d]"; cth=ComputeTheta ;
ComputeThetaValues::usage = "[b,d]"; ctv=ComputeThetaValues;
ComputeKm::usage = "[m,d]"; ckm=ComputeKm ;
PsiSmallXEpsilon::usage = "[b,m,d]"; pse=PsiSmallXEpsilon;
PsiLargeXEpsilon::usage = "[x,R,d]"; ple=PsiLargeXEpsilon;
PsiLargeXEpsilon2::usage = "[x,R,d]"; ple2=PsiLargeXEpsilon2;
VerifyThetaX::usage = "[a,b,d]"; vth=VerifyThetaX;
Psi::usage = "[x,d]";

(* Chapter 4: Classical equivalences*)
M::usage = "[x]";

(* Chapter 5: Euler's totient function*)
NOverPhi::usage = "[n,d]"; nop=NOverPhi;
NicolasInequalityQ::usage = "[n,d]";niq=NicolasInequalityQ;
Computef::usage = "[x,d]"; cof=Computef;
NicolasTwoFailures::usage = "[k1,k0,d]"; n2f=NicolasTwoFailures;
ComputeH::usage = "[n,d]"; coh=ComputeH;
ComputeDeltas::usage = "[d]"; cod=ComputeDeltas;
FHalf::usage = "[z,d]";fha=FHalf;
PlotEulerPhiRatio::usage = "[b]";pep=PlotEulerPhiRatio;
CheckCNk::usage = "[k0,d]";cnk=CheckCNk;
PlotCNk::usage = "[k1,k0,d]"; pcnk=PlotCNk;
(* psitheta - worker *)
ComputeC::usage = "[n,d]";coc=ComputeC;
(* PlotAsymptoticNSimp - worker *)
(* cN66 - worker *)
(* cN1 - worker *)
(* CheckNicolasThmTwo::usage = "[b,d]";*)

(* Chapter 6: A variety of abundant numbers*)
HighlyCompositeQ::usage = "[n]";hcq=HighlyCompositeQ;
HardyRamanujanQ::usage = "[n]";hrq=HardyRamanujanQ;
SuperAbundantQ::usage = "[n]";saq=SuperAbundantQ;
ColossallyAbundantQ::usage = "[n]";caq=ColossallyAbundantQ;
ColossallyAbundantInteger::usage = "[e]"; cai=ColossallyAbundantInteger;
ColossallyAbundantF::usage = "[x,a,d]";caf=ColossallyAbundantF;
CriticalEpsilonValues::usage = "[b,jmax,d]";cev=CriticalEpsilonValues;

(* Chapter 7: Robin's theorem *)
RobinsInequalityQ::usage = "[n,d]"; riq=RobinsInequalityQ;
PlotRobinsInequality::usage = "[a,b,d]";pri=PlotRobinsInequality;
SigmaOverN::usage = "[n]"; son=SigmaOverN;
LagariasInequalityQ::usage = "[n,d]";liq=LagariasInequalityQ;
UnitaryDivisorSigma::usage = "[n]"; uds=UnitaryDivisorSigma;

(* Chapter 8: Numbers not satisfying Robin's inequality*)
SetA::usage = ""; 
SetB::usage = ""; 
SetC::usage = ""; 
SetD::usage = "";
NthPrimeBounds::usage = "[n,d]";prb=NthPrimeBounds;
KthPrimorial::usage = "[k]"; kpp=KthPrimorial;
LogKthPrimorialBounds::usage = "[n,d]";lppb=LogKthPrimorialBounds;
KthApproximatePrimorial::usage = "[n,d]"; app=KthApproximatePrimorial;

(* Chapter 9: Left, Right and Extremely abundant numbers*)
GronwallG::usage = "[n,d]";grg=GronwallG;
ExtraordinaryQ::usage = "[n]";exq=ExtraordinaryQ;
RightAbundantQ::usage = "[n,b,d]";raq=RightAbundantQ;
LeftAbundantQ::usage = "[n,d]";laq=LeftAbundantQ;
ExtremelyAbundantQ::usage = "[n,d]";eaq=ExtremelyAbundantQ;

(* Chapter 10: Other equivalences to RH *)
FareyFractions::usage = "[n]";ffr=FareyFractions;
RedhefferMatrix::usage = "[n]";rem=RedhefferMatrix;
DivisibilityGraph::usage = "[n]";dvg=DivisibilityGraph;
PlotDivisibilityGraph::usage = "[n]";pdg=PlotDivisibilityGraph;
Xi::usage = "[z,d]";
Lambda::usage = "[z,d]";lam=Lambda;
S::usage = "[T,d]";
IntegerHeight::usage = "[n]"; iht=IntegerHeight;
ZetaZeroCount::usage = "[T,d]"; zzc=ZetaZeroCount;
ZetaZeroHeight::usage="[n,d]"; zzh=ZetaZeroHeight;

(* Appendix A: Tables *)
GenerateHC::usage = "[b]";ghc=GenerateHC;
GenerateHR::usage = "[b]";ghr=GenerateHR;
GenerateSA::usage = "[b]";gsa=GenerateSA;
GenerateCA::usage = "[b]";gca=GenerateCA;
GenerateLandau::usage = "[n]";gla=GenerateLandau;
FindNextXA::usage = "[n,d]";fnx=FindNextXA;
CAPrimes::usage = "[b]";cap=CAPrimes;
InitialNR::usage = ""; inr=InitialNR;
 InitialHC::usage = ""; ihc=InitialHC;
 InitialSA::usage = "";isa=InitialSA;
 InitialCA::usage = "";ica=InitialCA;
 InitialXA::usage = "";ixa=InitialXA;
 InitialLandau::usage = "";ila=InitialLandau;
 InitialCAPrimes::usage = "";icap=InitialCAPrimes;
 InitialCA150::usage = ""; ca1=InitialCA150;
 InitialSA100::usage = ""; sa1=InitialSA100;
HardyRamanujanToPrimorialForm::usage = "[n]"; hr2pf=HardyRamanujanToPrimorialForm;
PrimorialFormQ::usage = "[p1,p2,...]"; prq=PrimorialFormQ;
PrimorialFormToInteger::usage = "[{p1,p2,...}]"; pf2i=PrimorialFormToInteger;

(* Public Utilities: *)
KFreeQ::usage = "[n,k]";kfq=KFreeQ;
NthPrime::usage = "[n]";npr=NthPrime;
MaxPrime::usage = "[x]";mxp=MaxPrime;
NthHarmonicNumber::usage = "[n,d]"; nha=NthHarmonicNumber;
NiceFactorInteger::usage = "[n]";nfi=NiceFactorInteger;
AllPrimesFactor::usage = "[n]";apf=AllPrimesFactor;
PrimePower::usage = "[p,n]"; prp=PrimePower;
(*CompositeQ::usage = "[n]";*)
PowerfulQ::usage = "[n]"; pfq=PowerfulQ;



(* Parameters: *)
$a1 = 122/1000;
$a2 = 278/1000;
$a3 = 2510/1000;
$H = (545/100)* 10^8;
$R = 8;

$RosserSchoenfeldR = 515/((Sqrt[546] - Sqrt[322])^2);

$FordR = 8;

Begin["`Private`"];

(* Chapter 2: Zeta =================================================*)

Clear[Psi0];
Psi0::argx = "Psi0 takes three arguments";
Psi0::arg1 = "the first argument is a positive real number";
Psi0::arg2 = "the second argument is a natural number";
Psi0::arg3 = "the third argument is a natural number";
Psi0[] := (Message[Psi0::argx]; Abort[]);
Psi0[argseq___] := Module[{args = List[argseq], x, nz, d, ans, j, gj},
  If[Length[args] != 3, (Message[Psi0::argx]; Abort[])];
  x = args[[1]]; nz = args[[2]]; d = args[[3]];
  If[Not[PositiveRealQ[x ]], (Message[Psi0::arg1]; Abort[])];
  If[Not[NaturalNumberQ[nz]], (Message[Psi0::arg2]; Abort[])];
  If[Not[NaturalNumberQ[d ]], (Message[Psi0::arg3]; Abort[])];
  Return[Psi0Worker[x, nz, d]]]
Psi0Worker[x_, nz_, d_] := Module[{ans, j, gj},
   ans = x - Log[1 - 1/x^2]/2 - Log[2 Pi] -
    Sqrt[x] Sum[gj = Im[N[ZetaZero[j], d]]; 
      N[(Cos[gj Log[x]] + 2 gj Sin[gj Log[x]])/(1/4 + gj^2), d], {j, 
       1, nz}];
  Return[ans]]

(*-----------------------------------------------------*)

Clear[PlotVonMangoldtPsi];
PlotVonMangoldtPsi::argx = "the function takes four arguments";
PlotVonMangoldtPsi::arg1 = 
  "the first argument is a positive real number";
PlotVonMangoldtPsi::arg2 = 
  "the second argument is a positive real number";
PlotVonMangoldtPsi::arg1 = "the third argument is a natural number";
PlotVonMangoldtPsi::arg2 = "the fourth argument is a natural number"; 
PlotVonMangoldtPsi[] := (Message[PlotVonMangoldtPsi::argx]; Abort[]);
PlotVonMangoldtPsi[argseq___] := 
 Module[{args = List[argseq], a, b, nx, d},
  If[Length[args] != 4, (Message[PlotVonMangoldtPsi::argx]; Abort[])];
  a = args[[1]]; b = args[[2]]; nx = args[[3]]; d = args[[4]];
  If[Not[PositiveRealQ[a ]], (Message[PlotVonMangoldtPsi::arg1]; 
    Abort[])];
  If[Not[PositiveRealQ[b]], (Message[PlotVonMangoldtPsi::arg2]; 
    Abort[])];
  If[a >= b, Print["the first argment must be less than the second"]; 
   Abort[]];
  If[Not[NaturalNumberQ[nx ]], (Message[PlotVonMangoldtPsi::arg3]; 
    Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[PlotVonMangoldtPsi::arg4]; 
    Abort[])]; 
  Print[Plot[Psi0Worker[x, nx, d], {x, a, b}, ImageSize -> 500, 
    AxesStyle -> Thick, AxesLabel -> {"x", "\[Psi]0(x)"}]];
  Return[Psi0Worker[b, nx, d]]]

(*-------------------------------------------------------*)

(* precision output will be the minimum of T and d *)
Clear[BacklundB];
BacklundB::argx = "the function takes two arguments";
BacklundB::arg1 = "the first argument is a positive real number";
BacklundB::arg2 = "the second argument is a natural number";
BacklundB[] := (Message[BacklundB::argx]; Abort[]);
BacklundB[argseq___] := Module[{args = List[argseq], T, d},
  If[Length[args] != 2, (Message[BacklundB::argx]; Abort[])];
  T = args[[1]]; d = args[[2]]; 
  If[Not[PositiveRealQ[T ]], (Message[BacklundB::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[BacklundB::arg2]; Abort[])];
  Return[BacklundBWorker[T, d]]]
BacklundBWorker[T_, d_] := 
 Module[{a1 = N[122/1000, d], a2 = N[278/1000, d], 
   a3 = N[2510/1000, d]},
  Return[If[T < N[Exp[1]], Abort[], a1 Log[T] + a2 Log[Log[T]] + a3]]]

(*=========================================================*)



SetAttributes[ZetaZero, NHoldFirst];

Clear[ZetaZeroHeight];
ZetaZeroHeight::argx = "Two arguments required";
ZetaZeroHeight::arg1 = "First argument is a natural number";
ZetaZeroHeight::"Second argument is a natural number";
ZetaZeroHeight[] := (Message[ZetaZeroHeight::argx]; Abort[]);
ZetaZeroHeight[argseq___] := Module[{args = List[argseq], n, d, zn},
  If[Length[args] != 2, (Message[ZetaZeroHeight::argx]; Abort[])];
  n = Part[args, 1]; d = Part[args, 2];
  If[Not[NaturalNumberQ[n]], (Message[ZetaZeroHeight::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[ZetaZeroHeight::arg2]; Abort[])];
  Return[ZetaZeroHeightWorker[n, d]]]
ZetaZeroHeightWorker[n_, d_] := Im[N[ZetaZero[n], d]]


(*-----------------------------------------------*)


Clear[LargestZetaZero];
LargestZetaZero::argx = "the function takes two arguments";
LargestZetaZero::arg1 = "the first argument is a positive real number";
LargestZetaZero::arg2 = "the second argument is a natural number";
LargestZetaZero[] := (Message[LargestZetaZero::argx]; Abort[]);
LargestZetaZero[argseq___] := Module[{args = List[argseq], h, d},
  If[Length[args] != 2, (Message[LargestZetaZero::argx]; Abort[])];
  h = args[[1]]; d = args[[2]]; 
  If[Not[PositiveRealQ[h]], (Message[LargestZetaZero::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[LargestZetaZero::arg2]; 
    Abort[])];
  Return[LargestZetaZeroWorker[h, d]]]
LargestZetaZeroWorker[h_, d_] := Module[{n = 1, zh, zzh}, 
  While[(zh = ZetaZeroHeightWorker[n, d]) <= h, zzh = zh; n++];
  Return[zzh]]

(*------------------------------------------------------------------*)

Clear[NextZetaZero];
NextZetaZero::argx = "the function takes three arguments";
NextZetaZero::arg1 = "the first argument is a positive real number";
NextZetaZero::arg2 = "the second argument is a natural number";
NextZetaZero[] := (Message[NextZetaZero::argx]; Abort[]);
NextZetaZero[argseq___] := Module[{args = List[argseq], h, d},
  If[Length[args] != 2, (Message[NextZetaZero::argx]; Abort[])];
  h = args[[1]]; d = args[[2]];
  If[Not[PositiveRealQ[h]], (Message[NextZetaZero::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[NextZetaZero::arg2]; Abort[])];
  Return[NextZetaZeroWorker[h, d]]]
NextZetaZeroWorker[h_, d_] := Module[{n = 1, zh},
  While[(zh = ZetaZeroHeightWorker[n, d]) <=  h, n++];
  Return[zh]]

(* Chapter 3: Estimates \
=================================================*)

k[1] = 463/10000;
k[2] = 167/10^5;
k[3] = 744/10^7;
Clear[ComputeKm];
ComputeKm::argx = "the function takes two arguments";
ComputeKm::arg1 = "the first argument is a positive real number";
ComputeKm::arg2 = "the second argument is a natural number";
ComputeKm[] := (Message[ComputeKm::argx]; Abort[]);
ComputeKm[argseq___] := Module[{args = List[argseq], m, d, k},
  If[Length[args] != 2, (Message[ComputeKm::argx]; Abort[])];
  m = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[m ]], (Message[ComputeKm::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[ComputeKm::arg2]; Abort[])];
  k[1] = N[463/10000, d];
  k[2] = N[167/10^5, d];
  k[3] = N[744/10^7, d];
  ComputeKmWorker[m, k, d];
  Return[Table[k[j], {j, 1, m}]]]
ComputeKmWorker[m_, k_, d_] := 
 Module[{n, sum = 0, start = 0, startm = 0},
  If[m <= 3, Return[k[m]]];
  Do[start = start + 2/Im[N[ZetaZero[j], d]]^(m + 1);
   startm = startm + 2/Im[N[ZetaZero[j], d]]^(m), {j, 1, 10^2}];
  k[m] = start + (2/
       Im[N[ZetaZero[10^2 + 1], d]]) (ComputeKmWorker[m - 1, k, d] - 
       startm);
  Return[k[m]]]

(*--------------------------------------------------------------*)

Clear[PsiSmallXEpsilon];
PsiSmallXEpsilon::argx = "the function takes three arguments";
PsiSmallXEpsilon::arg1 = 
  "the first argument is a positive real number";
PsiSmallXEpsilon::arg2 = "the second argument is a natural number";
PsiSmallXEpsilon::arg3 = "the third argument is a natural number";
PsiSmallXEpsilon[] := (Message[PsiSmallXEpsilon::argx]; Abort[]);
PsiSmallXEpsilon[argseq___] := Module[{args = List[argseq], b, m, d},
  If[Length[args] != 3, (Message[PsiSmallXEpsilon::argx]; Abort[])];
  b = args[[1]]; m = args[[2]]; d = args[[3]];
  If[Not[PositiveRealQ[b]], (Message[PsiSmallXEpsilon::arg1]; 
    Abort[])];
  If[Not[NaturalNumberQ[m]], (Message[PsiSmallXEpsilon::arg2]; 
    Abort[])];
  If[Not[NaturalNumberQ[d ]], (Message[PsiSmallXEpsilon::arg3]; 
    Abort[])];
  Return[PsiSmallXEpsilonWorker[b, m, d]]]

PsiSmallXEpsilonWorker[b_, m_, d_] := 
 Module[{km, del, B, x, R, LogHO2Pi, num, den, xORLogH, eps, epsdash, 
   P, Q, a1, a2, a3, j, H, h},
  km = Last[ComputeKm[m, d]];
  R = $FordR;
  H = $H;
  LogHO2Pi = N[Log[H/(2 Pi)], d];
  x = N[Exp[b], d];
  xORLogH = N[x^(1/(R Log[H])), d];
  a1 = N[112/1000, d];
  a2 = 0 N[278/1000, d];
  a3 = N[2510/1000, d];
  P[h_] := (a1 Log[h] + a2 Log[Log[h]] + a3) 2/h^(m + 1);
  Q[h_] := 1/(2 Pi) + (a1 Log[h] + a2)/(h Log[h] Log[h/(2 Pi)]);
  num = (1 + m LogHO2Pi);
  den = m^2 H^m xORLogH (1 - num*b/(R m^2 Log[H]^2 LogHO2Pi));
  B = km/Sqrt[x] + P[H]/xORLogH + Q[H] num/den;
  If[B <= 0, Return[0]];
  del = 2 B^(1/(m + 1));
  eps = (del/
      2) ((1/2^m) Sum[Binomial[m, j] (1 + j del)^(m + 1), {j, 0, m}] +
       m);
  If[Not[2 b < R m Log[H]^2 &&  x > 1 + m del x], Print[False]];
  epsdash = 
   Max[eps + Log[2 Pi] /x - 1/(2 x^2), eps - Log[2 Pi]/x + 1/x^2];
  If[2 b < R m Log[H]^2 &&  x > 1 + m del x, Return[epsdash], 
   Return[0]]]

(*------------------------------------------------------------------*)

Clear[PsiLargeXEpsilon];
PsiLargeXEpsilon::argx = "the function takes three arguments";
PsiLargeXEpsilon::arg1 = 
  "the first argument is a positive real number";
PsiLargeXEpsilon::arg2 = 
  "the second argument is a positive real number not less than 16";
PsiLargeXEpsilon::arg3 = "the third argument is a natural number"; 
PsiLargeXEpsilon[] := (Message[PsiLargeXEpsilon::argx]; Abort[]);
PsiLargeXEpsilon[argseq___] := Module[{args = List[argseq], x, R, d},
  If[Length[args] != 3, (Message[PsiLargeXEpsilon::argx]; Abort[])];
  x = args[[1]]; R = args[[2]]; d = args[[3]];
  If[Not[PositiveRealQ[x ]], (Message[PsiLargeXEpsilon::arg1]; 
    Abort[])];
  If[Not[PositiveRealQ[R] && R >= 16], (Message[
     PsiLargeXEpsilon::arg2]; Abort[])];
  If[Not[NaturalNumberQ[d ]], (Message[PsiLargeXEpsilon::arg3]; 
    Abort[])];
  Return[PsiLargeXEpsilonWorker[x, R, d]]]
PsiLargeXEpsilonWorker[x_, R_, d_] := Module[{lx}, (* R>= 16 *)
  lx = N[Log[x], d];
  Return[N[Sqrt[lx] Exp[-Sqrt[lx/R]], d]]]

(*------------------------------------------------------------------*)

Clear[PsiLargeXEpsilon2];
PsiLargeXEpsilon2::argx = "the function takes three arguments";
PsiLargeXEpsilon2::arg1 = 
  "the first argument is a positive real number";
PsiLargeXEpsilon2::arg2 = 
  "the second argument is a positive real number in [8,18]";
PsiLargeXEpsilon2::arg3 = "the third argument is a natural number"; 
PsiLargeXEpsilon2[] := (Message[PsiLargeXEpsilon2::argx]; Abort[]);
PsiLargeXEpsilon2[argseq___] := Module[{args = List[argseq], x, R, d},
  If[Length[args] != 3, (Message[PsiLargeXEpsilon2::argx]; Abort[])];
  x = args[[1]]; R = args[[2]]; d = args[[3]];
  If[Not[PositiveRealQ[x ]], (Message[PsiLargeXEpsilon2::arg1]; 
    Abort[])];
  If[Not[PositiveRealQ[R] && 8 <= R <= 18], (Message[
     PsiLargeXEpsilon2::arg2]; Abort[])];
  If[Not[NaturalNumberQ[d ]], (Message[PsiLargeXEpsilon2::arg3]; 
    Abort[])];
  Return[PsiLargeXEpsilon2Worker[x, R, d]]]
PsiLargeXEpsilon2Worker[x_, R_, d_] := Module[{lx}, (* R>= 16 *)
  lx = N[Log[x], d];
  Return[N[Sqrt[2 *lx] Exp[-Sqrt[lx/R]], d]]]

(*------------------------------------------------------------------*)

Clear[VerifyThetaX];
VerifyThetaX::argx = "the function takes three arguments";
VerifyThetaX::arg1 = "the first argument is a positive real number";
VerifyThetaX::arg2 = "the second argument is a positive real number";
VerifyThetaX::arg2 = "the third argument is a natural number"; 
VerifyThetaX[] := (Message[VerifyThetaX::argx]; Abort[]);
VerifyThetaX[argseq___] := Module[{args = List[argseq], a, b, d},
  If[Length[args] != 3, (Message[VerifyThetaX::argx]; Abort[])];
  a = args[[1]]; b = args[[2]]; d = args[[3]];
  If[Not[PositiveRealQ[a ]], (Message[VerifyThetaX::arg1]; Abort[])];
  If[Not[PositiveRealQ[b]], (Message[VerifyThetaX::arg2]; Abort[])];
  If[Not[NaturalNumberQ[d ]], (Message[VerifyThetaX::arg3]; Abort[])];
  Return[VerifyThetaXWorker[a, b, d]]]
VerifyThetaXWorker[a_, b_, d_] := 
 Module[{t = N[Log[2], d], n, max, bad = {}, diff},
  max = t/2;
  n = 3;
  While[n <= b,
   If[Mod[n, 10^9] == 0, Print[{n, max, diff}]];
   If[PrimeQ[n], t = t + N[Log[n], d]];
   max = Max[max, t/n];
   diff = n - t;
   If[diff < N[10^(-8), IntegerPart[d/2]], Print[{n, diff}]; 
    bad = Append[bad, n]];
   n = n + 1];
  Return[bad]]

(*------------------------------------------------------------------*)

Clear[ComputeTheta];
ComputeTheta::argx = "the function takes two arguments";
ComputeTheta::arg1 = "the first argument is a positive real number";
ComputeTheta::arg2 = "the second argument is a natural number"; 
ComputeTheta[] := (Message[ComputeTheta::argx]; Abort[]);
ComputeTheta[argseq___] := Module[{args = List[argseq], x, d},
  If[Length[args] != 2, (Message[ComputeTheta::argx]; Abort[])];
  x = args[[1]]; d = args[[2]]; 
  If[Not[PositiveRealQ[x ]], (Message[ComputeTheta::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d ]], (Message[ComputeTheta::arg2]; Abort[])];
  Return[ComputeThetaWorker[x, d]]]
ComputeThetaWorker[x_, d_] := Module[{sum = N[0, d], n},
  Do[If[PrimeQ[n], sum = sum + N[Log[n], d]], {n, 2, x}];
  Return[sum]]

(*------------------------------------------------------------------*)

Clear[ComputeThetaValues];
ComputeThetaValues::argx = "the function takes two arguments";
ComputeThetaValues::arg1 = "the first argument is a natural number";
ComputeThetaValues::arg2 = "the second argument is a natural number"; 
ComputeThetaValues[] := (Message[ComputeThetaValues::argx]; Abort[]);
ComputeThetaValues[argseq___] := Module[{args = List[argseq], nx, d},
  If[Length[args] != 2, (Message[ComputeThetaValues::argx]; Abort[])];
  nx = args[[1]]; d = args[[2]]; 
  If[Not[NaturalNumberQ[nx ]], (Message[ComputeThetaValues::arg1]; 
    Abort[])];
  If[Not[NaturalNumberQ[d ]], (Message[ComputeThetaValues::arg2]; 
    Abort[])];
  Return[ComputeThetaValuesWorker[nx, d]]]
ComputeThetaValuesWorker[nx_, d_] :=  
 Module[{sum = N[0, d], n, ans = {}},
  Do[If[PrimeQ[n], sum = sum + N[Log[n], d]]; 
   ans = Append[ans, sum], {n, 2, nx}];
  Return[ans]]

(*------------------------------------------------------------------*)

Clear[Psi];
Psi::argx = "the function takes two arguments";
Psi::arg1 = "the first argument is a positive real number";
Psi::arg2 = "the second argument is a natural number";
Psi[] := (Message[Psi::argx]; Abort[]);
Psi[argseq___] := Module[{args = List[argseq], x, d},
  If[Length[args] != 2, (Message[Psi::argx]; Abort[])];
  x = args[[1]]; d = args[[2]]; 
  If[Not[PositiveRealQ[x ]], (Message[Psi::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[Psi::arg2]; Abort[])];
  Return[PsiWorker[x, d]]]
PsiWorker[x_, d_] := Module[{b, n, sum = N[0, d]},
  b = IntegerPart[Log[x]/N[Log[2], d]];
  Do[sum = sum + ComputeThetaWorker[x^(1/n), d], {n, 1, b}];
  Return[sum]]

(* Chapter 4: Classical \
equivalences================================== *)

(*------------------------------------------------------------------*)

Clear[M];
M::argx = "the function takes one argument";
M::arg1 = "the argument is a positive real number";
M[] := (Message[M::argx]; Abort[]);
M[argseq___] := Module[{args = List[argseq], x, d},
  If[Length[args] != 1, (Message[M::argx]; Abort[])];
  x = args[[1]]; 
  If[Not[PositiveRealQ[x ]], (Message[M::arg1]; Abort[])];
  Return[MWorker[x]]]
MWorker[x_] := Module[{n}, Return[Sum[MoebiusMu[n], {n, 1, x}]]]

(* Chapter 5: Euler's totient function=================================*)


Clear[NOverPhi];
NOverPhi::argx = "the function takes two arguments";
NOverPhi::arg1 = "the first argument is a natural number";
NOverPhi::arg2 = "the second argument is a natural number";
NOverPhi[] := (Message[NOverPhi::argx]; Abort[]);
NOverPhi[argseq___] := Module[{args = List[argseq], n, d},
  If[Length[args] != 2, (Message[NOverPhi::argx]; Abort[])];
  n = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[n]], (Message[NOverPhi::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[NOverPhi::arg2]; Abort[])];
  Return[NOverPhiWorker[n, d]]]
NOverPhiWorker[n_, d_] := Module[{iprod = 1, prod = N[1, d], j, ps, p},
  If[n == 1, Return[N[1, d]]];
  ps = Map[First, FactorInteger[n]];
  Do[p = ps[[j]];
   iprod = iprod* p/(p - 1), {j, 1, Length[ps]}];
  Return[N[iprod, d]]]

(*------------------------------------------------------------------*)

Clear[NicolasInequalityQ];
NicolasInequalityQ::argx = "the function takes two arguments";
NicolasInequalityQ::arg1 = "the first argument is a natural number";
NicolasInequalityQ::arg2 = "the second argument is a natural number";
NicolasInequalityQ[] := (Message[NicolasInequalityQ::argx]; Abort[]);
NicolasInequalityQ[argseq___] := Module[{args = List[argseq], n, d},
  If[Length[args] != 2, (Message[NicolasInequalityQ::argx]; Abort[])];
  n = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[n]], (Message[NicolasInequalityQ::arg1]; 
    Abort[])];
  If[Not[NaturalNumberQ[d ]], (Message[NicolasInequalityQ::arg2]; 
    Abort[])];
  Return[NicolasInequalityQWorker[n, d]]]
NicolasInequalityQWorker[n_, d_] := Module[{lhs, rhs, flag, eg, e},
  eg = Exp[EulerGamma];
  e = 2 d;
  lhs = NOverPhiWorker[n, d];
  rhs = N[eg Log[Log[n]], d];
  flag = If[lhs + N[10^(-e) , d] < rhs, True, False];
  Return[flag]]

(*------------------------------------------------------------------*)



Clear[PlotEulerPhiRatio];
PlotEulerPhiRatio::argx = "the function takes one argument";
PlotEulerPhiRatio::arg1 = 
  "the first argument is a positive real number";
PlotEulerPhiRatio[] := (Message[PlotEulerPhiRatio::argx]; Abort[]);
PlotEulerPhiRatio[argseq___] := Module[{args = List[argseq], b},
  If[Length[args] != 1, (Message[PlotEulerPhiRatio::argx]; Abort[])];
  b = args[[1]]; 
  If[Not[PositiveRealQ[b ]], (Message[PlotEulerPhiRatio::arg1]; 
    Abort[])];
  Return[PlotEulerPhiRatioWorker[b]]]
PlotEulerPhiRatioWorker[b_] := Module[{ans = {}, n, up = 0}, 
  Do[If[Exp[EulerGamma ] Log[Log[n]] < n/EulerPhi[n], up = up + 1]; 
   ans = Append[ans, up /n], {n, 2, b}];
  Print[ListLinePlot[ans, ImageSize -> 400, PlotRange -> {0, 0.2}, 
    AxesStyle -> Thick, AxesLabel -> {"n", "count"}]];
  Return[N[up/b]]]


(*------------------------------------------------------------------*)

Clear[Computef];
Computef::argx = "the function takes two arguments";
Computef::arg1 = "the first argument is a positive real number";
Computef::arg2 = "the second argument is a natural number";
Computef[] := (Message[Computef::argx]; Abort[]);
Computef[argseq___] := Module[{args = List[argseq], x, d},
  If[Length[args] != 2, (Message[Computef::argx]; Abort[])];
  x = args[[1]]; d = args[[2]]; 
  If[Not[PositiveRealQ[x ]], (Message[Computef::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[Computef::arg2]; Abort[])];
  Return[ComputefWorker[x, d]]]
ComputefWorker[x_, d_] := Module[{eg, th, prprod = N[1, d], j, p, n},
  eg = N[Exp[EulerGamma], d];
  th = ComputeThetaWorker[x, d];
  Do[If[PrimeQ[n], p = n; prprod = prprod* N[(p - 1)/p, d]], {n, 2, 
    x}];
  Return[eg* Log[th] * prprod]]

(*------------------------------------------------------------------*)

Clear[NicolasTwoFailures];
NicolasTwoFailures::argx = "the function takes three arguments";
NicolasTwoFailures::arg1 = "the first argument is a natural number";
NicolasTwoFailures::arg2 = "the second argument is a natural number";
NicolasTwoFailures::arg3 = "the third argument is a natural number";
NicolasTwoFailures[] := (Message[NicolasTwoFailures::argx]; Abort[]);
NicolasTwoFailures[argseq___] := 
 Module[{args = List[argseq], k0, k1, d},
  If[Length[args] != 3, (Message[NicolasTwoFailures::argx]; Abort[])];
  k1 = args[[1]]; k0 = args[[2]]; d = args[[3]];
  If[Not[NaturalNumberQ[k1 ]], (Message[NicolasTwoFailures::arg1]; 
    Abort[])];
  If[Not[NaturalNumberQ[k0]], (Message[NicolasTwoFailures::arg2]; 
    Abort[])];
  If[Not[NaturalNumberQ[d ]], (Message[NicolasTwoFailures::arg3]; 
    Abort[])];
  Return[NicolasTwoFailuresWorker[k1, k0, d]]]
NicolasTwoFailuresWorker[k1_, k0_, d_] := 
 Module[{k, rk, Nk, LogNk, cNk, failures = 0, eg, c},
  eg = N[Exp[EulerGamma], d];
  c = N[eg (4 + EulerGamma - Log[Pi] - 2 Log[2]), d];
  rk = Product[N[Prime[k]/(Prime[k] - 1), d], {k, 1, k1 - 1}];
  LogNk = Sum[N[Log[Prime[k]], d], {k, 1, k1 - 1}];
  Do[rk = rk* Prime[k]/(Prime[k] - 1);
   LogNk = LogNk + N[Log[Prime[k]], d];
   cNk = (rk - eg Log[LogNk]) Sqrt[LogNk];
   If[cNk > c, failures = failures + 1], {k, k1, k0}];
  Return[failures]]

(*------------------------------------------------------------------*)

Clear[PlotCNk];
PlotCNk::argx = "the function takes three arguments";
PlotCNk::arg1 = "the first argument is a natural number";
PlotCNk::arg2 = "the second argument is a natural number";
PlotCNk::arg3 = "the third argument is a natural number";
PlotCNk[] := (Message[PlotCNk::argx]; Abort[]);
PlotCNk[argseq___] := Module[{args = List[argseq], k1, k0, d},
  If[Length[args] != 3, (Message[PlotCNk::argx]; Abort[])];
  k1 = args[[1]]; k0 = args[[2]]; d = args[[3]];
  If[Not[NaturalNumberQ[k1]], (Message[PlotCNk::arg1]; Abort[])];
  If[Not[NaturalNumberQ[k0]], (Message[PlotCNk::arg2]; Abort[])];
  If[k0 < k1, 
   Print["the second argument should be greater than the first"]; 
   Abort[]];
  If[Not[NaturalNumberQ[d ]], (Message[PlotCNk::arg3]; Abort[])];
  Return[PlotCNkWorker[k1, k0, d]]]
PlotCNkWorker[k1_, k0_, d_] := 
 Module[{ans = {}, k, rk, Nk, LogNk, cNk, failures = 0, eg, c},
  eg = N[Exp[EulerGamma], d];
  c = N[eg (4 + EulerGamma - Log[Pi] - 2 Log[2]), d];
  rk = Product[N[Prime[k]/(Prime[k] - 1), d], {k, 1, k1 - 1}];
  LogNk = Sum[N[Log[Prime[k]], d], {k, 1, k1 - 1}];
  Do[rk = rk* Prime[k]/(Prime[k] - 1);
   LogNk = LogNk + N[Log[Prime[k]], d];
   cNk = N[(rk - eg Log[LogNk]) Sqrt[LogNk], d];
   ans = Append[ans, cNk],
    {k, k1, k0}];
  Print[ListLinePlot[ans, ImageSize -> 500, AxesStyle -> Thick ,   
    DataRange -> {k1, k0}, AxesLabel -> {"k", "c(Nk)"}]];
  Return[Last[ans]]]

(*------------------------------------------------------------------*)

Clear[CheckCNk];
CheckCNk::argx = "the function takes two arguments";
CheckCNk::arg1 = "the first argument is a natural number";
CheckCNk::arg2 = "the second argument is a natural number";
CheckCNk[] := (Message[CheckCNk::argx]; Abort[]);
CheckCNk[argseq___] := Module[{args = List[argseq], k0, d},
  If[Length[args] != 2, (Message[CheckCNk::argx]; Abort[])];
  k0 = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[k0 ]], (Message[CheckCNk::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[CheckCNk::arg2]; Abort[])];
  Return[CheckCNkWorker[k0, d]]]
CheckCNkWorker[k0_, d_] := 
 Module[{k, rk, Nk, LogNk, cNk, failures = {}, eg, cN1n, cN66n},
  eg = N[Exp[EulerGamma], d];
  c = N[eg (4 + EulerGamma - Log[Pi] - 2 Log[2]), d];
  rk = N[Prime[1]/(Prime[1] - 1), d];
  LogNk = N[Log[Prime[1]], d];
  cN1n = cN1[d];
  cN66n = cN66[d];
  Do[If[Mod[k, 10^6] == 0, Print["->", k]];
   rk = rk* Prime[k]/(Prime[k] - 1);
   LogNk = LogNk + N[Log[Prime[k]], d];
   cNk = (rk - eg Log[LogNk]) Sqrt[LogNk];
   If[cNk > cN66n || cNk < cN1n, Print[k]; 
    failures = Append[failures, k]], {k, 2, k0}];
  Return[failures]]

(*-------------------------------------------------------------*)

psitheta[x_] := Module[{sum = 0.0, n},(* worker *)
  Do[sum = 
    sum + If[
      PrimeQ[n], (IntegerPart[N[Log[x]/Log[n]]] - 1.0) N[Log[n]], 
      0], {n, 2, x}];
  Return[sum]]

(*------------------------------------------------------------------*)

Clear[ComputeDeltas];
ComputeDeltas::argx = "the function takes one argument";
ComputeDeltas::arg1 = "the first argument is a natural number";
ComputeDeltas[] := (Message[ComputeDeltas::argx]; Abort[]);
ComputeDeltas[argseq___] := Module[{args = List[argseq], x, d},
  If[Length[args] != 1, (Message[ComputeDeltas::argx]; Abort[])];
  d = args[[1]]; 
  If[Not[NaturalNumberQ[d]], (Message[ComputeDeltas::arg1]; Abort[])];
  Return[ComputeDeltasWorker[d]]]
ComputeDeltasWorker[d_] := 
 Module[{ans = {}, i, q, j, Delta, delta, min = N[100, d]},
  Do[q = Prime[j]^i;
   If[q <= 599^3, ans = Append[ans, {q, Prime[j], i}]], {i, 2, 
    80}, {j, 1, 2000}];
  ans = Sort[ans];
  Print[Table[ans[[i]], {i, 1, 11}]]; Delta[1] = N[2 Log[2], d];
  delta[1] = Delta[1]/Sqrt[N[8, d]];
  Print[Length[ans]];
  Do[Delta[i] = Delta[i - 1] + N[Log[(ans[[i]])[[2]]], d];
   delta[i] = Delta[i]/Sqrt[First[ans[[i + 1]]]];
   If[i > 10, min = Min[min, delta[i]]], {i, 2, Length[ans] - 1}];
  Print[Table[delta[i], {i, 1, 11}]];
  Return[min]]

(*------------------------------------------------------------------*)


Clear[FHalf];
FHalf::argx = "the function takes two arguments";
FHalf::arg1 = "the first argument is a positive real number";
FHalf::arg2 = "the second argument is a natural number";
FHalf[] := (Message[FHalf::argx]; Abort[]);
FHalf[argseq___] := Module[{args = List[argseq], x, d},
  If[Length[args] != 2, (Message[FHalf::argx]; Abort[])];
  x = args[[1]]; d = args[[2]];
  If[Not[PositiveRealQ[x ]], (Message[FHalf::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[FHalf::arg2]; Abort[])];
  Return[FHalfWorker[x, d]]]
FHalfWorker[x_, d_] := 
 N[2/(Sqrt[x] Log[x]) - ExpIntegralEi[-Log[x]/2 ], d]

(*------------------------------------------------------------------*)

Clear[ComputeC];
ComputeC::argx = "the function takes two arguments";
ComputeC::arg1 = "the first argument is a natural number";
ComputeC::arg2 = "the second argument is a natural number"; 
ComputeC[] := (Message[ComputeC::argx]; Abort[]);
ComputeC[argseq___] := Module[{args = List[argseq], n, d},
  If[Length[args] != 2, (Message[ComputeC::argx]; Abort[])];
  n = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[n]], (Message[ComputeC::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[ComputeC::arg2]; Abort[])]; 
  Return[ComputeCWorker[n, d]]]
ComputeCWorker[n_, d_] := 
 N[(n/EulerPhi[n] - Exp[EulerGamma] Log[Log[n]]) Sqrt[Log[n]], d]

(*------------------------------------------------------------------*)


PlotAsymptoticNsimp[sup_] := Module[{ans = {}, p, n}, (* worker *)
  Do[p = First[Last[FactorInteger[sup[[n]]]]];
   ans = Append[ans , N[Log[sup[[n]]]/p]], {n, 2, Length[sup]}];
  Print[ListLinePlot[ans, ImageSize -> 400]];
  Return[Null]]

cN66[d_] := Module[{k, rk, Nk, LogNk, cNk, eg, c},(* worker *)
  eg = N[Exp[EulerGamma], d];
  rk = N[Prime[1]/(Prime[1] - 1), d];
  LogNk = N[Log[Prime[1]], d];
  Do[rk = rk* Prime[k]/(Prime[k] - 1);
   LogNk = LogNk + N[Log[Prime[k]], d];
   cNk = (rk - eg Log[LogNk]) Sqrt[LogNk];
   , {k, 2, 66}];
  Return[cNk]]

cN1[d_] := Module[{k, rk, Nk, LogNk, cNk, eg}, (* worker *)
  eg = N[Exp[EulerGamma], d];
  rk = N[Prime[1]/(Prime[1] - 1), d];
  LogNk = N[Log[Prime[1]], d];
  cNk = (rk - eg Log[LogNk]) Sqrt[LogNk];
  Return[cNk]]


(* Chapter 6: A variety of abundant \
numbers===============================*)


(*------------------------------------------------------------------*)

Clear[HighlyCompositeQ];
HighlyCompositeQ::argx = "the function takes one argument";
HighlyCompositeQ::arg1 = "the argument is a natural number";
HighlyCompositeQ[] := (Message[HighlyCompositeQ::argx]; Abort[]);
HighlyCompositeQ[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[HighlyCompositeQ::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[HighlyCompositeQ::arg1]; 
    Abort[])];
  Return[HighlyCompositeQWorker[n]]]
HighlyCompositeQWorker[n_] := Module[{big, j, flag = False},
  big = Last[InitialHC];
  If[n <= big , Return[If[MemberQ[InitialHC, n], True, False]]];
  Do[If[DivisorSigma[0, n] <= DivisorSigma[0, j], flag = False; 
    Break[]], {j, big, n - 1}];
  Return[flag]]

(*------------------------------------------------------------------*)

Clear[HardyRamanujanQ];
HardyRamanujanQ::argx = "the function takes one argument";
HardyRamanujanQ::arg1 = "the argument is a natural number";
HardyRamanujanQ[] := (Message[HardyRamanujanQ::argx]; Abort[]);
HardyRamanujanQ[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[HardyRamanujanQ::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[HardyRamanujanQ::arg1]; 
    Abort[])];
  Return[HardyRamanujanQWorker[n]]]
HardyRamanujanQWorker[n_] := 
 Module[{fs, p, e, nextp, nexte, flag = True},
  If[n == 1, Return[False]];
  fs = FactorInteger[n];
  {p, e} = fs[[1]];
  If[p > 2, Return[False]];
  Do[{nextp, nexte} = fs[[j]];
   If[nextp > Prime[j], flag = False; Break[]];
   If[nexte > e, flag = False; Break[]]; e = nexte; 
   p = nextp, {j, 2, Length[fs]}];
  Return[flag]]

(*----------------------------------------------------------------------------*)

Clear[GenerateHR];
GenerateHR::argx = "the function takes one argument";
GenerateHR::arg1 = "the argument is a natural number";
GenerteHR[] := (Message[GenerateHR::argx]; Abort[]);
GenerateHR[argseq___] := Module[{args = List[argseq], b},
  If[Length[args] != 1, (Message[GenerateHR::argx]; Abort[])];
  b = args[[1]]; 
  If[Not[NaturalNumberQ[b]], (Message[GenerateHR::arg1]; Abort[])];
  Return[GenerateHRWorker[b]]]
GenerateHRWorker[b_] := Module[{ans = {}, n},
  Do[If[HardyRamanujanQWorker[n], ans = Append[ans, n]], {n, 1, b}];
  Return[ans]]
(*-------------------------------------------------*)



Clear[HardyRamanujanToPrimorialForm];
HardyRamanujanToPrimorialForm::argx = "the function takes one argument";
HardyRamanujanToPrimorialForm::arg1 = "the argument is a Hardy Ramanujan number";
HardyRamanujanToPrimorialForm[] := (Message[HardyRamanujanToPrimorialForm::argx]; Abort[]);
HardyRamanujanToPrimorialForm[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[HardyRamanujanToPrimorialForm::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[HardyRamanujanQ[n]], (Message[HardyRamanujanToPrimorialForm::arg1]; Abort[])];
  Return[HR2PrlWorker[n]]];

H2PStep[n_, prs_] := Module[{fs, p, e, ppi, pri, m},
  fs = FactorInteger[n];
  {p, e} = Last[fs];
  ppi = PrimePi[p];
  pri = Product[Prime[m], {m, 1, ppi}];
  Return[{n/pri^e, Append[prs, {p, e}]}]];

HR2PrlWorker[n_] := Module[{m, ans = {}},
  m = n;
  While[m > 1, {m, ans} = H2PStep[m , ans]];
  Return[HR2PrlWorker2[ans]]];
  
HR2PrlWorker2[lis_]:= Module[{ans={},p,m,i,j}, 
   Do[{p,m}=lis[[j]];
      ans=Join[ans, Table[p,{i,1,m}]],{j,1,Length[lis]}];
   Return[ans]];
(*-------------------------------------------------*)
Clear[PrimorialFormToInteger]

PrimorialFormToInteger::argx = "the function takes one argument";
PrimorialFormToInteger::arg1 = "the argument is a list of primes in non-increasing order";
PrimorialFormToInteger[argseq___]:= Module[{args=List[argseq],lis},
 If[Length[args] != 1, (Message[PrimorialFormToInteger::argx]; Abort[])];
  lis = args[[1]]; 
  If[Not[PrimorialFormQ[lis]], (Message[PrimorialFormToInteger::arg1]; Abort[])];
Return[PrimorialFormToIntegerWorker[lis]]];

PrimorialFormToIntegerWorker[lis_]:= Apply[Times, Map[prime2primorial,lis]];

prime2primorial[p_]:= Module[{prod=1,n}, 
   Do[If[PrimeQ[n], prod=prod*n], {n,2,p}];
   Return[prod]];
(*-------------------------------------------------*)

Clear[PrimorialFormQ]
PrimorialFormQ::argx = "the function takes one argument";
PrimorialFormQ::arg1 = "the argument is a list of primes in non-increasing order";
PrimorialFormQ[argseq___]:= Module[{args=List[argseq],lis},
 If[Length[args] != 1, (Message[PrimorialFormQ::argx]; Abort[])];
  lis = args[[1]]; 
  If[Not[NaturalNumberListQ[lis] && NonIncreasingQ[lis]], (Message[PrimorialFormQ::arg1]; Abort[])];
Return[PrimorialFormQWorker[lis]]]

PrimorialFormQWorker[lis_]:=Module[{flag=True,j},
If[Not[ListQ[lis] && NonIncreasingQ[lis]], Return[False]];
Do[If[Not[PrimeQ[lis[[j]]]], flag=False; Break[]], {j,1,Length[lis]}];
Return[flag]];

NonIncreasingQ[lis_]:= Sort[Reverse[lis]]== Reverse[lis];
(*--------------------------------------------------*)

Clear[SuperAbundantQ];
SuperAbundantQ::argx = "the function takes one argument";
SuperAbundantQ::arg1 = "the argument is a natural number";
SuperAbundantQ[] := (Message[SuperAbundantQ::argx]; Abort[]);
SuperAbundantQ[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[SuperAbundantQ::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[SuperAbundantQ::arg1]; Abort[])];
  Return[SuperAbundantQWorker[n]]]
SuperAbundantQWorker[n_] := Module[{big, j, flag = False},
  (* check initial *)
  big = Last[InitialSA];
  If[n <= big , Return[If[MemberQ[InitialSA, n], True, False]]];
  (* check Hardy-Ramanujan *)
  If[Not[HardyRamanujanQ[n]], Return[False]];
  (* exhaustive check *)
  Do[If[j*DivisorSigma[1, n] <= n* DivisorSigma[1, j], flag = False; 
    Break[]], {j, big, n - 1}];
  Return[flag]]

(*------------------------------------------------------------------*)

Clear[ColossallyAbundantQ];
ColossallyAbundantQ::argx = "the function takes one argument";
ColossallyAbundantQ::arg1 = "the argument is a natural number";
ColossallyAbundantQ[] := (Message[ColossallyAbundantQ::argx]; Abort[]);
ColossallyAbundantQ[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[ColossallyAbundantQ::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[ColossallyAbundantQ::arg1]; 
    Abort[])];
  Return[ColossallyAbundantQWorker[n]]]
ColossallyAbundantQWorker[n_] := 
 Module[{fs, j, flag = True, p, a, ints = {}, F, u, v, big, nu},
  (* check initial *)
  big = Last[InitialCA];
  If[n <= big  , Return[If[MemberQ[InitialCA, n], True, False]]];
  (* check superabundant *)
  If[Not[SuperAbundantQ[n]], Return[False]];
  (* prime power check *)
  F[p_, a_] := N[Log[(p^(a + 1) - 1)/(p (p^a - 1))]/Log[p]];
  fs = FactorInteger[n];
  Do[{p, a} = fs[[j]];
   ints = Append[ints, {F[p, a + 1], F[p, a]}]
   , {j, 1, Length[fs]}];
  u = Max[Map[First, ints]];
  Print[ints];
  v = Min[Map[#[[2]] &, ints]];
  If[u > v, Return[False]];
  nu = ColossallyAbundantInteger[(u + v)/2];
  If[nu != n, Return[False]];
  Return[True]]

(*------------------------------------------------------------------*)

Clear[ColossallyAbundantInteger];
ColossallyAbundantInteger::argx = "the function takes one argument";
ColossallyAbundantInteger::arg1 = 
  "the argument is a positive real number";
ColossallyAbundantInteger[] := (Message[
    ColossallyAbundantInteger::argx]; Abort[]);
ColossallyAbundantInteger[argseq___] := 
 Module[{args = List[argseq], x, d},
  If[Length[args] != 1, (Message[ColossallyAbundantInteger::argx]; 
    Abort[])];
  e = args[[1]];
  If[Not[PositiveRealQ[e ]], (Message[
     ColossallyAbundantInteger::arg1]; Abort[])];
  Return[ColossallyAbundantIntegerWorker[e]]]
ColossallyAbundantIntegerWorker[e_] := 
 Module[{p, j, a, exps = {}, n = 1},
  Do[p = Prime[j];
   a = IntegerPart[Log[(p^(e + 1) - 1)/(p^e - 1)]/Log[p]] - 1; 
   exps = Append[exps, a];
   n = n*p^a;
   If[a == 0, Break[]], {j, 1, 100}];
  Return[{n, Take[exps, Length[exps] - 1]}]]

(*------------------------------------------------------------------*)

Clear[ColossallyAbundantF];
ColossallyAbundantF::argx = "the function takes three arguments";
ColossallyAbundantF::arg1 = 
  "the first argument is a real number greater than 1";
ColossallyAbundantF::arg2 = "the second argument is a natural number";
ColossallyAbundantF::arg2 = "the third argument is a natural number not less than $MinPrecision"; 
ColossallyAbundantF[] := (Message[ColossallyAbundantF::argx]; Abort[]);
ColossallyAbundantF[argseq___] := Module[{args = List[argseq], x, a},
  If[Length[args] != 3, (Message[ColossallyAbundantF::argx]; Abort[])];
  x = args[[1]]; a = args[[2]]; d = args[[3]];
  If[Not[RealNumberQ[x ] && x > 1], (Message[
     ColossallyAbundantF::arg1]; Abort[])];
  If[Not[NaturalNumberQ[a]], (Message[ColossallyAbundantF::arg2]; 
    Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[ColossallyAbundantF::arg3]; 
    Abort[])]; Return[ColossallyAbundantFWorker[x, a, d]]]
ColossallyAbundantFWorker[x_, a_, d_] := 
 N[Log[((x^(a + 1) - 1)/(x (x^a - 1)))]/Log[x], d]

(*------------------------------------------------------------------*)


Clear[CriticalEpsilonValues]; (* fix to get a contiguous group of \
values *)
CriticalEpsilonValues::argx = "the function takes three arguments";
CriticalEpsilonValues::arg1 = "the first argument is a natural number";
CriticalEpsilonValues::arg2 = 
  "the second argument is a natural number";
CriticalEpsilonValues[] := (Message[CriticalEpsilonValues::argx]; 
   Abort[]);
CriticalEpsilonValues[argseq___] := 
 Module[{args = List[argseq], b, jmax},
  If[Length[args] != 2, (Message[CriticalEpsilonValues::argx]; 
    Abort[])];
  b = args[[1]]; jmax = args[[2]]; 
  If[Not[PositiveRealQ[b ]], (Message[CriticalEpsilonValues::arg1]; 
    Abort[])];
  If[Not[NaturalNumberQ[jmax]], (Message[
     CriticalEpsilonValues::arg2]; Abort[])];
  Return[CriticalEpsilonValuesWorker[b, jmax]]]
CriticalEpsilonValuesWorker[b_, jmax_] := Module[{EE, F, p, a},
  F[p_, a_] := Log[((p^(a + 1) - 1)/(p (p^a - 1)))]/Log[p];
  EE[p_] := Table[F[p, a], {a, 1, b}];
  Return[Sort[
    Apply[ Join, Table[EE[Prime[j]], {j, 1, jmax}]], #1 > #2 &]]]



(* Chapter 7: Robin's theorem=====================================*)

Clear[RobinsInequalityQ];
RobinsInequalityQ::argx = "the function takes two arguments";
RobinsInequalityQ::arg1 = "the first argument is a natural number";
RobinsInequalityQ::arg2 = "the second argument is a natural number"; 
RobinsInequalityQ[] := (Message[RobinsInequalityQ::argx]; Abort[]);
RobinsInequalityQ[argseq___] := Module[{args = List[argseq], n, d},
  If[Length[args] != 2, (Message[RobinsInequalityQ::argx]; Abort[])];
  n = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[n]], (Message[RobinsInequalityQ::arg1]; 
    Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[RobinsInequalityQ::arg2]; 
    Abort[])]; Return[RobinsInequalityQWorker[n, d]]]
RobinsInequalityQWorker[n_, d_] := Module[{sn, snon},
  If[n == 1, Return[True]];
  sn = DivisorSigma[1, n];
  snon = N[sn/n, d];
  Return[snon < N[Exp[EulerGamma] Log[Log[n]], d]]]

(*------------------------------------------------------------------*)

Clear[PlotRobinsInequality];
PlotRobinsInequality::argx = "the function takes three arguments";
PlotRobinsInequality::arg1 = "the first argument is a natural number";
PlotRobinsInequality::arg2 = "the second argument is a natural number";
PlotRobinsInequality::arg3 = "the third argument is a natural number";
PlotRobinsInequality[] := (Message[PlotRobinsInequality::argx]; 
   Abort[]);
PlotRobinsInequality[argseq___] := 
 Module[{args = List[argseq], a, b, d},
  If[Length[args] != 3, (Message[PlotRobinsInequality::argx]; 
    Abort[])];
  a = args[[1]]; b = args[[2]]; d = args[[3]];
  If[Not[NaturalNumberQ[a]], (Message[PlotRobinsInequality::arg1]; 
    Abort[])];
  If[Not[NaturalNumberQ[b]], (Message[PlotRobinsInequality::arg2]; 
    Abort[])];
  If[Not[NaturalNumberQ[d ]], (Message[PlotRobinsInequality::arg3]; 
    Abort[])];
  If[a > b, 
   Print["the first argument should be less than the second"]; 
   Abort[]];
  Return[PlotRobinsInequalityWorker[a, b, d]]]
PlotRobinsInequalityWorker[a_, b_, d_] := 
 Module[{snon = {}, rhs = {}, n},
  Do[snon = Append[snon, {n, DivisorSigma[1, n]/n}];
   rhs = Append[rhs, {n, N[Exp[EulerGamma] Log[Log[n]], d]}], {n, a, 
    b}];
  Print[ListLinePlot[{snon, rhs}, AxesOrigin -> {a, 0}, 
    AxesStyle -> Thick, 
    AxesLabel -> {"n", "\[Sigma](n)/n,e^\[Gamma] loglog(n)"}, 
    ImageSize -> 500]];
  Return[{N[Last[snon][[2]], d], Last[rhs][[2]]}]]

(*------------------------------------------------------------------*)

Clear[SigmaOverN];
SigmaOverN::argx = "the function takes one argument";
SigmaOverN::arg1 = "the argument is a natural number";
SigmaOverN[] := (Message[SigmaOverN::argx]; Abort[]);
SigmaOverN[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[SigmaOverN::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[SigmaOverN::arg1]; Abort[])];
  Return[SigmaOverNWorker[n]]]
SigmaOverNWorker[n_] := Module[{prod = 1, fs, j},
  If[n == 1, Return[1]];
  If[PrimeQ[n], Return[(n + 1)/n]];
  fs = FactorInteger[n];
  Do[{p, e} = fs[[j]];
   prod = prod*((p^(e + 1) - 1)/(p^e (p - 1))), {j, 1, Length[fs]}];
  Return[prod]]
(*------------------------------------------------------------------*)


Clear[ComputeH];
ComputeH::argx = "the function takes two arguments";
ComputeH::arg1 = 
  "the first argument is a real number not less than 16";
ComputeH::arg2 = "the second argument is a natural number"; 
ComputeH[] := (Message[ComputeH::argx]; Abort[]);
ComputeH[argseq___] := Module[{args = List[argseq], x, d},
  If[Length[args] != 2, (Message[ComputeH::argx]; Abort[])];
  x = args[[1]]; d = args[[2]];
  If[Not[RealNumberQ[x] && x >= 16], (Message[ComputeH::arg1]; 
    Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[ComputeH::arg2]; Abort[])]; 
  Return[ComputeHWorker[x, d]]]
ComputeHWorker[x_, d_] := Module[{sum, j, k},
  sum = 1;
  k = IntegerPart[Log[x]/Log[2]];
  If[x <= 16, Return[N[0, d]]];
  Do[sum = sum + N[x^(1/j - 1/3), d], {j, 4, k}];
  Return[sum]]

(*------------------------------------------------------------------*)

Clear[LagariasInequalityQ];
LagariasInequalityQ::argx = "the function takes two arguments";
LagariasInequalityQ::arg1 = "the first argument is a natural number";
LagariasInequalityQ::arg2 = "the second argument is a natural number";
LagariasInequalityQ[] := (Message[LagariasInequalityQ::argx]; Abort[]);
LagariasInequalityQ[argseq___] := Module[{args = List[argseq], n, d},
  If[Length[args] != 2, (Message[LagariasInequalityQ::argx]; Abort[])];
  n = args[[1]]; d = args[[2]]; 
  If[Not[NaturalNumberQ[n ]], (Message[LagariasInequalityQ::arg1]; 
    Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[LagariasInequalityQ::arg2]; 
    Abort[])];
  Return[LagariasInequalityQWorker[n, d]]]
LagariasInequalityQWorker[n_, d_] := Module[{Hn, ds},
  Hn = NthHarmonicNumber[n, d];
  ds = DivisorSigma[1, n];
  Return[ds <=  N[Hn + Exp[Hn] Log[Hn], d]]]

(*------------------------------------------------------------------*)

Clear[UnitaryDivisorSigma];
UnitaryDivisorSigma::argx = "the function takes one argument";
UnitaryDivisorSigma::arg1 = "the argument is a natural number";
UnitaryDivisorSigma[] := (Message[UnitaryDivisorSigma::argx]; Abort[]);
UnitaryDivisorSigma[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[UnitaryDivisorSigma::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[UnitaryDivisorSigma::arg1]; 
    Abort[])];
  Return[UnitaryDivisorSigmaWorker[n]]]
UnitaryDivisorSigmaWorker[n_] := Module[{sum = 0, j},
  Do[If[Mod[n, j] == 0 && GCD[j, n/j] == 1, sum = sum + j], {j, 1, n}];
  Return[sum]]

(* Chapter 8: Numbers not satisfying Robin's \
inequality===================*)


SetA = {1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 18, 20, 24, 30, 36, 48, 
   60, 72, 84, 120, 180, 240, 360, 720, 840, 2520, 5040};
SetB = {2, 3, 5, 6, 10, 30};
SetC = {1, 4, 8, 9, 16, 36};
SetD = {4, 8, 9, 16, 36, 72, 108, 144, 216, 900, 1800, 2700, 3600, 
   44100, 88200};

(*------------------------------------------------------------------*)

Clear[KFreeQ];
KFreeQ::argx = "the function takes two arguments";
KFreeQ::arg1 = "the first argument is a natural number";
KFreeQ::arg2 = "the second argument is a natural number";
KFreeQ[] := (Message[KFreeQ::argx]; Abort[]);
KFreeQ[argseq___] := Module[{args = List[argseq], n, k},
  If[Length[args] != 2, (Message[KFreeQ::argx]; Abort[])];
  n = args[[1]]; k = args[[2]];
  If[Not[PositiveRealQ[n]], (Message[KFreeQ::arg1]; Abort[])];
  If[Not[NaturalNumberQ[k]], (Message[KFreeQ::arg2]; Abort[])];
  Return[KFreeQWorker[n, k]]]
KFreeQWorker[n_, k_] := Module[{flag = True, j, fs, p, e},
  If[n == 1, Return[True]];
  fs = FactorInteger[n];
  Do[{p, e} = fs[[j]]; 
   If[e >= k, flag = False; Break[]], {j, 1, Length[fs]}];
  Return[flag]]

(*------------------------------------------------------------------*)

Clear[NthPrimeBounds];
NthPrimeBounds::argx = "the function takes two arguments";
NthPrimeBounds::arg1 = 
  "the first argument is a natural number greater than 20";
NthPrimeBounds::arg2 = "the first argument is a natural number"; 
NthPrimeBounds[] := (Message[NthPrimeBounds::argx]; Abort[]);
NthPrimeBounds[argseq___] := Module[{args = List[argseq], n, d},
  If[Length[args] != 2, (Message[NthPrimeBounds::argx]; Abort[])];
  n = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[n] && n > 20 ], (Message[
     NthPrimeBounds::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[NthPrimeBounds::arg2]; 
    Abort[])]; Return[NthPrimeBoundsWorker[n, d]]]
NthPrimeBoundsWorker[n_, d_] := 
 N[{n (Log[n] + Log[Log[n]] - 3/2), n (Log[n] + Log[Log[n]] - 1/2)}, d]

(*------------------------------------------------------------------*)

Clear[LogKthPrimorialBounds];
LogKthPrimorialBounds::argx = "the function takes two arguments";
LogKthPrimorialBounds::arg1 = 
  "the first argument is a natural number greater than 197";
LogKthPrimorialBounds::arg2 = "the second argument is a natural \
number"; LogKthPrimorialBounds[] := (Message[
   LogKthPrimorialBounds::argx]; Abort[]);
LogKthPrimorialBounds[argseq___] := Module[{args = List[argseq], n, d},
  If[Length[args] != 2, (Message[LogKthPrimorialBounds::argx]; 
    Abort[])];
  n = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[n]], (Message[LogKthPrimorialBounds::arg1]; 
    Abort[])];
  If[Not[NaturalNumberQ[n] && n >= 198], (Message[
     LogKthPrimorialBounds::arg1]; Abort[])]; 
  Return[LogKthPrimorialBoundsWorker[n, d]]]
LogKthPrimorialBoundsWorker[n_, d_] := 
 N[{n (Log[n] + Log[Log[n]] - 1 + (Log[Log[n]] - 21454/10^4)/Log[n]), 
   n (Log[n] + Log[Log[n]] - 1 + (Log[Log[n]] - 2)/Log[n])}, d]

(*------------------------------------------------------------------*)
Clear[KthPrimorial];
KthPrimorial::argx = "the function takes one argument";
KthPrimorial::arg1 = "the argument is a natural number";
KthPrimorial[] := (Message[KthPrimorial::argx]; Abort[]);
KthPrimorial[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[KthPrimorial::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[KthPrimorial::arg1]; Abort[])];
  Return[KthPrimorialWorker[n]]]
KthPrimorialWorker[n_] := Module[{prod = 1, j},
  Do[prod = prod*Prime[j], {j, 1, n}];
  Return[prod]]

(*------------------------------------------------------------------*)
Clear[KthApproximatePrimorial];
KthApproximatePrimorial::argx = "the function takes two arguments";
KthApproximatePrimorial::arg1 = "the first argument is a natural number";
KthApproximatePrimorial::arg2 = "the second argument is a natural number";KthApproximatePrimorial[] := (Message[KthApproximatePrimorial::argx]; Abort[]);
KthApproximatePrimorial[argseq___] := Module[{args = List[argseq], n,d},
  If[Length[args] != 2, (Message[KthApproximatePrimorial::argx]; Abort[])];
  n = args[[1]]; d= args[[2]];
  If[Not[NaturalNumberQ[n]], (Message[KthApproximatePrimorial::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[KthApproximatePrimorial::arg2]; Abort[])]; 
  Return[KthApproximatePrimorialWorker[n,d]]]
KthApproximatePrimorialWorker[n_,d_] := Module[{prod = N[1,d], j},
  Do[prod = prod* N[Prime[j],d], {j, 1, n}];
  Return[prod]]

(* Chapter 9: Left, Right and Extremely abundant numbers==================*)


Clear[GronwallG];
GronwallG::argx = "the function takes two arguments";
GronwallG::arg1 = "the first argument is a natural number";
GronwallG::arg1 = "the second argument is a natural number"; 
GronwallG[] := (Message[GronwallG::argx]; Abort[]);
GronwallG[argseq___] := Module[{args = List[argseq], n, d},
  If[Length[args] != 2, (Message[GronwallG::argx]; Abort[])];
  n = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[n]], (Message[GronwallG::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[GronwallG::arg2]; Abort[])]; 
  Return[GronwallGWorker[n, d]]]
GronwallGWorker[n_, d_] :=  DivisorSigma[1, n]/N[n Log[Log[n]], d]

(*------------------------------------------------------------------*)
Clear[RightAbundantQ];
RightAbundantQ::argx = "the function takes three arguments";
RightAbundantQ::arg1 = 
  "the first argument is a composite natural number";
RightAbundantQ::arg2 = "the second argument is a natural number";
RightAbundantQ::arg3 = "the third argument is a natural number"; 
RightAbundantQ[] := (Message[RightAbundantQ::argx]; Abort[]);
RightAbundantQ[argseq___] := Module[{args = List[argseq], n, b},
  If[Length[args] != 3, (Message[RightAbundantQ::argx]; Abort[])];
  n = args[[1]]; b = args[[2]]; d = args[[3]];
  If[Not[NaturalNumberQ[n] && CompositeQ[n]], (Message[
     RightAbundantQ::arg1]; Abort[])];
  If[Not[NaturalNumberQ[b]], (Message[RightAbundantQ::arg2]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[RightAbundantQ::arg3]; 
    Abort[])]; Return[RightAbundantQWorker[n, b, d]]]
RightAbundantQWorker[n_, b_, d_] := Module[{flag = True, j},
  Do[If[N[DivisorSigma[1, n j] Log[Log[n]], d] > 
     N[j DivisorSigma[1, n] Log[Log[n j]], d], flag = False; 
    Break[]], {j, 2, b}];
  Return[flag]]

(*------------------------------------------------------------------*)
Clear[LeftAbundantQ];
LeftAbundantQ::argx = "the function takes two arguments";
LeftAbundantQ::arg1 = 
  "the first argument is a composite natural number";
LeftAbundantQ::arg2 = "the second argument is a natural number"; 
LeftAbundantQ[] := (Message[LeftAbundantQ::argx]; Abort[]);
LeftAbundantQ[argseq___] := Module[{args = List[argseq], n, d},
  If[Length[args] != 2, (Message[LeftAbundantQ::argx]; Abort[])];
  n = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[n] && CompositeQ[n]], (Message[
     LeftAbundantQ::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[LeftAbundantQ::arg2]; Abort[])];
   Return[LeftAbundantQWorker[n, d]]]
LeftAbundantQWorker[n_, d_] := Module[{ps, p, j, flag = True},
  If[Not[CompositeQWorker[n]], Return[False]];
  ps = Map[First, FactorInteger[n]];
  Do[p = ps[[j]];
   If[GronwallGWorker[n/p, d] > GronwallGWorker[n, d], flag = False; 
    Break[]], {j, 1, Length[ps]}];
  Return[flag]]

(*------------------------------------------------------------------*)
Clear[ExtraordinaryQ];
ExtraordinaryQ::argx = "the function takes three arguments";
ExtraordinaryQ::arg1 = 
  "the first argument is a composite natural number";
ExtraordinaryQ::arg2 = "the second argument is a natural number";
ExtraordinaryQ::arg3 = "the third argument is a natural number";
ExtraordinaryQ[] := (Message[ExtraordinaryQ::argx]; Abort[]);
ExtraordinaryQ[argseq___] := Module[{args = List[argseq], n, b, d},
  If[Length[args] != 3, (Message[ExtraordinaryQ::argx]; Abort[])];
  n = args[[1]]; b = args[[2]]; d = args[[3]];
  If[Not[NaturalNumberQ[n] && CompositeQ[n]], (Message[
     ExtraordinaryQ::arg1]; Abort[])];
  If[Not[NaturalNumberQ[b]], (Message[ExtraordinaryQ::arg2]; Abort[])];
  If[Not[NaturalNumberQ[d ]], (Message[ExtraordinaryQ::arg3]; 
    Abort[])];
  Return[ExtraordinaryQWorker[n, b, d]]]
ExtraordinaryQWorker[n_, b_, d_] := 
 LeftAbundantQWorker[n, d] && RightAbundantQWorker[n, b, d]

(*------------------------------------------------------------------*)
Clear[ExtremelyAbundantQ];
ExtremelyAbundantQ::argx = "the function takes two arguments";
ExtremelyAbundantQ::arg1 = 
  "the first argument is a natural number not less than 10080";
ExtremelyAbundantQ::arg2 = "the second argument is a natural number"; 
ExtremelyAbundantQ[] := (Message[ExtremelyAbundantQ::argx]; Abort[]);
ExtremelyAbundantQ[argseq___] := Module[{args = List[argseq], n, d},
  If[Length[args] != 2, (Message[ExtremelyAbundantQ::argx]; Abort[])];
  n = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[n] && n >= 10080], (Message[
     ExtremelyAbundantQ::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[ExtremelyAbundantQ::arg2]; 
    Abort[])]; Return[ExtremelyAbundantQWorker[n, d]]]
ExtremelyAbundantQWorker[n_, d_] := Module[{flag = True, j},
  Do[If[GronwallGWorker[j, d] > GronwallGWorker[n, d], flag = False; 
    Break[]], {j, 2, n - 1}];
  Return[flag]]

(* Chapter 10: Other equivalences to \
RH==================================*)

Clear[FareyFractions];
FareyFractions::argx = "the function takes one argument";
FareyFractions::arg1 = "the argument is a natural number";
FareyFractions[] := (Message[FareyFractions::argx]; Abort[]);
FareyFractions[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[FareyFractions::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[FareyFractions::arg1]; Abort[])];
  Return[FareyFractionsWorker[n]]]
FareyFractionsWorker[n_] := Union[Range[0, n]/n]

(*------------------------------------------------------------------*)
Clear[RedhefferMatrix];
RedhefferMatrix::argx = "the function takes one argument";
RedhefferMatrix::arg1 = "the argument is a natural number";
RedhefferMatrix[] := (Message[RedhefferMatrix::argx]; Abort[]);
RedhefferMatrix[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[RedhefferMatrix::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[RedhefferMatrix::arg1]; 
    Abort[])];
  Return[RedhefferMatrixWorker[n]]]
RedhefferMatrixWorker[n_] := Module[{ans, i, j},
  Do[ans[i, j] = If[j == 1 || Mod[j, i] == 0, 1, 0], {i, 1, n}, {j, 1,
     n}];
  Return[Table[ans[i, j], {i, 1, n}, {j, 1, n}]]]

(*------------------------------------------------------------------*)
Clear[DivisibilityGraph];
DivisibilityGraph::argx = "the function takes one argument";
DivisibilityGraph::arg1 = "the argument is a natural number";
DivisibilityGraph[] := (Message[DivisibilityGraph::argx]; Abort[]);
DivisibilityGraph[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[DivisibilityGraph::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[DivisibilityGraph::arg1]; 
    Abort[])];
  Return[DivisibilityGraphWorker[n]]]
DivisibilityGraphWorker[n_] := Module[{ans, i, j},
  Do[ans[i, j] = 
    If[(j == 1 && i != 1) || ( Mod[j, i] == 0 && i != j), i -> j, 
     0], {i, 1, n}, {j, 1, n}];
  Return[DeleteCases[
    Apply[Join, Table[ans[i, j], {i, 1, n}, {j, 1, n}]], 0]]]

(*------------------------------------------------------------------*)
Clear[PlotDivisibilityGraph];
PlotDivisibilityGraph::argx = "the function takes one argument";
PlotDivisibilityGraph::arg1 = "the argument is a natural number";
PlotDivisibilityGraph[] := (Message[PlotDivisibilityGraph::argx]; 
   Abort[]);
PlotDivisibilityGraph[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[PlotDivisibilityGraph::argx]; 
    Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[PlotDivisibilityGraph::arg1]; 
    Abort[])];
  Return[PlotDivisibilityGraphWorker[n]]]
PlotDivisibilityGraphWorker[n_] := 
 GraphPlot[DivisibilityGraphWorker[n], DirectedEdges -> True, 
  VertexLabeling -> True]

(*------------------------------------------------------------------*)
Clear[Xi];
Xi::argx = "the function takes two arguments";
Xi::arg1 = "the first argument is a real or complex number";
Xi::arg2 = "the second argument is a natural number";
Xi[] := (Message[Xi::argx]; Abort[]);
Xi[argseq___] := Module[{args = List[argseq], z, d},
  If[Length[args] != 2, (Message[Xi::argx]; Abort[])];
  z = args[[1]]; d = args[[2]]; 
  If[Not[NumberQ[z]], (Message[Xi::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[Xi::arg2]; Abort[])];
  Return[XiWorker[z, d]]]
XiWorker[s_, d_] := N[s (s - 1) Pi^(-s/2) Gamma[s/2] Zeta[s]/2, d]

(*------------------------------------------------------------------*)
Clear[Lambda];
Lambda::argx = "the function takes two arguments";
Lambda::arg1 = "the first argument is a real or complex number";
Lambda::arg2 = "the second argument is a natural number";
Lambda[] := (Message[Lambda::argx]; Abort[]);
Lambda[argseq___] := Module[{args = List[argseq], z, d},
  If[Length[args] != 2, (Message[Lambda::argx]; Abort[])];
  z = args[[1]]; d = args[[2]];
  If[Not[NumberQ[z]], (Message[Lambda::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[Lambda::arg2]; Abort[])];
  Return[LambdaWorker[z, d]]]
LambdaWorker[s_, d_] := N[2 (2 Pi)^(-s) Gamma[s] Cos[Pi s/2], d]

(*------------------------------------------------------------------*)
Clear[IntegerHeight];
IntegerHeight::argx = "the function takes one argument";
IntegerHeight::arg1 = "the argument is a natural number";
IntegerHeight[] := (Message[IntegerHeight::argx]; Abort[]);
IntegerHeight[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[IntegerHeight::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[IntegerHeight::arg1]; Abort[])];
  Return[IntegerHeightWorker[n]]]
IntegerHeightWorker[n_] := Module[{fs, j, sum = 0, p, e},
  If[n == 1, Return[0]];
  fs = FactorInteger[n];
  Do[{p, e} = fs[[j]];
   sum = sum + p^e, {j, 1, Length[fs]}];
  Return[sum]]

(*------------------------------------------------------------------*)

InitialLandau = {1, 1, 2, 3, 4, 6, 6, 12, 15, 20, 30, 30, 60, 60, 84, 
   105, 140, 210, 210, 420, 420, 420, 420, 840, 840, 1260, 1260, 1540,
    2310, 2520, 4620, 4620, 5460, 5460, 9240, 9240, 13860, 13860, 
   16380, 16380, 27720, 30030, 32760, 60060, 60060, 60060, 60060, 
   120120, 120120};


Clear[GenerateLandau];
GenerateLandau::argx = "the function takes one argument";
GenerateLandau::arg1 = "the argument is a natural number";
GenerateLandau[] := (Message[GenerateLandau::argx]; Abort[]);
GenerateLandau[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[GenerateLandau::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[GenerateLandau::arg1]; Abort[])];
  Return[GenerateLandauWorker[n]]]
GenerateLandauWorker[
  B_] := (*Jean-François Alcover,Mar 07 2014,after Alois P.Heinz*)
 Module[{a, b, n, i}, 
  b[n_, i_] := 
   b[n, i] = 
    Module[{p}, p = If[i < 1, 1, Prime[i]]; 
     If[n == 0 || i < 1, 1, 
      Max[b[n, i - 1], 
       Table[p^j*b[n - p^j, i - 1], {j, 1, Log[p, n] // Floor}]]]]; 
  a[n_] := b[n, 
    If[n < 8, 3, PrimePi[Ceiling[1.328*Sqrt[n*Log[n] // Floor]]]]]; 
  Return[Table[{n, a[n]}, {n, 0, B}] ]]

(*------------------------------------------------------------------*)
Clear[ZetaZeroCount];
ZetaZeroCount::argx = "the function takes two arguments";
ZetaZeroCount::arg1 = "the first argument is a positive real number";
ZetaZeroCount::arg2 = "the second argument is a natural number";
ZetaZeroCount[] := (Message[ZetaZeroCount::argx]; Abort[]);
ZetaZeroCount[argseq___] := Module[{args = List[argseq], T, d},
  If[Length[args] != 2, (Message[ZetaZeroCount::argx]; Abort[])];
  T = args[[1]]; d = args[[2]];
  If[Not[PositiveRealQ[T ]], (Message[ZetaZeroCount::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[ZetaZeroCount::arg2]; Abort[])];
  Return[ZetaZeroCountWorker[T, d]]]
ZetaZeroCountWorker[T_, d_] := Module[{cnt = 0, j = 1},
  While[Im[N[ZetaZero[j], d]] <= T, cnt++; j++];
  Return[cnt]]

(*------------------------------------------------------------------*)
Clear[S];
S::argx = "the function takes two arguments";
S::arg1 = "the first argument is a positive real number";
S::arg2 = "the second argument is a natural number";
S[] := (Message[S::argx]; Abort[]);
S[argseq___] := Module[{args = List[argseq], t, d},
  If[Length[args] != 2, (Message[S::argx]; Abort[])];
  t = args[[1]]; d = args[[2]];
  If[Not[PositiveRealQ[t ]], (Message[S::arg1]; Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[S::arg2]; Abort[])];
  Return[SWorker[t, d]]]
SWorker[t_, d_] := 
 ZetaZeroCountWorker[t, d] - N[RiemannSiegelTheta[t]/Pi - 1, d]

(* Appendix A: 
Tables================================================= *)

(*------------------------------------------------------------------*)

Clear[GenerateHC];
GenerateHC::argx = "the function takes one argument";
GenerateHC::arg1 = "the argument is a natural number";
GenerateHC[] := (Message[GenerateHC::argx]; Abort[]);
GenerateHC[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[GenerateHC::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[GenerateHC::arg1]; Abort[])];
  Return[GenerateHCWorker[n]]]
GenerateHCWorker[b_] := Module[{a = 0, d, ans = {}},
  Do[d = DivisorSigma[0, n]; 
   If[d > a, a = d; ans = Append[ans, n]], {n, 1, b}];
  Return[ans]]
InitialHC = {1, 2, 4, 6, 12, 24, 36, 48, 60, 120, 180, 240, 360, 720, 
   840, 1260, 1680, 2520, 5040, 7560, 10080, 15120, 20160, 25200, 
   27720, 45360, 50400, 55440, 83160, 110880, 166320, 221760, 277200, 
   332640, 498960, 554400, 665280, 720720, 1081080, 1441440, 2162160};


(*-------------------------------------------------------------------*)
InitialSA = {1, 2, 4, 6, 12, 24, 36, 48, 60, 120, 180, 240, 360, 720, 
   840, 1260, 1680, 2520, 5040, 10080, 15120, 25200, 27720, 55440, 
   110880, 166320, 277200, 332640, 554400, 665280, 720720, 1441440, 
   2162160,
   3603600, 4324320, 7207200, 8648640, 10810800, 21621600, 36756720,
   61261200, 73513440, 122522400, 147026880, 183783600, 367567200,
   698377680, 735134400, 1102701600, 1163962800, 1396755360, 2327925600,
   2793510720, 3491888400, 6983776800, 13967553600, 20951330400,
   27935107200, 41902660800, 48886437600, 80313433200, 160626866400,
   321253732800, 481880599200, 642507465600, 963761198400, 
   1124388064800,
   1927522396800, 2248776129600, 3373164194400, 4497552259200, 
   4658179125600,
   6746328388800, 9316358251200, 13974537376800, 18632716502400,
   27949074753600, 32607253879200, 55898149507200, 65214507758400,
   97821761637600, 130429015516800, 144403552893600, 195643523275200,
   288807105787200, 433210658680800, 577614211574400, 866421317361600,
   1010824870255200, 1732842634723200, 2021649740510400, 
   3032474610765600,
   4043299481020800, 6064949221531200, 10685862914126400, 
   12129898443062400,
   21371725828252800, 24259796886124800, 30324746107656000,
   32057588742379200, 37400520199442400, 64115177484758400,
   74801040398884800, 112201560598327200, 149602080797769600,
    224403121196654400};

Clear[GenerateSA];
GenerateSA::argx = "the function takes one argument";
GenerateSA::arg1 = "the argument is a natural number";
GenerateSA[] := (Message[GenerateSA::argx]; Abort[]);
GenerateSA[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[GenerateSA::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[GenerateSA::arg1]; Abort[])];
  Return[GenerateSAWorker[n]]]
GenerateSAWorker[b_] := 
 Module[{supers = {1}, ratios, n, next = 1, previous = 1},
  Do[ratios[n] = DivisorSigma[1, n]/n, {n, 1, b}];
  Do[next = ratios[n];
   If[next > previous, supers = Append[supers, n]; 
    previous = next], {n, 2, b}];
  Return[supers]]

(*------------------------------------------------------------------*)
InitialSA100={{2}, {2, 2}, {3}, {3, 2}, {3, 2, 2}, {3, 3}, {3, 2, 2, 2}, {5, 
  2}, {5, 2, 2}, {5, 3}, {5, 2, 2, 2}, {5, 3, 2}, {5, 3, 2, 2}, {7, 2,
   2}, {7, 3}, {7, 2, 2, 2}, {7, 3, 2}, {7, 3, 2, 2}, {7, 3, 2, 2, 
  2}, {7, 3, 3, 2}, {7, 5, 2, 2}, {11, 3, 2}, {11, 3, 2, 2}, {11, 3, 
  2, 2, 2}, {11, 3, 3, 2}, {11, 5, 2, 2}, {11, 3, 3, 2, 2}, {11, 5, 2,
   2, 2}, {11, 3, 3, 2, 2, 2}, {13, 3, 2, 2}, {13, 3, 2, 2, 2}, {13, 
  3, 3, 2}, {13, 5, 2, 2}, {13, 3, 3, 2, 2}, {13, 5, 2, 2, 2}, {13, 3,
   3, 2, 2, 2}, {13, 5, 3, 2}, {13, 5, 3, 2, 2}, {17, 3, 3, 2}, {17, 
  5, 2, 2}, {17, 3, 3, 2, 2}, {17, 5, 2, 2, 2}, {17, 3, 3, 2, 2, 
  2}, {17, 5, 3, 2}, {17, 5, 3, 2, 2}, {19, 3, 3, 2}, {17, 5, 3, 2, 2,
   2}, {17, 5, 3, 3, 2}, {19, 5, 2, 2}, {19, 3, 3, 2, 2}, {19, 5, 2, 
  2, 2}, {19, 3, 3, 2, 2, 2}, {19, 5, 3, 2}, {19, 5, 3, 2, 2}, {19, 5,
   3, 2, 2, 2}, {19, 5, 3, 3, 2}, {19, 5, 3, 2, 2, 2, 2}, {19, 5, 3, 
  3, 2, 2}, {19, 7, 3, 2, 2}, {23, 5, 3, 2}, {23, 5, 3, 2, 2}, {23, 5,
   3, 2, 2, 2}, {23, 5, 3, 3, 2}, {23, 5, 3, 2, 2, 2, 2}, {23, 5, 3, 
  3, 2, 2}, {23, 7, 3, 2, 2}, {23, 5, 3, 3, 2, 2, 2}, {23, 7, 3, 2, 2,
   2}, {23, 7, 3, 3, 2}, {23, 7, 3, 2, 2, 2, 2}, {29, 5, 3, 2, 
  2}, {23, 7, 3, 3, 2, 2}, {29, 5, 3, 2, 2, 2}, {29, 5, 3, 3, 2}, {29,
   5, 3, 2, 2, 2, 2}, {29, 5, 3, 3, 2, 2}, {29, 7, 3, 2, 2}, {29, 5, 
  3, 3, 2, 2, 2}, {29, 7, 3, 2, 2, 2}, {29, 7, 3, 3, 2}, {29, 7, 3, 2,
   2, 2, 2}, {31, 5, 3, 2, 2}, {29, 7, 3, 3, 2, 2}, {31, 5, 3, 2, 2, 
  2}, {31, 5, 3, 3, 2}, {31, 5, 3, 2, 2, 2, 2}, {31, 5, 3, 3, 2, 
  2}, {31, 7, 3, 2, 2}, {31, 5, 3, 3, 2, 2, 2}, {31, 7, 3, 2, 2, 
  2}, {31, 7, 3, 3, 2}, {31, 7, 3, 2, 2, 2, 2}, {31, 7, 3, 3, 2, 
  2}, {37, 5, 3, 2, 2, 2}, {31, 7, 3, 3, 2, 2, 2}, {37, 5, 3, 2, 2, 2,
   2}, {31, 7, 3, 3, 2, 2, 2, 2}, {31, 7, 5, 3, 2, 2}, {37, 5, 3, 3, 
  2, 2}, {37, 7, 3, 2, 2}};
  (*--------------------------------------------------------------*)
  
InitialCA = {2, 6, 12, 60, 120, 360, 2520, 5040, 55440, 720720, 
   1441440, 4324320, 21621600, 367567200, 6983776800, 160626866400, 
   321253732800, 9316358251200, 288807105787200, 2021649740510400, 
   6064949221531200, 224403121196654400, 9200527969062830400, 
   395622702669701707200, 791245405339403414400, 
   37188534050951960476800, 1970992304700453905270400, 
   116288545977326780410953600, 581442729886633902054768000, 
   35468006523084668025340848000, 2376356437046672757697836816000, 
   168721307030313765796546413936000, 
   12316655413212904903147888217328000, 
   135483209545341953934626770390608000, 
   10703173554082014360835514860858032000, 
   21406347108164028721671029721716064000, 
   1776726809977614383898695466902433312000, 
   5330180429932843151696086400707299936000, 
   474386058264023040500951689662949694304000, 
   46015447651610234928592313897306120347488000, 
   598200819470933054071700080664979564517344000, 
   60418282766564238461241708147162936016251744000, 
   6223083124956116561507895939157782409673929632000, 
   665869894370304472081344865489882717835110470624000, 
   72579818486363187456866590338397216244027041298016000, 
   8201519488959040182625924708238885435575055666675808000, 
   1041592975097798103193492437946338450318032069667827616000, 
   136448679737811551518347509370970336991662201126485417696000, 
   18693469124080182558013608783822936167857721554328502224352000, 
   2598392208247145375563891620951388127332223296051661809184928000};

(*---------------------------------------------------------------------------\
*)

InitialCAPrimes = {2, 3, 2, 5, 2, 3, 7, 2, 11, 13, 2, 3, 5, 17, 19, 
   23, 2, 29, 31, 7, 3, 37, 41, 43, 2, 47, 53, 59, 5, 61, 67, 71, 73, 
   11, 79, 2, 83, 3, 89, 97, 13, 101, 103, 107, 109, 113, 127, 131, 
   137, 139, 2, 149, 151, 7, 157, 163, 167, 17, 173, 179, 181, 191, 
   193, 197, 199, 19, 211, 3, 223, 227, 229, 5, 233, 239, 241, 251, 2,
    257, 263, 269, 271, 277, 281, 283, 293, 23, 307, 311, 313, 317, 
   331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 
   409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 2, 463, 467, 29, 
   479, 487, 491, 499, 503, 509, 521, 523, 541, 31, 547, 11, 557, 563,
    3, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 
   641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 
   727, 733, 739, 743, 751, 757, 761, 37, 769, 773, 787, 797, 809, 
   811, 7, 821, 823, 827, 829, 839, 2, 853, 857, 859, 863, 877, 881, 
   883, 887, 13, 907, 911, 919, 5, 929, 41, 937, 941, 947, 953, 967, 
   971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 43, 1031, 1033, 
   1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 
   1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 
   1201, 1213, 1217, 47, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 
   1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 
   1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 
   1447, 1451, 1453, 1459, 1471, 3, 1481, 1483, 1487, 1489, 1493, 
   1499, 1511, 1523, 1531, 1543, 2, 53, 1549, 1553, 1559, 1567, 1571, 
   1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 
   1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 
   1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 
   1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 
   59, 1913, 1931, 1933, 1949, 1951, 17, 1973, 1979, 1987, 1993, 1997,
    1999, 2003, 2011, 2017, 2027, 2029, 2039, 61, 2053, 2063, 2069, 
   2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 
   2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 
   2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 
   2339, 2341
    2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399,
   2411, 2417, 2423, 2437, 2441, 2447, 67, 2459, 2467, 2473, 2477, 
   2503,
   2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609,
   2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687,
   2689, 2693, 19, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 
   2749} ;

Clear[CAPrimes];
CAPrimes::argx = "the function takes one argument";
CAPrimes::arg1 = "the argument is a natural number";
CAPrimes[] := (Message[CAPrimes::argx]; Abort[]);
CAPrimes[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[CAPrimes::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[CAPrimes::arg1]; Abort[])];
  Return[CAPrimesWorker[n]]]
CAPrimesWorker[b_] := 
 Module[{pFactor, f, x, lst, maxN, primes, n, i}, 
  pFactor[f_List] := 
   Module[{p = f[[1]], k = f[[2]]}, 
    N[Log[(p^(k + 2) - 1)/(p^(k + 1) - 1)]/Log[p]] - 1]; maxN = b; 
  f = {{2, 1}, {3, 0}}; primes = 1; lst = {2}; 
  x = Table[pFactor[f[[i]]], {i, primes + 1}]; 
  For[n = 2, n <= maxN, n++, i = Position[x, Max[x]][[1, 1]]; 
   AppendTo[lst, f[[i, 1]]]; f[[i, 2]]++; 
   If[i > primes, primes++; AppendTo[f, {Prime[i + 1], 0}]; 
    AppendTo[x, pFactor[f[[-1]]]]]; x[[i]] = pFactor[f[[i]]]];
  Return[lst ]]

(*------------------------------------------------------------------*)

Clear[GenerateCA];
GenerateCA::argx = "the function takes one argument";
GenerateCA::arg1 = "the argument is a natural number";
GenerateCA[] := (Message[GenerateCA::argx]; Abort[]);
GenerateCA[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[GenerateCA::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[GenerateCA::arg1]; Abort[])];
  Return[GenerateCAWorker[n]]]
GenerateCAWorker[b_] := Module[{prs, n, ans = {}, next = 1, p},
  prs = CAPrimes[b];
  Do[p = prs[[n]];
   next = next*p;
   ans = Append[ans, next], {n, 1, b}];
  Return[ans]]

(*--------------------------------------------------------------------------*)

InitialCA150={{2}, {3}, {3, 2}, {5, 2}, {5, 2, 2}, {5, 3, 2}, {7, 3, 2}, {7, 3, 2, 
  2}, {11, 3, 2, 2}, {13, 3, 2, 2}, {13, 3, 2, 2, 2}, {13, 3, 3, 2, 
  2}, {13, 5, 3, 2, 2}, {17, 5, 3, 2, 2}, {19, 5, 3, 2, 2}, {23, 5, 3,
   2, 2}, {23, 5, 3, 2, 2, 2}, {29, 5, 3, 2, 2, 2}, {31, 5, 3, 2, 2, 
  2}, {31, 7, 3, 2, 2, 2}, {31, 7, 3, 3, 2, 2}, {37, 7, 3, 3, 2, 
  2}, {41, 7, 3, 3, 2, 2}, {43, 7, 3, 3, 2, 2}, {43, 7, 3, 3, 2, 2, 
  2}, {47, 7, 3, 3, 2, 2, 2}, {53, 7, 3, 3, 2, 2, 2}, {59, 7, 3, 3, 2,
   2, 2}, {59, 7, 5, 3, 2, 2, 2}, {61, 7, 5, 3, 2, 2, 2}, {67, 7, 5, 
  3, 2, 2, 2}, {71, 7, 5, 3, 2, 2, 2}, {73, 7, 5, 3, 2, 2, 2}, {73, 
  11, 5, 3, 2, 2, 2}, {79, 11, 5, 3, 2, 2, 2}, {79, 11, 5, 3, 2, 2, 2,
   2}, {83, 11, 5, 3, 2, 2, 2, 2}, {83, 11, 5, 3, 3, 2, 2, 2}, {89, 
  11, 5, 3, 3, 2, 2, 2}, {97, 11, 5, 3, 3, 2, 2, 2}, {97, 13, 5, 3, 3,
   2, 2, 2}, {101, 13, 5, 3, 3, 2, 2, 2}, {103, 13, 5, 3, 3, 2, 2, 
  2}, {107, 13, 5, 3, 3, 2, 2, 2}, {109, 13, 5, 3, 3, 2, 2, 2}, {113, 
  13, 5, 3, 3, 2, 2, 2}, {127, 13, 5, 3, 3, 2, 2, 2}, {131, 13, 5, 3, 
  3, 2, 2, 2}, {137, 13, 5, 3, 3, 2, 2, 2}, {139, 13, 5, 3, 3, 2, 2, 
  2}, {139, 13, 5, 3, 3, 2, 2, 2, 2}, {149, 13, 5, 3, 3, 2, 2, 2, 
  2}, {151, 13, 5, 3, 3, 2, 2, 2, 2}, {151, 13, 7, 3, 3, 2, 2, 2, 
  2}, {157, 13, 7, 3, 3, 2, 2, 2, 2}, {163, 13, 7, 3, 3, 2, 2, 2, 
  2}, {167, 13, 7, 3, 3, 2, 2, 2, 2}, {167, 17, 7, 3, 3, 2, 2, 2, 
  2}, {173, 17, 7, 3, 3, 2, 2, 2, 2}, {179, 17, 7, 3, 3, 2, 2, 2, 
  2}, {181, 17, 7, 3, 3, 2, 2, 2, 2}, {191, 17, 7, 3, 3, 2, 2, 2, 
  2}, {193, 17, 7, 3, 3, 2, 2, 2, 2}, {197, 17, 7, 3, 3, 2, 2, 2, 
  2}, {199, 17, 7, 3, 3, 2, 2, 2, 2}, {199, 19, 7, 3, 3, 2, 2, 2, 
  2}, {211, 19, 7, 3, 3, 2, 2, 2, 2}, {211, 19, 7, 3, 3, 3, 2, 2, 
  2}, {223, 19, 7, 3, 3, 3, 2, 2, 2}, {227, 19, 7, 3, 3, 3, 2, 2, 
  2}, {229, 19, 7, 3, 3, 3, 2, 2, 2}, {229, 19, 7, 5, 3, 3, 2, 2, 
  2}, {233, 19, 7, 5, 3, 3, 2, 2, 2}, {239, 19, 7, 5, 3, 3, 2, 2, 
  2}, {241, 19, 7, 5, 3, 3, 2, 2, 2}, {251, 19, 7, 5, 3, 3, 2, 2, 
  2}, {251, 19, 7, 5, 3, 3, 2, 2, 2, 2}, {257, 19, 7, 5, 3, 3, 2, 2, 
  2, 2}, {263, 19, 7, 5, 3, 3, 2, 2, 2, 2}, {269, 19, 7, 5, 3, 3, 2, 
  2, 2, 2}, {271, 19, 7, 5, 3, 3, 2, 2, 2, 2}, {277, 19, 7, 5, 3, 3, 
  2, 2, 2, 2}, {281, 19, 7, 5, 3, 3, 2, 2, 2, 2}, {283, 19, 7, 5, 3, 
  3, 2, 2, 2, 2}, {293, 19, 7, 5, 3, 3, 2, 2, 2, 2}, {293, 23, 7, 5, 
  3, 3, 2, 2, 2, 2}, {307, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {311, 23, 7, 
  5, 3, 3, 2, 2, 2, 2}, {313, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {317, 23, 
  7, 5, 3, 3, 2, 2, 2, 2}, {331, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {337, 
  23, 7, 5, 3, 3, 2, 2, 2, 2}, {347, 23, 7, 5, 3, 3, 2, 2, 2, 
  2}, {349, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {353, 23, 7, 5, 3, 3, 2, 2, 
  2, 2}, {359, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {367, 23, 7, 5, 3, 3, 2, 
  2, 2, 2}, {373, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {379, 23, 7, 5, 3, 3, 
  2, 2, 2, 2}, {383, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {389, 23, 7, 5, 3, 
  3, 2, 2, 2, 2}, {397, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {401, 23, 7, 5, 
  3, 3, 2, 2, 2, 2}, {409, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {419, 23, 7, 
  5, 3, 3, 2, 2, 2, 2}, {421, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {431, 23, 
  7, 5, 3, 3, 2, 2, 2, 2}, {433, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {439, 
  23, 7, 5, 3, 3, 2, 2, 2, 2}, {443, 23, 7, 5, 3, 3, 2, 2, 2, 
  2}, {449, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {457, 23, 7, 5, 3, 3, 2, 2, 
  2, 2}, {461, 23, 7, 5, 3, 3, 2, 2, 2, 2}, {461, 23, 7, 5, 3, 3, 2, 
  2, 2, 2, 2}, {463, 23, 7, 5, 3, 3, 2, 2, 2, 2, 2}, {467, 23, 7, 5, 
  3, 3, 2, 2, 2, 2, 2}, {467, 29, 7, 5, 3, 3, 2, 2, 2, 2, 2}, {479, 
  29, 7, 5, 3, 3, 2, 2, 2, 2, 2}, {487, 29, 7, 5, 3, 3, 2, 2, 2, 2, 
  2}, {491, 29, 7, 5, 3, 3, 2, 2, 2, 2, 2}, {499, 29, 7, 5, 3, 3, 2, 
  2, 2, 2, 2}, {503, 29, 7, 5, 3, 3, 2, 2, 2, 2, 2}, {509, 29, 7, 5, 
  3, 3, 2, 2, 2, 2, 2}, {521, 29, 7, 5, 3, 3, 2, 2, 2, 2, 2}, {523, 
  29, 7, 5, 3, 3, 2, 2, 2, 2, 2}, {541, 29, 7, 5, 3, 3, 2, 2, 2, 2, 
  2}, {541, 31, 7, 5, 3, 3, 2, 2, 2, 2, 2}, {547, 31, 7, 5, 3, 3, 2, 
  2, 2, 2, 2}, {547, 31, 11, 5, 3, 3, 2, 2, 2, 2, 2}, {557, 31, 11, 5,
   3, 3, 2, 2, 2, 2, 2}, {563, 31, 11, 5, 3, 3, 2, 2, 2, 2, 2}, {563, 
  31, 11, 5, 3, 3, 3, 2, 2, 2, 2}, {569, 31, 11, 5, 3, 3, 3, 2, 2, 2, 
  2}, {571, 31, 11, 5, 3, 3, 3, 2, 2, 2, 2}, {577, 31, 11, 5, 3, 3, 3,
   2, 2, 2, 2}, {587, 31, 11, 5, 3, 3, 3, 2, 2, 2, 2}, {593, 31, 11, 
  5, 3, 3, 3, 2, 2, 2, 2}, {599, 31, 11, 5, 3, 3, 3, 2, 2, 2, 
  2}, {601, 31, 11, 5, 3, 3, 3, 2, 2, 2, 2}, {607, 31, 11, 5, 3, 3, 3,
   2, 2, 2, 2}, {613, 31, 11, 5, 3, 3, 3, 2, 2, 2, 2}, {617, 31, 11, 
  5, 3, 3, 3, 2, 2, 2, 2}, {619, 31, 11, 5, 3, 3, 3, 2, 2, 2, 
  2}, {631, 31, 11, 5, 3, 3, 3, 2, 2, 2, 2}, {641, 31, 11, 5, 3, 3, 3,
   2, 2, 2, 2}, {643, 31, 11, 5, 3, 3, 3, 2, 2, 2, 2}, {647, 31, 11, 
  5, 3, 3, 3, 2, 2, 2, 2}, {653, 31, 11, 5, 3, 3, 3, 2, 2, 2, 
  2}, {659, 31, 11, 5, 3, 3, 3, 2, 2, 2, 2}, {661, 31, 11, 5, 3, 3, 3,
   2, 2, 2, 2}};
   
   (*--------------------------------------------------------------------------*)
InitialXA = {10080, 
   8201519488959040182625924708238885435575055666675808000, 
   1041592975097798103193492437946338450318032069667827616000, 
   136448679737811551518347509370970336991662201126485417696000, 
   18693469124080182558013608783822936167857721554328502224352000, 
   2598392208247145375563891620951388127332223296051661809184928000, 
   5196784416494290751127783241902776254664446592103323618369856000};
Clear[FindNextXA];
 
FindNextXA::argx = "the function takes two argument";
FindNextXA::arg1 = 
  "the first argument is a natural number which must be extremely \
abundant";
FindNextXA::arg1 = "the second argument is a natural number"; 
FindNextXA[] := (Message[FindNextXA::argx]; Abort[]);
FindNextXA[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 2, (Message[FindNextXA::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[FindNextXA::arg1]; Abort[])];
  If[Not[NaturalNumberQ[n]], (Message[FindNextXA::arg1]; Abort[])]; 
  Return[FindNextXAWorker[n, d]]]
FindNextXAWorker[n_, d_] := Module[{next, gn},
  next = n + 1;
  gn = GronwallGWorker[n, d];
  While[GronwallGWorker[next, d] < gn, next++];
  Return[next]]

(*--------Nicolas reversed \
---------------------------------------------------------*)

InitialNR = {1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 22, 
   24, 26,
   28, 30, 36, 40, 42, 48, 50, 54, 60, 66, 70, 72, 78, 84, 90, 96, 
   102, 108,
   114, 120, 126, 132, 138, 140, 144, 150, 156, 162, 168, 174, 180, 
   186,
   192, 198, 204, 210, 216, 222, 228, 234, 240, 246, 252, 258, 264, 
   270,
   276, 294, 300, 306, 312, 330, 336, 342, 360, 378, 390, 396, 420, 
   450,
   462, 468, 480, 504, 510, 528, 540, 546, 570, 588, 600, 630, 660, 
   672,
   690, 714, 720, 750, 756, 780, 798, 810, 840, 858, 870, 882, 900, 
   924,
   930, 960, 966, 990, 1008, 1020, 1050, 1080, 1092, 1110, 1122, 1134,
   1140, 1170, 1176, 1200, 1218, 1230, 1260, 1290, 1302, 1320, 1350,
   1380, 1386, 1410, 1428, 1440, 1470, 1500, 1530, 1554, 1560, 1590,
   1596, 1620, 1638, 1650, 1680, 1710, 1722, 1740, 1770, 1800, 1830,
   1848, 1860, 1890, 1920, 1932, 1950, 1980, 2010, 2040, 2070, 2100,
   2130, 2142, 2160, 2184, 2190, 2220, 2250, 2280, 2310, 2340, 2370,
   2394, 2400, 2430, 2460, 2490, 2520, 2550, 2580, 2610, 2640, 2670,
   2700, 2730, 2760, 2772, 2790, 2820, 2850, 2856, 2880, 2910, 2940,
   2970, 3000, 3030, 3060, 3090, 3120, 3150, 3180, 3210, 3234, 3240,
   3270, 3276, 3300, 3330, 3360, 3390, 3420, 3450, 3480, 3510, 3540,
   3570, 3600, 3630, 3660, 3690, 3696, 3720, 3780, 3810, 3822, 3870,
   3900, 3930, 3960, 3990, 4020, 4080, 4110, 4140, 4158, 4170, 4200,
   4230, 4260, 4290, 4350, 4368, 4380, 4410, 4440, 4560, 4590, 4620,
   4650, 4680, 4770, 4830, 4920, 4950, 5040, 5070, 5082, 5100, 5130,
   5160, 5220, 5250, 5280, 5460, 5520, 5544, 5550, 5580, 5610, 5670,
   5700, 5850, 5880, 5940, 6006, 6090, 6120, 6210, 6240, 6270, 6300,
   6510, 6600, 6630, 6720, 6840, 6900, 6930, 6960, 7020, 7140, 7260,
   7350, 7410, 7560, 7590, 7650, 7770, 7800, 7854, 7920, 7980, 8160,
   8190, 8250, 8280, 8400, 8550, 8580, 8610, 8670, 8778, 8820, 8910,
   8970, 9030, 9120, 9180, 9240, 9282, 9360, 9450, 9570, 9660, 9690,
   9750, 9870, 9900};
(*=====Utilities=======================================================\
=*)

Clear[NthPrime];
NthPrime::argx = "the function takes one argument";
NthPrime::arg1 = "the argument is a natural number";
NthPrime[] := (Message[NthPrime::argx]; Abort[]);
NthPrime[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[NthPrime::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[NthPrime::arg1]; Abort[])];
  Return[NthPrimeWorker[n]]]
NthPrimeWorker[n_] := Prime[n]

(*------------------------------------------------------------------*)

Clear[MaxPrime];
MaxPrime::argx = "the function takes one argument";
MaxPrime::arg1 = "the first argument is a positive real number";
MaxPrime[] := (Message[MaxPrime::argx]; Abort[]);
MaxPrime[argseq___] := 
 Module[{args = List[argseq], x, nz, d, ans, j, gj},
  If[Length[args] != 1, (Message[MaxPrime::argx]; Abort[])];
  x = args[[1]];
  If[Not[PositiveRealQ[x ]], (Message[MaxPrime::arg1]; Abort[])];
  Return[MaxPrimeWorker[x]]]
MaxPrimeWorker[x_] := Prime[ PrimePi[IntegerPart[x]]]

(*------------------------------------------------------------------*)


(* Clear[NthPrimorial];
NthPrimorial::argx = "the function takes one argument";
NthPrimorial::arg1 = "the argument is a natural number";
NthPrimorial[] := (Message[NthPrimorial::argx]; Abort[]);
NthPrimorial[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[NthPrimorial::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[NthPrimorial::arg1]; Abort[])];
  Return[NthPrimorialWorker[n]]]
NthPrimorialWorker[n_] := Product[Prime[j], {j, 1, n}] *)

(*------------------------------------------------------------------*)
Clear[NthHarmonicNumber];
NthHarmonicNumber::argx = "the function takes two arguments";
NthHarmonicNumber::arg1 = "the first argument is a natural number";
NthHarmonicNumber::arg2 = "the second argument is a natural number";
NthHarmonicNumber[] := (Message[NthHarmonicNumber::argx]; Abort[]);
NthHarmonicNumber[argseq___] := Module[{args = List[argseq], n, d},
  If[Length[args] != 2, (Message[NthHarmonicNumber::argx]; Abort[])];
  n = args[[1]]; d = args[[2]];
  If[Not[NaturalNumberQ[n ]], (Message[NthHarmonicNumber::arg1]; 
    Abort[])];
  If[Not[NaturalNumberQ[d]], (Message[NthHarmonicNumber::arg2]; 
    Abort[])];
  Return[NthHarmonicNumberWorker[n, d]]]
NthHarmonicNumberWorker[n_, d_] := Module[{sum = N[0, d], j},
  Do[sum = sum + N[1/j, d], {j, 1, n}];
  Return[sum]]


(*------------------------------------------------------------------*)


Clear[NiceFactorInteger];
NiceFactorInteger::argx = "the function takes one argument";
NiceFactorInteger::arg1 = "the argument is a natural number";
NiceFactorInteger[] := (Message[NiceFactorInteger::argx]; Abort[]);
NiceFactorInteger[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[NiceFactorInteger::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[NiceFactorInteger::arg1]; 
    Abort[])];
  Return[NiceFactorIntegerWorker[n]]]
NiceFactorIntegerWorker[n_] := 
 If[n == 1, {1}, Map[FixFactor, FactorInteger[n]]]

(*------------------------------------------------------------------*)


Clear[AllPrimesFactor];
AllPrimesFactor::argx = "the function takes one argument";
AllPrimesFactor::arg1 = "the argument is a natural number";
AllPrimesFactor[] := (Message[AllPrimesFactor::argx]; Abort[]);
AllPrimesFactor[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[AllPrimesFactor::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[AllPrimesFactor::arg1]; 
    Abort[])];
  Return[AllPrimesFactorWorker[n]]]
AllPrimesFactorWorker[n_] := Module[{fs, P, J, ans = {}, primes, pj},
  If[n == 1, Return[{1}]];
  fs = FactorInteger[n];
  P = First[Last[fs]];
  J = PrimePi[P];
  primes = Map[First, fs];
  Do[pj = Prime[j]; 
   ans = Append[ans, 
     If[MemberQ[primes, pj], {pj, PrimePower[pj, n]}, {pj, 0}]], {j, 
    1, J}];
  Return[ans]]

(*------------------------------------------------------------------*)
Clear[PrimePower];
PrimePower::argx = "the function takes two arguments";
PrimePower::arg1 = "the first argument is a natural number greater than 1";
PrimePower::arg2 = "the second argument is a natural number";
PrimePower::args = "the first argument must divide the second";
PrimePower[] := (Message[PrimePower::argx]; Abort[]);
PrimePower[argseq___] := Module[{args = List[argseq], p, n},
  If[Length[args] != 2, (Message[PrimePower::argx]; Abort[])];
  p = args[[1]]; n = args[[2]]; 
  If[Not[NaturalNumberQ[p] && p>1], (Message[PrimePower::arg1]; Abort[])];
  If[Not[NaturalNumberQ[n ]], (Message[PrimePower::arg2]; Abort[])];
  If[Not[Mod[n,p]==0], (Message[PrimePower::args]; Abort[])];
  Return[PrimePowerWorker[p, n]]]
PrimePowerWorker[p_, n_] := Module[{},
  If[Mod[n,p]!= 0,Return[0]];
  If[p==n, 1, Return[1+PrimePowerWorker[p,n/p]]]]


(*------------------------------------------------------------------*)


(* Clear[CompositeQ];
CompositeQ::argx = "the function takes one argument";
CompositeQ::arg1 = "the argument is a natural number";
CompositeQ[] := (Message[CompositeQ::argx]; Abort[]);
CompositeQ[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[CompositeQ::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[CompositeQ::arg1]; Abort[])];
  Return[CompositeQWorker[n]]] *)
CompositeQWorker[n_] := IntegerQ[n] && n > 3 && Not[PrimeQ[n]]

(*------------------------------------------------------------------*)
Clear[PowerfulQ];
PowerfulQ::argx = "the function takes one argument";
PowerfulQ::arg1 = "the argument is a natural number";
PowerfulQ[] := (Message[PowerfulQ::argx]; Abort[]);
PowerfulQ[argseq___] := Module[{args = List[argseq], n},
  If[Length[args] != 1, (Message[PowerfulQ::argx]; Abort[])];
  n = args[[1]]; 
  If[Not[NaturalNumberQ[n]], (Message[PowerfulQ::arg1]; Abort[])];
  Return[PowerfulQWorker[n]]]
PowerfulQWorker[n_] := Module[{pows, j, ans = True},
  If[n == 1 || SquareFreeQ[n], Return[False]];
  pows = Map[#[[2]] &, FactorInteger[n]];
  Do[If[pows[[j]] < 2, ans = False; Break[]], {j, 1, Length[pows]}];
  Return[ans]]
 (* -------------------------------------------------------------*)



(* internal utilities *)

NaturalNumberQ[n_] := IntegerQ[n] && n > 0

Omega[n_] := If[n == 1, 0, Apply[Plus, Map[First, FactorInteger[n]]]]

FixFactor[{p_, e_}] := If[e == 1, p, {p, e}]

RealNumberQ[x_] := Element[x, Reals]

PositiveRealQ[x_] := NumberQ[x] && RealNumberQ[x] && x > 0 

NaturalNumberListQ[x_]:= ListQ[x] && Union[Map[NaturalNumberQ, x]]=={True}
(* ----------------end of RHpack source------------------------*)

End[];
Print[StyleForm["RHpack loaded.",FontWeight->"Bold"]];
EndPackage[];