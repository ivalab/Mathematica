(*  configuration.m

  This package provides the necessary function andmanipulations required
  to implement the dynamics of a configuration manifold with a principle
  connection.  The configuration manifold in this case is equivalent to
  Q=MxG where the base space M is a manifold, and the position space G is 
  a Lie group.

  Patricio Vela
*)

BeginPackage["BKMM`"];

Needs["Objects`","mathematica/libs/objects.m"];
Needs["Euclidean`","mathematica/libs/euclidean.m"];
Needs["Vectors`","mathematica/libs/vectors.m"];
Needs["Manifolds`","mathematica/libs/manifolds.m"];
Needs["LieGroups`","mathematica/libs/groups.m"];
Needs["LagrangianMechanics`","mathematica/libs/lagrangian.m"];
Needs["Equations`","mathematica/libs/equations.m"];
Needs["PrincipalBundle`","mathematica/libs/principal.m"];

pInit::usage="--";
fDimBase::usage="--";
fDimGroup::usage="--";
fDimMom::usage="--";
pDefBase::usage="--";
pDefGroup::usage="--";
pDefTotal::usage="--";
pDefDimMom::usage="--";
pDefMomentum::usage="--";

BKMM::usage="Provides functionality of principal bundle Q = R^n x SE(2)";
oInit::usage="oInit[Q, name_, params ... ] \n Initializes the principal
bundle object given by Q and sets it's name plus additional optional parameters.";
dBase::usage="--";
dGroup::usage="--";
dTotal::usage="--";
dMomentum::usage="--";

gBase::usage="--";
gGroup::usage="--";
gTotal::usage="--";
gMomentum::usage="--";

dDimMom::usage="--";
gDimBase::usage="--";
gDimGroup::usage="--";
gDimMom::usage="--";

dBasis::usage="--";
dT::usage="--";
dGamma::usage="--";
dSC::usage="--";
dLagrangian::usage="--";
dLI::usage="--";

dConnection::usage="--";
dAffine::usage="--";
dConnectionNoSimp::usage="--";

dl1::usage="--";
dl2::usage="--";

gBase::usage="--";
gGroup::usage="--";
gMomentum::usage="--";
gBasis::usage="--";
gT::usage="--";
gGamma::usage="--";
gSC::usage="--";
gConnection::usage="--";
gAffine::usage="--";
gCurvature::usage="--";
gLagrangian::usage="--";
gLI::usage="--";
gLIinv::usage="--";
gl1::usage="--";
gl2::usage="--";
gMomCoefficients::usage="--";
gBaseCoefficients::usage="--";
gBaseVariation::usage="--";
gDynamics::usage="--";
cReduceLagrangian::usage="--";
cTransformLagrangian::usage="--";
cCoefficients::usage="--";
cMomCoefficients::usage="--";
cBaseCoefficients::usage="--";
cVariation::usage="--";
cDynamics::usage="--";
gMomDynamics::usage="--";
gLDP::usage="--";
gEOM::usage="--";

fEOM::usage="--";
fDSolve::usage="--";

Diff::usage="Take the differential of a function with respect to the configuration space.";
fLie::usage="fLie[X, Y] \n Take the Lie bracket of the two reduced
vector fields.";

CalcGamma::usage="--";

Begin["`Private`"];

(* R_SE2 is Q=R(m)xSE(2).  *)

iNAME		= 1;    (* Name of the object for Output*)
iDIMBASE 	= 2;	(* dimension of the base space. *)
iDIMMOM 	= 3;	(* dimension of the momentum space. *)
iBASE 		= 4; 	(* base variables. *)
iGROUP		= 5;	(* group variables. *)
iTOTAL		= 6;	(* total configuration. *)
iMOMENTUM	= 7; 	(*the momentum variable(s). *)
iBASIS		= 8;	(* basis vectors for the Lie algebra. *)
iT		= 9;	(* transformation matrix for the Lie algebra. *)
iGAMMA		= 10;	(* partial of the basis with respect to base. *)
iSC		= 11;	(* structure constants. *)
iCONN		= 12; 	(* Principle connection. *)
iAFFINE		= 13; 	(* Affine part of Principle connection. *)
iCURV		= 14;	(* Curvature of the Principle Connection. *)
iLAGRANGIAN	= 15;	(* (Reduced) Lagrangian of the system. *)
iLI		= 16;   (* Locked Inertia Tensor. *)
iLIinv		= 17;	(* Inverse of the Locked Inertia Tensor. *)
il1		= 18;	(* Partials of the Lagrangian. *)
il2		= 19;	(* Partials of the Lagrangian. *)
iMOMDYN		= 20;	(* Dynamic equations of the momentum. *)
iBASEDYN	= 21;	(* Dynamic equations of the base space. *)
iBASEVAR	= 22;	(* Variation of Lagrangian w.r.t. base vars. *)
iDYNAMICS	= 23;	(* The dynamical equations of motion. *)
iSYSTEM		= 24;	(* This is the system type *)
iMAX 		= 24;	(* max index. *)

(*--------------------oInit--------------------*)
(*--  Initialize the configuration object.
*)

oInit[BKMM, name_, m_] ^:= Module[
  {S},
  S = Table[Null, {i, iMAX}];
  S = ReplacePart[S, name, iNAME];
  S = ReplacePart[S, m, iDIMBASE];
  BKMM[S]
];

oInit[BKMM, name_, m_, p_] ^:= Module[
  {S} , 
  S = Table[0 , {i, iMAX}];
  S = ReplacePart[S, name, iNAME];
  S = ReplacePart[S, m , iDIMBASE];
  S = ReplacePart[S, p, iDIMMOM];
  BKMM[S]
];

BKMM[pInit, name_, m_, p_] ^:= Module[
  {S} , 
  S = Table[0 , {i, iMAX}];
  S = ReplacePart[S, name, iNAME];
  S = ReplacePart[S, m , iDIMBASE];
  S = ReplacePart[S, p, iDIMMOM]
];

(*--------------------State Representations--------------------*)

(*--Base and group variables--*)
dBase[BKMM[obj_], r_] ^:= BKMM[ReplacePart[obj, r, iBASE]];
dGroup[BKMM[obj_], g_] ^:= BKMM[ReplacePart[obj, g, iGROUP]];
dTotal[BKMM[obj_] , q_] ^:= BKMM[ReplacePart[obj, q, iTOTAL]];

gDimBase[BKMM[obj_]] ^:= obj[[iDIMBASE]];
gDimGroup[BKMM[obj_]] ^:= gDim[SE2];
gDimMom[BKMM[obj_]] ^:= obj[[iDIMMOM]];

(*--Momentum variable for reduced equations--*)
dDimMom[BKMM[obj_], dimP_] ^:= BKMM[ReplacePart[obj, dimP, iDIMMOM]];
dMomentum[BKMM[obj_], p_] ^:= BKMM[ReplacePart[obj,
  Table[ p[i] , {i, 1, gDimMom[BKMM[obj]]} ] , iMOMENTUM] ];

(*--------------------Differential Geometry--------------------*)
dBasis[BKMM[S_], e_] ^:= BKMM[ReplacePart[S, e, iBASIS]];
dT[BKMM[S_], T_] := BKMM[ReplacePart[S, T, iT]];
dGamma[BKMM[S_], gamma_] ^:= BKMM[ReplacePart[S, gamma, iGAMMA]];
dSC[BKMM[S_], SC_] ^:= BKMM[ReplacePart[S, SC, iSC]];


(*--------------------Lagrangian Mechanics--------------------*)
dLagrangian[BKMM[obj_], L_] ^:= BKMM[ReplacePart[obj, L, iLAGRANGIAN]];

dLI[BKMM[S_], LI_] ^:= Module[
  {LIinv} , 
  
  LIinv = Simplify[Inverse[LI]];
  BKMM[ReplacePart[ReplacePart[S, LI, iLI], LIinv, iLIinv]]
];

(*--------------------Constraints--------------------*)
dConnection[BKMM[S_], A_] ^:= Module[
  {dimG, dimM, tS, r, fA, fSC, B1, B2, B3, gamma, curv, curvB, B} , 
  r = gLocal[gBase[BKMM[S]]];
  tS = ReplacePart[S, A, iCONN];
  fA[r_] := Evaluate[A];
  fSC[r_] := Evaluate[ S[[iSC]] ]; 
  gamma = gGamma[BKMM[S]];
  dimG = gDimGroup[BKMM[S]];
  dimM = gDimBase[BKMM[S]];

  B1 = Table[ Simplify[ 
	 D[fA[r][[b, alpha]], r[[beta]]] 
           - D[fA[r][[b, beta]], r[[alpha]]] ], 
       {b, 1, dimG}, {alpha, 1, dimM}, {beta, 1, dimM}];
  B2 =Table[ Simplify[ 
	 Sum[ fA[r][[d, alpha]] gamma[[b, d, beta]]
           - fA[r][[d, beta]] gamma[[b, d, alpha]] , {d, 1, dimG} ] ],
       {b, 1, dimG}, {alpha, 1, dimM}, {beta, 1, dimM}];
  B3 = Table[  Simplify[
         - Sum[fSC[r][[b, a, c]] fA[r][[a, alpha]] fA[r][[c, beta]], 
             {a, 1, dimG}, {c, 1, dimG}] ], 
       {b, 1, dimG}, {alpha, 1, dimM}, {beta, 1, dimM}];
  curvB = Simplify[B1+B2+B3];
  BKMM[ReplacePart[tS, curvB, iCURV]]
];

dConnectionNoSimp[BKMM[S_], A_] ^:= Module[
  {dimG, dimM, tS, r, fA, fSC, B1, B2, B3, gamma, curv, curvB, B} , 
  r = gLocal[gBase[BKMM[S]]];
  tS = ReplacePart[S, A, iCONN];
  fA[r_] := Evaluate[A];
  fSC[r_] := Evaluate[ S[[iSC]] ]; 
  gamma = gGamm[BKMM[S]];
  dimG = gDimGroup[BKMM[S]];
  dimM = gDimBase[BKMM[S]];

  B1 = Table[ D[fA[r][[b, alpha]], r[[beta]]] 
           - D[fA[r][[b, beta]], r[[alpha]]] , 
       {b, 1, dimG}, {alpha, 1, dimM}, {beta, 1, dimM}];
  B2 = Table[ 
	 Sum[ fA[r][[d, alpha]] gamma[[b, d, beta]]
           - fA[r][[d, beta]] gamma[[b, d, alpha]] , {d, 1, dimG} ] ,
       {b, 1, dimG}, {alpha, 1, dimM}, {beta, 1, dimM}];
  B3 = Table[  
         - Sum[fSC[r][[b, a, c]] fA[r][[a, alpha]] fA[r][[c, beta]], 
             {a, 1, dimG}, {c, 1, dimG}] , 
       {b, 1, dimG}, {alpha, 1, dimM}, {beta, 1, dimM}];
  curvB = B1+B2+B3;
  BKMM[ReplacePart[tS, curvB, iCURV]]
];


dAffine[BKMM[S_], Aff_] ^:= BKMM[ReplacePart[S, Aff, iAFFINE]];
dl1[BKMM[S_], l1_] ^:= BKMM[ReplacePart[S, l1, il1]];
dl2[BKMM[S_], l2_] ^:= BKMM[ReplacePart[S, l2, il2]];

gBase[BKMM[S_]] ^:= Part[S, iBASE];
gGroup[BKMM[S_]] ^:= Part[S, iGROUP];
gMomentum[BKMM[S_]] ^:= Part[S, iMOMENTUM];
gBasis[BKMM[S_]] ^:= Part[S, iBASIS];
gT[BKMM[S_]] ^:= Part[S, iT];
gGamma[BKMM[S_]] ^:= Part[S, iGAMMA];
gSC[BKMM[S_]] ^:= Part[S, iSC];
gConnection[BKMM[S_]] ^:= Part[S, iCONN];
gAffine[BKMM[S_]] ^:= Part[S, iAFFINE];
gCurvature[BKMM[S_]] ^:= Part[S, iCURV];
gLagrangian[BKMM[S_]] ^:= Part[S, iLAGRANGIAN];
gLI[BKMM[S_]] ^:= Part[S, iLI];
gLIinv[BKMM[S_]] ^:= Part[S, iLIinv];
gl1[BKMM[S_]] ^:= Part[S, il1];
gl2[BKMM[S_]] ^:= Part[S, il2];
gMomCoefficients[BKMM[S_]] ^:= Part[S, iMOMDYN];
gBaseCoefficients[BKMM[S_]] ^:= Part[S, iBASEDYN];
gBaseVariation[BKMM[S_]] ^:= Part[S, iBASEVAR];
gDynamics[BKMM[S_]] ^:= Part[S, iDYNAMICS];


cReduceLagrangian[BKMM[S_], xi_] ^:= Module[
  {tS, L, T, g, e} , 

  L = gLagrangian[BKMM[S]];
  g = gGroup[BKMM[S]];
  e = gGroupIdentity[SE2];
  tS = ReplacePart[S, fReduceLagrangian[ L, g, e, xi], iLAGRANGIAN];
  T = gT[BKMM[tS]];
  dT[BKMM[tS], fReduceTransformation[T, g, e]]
];

cTransformLagrangian[BKMM[S_], xi_] ^:= Module[
  {L, T, e} , 

  L = gLagrangian[BKMM[S]];
  T = gT[BKMM[S]];
  e = gLocal[gBasis[BKMM[S]]];
  L = Simplify[ fTransformLagrangian[ L , Inverse[T], xi, e] ];
  BKMM[ReplacePart[S, L, iLAGRANGIAN]]
];

cCoefficients[BKMM[S_]] ^:= Module[
  {L, l1, l2, LI, LIinv, e, r, tS, dimG, dimM} , 

  L = gLagrangian[BKMM[S]];
  e = gLocal[gBasis[BKMM[S]]];
  r = gLocal[gBase[BKMM[S]]];
  dimG = gDimGroup[BKMM[S]];
  dimM = gDimBase[BKMM[S]];

  l1 = Simplify[
    Table[ D[ D[L, e[[i]]],  r[[j]]' ] , {i, 1, dimG}, {j, 1, dimM}]];
  LI = Simplify[
    Table[ D[ D[L, e[[i]]],  e[[j]] ], {i, 1, dimG}, {j, 1, dimG}]];
  l2 = LI;
  tS = dLI[BKMM[S],LI];
  tS = dl1[tS, l1];
  dl2[tS, l2]
];

cMomCoefficients[BKMM[S_]] ^:= Module[
  { lambda, D1, D2, eq1, eq2 , pd , 
    l1, l2, SC, A, B, LIinv,  gamma, dimG, dimM , dimP} , 

  gamma = gGamma[BKMM[S]];
  SC = gSC[BKMM[S]];
  A = gConnection[BKMM[S]];
  B = gCurvature[BKMM[S]];
  LIinv = gLIinv[BKMM[S]];
  l1 = gl1[BKMM[S]];
  l2 = gl2[BKMM[S]];
  dimG = gDimGroup[BKMM[S]];
  dimM = gDimBase[BKMM[S]];
  dimP = gDimMom[BKMM[S]];

  lambda = Table[ l1[[aa , alpha]] 
    - Sum[ l2[[aa, bb]] A[[bb , alpha]] , {bb, 1, dimG}],
  {aa, 1, dimG}, {alpha, 1, dimM}];

  D1 = Table[ Sum[ -SC[[c, aa, b]] A[[aa, alpha]] , {aa, 1, dimG}] 
    + gamma[[c, b, alpha]] 
    + Sum[ lambda[[aa, alpha]] 
      * Sum[ SC[[aa, a, b]] LIinv[[a, c]] , {a, 1, 1}] , {aa, 1,dimG}] ,
  {c, 1, dimP} , {b, 1, dimP}, {alpha, 1, dimM}];

  D2 = Table[ Sum[ lambda[[aa, alpha]] (gamma[[aa, b, beta]] 
       - Sum[ SC[[aa, bb, b]] A[[bb, beta]] , {bb, 1, dimG}]), 
    {aa, 1+dimP, dimG}], 
  {alpha, 1, dimM}, {beta, 1, dimM}, {b , 1, dimP}];

 pd =  Table[0 , {i, 4}];
 pd = ReplacePart[ pd, 
	Table[ Sum[ SC[[c, a, b]] LIinv[[a , dcnt]] , {a , 1, dimP} ] , 
          {c, 1, dimP} , {dcnt , 1, dimP}, {b, 1, dimP} ], 1];
 pd = ReplacePart[ pd , D1, 2];
 pd = ReplacePart[ pd , D2, 3];
 pd = ReplacePart[ pd, lambda, 4];
 BKMM[ReplacePart[S, pd, iMOMDYN]]
];

cBaseCoefficients[BKMM[S_]] ^:= Module[
  { pd , coeff, A, l1, l2, LIinv, lambda, D1, D2, r, k, m, n, index} ,

  index = 1;
  r = gLocal[gBase[BKMM[S]]];
  k = gDimGroup[BKMM[S]];
  n = gDimMom[BKMM[S]];
  m = gDimBase[BKMM[S]];
  LIinv = gLIinv[BKMM[S]];
  A = gConnection[BKMM[S]];
  B = gCurvature[BKMM[S]];
  l1 = gl1[BKMM[S]];
  l2 = gl2[BKMM[S]];

  If[ n > 0 , 
  Module[ {}, 
    lambda = Table[ l1[[aa , alpha]] 
      - Sum[ l2[[aa, bb]] A[[bb , alpha]] , {bb, 1, k}],
    {aa, 1, k}, {alpha, 1, m}];
    pd = Table [0, {i, 5}];
    coeff = gMomCoefficients[BKMM[S]];
    D1 = coeff[[2]];
    D2 = coeff[[3]];

    coeff = Table[ D[ LIinv[[a,b]] , r[[alpha]] ] , {a, 1, n} , {b, 1, n},
      {alpha, 1, m}];
    pd = ReplacePart[ pd, coeff, index++];

    coeff = Table[ Sum[D1[[c, b, alpha]] LIinv[[b, a]] , {b, 1, n}] ,
      {c, 1,n} , {a, 1, n}, {alpha, 1, m}];
    pd = ReplacePart[ pd, coeff, index++];

    pd = ReplacePart[ pd, B , index++];

    coeff = Table[ Sum[ D2[[beta, alpha, b]] LIinv[[b, c]] , {b, 1, n}] ,
      {c, 1, n}, {beta, 1, m}, {alpha, 1, m}];
    pd = ReplacePart[ pd, coeff, index++];

    coeff = Table[ Sum[ lambda[[aa, alpha]] B[[aa, beta, gammacnt]] , 
      {aa, 1, k}] , {alpha, 1, m}, {beta, 1, m}, {gammacnt, 1, m} ];
    pd = ReplacePart[ pd, coeff, index++];
  ] ,  (* else *)
    pd = { {}, {} , B, {}, {} }; 
  ];

  BKMM[ReplacePart[S, pd, iBASEDYN]]
];

cVariation[BKMM[S_], p_, t_] ^:= Module[
  {m, n, k, r, rp, e, ec, R, z, Z, Zc, L, A, LIinv, pd, tvar} , 

  m = gDimBase[BKMM[S]];
  n = gDimMom[BKMM[S]];
  k = gDimGroup[BKMM[S]];
  r = gLocal[gBase[BKMM[S]]];
  e = gLocal[gBasis[BKMM[S]]];
  L = gLagrangian[BKMM[S]];

  R = Table[ r[[i]][t] , {i, 1, m} ];
  R' = Table[ r[[i]]'[t] , {i, 1, m}];
  r' = Table[ rp[i] , {i, 1, m}];
  z = Join[r, r', e];
  Z = Join[R, R', Table[ e[[i]][t], {i, 1, k}]];
  P = Table[ p[i][t] , {i, n}];

  For[ i=Length[r] , i>0 , i-- , L = L /. {r[[i]]' -> rp[i]}];
  For[ i=1 , i<=Length[Z] , i++ , L = L /. {z[[i]] -> Z[[i]]} ];

  A = gConnection[BKMM[S]];
  For[ i=Length[r] , i>0 , i-- , A = A /. {r[[i]] -> R[[i]]}];
  LIinv = BKMM[fGetLIinv, S];
  LIinv = LinearAlgebra`MatrixManipulation`TakeColumns[LIinv, 
    gDimMom[BKMM[S]] ];
  For[ i=Length[r] , i>0 , i-- , LIinv = LIinv /. {r[[i]] -> R[[i]]}];
  ec = -A.R' + LIinv.P;
  Zc = Join[R, R', ec];

  For[ i=1 , i<=Length[e] , i++ , L = L /. {e[[i]][t] -> ec[[i]]} ];
  pd = Table[Dt[D[L, Zc[[i + m]]], t], {i, 1, m}];

  tvar = pd;
  pd = Table[D[L, Zc[[i]] ] , {i, 1, m}];
  tvar = Append[{tvar}, pd];

  BKMM[ReplacePart[S, tvar, iBASEVAR]]
];

cDynamics[BKMM[S_], t_] ^:= Module[
  {m, n, k, r, R, p, P, lc, dLIinv, D1LIinv, B, D2LIinv, K, SCI, D1, D2, 
   var1, var2, pd},

  m = gDimBase[BKMM[S]];
  n = gDimMom[BKMM[S]];
  r = gLocal[gBase[BKMM[S]] ];
  R = Table[ r[[i]][t] , {i, 1, m} ];
  p = gMomentum[BKMM[S]];
  lc = gLagrangian[BKMM[S]];
  pd = gBaseCoefficients[BKMM[S]];
  For[ i=Length[r] , i>0 , i-- , pd = pd /. {r[[i]] -> R[[i]]}];
  dLIinv = Part[pd, 1];
  D1LIinv = Part[pd, 2];
  B = Part[pd, 3];
  D2LIinv = Part[pd, 4];
  K = Part[pd, 5];
  pd = gMomCoefficients[BKMM[S]];
  For[ i=Length[r] , i>0 , i-- , pd = pd /. {r[[i]] -> R[[i]]}];
  SCI = Part[pd,1];
  D1 = Part[pd,2];
  D2 = Part[pd,3];
  pd = gBaseVariation[BKMM[S]];
  var1 = Part[pd,1];
  var2 = Part[pd,2];
  
  pd = Table[ var1[[i]] - var2[[i]]  + 
    Sum[ (dLIinv[[a,b,i]] + D1LIinv[[a,b,i]]) p[[a]][t] p[[b]][t],
      {a,1,n}, {b,1,n} ] 
    + Sum[ B[[c, i, beta]] p[[c]][t] r[[beta]]'[t] + D2LIinv[[c, beta, i]] ,
	{c, 1, n}, {beta, 1, m}] == 0
    , {i, 1, m}];
  pd = Join[ pd, Table[ p[[a]]'[t]  
    - Sum[ SCI[[b, c, a]] p[[b]][t] p[[c]][t], {b, 1, n}, {c, 1, n}]
    - Sum[ D1[[b, a, alpha]] p[[b]][t] r[[alpha]]'[t] ,  {b,1,n}, {alpha,1,m}]
    - Sum[ D2[[alpha, beta, a]] r[[alpha]]'[t] r[[beta]]'[t] , 
      {alpha, 1, m}, {beta, 1, m} ] == 0 
    , {a , 1, n} ] ];
  BKMM[ReplacePart[S, pd, iDYNAMICS]]
];


gMomDynamics[BKMM[S_]] ^:= Module[ {tMomDyn, dimM, dimP} , 
  dimM = gDimBase[BKMM[S]];
  dimP = gDimMom[BKMM[S]];
  TMomDyn = S[[iDynamics]];
  Take[TMomDyn , {dimM+1, dimM+dimP}]
];


gLDP[BKMM[S_], t_] := Module[
  {m, n, k, r, R, p, P, lc, dLIinv, D1LIinv, B, D2LIinv, K, SCI, D1, D2, 
   var1, var2, pd},

  m = gDimBase[BKMM[S]];
  n = gDimMom[BKMM[S]];
  r = gLocal[gBase[BKMM[S]] ];
  R = Table[ r[[i]][t] , {i, 1, m} ];
  p = gMomentum[BKMM[S]];
  lc = gLagrangian[BKMM[S]];
  pd = gBaseCoefficients[BKMM[S]];
  For[ i=Length[r] , i>0 , i-- , pd = pd /. {r[[i]] -> R[[i]]}];
  dLIinv = Part[pd, 1];
  D1LIinv = Part[pd, 2];
  B = Part[pd, 3];
  D2LIinv = Part[pd, 4];
  K = Part[pd, 5];
  pd = gMomCoefficients[BKMM[S]];
  For[ i=Length[r] , i>0 , i-- , pd = pd /. {r[[i]] -> R[[i]]}];
  SCI = Part[pd,1];
  D1 = Part[pd,2];
  D2 = Part[pd,3];
  pd = gBaseVariation[BKMM[S]];
  var1 = Part[pd,1];
  var2 = Part[pd,2];
  
  pd = Table[ var1[[i]] - var2[[i]]  + 
    Sum[ (dLIinv[[a,b,i]] + D1LIinv[[a,b,i]]) p[[a]][t] p[[b]][t],
      {a,1,n}, {b,1,n} ] 
    + Sum[ B[[c, i, beta]] p[[c]][t] r[[beta]]'[t] + D2LIinv[[c, beta, i]] ,
	{c, 1, n}, {beta, 1, m}] 
    , {i, 1, m}];
  pd = Join[ pd, Table[ p[[a]]'[t]  
    - Sum[ SCI[[b, c, a]] p[[b]][t] p[[c]][t], {b, 1, n}, {c, 1, n}]
    - Sum[ D1[[b, a, alpha]] p[[b]][t] r[[alpha]]'[t] ,  {b,1,n}, {alpha,1,m}]
    - Sum[ D2[[alpha, beta, a]] r[[alpha]]'[t] r[[beta]]'[t] , 
      {alpha, 1, m}, {beta, 1, m} ] , {a , 1, n} ] ]
];


(*--------------------fGetEOM--------------------*)
(*  Once the configuration has been setup, this function will
    return the equations of motion.  They will depend on the type
    of system that we are dealing with.

    This is not complete yet and will depend on how the whole thing
    is implemented.  This is a challenge, since the implementation
    is determined by what I've needed up til then.

    For now, the following are possible for fGetEOM:

    1) Principal kinematic case
*)

gEOM[BKMM[S_]] := Module[
  {Pconn},
  Pconn = gConnection[BKMM[S]]
];
  

(*--------------------fEOM--------------------*)
(*
  Give equations of motion for first order nonlinear matrix ODE of the form
    r' = u(r,t)
    g' = g A(r) u(r,t)
  over the given interval. q = (r, g)
*)
gEOM[BKMM[S_], A_, u_, q_ , q0_ , int_] := Module [
  {h, hp, hh, Meqns, Geqns, dimM, dimG} , 
  dimG = SE2[fDim];
  dimM = Length[q] - dimG;
  hh = SE2[fInit,q[[dimM+1]],q[[dimM+2]],q[[dimM+3]]];
  h  = SE2[fInit,q[[dimM+1]][int[[1]]],q[[dimM+2]][int[[1]]],
    q[[dimM+3]][int[[1]]]];
  hp = SE2[fInit,q[[dimM+1]]'[int[[1]]],q[[dimM+2]]'[int[[1]]],
    q[[dimM+3]]'[int[[1]]]];
  Geqns = Join[ Table[hp[[i]] == ((SE2[fT,h]).A.u)[[i]] , {i,dimG}], 
	        Table[hh[[i]][int[[2]]] == q0[[dimM+i]] , {i,dimG}] ];
  Meqns = Join[ Table[ q'[[i]][int[[1]]] == u[[i]] , {i, dimM}],
    Table[ q[[i]][int[[2]]] == q0[[i]], {i,dimM}] ];

  Join[Meqns, Geqns]
]

(*--------------------fDSolve--------------------*)
(*
  Solve first order nonlinear matrix ODE of the form
    r' = u(r,t)
    g' = g A(r) u(r,t)
  over the given interval. q = (r, g)
*)
BKMM[fDSolve, A_, u_, q_ , q0_ , int_] := Module [
  {h, hp, hh, Meqns, Geqns, dimM, dimG} , 
  dimG = SE2[fDim];
  dimM = Length[q] - dimG;
  hh = SE2[fInit,q[[dimM+1]],q[[dimM+2]],q[[dimM+3]]];
  h  = SE2[fInit,q[[dimM+1]][int[[1]]],q[[dimM+2]][int[[1]]],
    q[[dimM+3]][int[[1]]]];
  hp = SE2[fInit,q[[dimM+1]]'[int[[1]]],q[[dimM+2]]'[int[[1]]],
    q[[dimM+3]]'[int[[1]]]];
  Geqns = Join[ Table[hp[[i]] == ((SE2[fT,h]).A.u)[[i]] , {i,dimG}], 
	        Table[hh[[i]][int[[2]]] == q0[[dimM+i]] , {i,dimG}] ];
  Meqns = Join[ Table[ q'[[i]][int[[1]]] == u[[i]] , {i, dimM}],
    Table[ q[[i]][int[[2]]] == q0[[i]], {i,dimM}] ];

  NDSolve[Join[Meqns, Geqns], q, int]
]

(*--------------------Diff--------------------*)

Diff[Q_ , S_, f_ ] := Module[ 
  {q} ,

  q = Join[Q[fGetBase, S], Q[fGetGroup, S]];
  Transpose[Table[ D[f,q[[i]]] , {i, Length[q]}]]
];

(*--Functions that support the configuration object.
    They are not bound to any instance in particular.  --*)

CalcGamma[S_, T_ , Q_, G_] := Module[
  {r, gamma, dimG, dimM, eb},

  dimG = G[LieGroups`fDim];
  dimM = Q[fDimBase, S];
  r = Q[fGetBase, S];
  eb[i_] := Transpose[T][[i]];
  gamma = Simplify[Table[ Inverse[T].D[eb[j], r[[alpha]] ], {j,1,dimG},
    {alpha, 1, dimM}]];
  gamma = Table[ gamma[[j, alpha, i]], {i, 1, dimG}, {j,1,dimG},
  {alpha,1,dimM}]
];

(*--------------------fLie--------------------*)

BKMM[fLie, S_, X_, Y_, q_] := Module[{dimM, dimG, dimP, 
  rX, rY, xX, xY, vX, vY, pX, pY, fX, fY, fZ, Z, ql},
  dimM = BKMM[fDimBase, S];
  dimG = SE2[fDim];
  dimP = BKMM[fDimMom, S];

  Z = Table[ 0 , {i, 1, 2 dimM + dimG + dimP}];
  rX = Table[ X[[i]], {i, 1, dimM}];
  xX = Table[ X[[i]], {i, 1+dimM, dimM + dimG}];
  vX = Table[ X[[i]], {i, dimM+dimG, 2 dimM + dimG}];
  pX = Table[ X[[i]], {i, 2 dimM + dimG, 2 dimM + dimG + dimP}];

  rY = Table[ Y[[i]], {i, 1, dimM}];
  xY = Table[ Y[[i]], {i, 1+dimM, dimM + dimG}];
  vY = Table[ Y[[i]], {i, dimM+dimG, 2 dimM + dimG}];
  pY = Table[ Y[[i]], {i, 2 dimM + dimG, 2 dimM + dimG + dimP}];

  fX = Join[ rX, vX, pX];
  fY = Join[ rX, vX, pX];
  ql = Join[ Table[ q[[i]], {1, dimM} ] , 
    Table[ q[[i]], {1+dimM + dimG, 2 dimM + dimG + dimP}] ];

  Z = Lie[X, Y, q];
  fZ = Join[ Table[ 0, {i, 1, dimM}] , se2[xX] ^ se2[xY], 
    Table[ 0, {i, 1+dimM+dimG, 2 dimM + dimG + dimP}] ];
  Z + fZ 
 ]


Format[BKMM[object_]] ^:= object[[iNAME]];

End[];



EndPackage[];
