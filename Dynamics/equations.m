(*--------------------Equations of Motion--------------------*)
(*

  This package implements an object that keeps track of the
  equations of motion for a dynamical system.  These are typically
  tied in to a particular manifold object and do not necessarily
  function alone.

  Author: Patricio Vela
  Date: June 7, 2002

  Notes: At some point, the different objects need to be separated.

*)
BeginPackage["Equations`"];

Needs["Objects`","mathematica/libs/objects.m"];
Needs["Euclidean`","mathematica/libs/euclidean.m"];
Needs["Manifolds`","mathematica/libs/manifolds.m"];
Needs["Vectors`","mathematica/libs/vectors.m"];
Needs["Bundles`","mathematica/libs/bundles.m"];
Needs["Tangents`","mathematica/libs/tangents.m"];
Needs["TangentManifolds`","mathematica/libs/tmanifolds.m"];
Needs["LieGroups`","mathematica/libs/liegroups.m"];
Needs["LagrangianMechanics`","mathematica/libs/lagrangian.m"];

EulerLagrange::usage=
  "Object which stores the Euler-Lagrange equations of motion.";
EulerPoincare::usage=
  "Object which stores the Euler-Poincare equations of motion.";
Hamel::usage=
  "Object which stores Hamel's equations of motion.";
Kinematic::usage=
  "Object which is used to describe kinematic constraints.";
PKinematic::usage=
  "Object which is used to describe G-invariant kinematic constraints.";

(*-----Lagrange-d'Alembert-Poincare' equations-----*)


(*-----Lagrangian functions-----*)

dLagrangian::usage="dLagrangian[Obj, Lag] \n Define the Lagrangian.";
gLagrangian::usage="gLagrangian[Obj] \n Returns the Lagrangian.";

oReduceLagrangian::usage="oReduceLagrangian[Onj, xi] \n Reduce the Lagrangian.";

(*-----Lie algebra basis-----*)

dBasis::usage="dBasis[Obj, basis]";
gBasis::usage="gBasis[Obj]";

(*-----General object functions-----*)

cEquations::usage="cEquations[Obj, Forces, etc.]";
gReconstruction::usage="gReconstruction[Obj]";
gEquations::usage="gEquations[Obj]";
rEquations::usage="rEquations[Obj]";
pEquations::usage="pEquations[Obj]";
cAcceleration::usage="cAcceleration[Obj]";
cFeedbackLin::usage="cFeedbackLin[Obj]";

(*-----External Forcing-----*)

dForces::usage="dForces[Obj, Forces]";
gForces::usage="gForces[Obj] \n Returns external forces on the system.";
dControls::usage="dControls[Obj, Controls]";
dAccels::usage="dAccels[Obj, accel]";

NDSolve::usage="NDSolve[Obj, q0, {t, t0, t1}, optional params]";
rNDSolve::usage="rNDSolve[Obj, q0, {t, t0, t1}, optional params]";

(*-----For use with Bundled Equations of Motion-----*)
(*--These have base and fiber spaces--*)

dBaseDynamics::usage="dBaseDynamics[Obj, r, t]";

Begin["`Private`"];

(*-----The different indices for the equation objects. -----*)

iQ		= 1;
iQP		= 2;
iLAG		= 3;
iMETRIC		= 4;
iFORCES		= 5;
iLEQ		= 6;
iREQ		= 7;

lenEulerLagrange = 7;

iG		= 1;
iGP		= 2;
iXI		= 2;
iMOM		= 2;

iBASIS		= 8;

lenEulerPoincare = lenEulerLagrange+1;

iMETRIC1	= iMETRIC;
iFORCES1	= iFORCES;
iLEQ1		= iLEQ;
iREQ1		= iREQ;

iQ2		= 9;
iQP2		= 10;
iMETRIC2	= 11;
iCOUPLE		= 12;
iFORCES2	= 13;
iLEQ2		= 14;
iREQ2		= 15;
iFULLMETRIC	= 16;

lenSplitEulerLagrange = lenEulerPoincare+8;

iR		= iQ2;
iRP		= iQP2;

lenHamel	= lenSplitEulerLagrange;

(*--------------------oInit--------------------*)

oInit[EulerLagrange, q_, qp_] ^:= Module[{obj},

  obj = Table[ 0, {i, lenEulerLagrange}];

  obj = ReplacePart[obj, q, iQ];
  obj = ReplacePart[obj, qp, iQP];
  obj = ReplacePart[obj, Table[0, {i, gDim[qp]}, {j, gDim[qp]}], iMETRIC];
  obj = ReplacePart[obj, oInit[Covector, Table[0, {i, gDim[qp]}]], iFORCES];

  EulerLagrange[obj]
];

oInit[EulerLagrange, qp_] ^:= 
  oInit[EulerLagrange, gBase[qp], gVector[qp]];

oInit[EulerPoincare, g_, v_] ^:= Module[{obj},

  obj = Join[oContents[oInit[EulerLagrange, g, v]], 
    Table[0, {i, lenEulerPoincare-lenEulerLagrange}] ];

  obj = ReplacePart[obj, 
    oInit[gDualLieAlgebra[g], Table[0, {i, gDim[v]}]], iFORCES];
  dBasis[EulerPoincare[obj], IdentityMatrix[gDim[g]]]
];

oInit[EulerPoincare, q_] ^:= 
  oInit[EulerPoincare, gBase[q], gFiber[q]];

oInit[Hamel, r_, rp_ , g_, v_] ^:= Module[{obj},

  obj = Join[oContents[oInit[EulerPoincare, g, v]],
    Table[Null, {i, lenHamel - lenEulerPoincare}] ];

  obj = ReplacePart[obj, r, iR];
  obj = ReplacePart[obj, rp, iRP];

  Hamel[obj]
];
  
(*--------------------Basis--------------------*)

dBasis[EulerPoincare[obj_], basis_] ^:= 
  EulerPoincare[ReplacePart[obj, basis, iBASIS]];

dBasis[Hamel[obj_], basis] := 
  Hamel[oContents[dBasis[EulerPoincare[obj, basis]]]];

gBasis[EulerPoincare[obj_]] ^:= Part[obj, iBASIS];

gBasis[Hamel[obj_]] := gBasis[EulerPoincare[obj]];

oReduceBasis[Hamel[obj_]] ^:= Module[{Grp, xi},
  xi = Part[obj, iGP];

  (*--Incomplete--*)
];

(*--------------------Lagrangian--------------------*)

dLagrangian[EulerLagrange[obj_], Lag_] ^:=
  EulerLagrange[ReplacePart[obj, Lag, iLAG]];

dLagrangian[EulerPoincare[obj_], Lag_] ^:= 
  EulerPoincare[ReplacePart[obj, Lag, iLAG]];

dLagrangian[Hamel[obj_], Lag_] ^:= 
  Hamel[ReplacePart[obj, Lag, iLAG]];

gLagrangian[EulerLagrange[obj_]] ^:= Part[obj, iLAG];
gLagrangian[EulerPoincare[obj_]] ^:= Part[obj, iLAG];
gLagrangian[Hamel[obj_]] ^:= Part[obj, iLAG];

oReduceLagrangian[EulerPoincare[obj_], rgp_] ^:= Module[
  {lag, gp, g, X, subs, nobj},

  g = Part[obj, iG];
  gp = Part[obj, iGP];
  X = oInit[gTangent[g], g, gp];
  lag = Part[obj, iLAG];

  (*--Might want to check reduce algebra versus given algebra--*)
  lag = cReduce[lag, X, rgp];

  nobj = ReplacePart[obj, rgp, iGP];
  nobj = ReplacePart[nobj, lag, iLAG];
  EulerPoincare[nobj]
];

oReduceLagrangian[Hamel[obj_], rgp_] ^:= Module[{nobj},
  
   nobj = oContents[oReduceLagrangian[EulerPoincare[obj], rgp]];
   Hamel[nobj]
];

(*--------------------cEquations--------------------*)

cEquations[EulerLagrange[obj_]] ^:= Module[
  {Lag, LM, LMinv, q, qp, v, qpp, nobj, leq, req, forces, dimQ},

  q = Part[obj, iQ];
  qp = Part[obj, iQP];
  v = gLocal[qp];
  dimQ = gDim[q];
  qpp = Table[ v[[i]]' , {i, dimQ} ];

  Lag = Part[obj, iLAG];

  LM = Simplify[gLagrangianMetric[Lag, v]];
  leq = LM.v;

  forces = Part[obj, iFORCES];

  (*--configuration space EOM--*)
  req = 
    (1/2) (Table[D[v.LM.v, gLocal[q][[i]] ], {i, dimQ}]);
 
  req -= 
    v.(Table[D[LM.v, gLocal[q][[i]] ], {i,dimQ}]);

  req += gLocal[forces];

  leq = LM.qpp;

  nobj = ReplacePart[obj, leq, iLEQ];
  nobj = ReplacePart[nobj, req, iREQ];
  nobj = ReplacePart[nobj, LM , iMETRIC];
  EulerLagrange[nobj]
];

cEquations[EulerPoincare[obj_]] ^:= Module[
  {Lag, MGG, MinvGG, g, v, vp, nobj, leq, req, SC, noforces, basis},

  g = Part[obj, iG];
  v = Part[obj, iGP];
  vp = Table[ gLocal[v][[i]]', {i, gDim[v]} ];
  basis = gBasis[EulerPoincare[obj]];
  SC = StructureConstants[ gLieAlgebra[gGroup[g]], basis]; 
  Lag = Part[obj, iLAG];

  MGG = Simplify[gLagrangianMetric[Lag, gLocal[v]]];
  leq = MGG.vp;

  If[ oCompare[ oType[v], gLieAlgebra[gGroup[g]] ] ,
    req = (MGG.gLocal[v]).(SC.gLocal[v]);
  , If[ oCompare[ oType[v], gDualLieAlgebra[gGroup[g]] ] ,
    req = (gLocal[v]).(SC.MGG.gLocal[v]);
  ] ];

  noforces = oInit[gDualLieAlgebra[gGroup[g]], Table[0 , {i, gDim[v]}] ];

  nobj = ReplacePart[obj, Lag, iLAG];
  nobj = ReplacePart[nobj, MGG, iMETRIC];
  nobj = ReplacePart[nobj, leq, iLEQ];
  nobj = ReplacePart[nobj, req, iREQ];
  nobj = ReplacePart[nobj, noforces, iFORCES];
  EulerPoincare[nobj]
]

cEquations[EulerPoincare[obj_], Forces_] ^:= Module[
  {nobj},

  nobj = cEquations[EulerPoincare[obj]];
  dForces[nobj, Forces]
];

cEquations[Hamel[obj_], Forces_] ^:= 
Module[
  {L, MSS, MSG, MGS, MGG, dimQ, dimG, basis, SC,  rp,
  rpp, xip, leq1, leq2, req1, req2, r, g, xi, G, dimM, nobj},

  nobj = oContents[cEquations[EulerPoincare[obj], Forces[[2]] ] ];

  r = Part[nobj, iR];
  rp = Part[nobj, iRP];
  g = Part[nobj, iG];
  gp = Part[nobj, iGP];

  dimG = gDim[g];
  dimM = gDim[r];
  dimQ = dimM + dimG;

  rpp = Table[ gLocal[rp][[i]]', {i, dimM}];
  xip = Table[ gLocal[gp][[i]]', {i, dimG}];
  L = Part[nobj, iLAG];
  lm = Simplify[gLagrangianMetric[L, Join[gLocal[rp], gLocal[gp]] ] ];

  basis = gBasis[Hamel[obj]];
  SC = StructureConstants[ gLieAlgebra[gGroup[g]], basis]; 

  MGG = Part[nobj, iMETRIC1];
  (*MGS = Table[ lm[[i, j]] , {i, dimM + 1, dimQ}, {j, dimM}];*)
  MSS = Simplify[Table[ lm[[i, j]] , {i, dimM}, {j, dimM}]];
  MSG = Simplify[Table[ lm[[i, j]] , {i, dimM}, {j, dimM+1, dimQ}]];

  req1 = Part[nobj, iREQ1];
  leq1 = Part[nobj, iLEQ1];

  (*--Shape space EOM--*)
  req2 = 
    (1/2) (Table[D[MSS, gLocal[r][[i]] ], {i, dimM}].gLocal[rp]).gLocal[rp]
     + (Table[D[MSG, gLocal[r][[i]] ], {i, dimM}].gLocal[gp]).gLocal[rp]
     + (1/2) (Table[D[MGG, gLocal[r][[i]] ], {i,dimM}].gLocal[gp]).gLocal[gp];
  (*   + gLocal[Forces[[1]]];*)
 
  req2 -= 
    (gLocal[rp].Table[D[MSG, gLocal[r][[i]] ], {i, dimM}]).gLocal[gp]
    + (gLocal[rp].Table[D[MSS, gLocal[r][[i]] ], {i,
    dimM}]).gLocal[rp];

  leq2 = MSG.xip + MSS.rpp;

  (*--Group space EOM--*)
  (*req1 -= (MGG.xi).(SC.xi) + (MGG.xi).(gamma.rp) + (rp.MSG).(SC.xi)
    + (rp.MSG).(gamma.rp) + Forces[[2]];*)

  req1 += (gLocal[rp].MSG).(SC.gLocal[gp]);
  req1 -= (gLocal[rp].Table[D[MGG, gLocal[r][[i]]], {i, dimM}]).gLocal[gp] 
    +  (gLocal[rp].(gLocal[rp].Table[D[MSG, gLocal[r][[i]]], {i, dimM}]));

  leq1 +=  rpp.MSG; 

  nobj = ReplacePart[nobj, leq1, iLEQ1];
  nobj = ReplacePart[nobj, leq2, iLEQ2];
  nobj = ReplacePart[nobj, req1, iREQ1];
  nobj = ReplacePart[nobj, req2, iREQ2];
  nobj = ReplacePart[nobj, MSS , iMETRIC2];
  nobj = ReplacePart[nobj, MSG , iCOUPLE]; 
  nobj = ReplacePart[nobj, lm , iFULLMETRIC]; 
  nobj = ReplacePart[nobj, Forces[[1]], iFORCES2]; 
  Hamel[nobj]
];


(*--------------------gReconstruction--------------------*)

gReconstruction[EulerPoincare[obj_]] ^:=  Module[{gp, g, lm} ,

  gp = Part[obj, iGP];
  g = Part[obj, iG];
  lm = Part[obj, iMETRIC];
  
  If[ oCompare[oType[gp], gDualLieAlgebra[g]] ,
    gp = oInit[gLieAlgebra[g] , Inverse[lm].gLocal[gp] ] ];
    
  gEOM[Part[obj, iG], gp]
]

gReconstruction[EulerPoincare[obj_], t_] ^:= Module[{eom, g, h, hp, lm},

  g = Part[obj, iG];
  lm = Part[obj, iMETRIC];
  hp = oInit[oType[Part[obj, iGP]] , Table[ gLocal[Part[obj, iGP]][[i]][t],
    {i, gDim[Part[obj, iGP]]} ] ];

  If[ oCompare[oType[hp], gDualLieAlgebra[g]] ,
    hp = oInit[gLieAlgebra[g] , Inverse[lm].gLocal[hp] ] ];

  eom = gEOM[g ,hp];
  eom /. Join[Table[gLocal[g][[i]]' -> gLocal[g][[i]]'[t],{i, gDim[g]}],
    Table[ gLocal[g][[i]] -> gLocal[g][[i]][t] , {i, gDim[g]}] ]
]

(*--------------------gEquations--------------------*)

gEquations[EulerLagrange[object_]] ^:= Module[
  {leq, req, force} ,

  leq = Part[object, iLEQ];
  req = Part[object, iREQ];
  forces = gLocal[Part[object, iFORCES]];

  Table[ leq[[i]] == req[[i]] + forces[[i]], {i, Length[leq]} ]
]

gEquations[EulerLagrange[object_], time_] ^:= Module[
  {eqns, subs, q, qp, v, dim},

  eqns = gEquations[EulerLagrange[object]];

  q = gLocal[Part[object, iQ]];
  qp = gLocal[Part[object, iQP]];
  dim = gDim[Part[object, iQ]];

  subs = Join[ Table[qp[[i]]' -> qp[[i]]'[time] , {i, dim}],
      Table[ qp[[i]] -> qp[[i]][time] , {i, dim}], 
      Table[ q[[i]] -> q[[i]][time] , {i, dim}] ];

  eqns /. subs
];

gEquations[EulerPoincare[object_]] ^:= Module[
  {leq, req, forces} , 

  leq = Part[object, iLEQ];
  req = Part[object, iREQ];
  forces = gLocal[Part[object, iFORCES]];

  Table[ leq[[i]] == req[[i]] + forces[[i]] , {i, Length[leq]} ]
]

gEquations[EulerPoincare[object_], time_] ^:= Module[
  {eqns, subs, g, gp, v, dim},

  eqns = gEquations[EulerPoincare[object]];

  g = gLocal[Part[object, iG]];
  gp = gLocal[Part[object, iGP]];
  dim = gDim[Part[object, iG]];

  subs = Join[ Table[gp[[i]]' -> gp[[i]]'[time] , {i, dim}],
      Table[ gp[[i]] -> gp[[i]][time] , {i, dim}] ];

  eqns /. subs
];

rEquations[EulerPoincare[obj_]] ^:=
  Join[gReconstruction[EulerPoincare[obj]], 
    gEquations[EulerPoincare[obj]] ];

rEquations[EulerPoincare[obj_], t_] ^:=
  Join[gReconstruction[EulerPoincare[obj],t], 
    gEquations[EulerPoincare[obj],t] ];

gEquations[Hamel[object_]] ^:= Module[
  {leq1, req1, leq2, req2, f1, f2, LHS, RHS},

  If[ Length[object] == lenHamel ,
    leq1 = Part[object, iLEQ1];
    leq2 = Part[object, iLEQ2];
    req1 = Part[object, iREQ1];
    req2 = Part[object, iREQ2];
    f1   = gLocal[Part[object, iFORCES1]];
    f2   = gLocal[Part[object, iFORCES2]];

    LHS = Join[leq1, leq2];
    RHS = Join[req1 + f1, req2 + f2];
    Table[LHS[[i]] == RHS[[i]], {i, Length[LHS]} ]
  ]
]

gEquations[Hamel[object_], time_] ^:= Module[
  {f1, f2, leq1, req1, leq2, req2, LHS, RHS, subs, r, xi, rp, rps, rpps, eta, 
  dimM, dimG},

  If[ Length[object] == lenHamel ,
    leq1 = Part[object, iLEQ1];
    leq2 = Part[object, iLEQ2];
    req1 = Part[object, iREQ1];
    req2 = Part[object, iREQ2];
    f1   = gLocal[Part[object, iFORCES1]];
    f2   = gLocal[Part[object, iFORCES2]];

    xi = gLocal[Part[object, iXI]];
    r = gLocal[Part[object, iR]];
    rp = gLocal[Part[object, iRP]];
    dimM = Length[r];
    dimG = Length[xi];

    subs = Join[ Table[rp[[i]]' -> rpps[i] , {i, dimM}],
      Table[ rp[[i]] -> rps[i] , {i, dimM}] ,
      Table[ xi[[i]]' -> eta[i] , {i, dimG}] ];
    
    LHS = Join[leq1, leq2] /. subs;
    RHS = Join[req1 + f1, req2 + f2] /. subs;
    
    subs = Join[ Table[ r[[i]] -> r[[i]][time] , {i, dimM}],
      Table[ rps[i] -> rp[[i]][time] , {i, dimM}],
      Table[ rpps[i] -> rp[[i]]'[time] , {i, dimM}],
      Table[ xi[[i]] -> xi[[i]][time] , {i, dimG}],
      Table[ eta[i] -> xi[[i]]'[time] , {i, dimG}] ];
    LHS = LHS /. subs;
    RHS = RHS /. subs;

    Table[LHS[[i]] == RHS[[i]] , {i, Length[LHS]} ]
  ]
]

pEquations[Hamel[object_]] ^:= Module[
  {leq1, req1, leq2, req2, LHS, RHS},

  If[ Length[object] == lenHamel,
    leq1 = Part[object, iLEQ1];
    leq2 = Part[object, iLEQ2];
    req1 = Part[object, iREQ1];
    req2 = Part[object, iREQ2];

    LHS = Join[leq1, leq2];
    RHS = Join[req1, req2];
    MatrixForm[ Table[LHS[[i]] == RHS[[i]] , {i, Length[LHS]} ] ]
  ]
]

(*--------------------dForces,dControls,dAccels--------------------*)

dForces[EulerLagrange[obj_], Forces_] ^:= 
  EulerLagrange[ReplacePart[obj, Forces, iFORCES]];

dForces[EulerPoincare[obj_], Forces_] ^:=
  EulerPoincare[oContents[dForces[EulerLagrange[obj], Forces] ] ];

dForces[Hamel[obj_], Forces_] ^:= Module[{nobj},
  nobj = oContents[dForces[EulerPoincare[obj], Forces[[2]] ] ];
  Hamel[ReplacePart[nobj, Forces[[1]] , iFORCES2] ]
];

gForces[EulerLagrange[obj_]] ^:= Part[obj, iFORCES];

gForces[EulerPoincare[obj_]] ^:= gForces[EulerLagrange[obj]];

gForces[Hamel[obj_], Forces_] ^:= 
  Join[gForces[EulerPoincare[obj]], Part[nobj,iFORCES2] ];


dControls[Obj_, u_] := dForces[Obj, u];

aControls[Hamel[obj_], u_] ^:= Module[{F1, F2, nobj},

  F1 = Part[obj, iFORCES] + u[[1]];
  F2 = Part[obj, iFORCES2] + u[[2]];

  nobj = ReplacePart[obj, F1, iFORCES];
  nobj = ReplacePart[nobj, F2, iFORCES2];

  Hamel[nobj];
];

dAccels[Obj_, u_] := dForces[Obj, u];


(*--------------------cAcceleration--------------------*)
(*
  Compute the accelerations given the determined equations
  of motion from the Euler Lagrange operator or an equivalent
  operator.
*)

cAcceleration[EulerLagrange[obj_]] ^:= Module[
  {v, vp, eqns, solved, nobj},

  v = Part[obj, iQP];
  vp = Table[ gLocal[v][[i]]' , {i, gDim[v]} ];

  eqns = gEquations[EulerLagrange[obj]];
  solved = Solve[ eqns, vp ];

  eqns = Flatten[ vp /. solved];

  nobj = ReplacePart[obj, vp, iLEQ];
  nobj = ReplacePart[nobj, eqns, iREQ];
  nobj = ReplacePart[nobj, oInit[Covector, Table[0, {i,gDim[v]}]], iFORCES];
  EulerLagrange[nobj]
];

cAcceleration[EulerPoincare[obj_]] ^:= Module[
  {v, vp, eqns, solved, nobj},

  v = Part[obj, iQP];
  vp = Table[ gLocal[v][[i]]', {i, gDim[v]} ];

  eqns = gEquations[EulerPoincare[obj]];
  solved = Solve[ eqns , vp ];

  eqns = Flatten[vp /. solved];

  nobj = ReplacePart[obj, vp, iLEQ];
  nobj = ReplacePart[nobj, eqns, iREQ];
  nobj = ReplacePart[nobj, oInit[gDual[v], Table[0, {i,gDim[v]}]], iFORCES];
  EulerPoincare[nobj]
];

cAcceleration[Hamel[obj_]] ^:= Module[
  {v1, v2, v1p, v2p, eqns, solved, nobj},

  v1 = Part[obj, iGP];
  v2 = Part[obj, iRP];
  v1p = Table[ gLocal[v1][[i]]' , {i, gDim[v1]}];
  v2p = Table[ gLocal[v2][[i]]' , {i, gDim[v2]}];

  eqns = gEquations[Hamel[obj]];
  solved = Solve[eqns, Join[v1p, v2p]];

  eqns = Flatten[v1p /. solved];
  nobj = ReplacePart[obj, v1p, iLEQ1];
  nobj = ReplacePart[nobj, eqns, iREQ1];

  eqns = Flatten[v2p /. solved];
  nobj = ReplacePart[nobj, v2p, iLEQ2];
  nobj = ReplacePart[nobj, eqns, iREQ2];

  nobj = ReplacePart[nobj, 
    oInit[gDual[v1], Table[0, {i,gDim[v1]}]], iFORCES1];
  noobj = ReplacePart[nobj, 
    oInit[Covector, Table[0, {i,gDim[v2]}]], iFORCES2];
  Hamel[nobj]
];

(*--------------------cFeedbackLin--------------------*)
(*
  Compute the feedback linearized equations of motion
  given computed accelerations.
*)

cFeedbackLin[Hamel[obj_], u_] ^:= Module[
  {req2, nobj},

  req2 = Part[obj, iREQ2];
  nobj = ReplacePart[obj, Table[0, {i, Length[req2]}], iREQ2];

  dControls[Hamel[nobj], u]
];

cPartialFeedbackLin[Hamel[obj_], x_ ] ^:= Module[
  {MinvGG, MinvSS, MGG, MSS, F1, F2, req2, xi, r, xip, rpp, G,
  nobj},


  If[ (Length[obj] != lenHamel) ,
    Hamel[obj]
  ,
    F1 = Part[obj, iREQ1];
    F2 = Part[obj, iREQ2];
    MGG  = Part[obj, iMGG];
    MSG  = Part[obj, iMSG];
    MSS  = Part[obj, iMSS];
    r    = Part[obj, iR];
    xi   = Part[obj, iXI];
    G    = Part[obj, iG];

    xip = Table[ xi[[i]]', {i,Length[xi]}];
    rpp = Table[ r[[i]]'', {i, Length[r]}];
    MinvSS = Inverse[MSS];
    MinvGG = Inverse[MGG];

    req2 = MinvGG.(F2 - (MinvSS.(F1 - MSG.xip)).MSG );

    nobj = ReplacePart[obj, rpp, iLEQ1];
    nobj = ReplacePart[nobj, Table[0, {i, Length[rpp]}], iREQ1];
    nobj = ReplacePart[nobj, xip, iLEQ2];
    noobj = ReplacePart[nobj, req2, iREQ2];
    Hamel[nobj]
  ]
];


(*--------------------NDSolve--------------------*)
(*
  Numerically integrate the equations of motion.
*)

NDSolve[EulerLagrange[obj_], vq0_, time_, params___] ^:= Module[
  {eom, q, qp, q0, qp0, dim},

  q = gLocal[Part[obj, iQ]];
  qp = gLocal[Part[obj, iQP]];
  dim  = gDim[Part[obj, iQP]];
  q0 = gLocal[gBase[vq0]];
  qp0 = gLocal[gVector[vq0]];

  eom = Join[gEquations[EulerLagrange[obj], time[[1]]] ,
    Table[q[[i]]'[time[[1]]] == qp[[i]][time[[1]]] , {i, dim}],
    Table[q[[i]][0] == q0[[i]] , {i, dim}], 
    Table[qp[[i]][0] == qp0[[i]] , {i, dim}] ];

  (*{eom , Join[q, qp], time, params};*)
  NDSolve[eom , Join[q, qp], time, params]
];

NDSolve[EulerPoincare[obj_], q0_, time_, params___] ^:= Module[
  {eom, gp, dim},

  gp = gLocal[Part[obj, iGP]];
  dim  = gDim[Part[obj, iGP]];

  eom = Join[gEquations[EulerPoincare[obj], time[[1]]] ,
    Table[gp[[i]][0] == q0[[i]] , {i, dim}] ];

  NDSolve[eom , gp, time, params]
];

rNDSolve[EulerPoincare[obj_], q0_, time_, params___] ^:= Module[
  {eom, g, gp},

  g = Part[obj, iG];
  gp = Part[obj, iGP];

  eom = Join[ rEquations[EulerPoincare[obj], time[[1]]] , 
    Table[gLocal[g][[i]][time[[2]]] == gLocal[q0][[i]] , {i, gDim[g]} ],
    Table[gLocal[gp][[i-gDim[g]]][time[[2]]] == gLocal[q0][[i]] , 
      {i, 1+gDim[g] , gDim[g]+gDim[gp]} ] ];

  NDSolve[eom, Join[ gLocal[g], gLocal[gp] ], time, params]
];
  
  

NDSolve[Hamel[obj_], q0_, time_] ^:= Module[
  {eom, eom0, eom1, eom2, q, r, g, G, LA, rp, xi, xip, dimM, dimG},

  eom2 = gEquations[Hamel[obj], time[[1]]];
  r = gLocal[Part[obj, iR]];
  xi = gLocal[Part[obj, iXI]];
  g = Part[obj, iG];
  G = oType[g];
  g = gLocal[g];
  LA = gLieAlgebra[G];
  dimM = Length[r];
  dimG = Length[xi];
  rp = gLocal[Part[obj, iRP]];
  q = Join[r, g, xi];
  eom0 = Join[ Table[ r[[i]][0] == q0[[i]] , {i, dimM}],
    (*Table[ g[[i]][0] == q0[[dimM+i]] , {i, dimG}],*)
    Table[ rp[[i]][0] == q0[[dimM+dimG+i]], {i, dimM}],
    Table[ xi[[i]][0] == q0[[2 dimM+dimG+i]], {i, dimG}] ];

  eom1 = gEOM[G[g], LA[Table[xi[[i]][time[[1]]], {i, dimG}]], 
    G[Table[q0[[i]],{i, 1+dimM, dimM+dimG}]], time[[1]]]; 

  eom = Join[ eom0, eom1, eom2 ];

  NDSolve[ eom , q,  time ]
];  

NDSolve[Hamel[obj_], q0_, time_, params_] ^:= Module[
  {eom, eom0, eom1, eom2, q, r, g, G, LA, rp, xi, xip, dimM, dimG},

  eom2 = gEquations[Hamel[obj], time[[1]]];
  r = gLocal[Part[obj, iR]];
  xi = gLocal[Part[obj, iXI]];
  g = Part[obj, iG];
  G = oType[g];
  g = gLocal[g];
  LA = gLieAlgebra[G];
  dimM = Length[r];
  dimG = Length[xi];
  rp = gLocal[Part[obj, iRP]];
  q = Join[r, g, xi];
  eom0 = Join[ Table[ r[[i]][0] == q0[[i]] , {i, dimM}],
    (*Table[ g[[i]][0] == q0[[dimM+i]] , {i, dimG}],*)
    Table[ rp[[i]][0] == q0[[dimM+dimG+i]], {i, dimM}],
    Table[ xi[[i]][0] == q0[[2 dimM+dimG+i]], {i, dimG}] ];

  eom1 = gEOM[LA[Table[xi[[i]][time[[1]]], {i, dimG}]], 
    G[Table[q0[[i]],{i, 1+dimM, dimM+dimG}]], G[g], time]; 

  eom = Join[ eom0, eom1, eom2 ];

  NDSolve[ eom , q,  time, params ]
];  


(*--------------------Format--------------------*)
(*
  Describe the output format of the object.
*)

Format[EulerLagrange[object_]] := "Euler Lagrange equations for " <>
  ToString[object[[iQ]]];

Format[EulerPoincare[object_]] := "Euler Poincare equations for " <>
  ToString[Join[gLocal[object[[iQ]]], gLocal[object[[iQP]]] ] ];

Format[Hamel[object_]] := "Hamel equations for " <> 
  ToString[Join[ gLocal[object[[iR]]], gLocal[object[[iG]]], 
    gLocal[object[[iXI]]] ]];


(*===============================================================*)
(*==============================Kinematic========================*)
(*===============================================================*)



End[];

EndPackage[];
