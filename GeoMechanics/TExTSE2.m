(*  TExTSE2.m

  This package provides the necessary function and manipulations required
  to implement the tangent space to the principal fiber bundle ExSE2.
  The cotangent space is also included, as are the equivalent reduced
  spaces for the case of symmetries.  

  Author:	Patricio Vela
  Date: 	08/08/2002
  MOdified:	08/08/2002
*)

TExTSE2::usage="Provides functionality of tangent to principal bundle Q = R^n x SE(2)";
dExdSE2::usage="Provides functionality of dual reduced tangent to principal bundle Q = R^n x SE(2)";

TExse2::usage="Provides functionality of reduced tangent to principal bundle Q = R^n x SE(2)";
TExdse2::usage="Provides functionality of dual reduced tangent to principal bundle Q = R^n x SE(2)";

(*--Declare functions--*)

gGroupVector::usage="gGroupVector[bvector]";


(*--Define functione--*)

Begin["`Private`"];

(* TExTSE2 is TQ=T(R(m)xSE(2)).  *)

iBASE 		= PrincipalBundle`Private`BASE; 	(* base variables. *)
iGROUP		= PrincipalBundle`Private`FIBER;	(* group variables. *)

lenExSE2	= PrincipalBundle`Private`lenExSE2;

iBASEV		= TBundles`Private`iBASEV;
iGROUPV		= TBundles`Private`iFIBERV;

lenTExTSE2	= TBundles`PRivate`lenTBUNDLE;

(*===================================================================*)
(*-----------------------------TExTSE2-------------------------------*)
(*===================================================================*)

(*--------------------oInit--------------------*)
(*  Initialize the configuration object.
*)

oInit[TExTSE2] ^:= TExTSE2[oContents[oInit[TBundle]]];

oInit[TExTSE2, Euclidean[r_], SE2[g_]] ^:=
  TExTSE2[oContents[oInit[TBundle, Euclidean[r], SE2[g]]]];

oInit[TExTSE2, ExSE2[q_]] ^:=
  TExTSE2[oContents[oInit[ TBundle, Bundle[q] ]]];

oInit[TExTSE2, ExSE2[q_], Vector[v_]] ^:=
  dVector[oInit[TExTSE2, ExSE2[q]], Vector[v]];

oInit[TExTSE2, Euclidean[r_], SE2[g_] , Vector[rp_] , Vector[gp_]] ^:= 
  oInit[TExTSE2, oInit[ExSE2, Euclidean[r], SE2[g]], 
    {Vector[rp], Vector[gp]}];

oInit[TExse2, Euclidean[r_] , Vector[rp_], se2[xi_]] ^:=
  oInit[TBundle];

(*---------------Define and Get components---------------*)

dBase[TExTSE2[v_], ExSE2[q_]] ^:=
  TExTSE2[oContents[dBase[ TBundle[v], Bundle[q] ]]];

dVector[TExTSE2[v_], Vector[rv_], Vector[gv_]] ^:=
  TExSE2[oContents[dVector[ TBundle[v], Vector[rv], Vector[gv] ]]];

dVector[TExTSE2[v_], Vector[vec_]] ^:=
  TExTSE2[oContents[dVector[ TBundle[v], Vector[vec] ]]];

gBase[TExTSE2[v_]] ^:= ExSE2[oContents[gBase[TBundle[v]]]];
gPoint[TExTSE2[v_]] ^:= gBase[TExTSE2[v]];
gVector[TExTSE2[v_]] ^:= gVector[TBundle[v]];
gBaseVector[TExTSE2[bv_]] ^:= gBaseVector[TBundle[bv]];
gFiberVector[TExTSE2[bv_]] ^:= gFiberVector[TBundle[bv]];
gGroupVector[TExTSE2[bv_]] ^:= gFiberVector[TExTSE2[bv]];

gDim[TExTSE2[v_]] ^:= gDim[TBundle[v]];

gLocal[TExTSE2[v]] ^:= gLocal[TBundle[v]];


(*--------------------Manifold Representations--------------------*)

(*
gTM[TExTSE2[obj_]] ^:= Euclidean`TEuclidean;
gTG[TExTSE2[obj_]] ^:= LieGroups`TSE2;

gLieAlgebra[ExSE2[obj_]] ^:= LieGroups`se2;
*)

(*--Base and group variables--*)
(*
dGroup[TExTSE2[obj_], g_] ^:= TExTSE2[ReplacePart[obj, g, iGROUP]];
dBase[TExTSE2[obj_], r_] ^:= TExTSE2[ReplacePart[obj, r, iBASE]];
dBaseVector[TExTSE2[obj_], rp_] ^:= 
  TExTSE2[ReplacePart[obj, rp, iBASEV]];
dGroupVector[TExTSE2[obj_], gp_] ^:= 
  TExTSE2[ReplacePart[obj, gp, iGROUPV]];

gBase[TExTSE2[obj_]] ^:= Part[obj, iBASE];
gGroup[TExTSE2[obj_]] ^:= Part[obj, iGROUP];

gBaseVector[TExTSE2[obj_]] ^:= Part[obj, iBASEV];
gGroupVector[TExTSE2[obj_]] ^:= Part[obj, iGROUPV];

gVector[TExTSE2[obj_]] ^:=
  {gBaseVector[TExTSE2[obj]] , gGroupVector[TExTSE2[obj]]};

gTBase[TExTSE2[obj_]] ^:= 
  oInit[gTM[TExTSE2[obj]], gBase[TExTSE2[obj]], gBaseVector[TExTSE2[obj]]];
gTGroup[TExTSE2[obj]] ^:= 
  oInit[gTG[TExTSE2[obj]], gGroup[TExTSE2[obj]], gGroupVector[TExTSE2[obj]]];
*)


(*--------------------fLie--------------------*)

(*
ExSE2[fLie, S_, X_, Y_, q_] := Module[{dimM, dimG, dimP, 
  rX, rY, xX, xY, vX, vY, pX, pY, fX, fY, fZ, Z, ql},
  dimM = ExSE2[fDimBase, S];
  dimG = SE2[fDim];
  dimP = ExSE2[fDimMom, S];

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
*)

(*===================================================================*)
(*-----------------------------dExdSE2-------------------------------*)
(*===================================================================*)

(*--------------------oInit--------------------*)
(*  Initialize the configuration object.
*)

oInit[dExdSE2] ^:= dExdSE2[oContents[oInit[dBundle]]];

oInit[dExdSE2, Euclidean[r_], SE2[g_]] ^:=
  dExdSE2[oContents[oInit[dBundle, Euclidean[r], SE2[g]]]];

oInit[dExdSE2, ExSE2[q_]] ^:=
  dExdSE2[oContents[oInit[ dBundle, Bundle[q] ]]];

oInit[dExdSE2, ExSE2[q_], Covector[v_]] ^:=
  dCovector[oInit[dExdSE2, ExSE2[q]], Covector[v]];

oInit[dExdSE2, Euclidean[r_], SE2[g_] , Covector[rp_] , Covector[gp_]] ^:= 
  oInit[dExdSE2, oInit[ExSE2, Euclidean[r], SE2[g]], 
    {Covector[rp], Covector[gp]}];

(*---------------Define and Get components---------------*)

dBase[dExdSE2[cv_], ExSE2[q_]] ^:=
  dExdSE2[oContents[dBase[ dBundle[cv], Bundle[q] ]]];

dCovector[dExdSE2[cv_], Covector[rcv_], Covector[gcv_]] ^:=
  dExdSE2[oContents[dVector[ dBundle[cv], Covector[rcv], Covector[gcv] ]]];

dCovector[dExdSE2[cv_], Covector[cvec_]] ^:=
  dExdSE2[oContents[dCovector[ dBundle[cv], Covector[cvec] ]]];

gBase[dExdSE2[cv_]] ^:= ExSE2[oContents[gBase[dBundle[cv]]]];
gPoint[dExdSE2[cv_]] ^:= gBase[dExdSE2[cv]];
gCovector[dExdSE2[cv_]] ^:= gCovector[dBundle[cv]];

gDim[dExdSE2[cv_]] ^:= gDim[dBundle[cv]];

gLocal[dExdSE2[cv]] ^:= gLocal[dBundle[cv]];


(*------------------------------------------------------------------*)
(*-----------------Mathematica Functions and Operations-------------*)
(*------------------------------------------------------------------*)


(*====================  Arithmetic  ====================*)

(*----------TExTSE2----------*)

alpha_ * TExTSE2[v_] ^:=
  oInit[TExTSE2, gBase[TExTSE2[v]], alpha * gVector[TExTSE2[v]] ];

TExTSE2[v1_] + TExTSE2[v2_] ^:=
  TExTSE2[oContents[TBundle[v1] + TBundle[v2] ]];

(*----------dExdSE2----------*)

alpha_ * dExdSE2[cv_] ^:=
  oInit[dExdSE2, gBase[dExdSE2[cv]], alpha * gCovector[dExdSE2[cv]] ];

dExdSE2[cv1_] + dExdSE2[cv2_] ^:=
  dExdSE2[oContents[dBundle[cv1] + dBundle[cv2] ]];

(*====================  Products  ====================*)

dExdSE2[cov_] * TExTSE2[v_] ^:= dBundle[cov] * TBundle[v];
  (*If[ gBase[dExdSE2[cov]] == gBase[TExTSE2[v]] ,
  ,
    "Not Possible", "Not Possible"];*)


(*====================  Differential  ====================*)

(*====================  D  ====================*)

D[f_, ExSE2[q_]] ^:= dExdSE2[oContents[D[f, Bundle[q]]]];


(*Diff[f_ , ExSE2[q_] ] ^:=  Module[{vec, cov},

   vec = oInit[Vector, 
     Join[gLocal[gBase[ExSE2[q]]], gLocal[gGroup[ExSE2[q]]] ] ];
   cov = oContents[Diff[f, vec]];
   dExdSE2[cov]
];*)

(*
D[f_ , ExSE2[q_] ] ^:= Module[{dg, dr, g, r},
  g = gGroup[ExSE2[q]];
  r = gBase[ExSE2[q]];

  dg = D[f, g];
  dr = D[f, r];

  {dr, dg}
];
*)

(*
D[SE2[f_] , ExSE2[q_] ] ^:= Module[{dg, dr, g, r},
  g = gGroup[ExSE2[q]];
  r = gBase[ExSE2[q]];

  dg = D[SE2[f], g];
  dr = Table[ D[gLocal[f][[i]], r], {i, gDim[SE2[f]]} ];

  {dr, dg}
];
*)

(*
D[f_ , ExSE2[q_] ] ^:=  Module[{vec, cov},

   vec = oInit[Vector, 
     Join[gLocal[gBase[ExSE2[q]]], gLocal[gGroup[ExSE2[q]]] ] ];
   cov = oContents[Diff[f, vec]];
   dExdSE2[cov]
];
*)


(*--------------------Graphical Representation--------------------*)

(*
dImage[ExSE2[obj_], grafx_] ^:= ExSE2[ReplacePart[obj, grafx, iGRAFIX]];
gImage[ExSE2[obj_]] ^:= Part[obj, iGRAFIX];
*)

(*--------------------Mathematica Representation--------------------*)

Format[TExTSE2[v_]] ^:=  "Tangent Vector " <>
  ToString[gLocal[gVector[TExTSE2[v]]]] <>
  " of principal fiber bundle at the point " <>
  ToString[{gLocal[gBase[gPoint[TExTSE2[v]]]],
    gLocal[gFiber[gPoint[TExTSE2[v]]]]}];

Format[dExdSE2[cv_]] ^:=  "Cotangent Vector " <>
  ToString[gLocal[gCovector[dExdSE2[cv]]]] <>
  " of principal fiber bundle at the point " <>
  ToString[{gLocal[gBase[gPoint[dExdSE2[cv]]]],
    gLocal[gFiber[gPoint[dExdSE2[cv]]]]}];



End[];

