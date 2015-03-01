(*  configuration.m

  This package provides the necessary function and manipulations required
  to implement the dynamics of a configuration manifold with a principle
  connection.  The configuration manifold in this case is equivalent to
  Q=MxG where the base space M is a manifold, and the position space G is 
  a Lie group.

  Author:   Patricio Vela
*)
BeginPackage["ExSE2`"];

ExSE2::usage="Provides functionality of principal bundle Q = R^n x SE(2)";

Begin["`Private`"];

(* ExSE2 is Q=E(m)xSE(2).  *)

lenExSE2	= Bundles`Private`lenBUNDLE;

(*--------------------oInit--------------------*)
(*  Initialize the configuration object.
*)

oInit[ExSE2] ^:= ExSE2[oContents[oInit[Bundle]]];

oInit[ExSE2, Euclidean[r_], SE2[g_]] ^:= 
  ExSE2[oContents[oInit[Bundle, Euclidean[r], SE2[g] ]]];


(*--------------------Manifold Representations--------------------*)

gQ[Obj_] := oType[Obj];
gM[ExSE2[obj_]] ^:= Euclidean`Euclidean;
gG[ExSE2[obj_]] ^:= LieGroups`SE2;
gLieAlgebra[ExSE2[obj_]] ^:= LieGroups`se2;

(*--Base and group variables--*)
dBase[ExSE2[obj_], r_] ^:= ExSE2[oContents[dBase[Bundle[obj], r]]];
dGroup[ExSE2[obj_], g_] ^:= ExSE2[oContents[dFiber[Bundle[obj], g]]];

gBase[ExSE2[obj_]] ^:= gBase[Bundle[obj]];
gFiber[ExSE2[obj_]] ^:= gFiber[Bundle[obj]];
gGroup[ExSE2[obj_]] ^:= gFiber[ExSE2[obj]];

gDim[ExSE2[obj_]] ^:= gDim[Bundle[obj]];
gDimBase[ExSE2[obj_]] ^:= gDimBase[Bundle[obj]];
gDimGroup[ExSE2[obj_]] ^:= gDimFiber[Bundle[obj]];

gLocal[ExSE2[obj_]] ^:= 
  Join[ gLocal[gBase[ExSE2[obj]]] , gLocal[gGroup[ExSE2[obj]]] ];


(*-----------------------------------------------------------*)
(*--------------------Equations of Motion--------------------*)
(*-----------------------------------------------------------*)

(*--------------------gEOM--------------------*)
(*
  Give equations of motion for first order nonlinear matrix ODE of the form
    r' = u(r,t)
    g' = g A(r) u(r,t)
  over the given interval. q = (r, g)
*)
gEOM[se2[A_], Vector[u_], ExSE2[q_] , ExSE2[q0_] , int_] := Module [
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

(*--------------------NDSolve--------------------*)
(*
  Solve first order nonlinear matrix ODE of the form
    r' = u(r,t)
    g' = g A(r) u(r,t)
  over the given interval. q = (r, g)
*)

NDSolve[se2[A_], Vector[u_], ExSE2[q_], ExSE2[q0_], int_] ^:= Module [
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

(*--------------------cNDSolve--------------------*)

(*
Should be part of another objeect.

cNDSolve[Obj_, q0_, time_, params___] := Module[{eom, type, cNDSolveParams},

  params = Flatten[{params}];
  eom = gEOM[Obj];

  type = oType[eom];

  Which[ type == Hamel ,
    "Testing";
  ];
]*)
  
(*---------------------------------------------------------------*)
(*-------------Mathematica Operators and Functions---------------*)
(*---------------------------------------------------------------*)

(*--------------------Graphical Representation--------------------*)

(*dImage[ExSE2[obj_], grafx_] ^:= ExSE2[ReplacePart[obj, grafx, iGRAFIX]];
gImage[ExSE2[obj_]] ^:= Part[obj, iGRAFIX];
*)

(*--------------------Mathematica Representation--------------------*)
Format[ExSE2[object_]] := 
   "Group " Format[Bundle[object]];

End[];

EndPackage[];
