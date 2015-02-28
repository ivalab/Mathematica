BeginPackage["Constraints`"];

Needs["LinearAlgebra`MatrixManipulation`"];

pDimConstraints::usage="fDimConstraints[W] \n Determines the dimension
of the constraint distribution given by W q' = 0.";

pDimMomentum::usage="fDimMomentum[W, dimM] \n Determines the difference \
in dimension between the configuration manifold and the the constraint \
distribution given by W q' = 0.";

pCalcConnection::usage="";

pReduceConnection::usage="";

Begin["`Private`"];

(*--------------------pDimConstraints--------------------*)

pDimConstraints[W_] := Length[NullSpace[W]];

(*--------------------pDimMomentum--------------------*)

pDimMomentum[W_, dimM_] := Length[NullSpace[W]] - dimM;

(*--------------------pCalcConnection--------------------*)
(*  Given the (affine) linear constraint forms, determines
    the local Ehresmann connection form.  The bundle is assumed 
    to be, in a trivialization, Q = M x S.
*)

pCalcConnection[W1_, dimQ_, dimM_] := Module[
  {TW1, E,F, rLen, cLen, tConn, tDet, tAff},

  TW1 = Transpose[W1];
  F = Transpose[Table[ TW1[[i]] , {i, dimM}]]; 
  E = Transpose[Table[ TW1[[i]] , {i, 1+dimM, dimQ}]];
  tDet = Det[E];
  tConn = Inverse[E/tDet].F;
  {tDet, tConn}
]

pCalcConnection[W1_, W2_, dimQ_, dimM_] := Module[
  {TW1, E,F, rLen, cLen, tConn, tDet, tAff},

  TW1 = Transpose[W1];
  F = Transpose[Table[ TW1[[i]] , {i, dimM}]]; 
  E = Transpose[Table[ TW1[[i]] , {i, 1+dimM, dimQ}]];
  tDet = Det[E];
  tConn = Inverse[E/tDet].F;
  tAff = Inverse[E/tDet].W2;
  {tDet, tConn, tAff}
]

(*--------------------pReduceConnection--------------------*)
(*  Given the local Ehresmann connection form for a group invariant 
    Ehresmann connection, returns a principal connection.
*)

pReduceConnection[Econn_, el_ ] := Module[
  {ident, Pconn},

  ident = gGroupIdentity[el];
  Pconn = Econn;
  For[i = 1, i <= group[LieGroups`fDim], i++,
    Pconn = Pconn /. {el[[i]] -> ident[[i]]} ];
  Pconn
];


End[];

EndPackage[];
