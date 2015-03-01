(*======================== Extra Math Functions ========================*)
(*

  This package implements an object that keeps track of the
  equations of motion for a dynamical system.  These are typically
  tied in to a particular manifold object and do not necessarily
  function alone.



  Author:       Patricio A. Vela
  Created:      2002/06/07
  Modified:     2015/03/01


*)
(*======================== Extra Math Functions ========================*)
BeginPackage["Extras`"];


(*=============== Function help declarations: ==============*)
doSimp::usage="doSimp[simpType], where simpType is one of the following \n 
  stNone (or None), stSimp, stFull, stTrigE, stTrigF, stTrigR. stChop \n";

fLinearize::usage="fLinearize[X, q, eval]";

fSign::usage="fSign[x] \n Returns +1 if non-negative and -1 if negative."

ArcTanQ::usage="ArcTanQ[x,y] \n Returns angle in  +/- \[Pi].";
ArcTanH::usage="ArcTanQ[x,y] \n Returns angle in  +/- \[Pi]/2.";

(*===== Numerical constants associated with the package. =====*)

(* Numerical constants define simplification type sought under doSimp. *)
stNone = 0;     (* No simplification. *)
stSimp = 1;     (* Simplification. *)
stFull = 2;     (* Full simplification. *)
stTrigE = 3;    (* Trigonometric expansion. *)
stTrigF = 4;    (* Trigonometric factoring. *)
stTrigR = 5;    (* Trigonometric reduction. *)
stChop = 6;     (* Numerical simplification by chopping. *)


(*========================== Package Functions =========================*)
Begin["`Private`"];

(*------------------------------ doSimp ------------------------------*)
doSimp[simpType_] := Which[
  simpType == stFull, FullSimplify,
  simpType == stTrigE, TrigExpand,
  simpType == stTrigF, TrigFactor,
  simpType == stTrigR, TrigReduce,
  simpType == stChop, Chop,
  simpType == True || simpType == stSimp, Simplify,
  simpType == False || simpType == stNone, Identity,
  True, Context[stSimp]];

(*---------------------------- fLinearize ----------------------------*)
fLinearize[X_, q_, p_] := Module[{reprule, ldim},
  ldim = Length[q];
  reprule = Table[ q[[i]] -> p[[i]] , {i, ldim}];
  Table[ (D[X, q[[i]]]) /. reprule , {i, ldim}]
];

(*------------------------------- fSign ------------------------------*)
fSign[el_] := If[ el >= 0 , 1 , -1];

(*------------------------------ ArcTanH -----------------------------*)
(*
    Arc tangent for the half plane.
*)
ArcTanH[x_, y_] := If[ x == 0 , Sign[y] Pi/2 , ArcTan[y/x] ]

(*------------------------------ ArcTanH -----------------------------*)
(*
    Arc tangent for the quarter plane, no matter where point lies.
*)
ArcTanQ[x_, y_] := Module[ {angle},
  If[ x == 0 , angle = Sign[y] Pi/2 , angle = ArcTan[y/x] ];
  Which[ x < 0 && angle >= 0 , -Pi + angle, x < 0 && angle < 0, Pi +
  angle , True , angle ]
];

End[];

EndPackage[];
(**)
(*======================== Extra Math Functions ========================*)
