(*--------------------Manifolds--------------------*)
(*

  tmanifolds.m - This package provides the necessary functions and
  manipulations required to implement the Tangent space to a
  Manifold, object functionality in Mathematica.  All basic manifold 
  operations are described here.

  Patricio Vela
*)
BeginPackage["TangentManifolds`"];

(*--Load needed libraries--*)
Needs["Objects`","ivamatica/Basic/objects.m"];
Needs["Euclidean`","ivamatica/DiffGeometry/euclidean.m"];
Needs["Bundles`","ivamatica/DiffGeometry/bundles.m"];
Needs["Vectors`","ivamatica/DiffGeometry/vectors.m"];
Needs["Tangents`","ivamatica/DiffGeometry/tangents.m"];
Needs["Manifolds`","ivamatica/DiffGeometry/manifolds.m"];


(*--Define objects--*)
TManifold::usage="Provides object functionality of a tangent manifold object."; 
dManifold::usage="Provides object functionality of a cotangent manifold object."; 


(*--Define additional functionality--*)





(*--Function definitions follow inthe private context.--*)

Begin["`Private`"];

(*--Private constants for implementation--*)
iBASE = Tangents`Private`iBASE;
iVECTOR = Tangents`Private`iFIBER;



(*-----Tangent Manifold space object functions-----*)

oInit[TManifold, Manifold[q_], Vector[v_]] ^:= 
  TManifold[oContents[oInit[Bundle, Manifold[q], Vector[v]]]];

Format[TManifold[q_]] := Format[TEuclidean[q]];

gDim[TManifold[qp_]] ^:= gDim[TEuclidean[qp]];

dBase[TManifold[qp_], Manifold[q_]] ^:=
  TManifold[oContents[dBase[TEuclidean[qp], Manifold[q]]]];
dVector[TManifold[qp_], Vector[v_]] ^:=
  TManifold[oContents[dVector[TEuclidean[qp], Vector[v]]]];

gBase[TManifold[qp_]] ^:= gBase[TEuclidean[qp]];
gVector[TManifold[qp_]] ^:= gVector[TEuclidean[qp]];

gLocal[TManifold[qp_]] ^:= gLocal[TEuclidean[qp]];

gDual[TManifold] ^:= dManifold;

(*-----Cotanget Manifold space object functions-----*)

oInit[dManifold, Manifold[q_], Covector[cov_]] ^:= 
  dManifold[oContents[oInit[Bundle, Manifold[q], Covector[cov] ]] ];

Format[dManifold[q_]] := Format[dEuclidean[q]];

gDim[dManifold[qp_]] ^:= gDim[dEuclidean[qp]];

dBase[dManifold[dq_], Manifold[q_]] ^:=
  dManifold[oContents[dBase[dEuclidean[dq], Manifold[q]]]];
dCovector[dManifold[qp_], Covector[cov_]] ^:=
  dManifold[oContents[dCovector[dEuclidean[dq], Covector[cov]]]];

gBase[dManifold[dq_]] ^:= gBase[dEuclidean[dq]];
gCovector[dManifold[dq_]] ^:= gCovector[dEuclidean[dq]];

gLocal[dManifold[dq_]] ^:= gLocal[dEuclidean[dq]];



End[];


(*----Overriding globally define Mathematica functions-----*)

(*--Define linear structure--*)
a_ * TManifold[qp_] ^:=
  oInit[TManifold, gBase[TManifold[qp]], a gVector[TManifold[qp]] ];
a_ * dManifold[dq_] ^:=
  oInit[dManifold, gBase[dManifold[dq]], a gCovector[dManifold[dq]]
];

TManifold[qp1_] + TManifold[qp2_] ^:=
  If[ gBase[TManifold[qp1]] == gBase[TManifold[qp2]] ,
    Module[{v1, v2},
    v1 = gVector[TManifold[qp1]];
    v2 = gVector[TManifold[qp2]];
    oInit[TManifold, gBase[TManifold[qp1]] , v1+v2] ]
  ,
    "Not possible." , "Not possible." ];

dManifold[dq1_] + dManifold[dq2_] ^:=
  If[ gBase[dManifold[dq1]] == gBase[dManifold[dq2]] ,
    Module[{v1, v2},
    v1 = gCovector[dManifold[dq1]];
    v2 = gCovector[dManifold[dq2]];
    oInit[dManifold, gBase[dManifold[dq1]] , v1+v2] ]
  ,
    "Not possible." , "Not possible." ];


(*--Define the natural pairing between covector and vector.--*)

dManifold[dq_] * TManifold[qp_] ^:=
  If[ gBase[dManifold[dq]] == gBase[TManifold[qp]] ,
    Module[{cv, v},
    cv = gCovector[dManifold[dq]];
    v  = gVector[TManifold[qp]];
    cv.v ]
  ,
    "Not possible." , "Not possible." ];

dManifold[dq_] . TManifold[qp_] ^:= dManifold[dq] * TManifold[qp];

D[f_ , Manifold[q_] ] ^:=
  oInit[ dManifold, Manifold[q],
    Covector[
      Table[ D[f,gLocal[Manifold[q]][[i]]] , 
        {i, Length[gLocal[Manifold[q]]]}] ] ];


EndPackage[];
