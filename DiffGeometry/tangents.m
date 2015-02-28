(*--------------------Tangents--------------------*)
(*

  tangents.m - This package provides the necessary functions and
  manipulations required to implement Euclidean Tangent and Cotangent 
  spaces as fiber bundle objects in Mathematica.  

  Author: Patricio Vela
  Date: 6/17/2002
*)
BeginPackage["Tangents`"];

(*--Load needed libraries--*)
Needs["Objects`","mathematica/libs/objects.m"];
Needs["Euclidean`","mathematica/libs/euclidean.m"];
Needs["Bundles`","mathematica/libs/bundles.m"];
Needs["Vectors`","mathematica/libs/vectors.m"];


(*--Define tangent and cotangent spaces--*)
TEuclidean::usage="The Tangent Euclidean space.";
dEuclidean::usage="The dual to the Tangent Euclidean space.";


(*--Define additional functionality--*)

gPoint::usage="Get base point of the tangent space.";

dVector::usage="dVector[qp, v] \n Define the vector.";
dCovector::usage="dCovector[dq, cov] \n Define the covector.";

gVector::usage="gVector[qp] \n Get the vector.";
gCovector::usage="gCovector[dq] \n Get the covector.";

gDual::usage="gDual[TQ] \n Returns dual space to the tangent space TQ.";

Lie::usage="Lie[X,Y] \n Takes the Jacobi-Lie bracket of X and Y.";


(*--Function definitions follow in Private context.--*)

Begin["`Private`"];

(*--Private constants for implementation--*)
iBASE = Bundles`Private`iBASE;
iVECTOR = Bundles`Private`iFIBER;



(*----------Tangent Euclidean Space object functions-----*)

oInit[TEuclidean, Euclidean[q_], Vector[qp_]] ^:= 
  TEuclidean[oContents[oInit[Bundle, Euclidean[q], Vector[qp]]] ];

Format[TEuclidean[qp_]] := 
  Subscript[Part[qp,iVECTOR], Part[qp, iBASE]];


gDim[TEuclidean[qp_]] ^:= gDim[Bundle[qp]];

dBase[TEuclidean[qp_], q_] ^:= 
  TEuclidean[oContents[dBase[Bundle[qp], q]]];
dVector[TEuclidean[qp_], v_] ^:= 
  TEuclidean[oContents[dFiber[Bundle[qp], v]]];

gBase[TEuclidean[qp_]] ^:= gBase[Bundle[qp]];
gPoint[TEuclidean[qp_]] ^:= gBase[TEuclidean[qp]];
gVector[TEuclidean[qp_]] ^:= gFiber[Bundle[qp]];

gLocal[TEuclidean[qp_]] ^:= 
  Join[ gLocal[gVector[TEuclidean[qp]]] , gLocal[gBase[TEuclidean[qp]]] ];


(*-----Equations of Motion-----*)

gEOM[ TEuclidean[X_] ] :=
  Module[ {q, v},

  q = gBase[TEuclidean[X]];
  v = gVector[TEuclidean[X]];

  gEOM[v,q]
]

gEOM[ TEuclidean[X_], time_ ] :=
  Module[ {q, v},

  q = gBase[TEuclidean[X]];
  v = gVector[TEuclidean[X]];

  gEOM[v,q, time]
]

gEOM[ TEuclidean[X_], Euclidean[q0_], time_ ] :=
  Module[ {q, v},

  q = gBase[TEuclidean[X]];
  v = gVector[TEuclidean[X]];

  gEOM[v, Euclidean[q0], q, time]
];

(*-----Jacobi Lie Bracket/Lie Derivative-----*)

Lie[TEuclidean[X_], TEuclidean[Y_]] ^:=
  Module[ {v1, v2, q1, q2},

  q1 = gBase[TEuclidean[X]];
  q2 = gBase[TEuclidean[Y]];

  If[ q1 == q2 ,
    v1 = gVector[TEuclidean[X]];
    v2 = gVector[TEuclidean[Y]];
    oInit[TEuclidean, q1, oInit[Vector, D[v2, q1] TEuclidean[X] 
      - D[v1, q1] TEuclidean[Y] ]]
  ,
    oInit[Vector, Table[0, {i, Length[gDim[q1]]}]]
  ]
 ];

(*----------Cotangent Euclidean Space object functions----*)

oInit[dEuclidean, Euclidean[q_], Covector[qp_]] ^:= 
  dEuclidean[oContents[oInit[Bundle, Euclidean[q], Covector[qp]]] ];

Format[dEuclidean[qp_]] := 
  Subscript[Part[qp,iVECTOR], Part[qp, iBASE]];
(*  ToString[Part[qp,iVECTOR]] <>
  " at the point " <> ToString[Part[qp, iBASE]]; *)


gDim[TEuclidean[qd_]] ^:= gDim[Bundle[qd]];

dBase[dEuclidean[qd_], q_] ^:= 
  dEuclidean[oContents[dBase[Bundle[qd], q]]];
dCovector[dEuclidean[qd_], cov_] ^:= 
  dEuclidean[oContents[dFiber[Bundle[qd], cov]]];

gBase[dEuclidean[qd_]] ^:= gBase[Bundle[qd]];
gCovector[dEuclidean[qd_]] ^:= gFiber[Bundle[qd]];

gLocal[dEuclidean[dq_]] ^:= 
  Join[ gLocal[gBase[dEuclidean[dq]]] , gLocal[gVector[dEuclidean[dq]]] ];

(*----Overriding globally define Mathematica functions-----*)


(*--------------------D--------------------*)

D[f_ , Euclidean[q_] ] ^:= 
  oInit[ dEuclidean, Euclidean[q], 
    Covector[Table[ D[f,q[[i]]] , {i, Length[q]}]] ];

D[Vector[vq_], Euclidean[q_]] ^:=
  oInit[ dEuclidean, Euclidean[q], 
    Covector[ Transpose[Outer[D, vq, q ]] ] ];

D[TEuclidean[vq_], Euclidean[q_]] ^:=
  If[ gBase[TEuclidean[vq]] == Euclidean[q] ,
    D[ gVector[TEuclidean[vq]], Euclidean[q] ]
  ,
    "Not possible"
  ];
  
(*----Equations of Motion-----*)


gEOM[TEuclidean[vq_]] ^:=  
  gEOM[ gVector[TEuclidean[vq]], gBase[TEuclidean[vq]] ];

gEOM[TEuclidean[vq_], time_] ^:= Module[ {tvq, tvec, tq},
  tvq = TEuclidean[vq];
  tq = gBase[tvq];
  tvec = gVector[tvq] /. 
      Table[gLocal[tq][[i]] -> gLocal[tq][[i]][time], {i, gDim[tq]}];
  tvq = dVector[tvq, tvec];
  gEOM[ gVector[tvq], gBase[tvq], time]
];

NDSolve[TEuclidean[vq_], Euclidean[q0_], time_, params___] ^:= 
  Module[ {tvq, tvec, tq},
    tvq = TEuclidean[vq];
    tq = gBase[tvq];
    tvec = gVector[tvq] /. 
        Table[gLocal[tq][[i]] -> gLocal[tq][[i]][time[[1]]], {i, gDim[tq]}];
    tvq = dVector[tvq, tvec];
    NDSolve[ gVector[tvq], Euclidean[q0], gBase[tvq], time, params]
];

End[];


(*--Define linear structure--*)
a_ * TEuclidean[qp_] ^:= 
  oInit[TEuclidean, gBase[TEuclidean[qp]], a gVector[TEuclidean[qp]] ];
a_ * dEuclidean[dq_] ^:= 
  oInit[dEuclidean, gBase[dEuclidean[dq]], a gCovector[dEuclidean[dq]] ];

TEuclidean[qp1_] + TEuclidean[qp2_] ^:=
  If[ gBase[TEuclidean[qp1]] == gBase[TEuclidean[qp2]] ,
    Module[{v1, v2},
    v1 = gVector[TEuclidean[qp1]];
    v2 = gVector[TEuclidean[qp2]];
    oInit[TEuclidean, gBase[TEuclidean[qp1]] , v1+v2] ]
  ,
    "Not possible." , "Not possible." ];

TEuclidean[qp1_] . TEuclidean[qp2_] ^:=
    oInit[ Vector, 
      D[ TEuclidean[qp1] , gBase[TEuclidean[qp2]] ] . TEuclidean[qp2]];

TEuclidean[qp1_] * TEuclidean[qp2_] ^:= "Not possible";

dEuclidean[dq1_] + dEuclidean[dq2_] ^:=
  If[ gBase[dEuclidean[dq1]] == gBase[dEuclidean[dq2]] ,
    Module[{v1, v2},
    v1 = gCovector[dEuclidean[dq1]];
    v2 = gCovector[dEuclidean[dq2]];
    oInit[dEuclidean, gBase[dEuclidean[dq1]] , v1+v2] ]
  ,
    "Not possible." , "Not possible." ];


(*--Define the natural pairing between covector and vector.--*)

dEuclidean[dq_] * TEuclidean[qp_] ^:=
  If[ gBase[dEuclidean[dq]] == gBase[TEuclidean[qp]] ,
    Module[{cv, v},
    cv = gCovector[dEuclidean[dq]];
    v  = gVector[TEuclidean[qp]];
    cv.v ]
  ,
    "Not possible." , "Not possible." ];

dEuclidean[dq_] . TEuclidean[qp_] ^:= dEuclidean[dq] * TEuclidean[qp];

(*--Define pairings between vectors and points--*)

Euclidean[q_] + TEuclidean[vq_] ^:=
  If[ Euclidean[q] == gBase[TEuclidean[vq]] ,
    Euclidean[q] + gVector[TEuclidean[vq]]
  ,
    "Not possible"
  ];



EndPackage[];
