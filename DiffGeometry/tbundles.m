(*--------------------Tangent Fiber Bundles--------------------*)
(*

  tbundles.m - This package provides the necessary functions and
  manipulations required to implement a Fiber Bundle object functionality
  in Mathematica.  All basic fiber bundle operations are described here.

  Patricio Vela
*)
BeginPackage["TBundles`"];

(*--Load libraries--*)
Needs["Objects`","ivamatica/Basic/objects.m"];
Needs["Euclidean`","ivamatica/DiffGeometry/euclidean.m"];
Needs["Vectors`","ivamatica/DiffGeometry/vectors.m"];
Needs["Tangents`","ivamatica/DiffGeometry/tangents.m"];
Needs["Bundles`","ivamatica/DiffGeometry/bundles.m"];
Needs["Manifolds`","ivamatica/DiffGeometry/manifolds.m"];
Needs["TangentManifolds`","ivamatica/DiffGeometry/tmanifolds.m"];


(*--Define objects--*)
TBundle::usage="Provides object functionality of the tangent space to a fiber bundle object."; 
dBundle::usage="Provides object functionality of the cotangent space to a fiber bundle object."; 


(*--Declare functions--*)
gBaseVector::usage="gBaseVector[vector]";
gFiberVector::usage="gFiberVector[vector]";

(*--Function definitions follow in Private context.--*)

Begin["`Private`"];

(*--Private constants for implementation.--*)
iBASE = Bundles`Private`iBASE;
iFIBER = Bundles`Private`iFIBER;

lenBUNDLE = Bundles`Private`lenBUNDLE;

iBASEV = 3;
iFIBERV = 4;

lenTBUNDLE = 4;

(*----------Tangent Space to Fiber Bundle functionality---------*)

oInit[TBundle] ^:= TBundle[Table[{}, {i, lenTBUNDLE}]];

oInit[TBundle, Bundle[q_]] ^:=
  dBase[oInit[TBundle], Bundle[q]];

oInit[TBundle, r_, y_] ^:=
  oInit[TBundle, oInit[Bundle, r, y]];

oInit[TBundle, Bundle[q_], Vector[v_] ] ^:= 
  dVector[oInit[TBundle, Bundle[q]] , Vector[v] ];

oInit[TBundle, r_, y_, Vector[rv_], Vector[yv_]] ^:=
  oInit[TBundle, oInit[Bundle, r, y], {Vector[rv], Vector[yv]}];

(*--------------------Define and Get components--------------------*)

dBase[TBundle[v_] , Bundle[q_]] ^:=
  TBundle[oContents[dFiber[dBase[Bundle[v],gBase[Bundle[q]] ],
    gFiber[Bundle[q]] ] ]];

dVector[TBundle[v_], Vector[bv_], Vector[fv_]] ^:=  
  TBundle[ReplacePart[ReplacePart[v,Vector[bv],iBASEV], 
    Vector[fv], iFIBERV] ];

dVector[TBundle[v_], Vector[vec_]] ^:= 
  TBundle[ReplacePart[ReplacePart[v, oInit[Vector,vec[[1]]], iBASEV],
    oInit[Vector,vec[[2]]], iFIBERV] ];

gBase[TBundle[bv_]] ^:= Bundle[ Take[bv, lenBUNDLE] ];
gPoint[TBundle[bv_]] ^:= gBase[TBundle[bv]];
gVector[TBundle[bv_]] ^:= Take[bv, {1+lenBUNDLE, lenTBUNDLE}];
gBaseVector[TBundle[bv_]] ^:= Part[bv, iBASEV];
gFiberVector[TBundle[bv_]] ^:= Part[bv, iFIBERV];

gDim[TBundle[bv_]] ^:= gDim[gBase[TBundle[bv]]] + 
  gDim[Part[bv, iBASEV]] + gDim[Part[bv, iFIBERV]];


gLocal[TBundle[bv_]] ^:= 
  Join[gLocal[gBase[TBundle[bv]]] , gLocal[Part[bv, iBASEV]] ,
    gLocal[Part[bv, iFIBERV]] ];

Format[TBundle[bv_]] := "Tangent Vector " <> 
  ToString[gLocal[gVector[TBundle[bv]]]] <> 
  " of fiber bundle at the point " <> 
  ToString[{gLocal[gBase[gPoint[TBundle[bv]]]],
    gLocal[gFiber[gPoint[TBundle[bv]]]]}];
 

(*----------Cotangent Space to Fiber Bundle functionality---------*)

oInit[dBundle] ^:= dBundle[Table[{}, {i, lenTBUNDLE}]];

oInit[dBundle, Bundle[q_]] ^:=
  dBase[oInit[dBundle], Bundle[q]];

oInit[dBundle, r_, y_] ^:= 
  oInit[dBundle, oInit[Bundle, r, y]];

oInit[dBundle, Bundle[q_], Covector[cv_] ] ^:= 
  dCovector[oInit[dBundle,Bundle[q]], Covector[cv]];

oInit[dBundle, r_, y_, Covector[dr_], Covector[dy_]] ^:=
  oInit[dBundle, oInit[Bundle[r, y]] , {Covector[dr], Covector[dy]} ];


(*--------------------Define and Get components--------------------*)

dBase[dBundle[cv_] , Bundle[q_]] ^:=
  dBundle[oContents[dFiber[dBase[Bundle[cv],gBase[Bundle[q]] ],
    gFiber[Bundle[q]] ] ]];

dCovector[dBundle[cv_], Covector[bcv_], Covector[fcv_]] ^:= 
  dBundle[ReplacePart[ReplacePart[cv,Covector[bcv],iBASEV], 
    Covector[fcv], iFIBERV] ];

dCovector[dBundle[cv_], Covector[covec_]] ^:=
  dBundle[ReplacePart[ReplacePart[cv, oInit[Covector,covec[[1]]], iBASEV],
    oInit[Covector,covec[[2]]], iFIBERV] ];

gBase[dBundle[bcv_]] ^:= Bundle[ Take[bcv, lenBUNDLE] ];
gPoint[dBundle[bcv_]] ^:= gBase[dBundle[bcv]];
gCovector[dBundle[bcv_]] ^:= Take[bcv, {1+lenBUNDLE, lenTBUNDLE}];

gDim[dBundle[bcv_]] ^:= gDim[gBase[dBundle[bcv]]] + 
  gDim[Part[bcv, iBASEV]] + gDim[Part[bcv, iFIBERV]];

gLocal[dBundle[bcv_]] ^:= 
  Join[gLocal[gBase[dBundle[bcv]]] , gLocal[Part[bcv, iBASEV]] ,
    gLocal[Part[bcv, iFIBERV]] ];

Format[dBundle[bcv_]] := "Cotangent Vector " <>
  ToString[gLocal[gCovector[dBundle[bcv]]]] <> 
  " of fiber bundle at the point " <> 
  ToString[{gLocal[gBase[gPoint[dBundle[bcv]]]],
    gLocal[gFiber[gPoint[dBundle[bcv]]]]}]; 

(*---------------------------------------------------------------------*)
(*-----------------Mathematica Functions and Operators-----------------*)
(*---------------------------------------------------------------------*)


(*====================  Arithmetic  ====================*)

(*----------TBundle----------*)

alpha_ * TBundle[v_] ^:=
  oInit[TBundle, gBase[TBundle[v]], alpha * gVector[TBundle[v]] ];

TBundle[bv1_] + TBundle[bv2_] ^:=
  If[ gBase[TBundle[bv1]] == gBase[TBundle[bv2]] ,
    Module[{v1, v2},

    v1 = gVector[TBundle[bv1]];
    v2 = gVector[TBundle[bv2]];

    oInit[TBundle, gBase[TBundle[bv1]] , v1 + v2 ]
    ]
  ,
    "Not Possible." , "Not possible." ];

(*----------dBundle----------*)

alpha_ * dBundle[v_] ^:=
  oInit[dBundle, gBase[dBundle[v]], alpha * gCovector[dBundle[v]] ];

dBundle[bcv1_] + dBundle[bcv2_] ^:=
  If[ gBase[dBundle[bcv1]] == gBase[dBundle[bcv2]] ,
    Module[{cv1, cv2},

    cv1 = gCovector[dBundle[bcv1]];
    cv2 = gCovector[dBundle[bcv2]];

    oInit[dBundle, gBase[dBundle[bcv1]] , cv1 + cv2 ]
    ]
  ,
    "Not Possible." , "Not possible." ];

(*====================  Products  ====================*)

dBundle[dq_] * TBundle[qv_] ^:=
  If[ gBase[dBundle[dq]] == gBase[TBundle[qv]] ,
    Module[{cv, v},

    cv = gLocal[gCovector[dBundle[dq]]];
    v = gLocal[gVector[TBundle[qv]]];

    Sum[ oInit[Covector, cv[[i]]]*oInit[Vector, v[[i]]] , {i, Length[cv]} ]
    ]
  ,
    "Not Possible." , "Not possible." ];



(*====================  Differential  ====================*)

D[f_, Bundle[q_]] ^:= oInit[dBundle, Bundle[q] , 
  {gCovector[D[f, gBase[Bundle[q]]]], gCovector[D[f, gFiber[Bundle[q]]]]}];

End[];


EndPackage[];
