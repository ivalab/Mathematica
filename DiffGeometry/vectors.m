(*--------------------Manifolds--------------------*)
(*

  vectors.m - This package provides the necessary functions and
  manipulations required to implement a Euclidean vector space 
  in Mathematica.  

  Author: Patricio Vela
  Date: 6/17/2002
*)
BeginPackage["Vectors`"];

Needs["Objects`","mathematica/libs/objects.m"];
(*Needs["Euclidean`","mathematica/libs/euclidean.m"];*)

Vector::usage="Provides object functionality of a vector space."; 
Covector::usage="Provides object functionality of a covector space.";

gLocal::usage="gLocal[Obj] \n Gives the local representation.";

Begin["`Private`"];

oInit[Vector, q_] ^:= Vector[q];
oInit[Covector, q_] ^:= Covector[q];

gDim[Vector[q_]] ^:= Length[q];
gDim[Covector[q_]] ^:= Length[q];

gLocal[Vector[q_]] ^:= q;
gLocal[Covector[q_]] ^:= q;

Format[Vector[q_]] := MatrixForm[q];
Format[Covector[q_]] :=  MatrixForm[{q}];

(*--------------------Diff--------------------*)

D[f_ , Vector[q_] ] ^:=  
  Covector[Diff[f, Euclidean[q]]];

End[];


(*-----Overriding globally defined Mathematica functions-----*)

Join[Vector[v1_], Vector[v2_]] ^:= Vector[Join[v1, v2]];
Join[Covector[v1_], Covector[v2_]] ^:= Covector[Join[v1, v2]];

List[Vector[v1_], Vector[v2_]] ^:= Vector[List[v1, v2]];
List[Covector[v1_], Covector[v2_]] ^:= Covector[List[v1, v2]];


(*--Define linear structure--*)
a_ * Vector[v_] ^:= Vector[a v];
a_ * Covector[v_] ^:= Covector[a v];

Vector[v1_] + Vector[v2_] ^:= 
  If[ gDim[Vector[v1]] == gDim[Vector[v2]] ,
    oInit[Vector, v1+v2] ];

Covector[v1_] + Covector[v2_] ^:= 
  If[ gDim[Covector[v1]] == gDim[Covector[v2]] ,
    oInit[Covector, v1+v2] ];

(*--Define the natural pairing between covector and vector.--*)
Covector[cov_] * Vector[v_] ^:= Module[{cvlen, vlen},
  cvlen = Length[Dimensions[cov]];
  vlen = Length[Dimensions[v]];
  
  Which[ (cvlen > 1) &&  (vlen == 1), Transpose[cov].v ,
    (cvlen == 1) && (vlen > 1),  cov.Transpose[v],
    (cvlen > 1 ) && (vlen > 1) , Transpose[cov].Transpose[v],
    (cvlen == 1) && (vlen == 1) , cov.v
  ]
 ];

Covector[cov_] . Vector[v_] ^:= Covector[cov]*Vector[v]; 

EndPackage[];
