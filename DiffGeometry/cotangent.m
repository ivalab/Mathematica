(*--------------------Manifolds--------------------*)
(*

  euclidean.m - This package provides the necessary functions and
  manipulations required to implement Euclidean space as a manifold object 
  in Mathematica.  

  Author: Patricio Vela
  Date: 6/17/2002
*)
BeginPackage["Euclidean`"];

Needs["Objects`","ivamatica/Basic/objects.m"];

Euclidean::usage="Provides object functionality of Euclidean space."; 
TEuclidean::usage="The Tangent Euclidean space.";
dEuclidean::usage="The dual to the Tangent Euclidean space.";

gDim::usage="gDim[Obj]; \n Returns the dimension of the object";
gLocal::usage="gLocal[Obj] \n Gives the local representation.";

(*-----Where does this go??-----*)
Diff::usage="Take the differential of a function with respect to the configuration space.";


Begin["`Private`"];

iBASE = 1;
iVECTOR = 2;

(*-----Euclidean object functions-----*)

oInit[Euclidean, q_] ^:= Euclidean[q];

gDim[Euclidean[q_]] ^:= Length[q];
gLocal[Euclidean[q_]] ^:= q;

Format[Euclidean[q_]] := q;


(*-----Tangent Euclidean Space object functions-----*)

oInit[TEuclidean, Euclidean[q_], Vector[qp_]] ^:= 
  TEuclidean[{Euclidean[q], Vector[qp]}];

Format[TEuclidean[qp_]] := "Tangent vector " <> ToString[Part[qp,iVECTOR]] <>
  " at the point " <> ToString[Part[qp, iBASE]];

(*-----Cotangent Euclidean Space object functions----*)


(*--------------------Diff--------------------*)

Diff[f_ , Euclidean[q_] ] ^:= 
  oInit[ dEuclidean, Transpose[Table[ D[f,q[[i]]] , {i, Length[q]}]] ];


End[];

EndPackage[];
