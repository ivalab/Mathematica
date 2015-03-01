(*--------------------Manifolds--------------------*)
(*

  manifolds.m - This package provides the necessary functions and
  manipulations required to implement a Manifold object functionality
  in Mathematica.  All basic manifold operations are described here.

  Patricio Vela
*)
BeginPackage["Manifolds`"];

(*--Load needed libraries--*)
Needs["Objects`","ivamatica/Basic/objects.m"];
Needs["Euclidean`","ivamatica/DiffGeometry/euclidean.m"];


(*--Define objects--*)
Manifold::usage="Provides object functionality of a manifold object."; 


(*--Define additional object functionality.--*)





(*--Function definitions follow in the Private context.--*)

Begin["`Private`"];

(*--Private constants for implementation.--*)
iDIM = 1;
iLOCAL = 2;


(*-----Manifold object functions-----*)

oInit[Manifold, dimQ_ , q_] ^:= Manifold[{dimQ, q}];

gDim[Manifold[q_]] ^:= Part[q, iDIM];
gLocal[Manifold[q_]] ^:= Part[q, iLOCAL];

Format[Manifold[q_]] := gLocal[Manifold[q]];




End[];

EndPackage[];
