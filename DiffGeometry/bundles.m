(*--------------------Manifolds--------------------*)
(*

  bundles.m - This package provides the necessary functions and
  manipulations required to implement a Fiber Bundle object functionality
  in Mathematica.  All basic fiber bundle operations are described here.

  Patricio Vela
*)
BeginPackage["Bundles`"];

(*--Load libraries--*)
Needs["Objects`","mathematica/libs/objects.m"];
Needs["Vectors`","mathematica/libs/vectors.m"];
Needs["Euclidean`","mathematica/libs/euclidean.m"];


(*--Define objects--*)
Bundle::usage="Provides object functionality of a fiber bundle object."; 



(*--Define functionality--*)
dBase::usage="Define base space variables.";
dFiber::usage="Define the fiber space variables.";

gBase::usage="Get base space variables.";
gFiber::usage="Get fiber space variables.";

gDimBase::usage="Get dimension of the base space.";
gDimFiber::usage="Get dimension of the fiber space.";



(*--Function definitions follow in Private context.--*)

Begin["`Private`"];

(*--Private constants for implementation.--*)
iBASE = 1;
iFIBER = 2;

lenBUNDLE = 2;

(*----------Bundle functionality---------*)

oInit[Bundle] ^:= Table[ {}, {i, lenBUNDLE}];
oInit[Bundle, r_, s_ ] ^:= Bundle[{r, s}];

dBase[Bundle[q_], r_] ^:= Bundle[ReplacePart[q,r,iBASE]];
dFiber[Bundle[q_], s_] ^:= Bundle[ReplacePart[q,s,iFIBER]];

gBase[Bundle[q_]] ^:= Part[q,iBASE];
gFiber[Bundle[q_]] ^:= Part[q,iFIBER];

gDim[Bundle[q_]] ^:= gDim[gBase[Bundle[q]]] + gDim[gFiber[Bundle[q]]] ;
gDimBase[Bundle[q_]] ^:= gDim[gBase[Bundle[q]]];
gDimFiber[Bundle[q_]] ^:= gDim[gFiber[Bundle[q]]];

gLocal[Bundle[q_]] ^:= 
  Join[gLocal[gBase[Bundle[q]]],gLocal[gFiber[Bundle[q]]]];

Format[Bundle[q_]] := "Fiber " <> ToString[gLocal[gFiber[Bundle[q]]]] <>
  " with base " <> ToString[gLocal[gBase[Bundle[q]]]];

(*---------------------------------------------------------------------*)
(*-----------------Mathematica Functions and Operators-----------------*)
(*---------------------------------------------------------------------*)

Bundle[q1_] == Bundle[q2_] ^:= 
  And[gBase[Bundle[q1]] == gBase[Bundle[q2]], 
      gFiber[Bundle[q1]] == gFiber[Bundle[q2]] ];

End[];

EndPackage[];
