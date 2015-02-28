(*--------------------TangentLieGroups--------------------

  tliegroups.m - This package provides the necessary functions and
  manipulations required to implement Tangent spaces to Lie groups 
  in Mathematica.

  Patricio Vela
*)
BeginPackage["TangentLieGroups`"];

Needs["Objects`","mathematica/libs/objects.m"];
Needs["Euclidean`","mathematica/libs/euclidean.m"];
Needs["Bundles`","mathematica/libs/bundles.m"];
Needs["Vectors`","mathematica/libs/vectors.m"];
Needs["Tangents`","mathematica/libs/tangents.m"];
Needs["Manifolds`","mathematica/libs/manifolds.m"];
Needs["TangentManifolds`","mathematica/libs/tmanifolds.m"];
Needs["LieGroups`","mathematica/libs/liegroups.m"];

(*--This package reads in from several files since there are many
  Lie groups that can be defined.  This file contains the basic
  declarations that all Tangent Lie group spaces should have.  
--*)


(*-----Additional functionality declarations-----*)


Begin["`Private`"];

End[];

<<mathematica/libs/TSE.m;

EndPackage[];
