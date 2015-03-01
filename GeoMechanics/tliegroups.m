(*--------------------TangentLieGroups--------------------

  tliegroups.m - This package provides the necessary functions and
  manipulations required to implement Tangent spaces to Lie groups 
  in Mathematica.

  Patricio Vela
*)
BeginPackage["TangentLieGroups`"];

Needs["Objects`","ivamatica/Basic/objects.m"];
Needs["Euclidean`","ivamatica/DiffGeometry/euclidean.m"];
Needs["Bundles`","ivamatica/DiffGeometry/bundles.m"];
Needs["Vectors`","ivamatica/DiffGeometry/vectors.m"];
Needs["Tangents`","ivamatica/DiffGeometry/tangents.m"];
Needs["Manifolds`","ivamatica/DiffGeometry/manifolds.m"];
Needs["TangentManifolds`","ivamatica/DiffGeometry/tmanifolds.m"];
Needs["LieGroups`","ivamatica/GeoMechanics/liegroups.m"];

(*--This package reads in from several files since there are many
  Lie groups that can be defined.  This file contains the basic
  declarations that all Tangent Lie group spaces should have.  
--*)


(*-----Additional functionality declarations-----*)


Begin["`Private`"];

End[];

<<mathematica/libs/TSE.m;

EndPackage[];
