(*

  groups.m - This package provides the necessary functions and
  manipulations required to implement Lie groups in Mathematica.

  Patricio Vela
*)
BeginPackage["LieGroups`"];

Needs["Objects`","ivamatica/Basic/objects.m"];
Needs["Euclidean`","ivamatica/DiffGeomoetry/euclidean.m"];
Needs["Vectors`","ivamatica/DiffGeomoetry/vectors.m"];
Needs["Manifolds`","ivamatica/DiffGeomoetry/manifolds.m"];
Needs["Tangents`","ivamatica/DiffGeomoetry/tangents.m"];
Unprotect[Inverse];
Unprotect[Left];
Unprotect[Right];

(*--This package reads in from several files since there are many
  Lie groups that can be defined.  This file contains the basic
  declarations that all Lie groups should have.
--*)


(*-----Additional functionality declarations.----*)

(*--Object initialization/implementation--*)

oName::usage"oName[Group, name] \n Define variable name for the group";

(*--Group data/information--*)
gGroup::usage="gGroup[Obj] \n Returns the Lie group associated to the object.";
gTangent::usage"gTangent[Obj] \n Returns the Tangent bundle of the object.";
gReduce::usage"gReduce[Obj] \n Gets the reduced structure of the object.";
gGroupIdentity::usage="Returns the group identity.";
gLieAlgebra::usage="Returns the Lie Algebra.";
gDualLieAlgebra::usage="Returns the dual to the Lie Algebra.";

(*--Groups actions--*)
(*Left::usage="Left[g, X] \n The left action of g on X.";*)
(*Right::usage="Right[g, X] \n The right action of g on X.";*)

Action::usage="Action[g, X] \n Used to select which action to use. May
be set to Left, Right, Ad, etc. (default: Left)";

Action = Left;

(*--Lie algebra actions--*)
InfGen::usage="--";
Ad::usage="Ad[element, g] \n Takes the Adjoint at g of the element. \ 
  May be in the group, lie algebra, or dual lie algebra. For now, \
  though, only valid for dual.";
ad::usage="ad[xi1, xi2] \n Computes the adhjoint/Lie bracket [xi1,xi2]";

(*--Reduction--*)
gReduce::usage="oReduce[func, g] \n
Reduces the G-invariant function of G by the group G.\n Reduces the G-invariant function of TG by the group G.  Returns a function of the Lie algebra of G";
cReduce::usage="cReduce[func, ... ] \n Reduces function.  Different for
each Lie group. This is an overloading operator.";

StructureConstants::usage="StructureConstants[LA, T] \n Obtains the
structure constants of the Lie Algebra LA with respect to the
transformation of basis T.";

Begin["`Private`"];

gGroup[obj_] := oType[obj];
gReduce[obj_] := gGroupIdentity[obj];

gReduce[func_, g_] := Module[{e, h},
  e = gLocal[gReduce[g]];
  h = gLocal[g];
  func /. Table[ h[[i]] -> e[[i]], {i, gDim[g]} ]
];


(*-------------------StructureConstants--------------------*)

StructureConstants[LA_, T_] := Module[
  {eb, SC, dimG},

  dimG = gDim[LA];
  eb[i_] := Transpose[T][[i]];
  SC = Table[ Inverse[T].gLocal[ oad[LA[eb[j]], LA[eb[k]]] ],
      {j,1,dimG}, {k,1,dimG}];
  SC = Table[ -SC[[j,k,i]], {i, 1, dimG},{j,1,dimG},{k,1,dimG}]
];


End[];

<<ivamatica/GeoMechanics/SE.m;
<<ivamatica/GeoMechanics/SO.m;


Protect[Inverse];
Protect[Left];
Protect[Right];
EndPackage[];
