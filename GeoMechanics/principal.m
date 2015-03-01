(*  configuration.m

  This package provides the necessary function and manipulations required
  to implement the dynamics of a configuration manifold with a principle
  connection.  The configuration manifold in this case is equivalent to
  Q=MxG where the base space M is a manifold, and the position space G is 
  a Lie group.

  Patricio Vela
*)
BeginPackage["PrincipalBundle`"];

Needs["Objects`","ivamatica/Basic/objects.m"];
Needs["Euclidean`","ivamatica/DiffGeometry/euclidean.m"];
Needs["Vectors`","ivamatica/DiffGeometry/vectors.m"];
Needs["Manifolds`","ivamatica/DiffGeometry/manifolds.m"];
Needs["Bundles`","ivamatica/DiffGeometry/bundles.m"];
Needs["LieGroups`","ivamatica/GeoMechanics/liegroups.m"];

(*Needs["LagrangianMechanics`","ivamatica/libs/lagrangian.m"];*)
(*Needs["Equations`","ivamatica/libs/equations.m"];*)

gQ::usage="Returns the principal bundle manifold object type.";
gM::usage="Returns the base space object type.";
gG::usage="Returns the group space object type.";

gLieAlgebra::usage="Returns the Lie algebra object type.";

dGroup::usage="Define group space variables.";
dTotal::usage="Define total space variables.";
gGroup::usage="Get Group space variables.";
gTBase::usage="Get Base space variables.";
gTGroup::usage="Get Group space variables.";

gDimGroup::usage="--";

(*
dBasis::usage="--";
gBasis::usage="--";
*)

(*
dEOM::usage="Define the equations of motion.";
cEOM::usage="Calculate the equations of motion.";
gEOM::usage="Get the equations of motion.";
pEOM::usage="Print out the equations of motion.";
*)

(*cNDSolve::usage="cNDSolve[Obj, q0, {t, t0, t1}, Params] \n Numerically
integrate equations of motion.";*)

(*dImage::usage="dImage[Obj, Grafix]";
gImage::usage="gImage[Obj]";*)


(*CalcGamma::usage="--";*)

(*-----Where does this go??-----*)
(*fLie::usage="fLie[X, Y] \n Take the Lie bracket of the two reduced
vector fields.";*)

Begin["Private`"];

iBASE 		= Bundles`Private`iBASE;	(* base variables. *)
iGROUP		= Bundles`Private`iFIBER;;	(* group variables. *)

End[];

<<mathematica/libs/ExSE2.m


EndPackage[];
