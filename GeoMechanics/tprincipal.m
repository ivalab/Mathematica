(*  configuration.m

  This package provides the necessary function and manipulations required
  to implement the dynamics of a configuration manifold with a principle
  connection.  The configuration manifold in this case is equivalent to
  Q=MxG where the base space M is a manifold, and the position space G is 
  a Lie group.

  Patricio Vela
*)
BeginPackage["TPrincipalBundle`"];

Needs["Objects`","mathematica/libs/objects.m"];
Needs["Euclidean`","mathematica/libs/euclidean.m"];
Needs["Vectors`","mathematica/libs/vectors.m"];
Needs["Tangents`","mathematica/libs/tangents.m"];
Needs["Bundles`","mathematica/libs/bundles.m"];
Needs["Manifolds`","mathematica/libs/manifolds.m"];
Needs["TangentManifolds`","mathematica/libs/tmanifolds.m"];
Needs["LieGroups`","mathematica/libs/liegroups.m"];
(*Needs["LagrangianMechanics`","mathematica/libs/lagrangian.m"];
Needs["Equations`","mathematica/libs/equations.m"];*)
Needs["PrincipalBundle`","mathematica/libs/principal.m"];
Needs["TBundles`","mathematica/libs/tbundles.m"];


gTM::usage="Returns the tangent base space object type.";
gTG::usage="Returns the tangent group space object type.";

(*dBaseVector::usage="Define Base space vector.";
dGroupVector::usage="Define Group space vector.";
gBaseVector::usage="Get Base space vector.";
gGroupVector::usage="Get Group space vector.";
gTBase::usage="Get Base space variables.";
gTGroup::usage="Get Group space variables.";*)

(*dBasis::usage="--";
gBasis::usage="--";
*)

(*dEOM::usage="Define the equations of motion.";
cEOM::usage="Calculate the equations of motion.";
gEOM::usage="Get the equations of motion.";
pEOM::usage="Print out the equations of motion.";

cNDSolve::usage="cNDSolve[Obj, q0, {t, t0, t1}, Params] \n Numerically integrate equations of motion.";

dImage::usage="dImage[Obj, Grafix]";
gImage::usage="gImage[Obj]";


*)
CalcGamma::usage="--";

(*-----Where does this go??-----*)
(*fLie::usage="fLie[X, Y] \n Take the Lie bracket of the two reduced
vector fields.";*)

Begin["`Private`"];

iBASE 		= PrincipalBundle`Private`iBASE;	(* base variables. *)
iGROUP		= PrincipalBundle`Private`iGROUP;	(* group variables. *)

iBASEV		= TBundles`Private`iBASEV;
iGROUPV		= TBundles`Private`iFIBERV;

End[];

<<mathematica/libs/TExTSE2.m

Begin["Private`"];

(*============================================================*)
(*--
	Functions that support the principal bundle object.
    	They are not bound to any instance in particular.  

	Might have to move elsewhere, we'll see.
--*)
(*============================================================*)

CalcGamma[T_, r_] := Module[
  {rl, gamma, dimG, dimM, eb},

  dimG = Length[T];
  dimM = gDim[Obj];
  rl = gLocal[Obj];
  eb[i_] := Transpose[T][[i]];
  gamma = Simplify[Table[ Inverse[T].D[eb[j], r[[alpha]] ], {j,1,dimG},
    {alpha, 1, dimM}]];
  gamma = Table[ gamma[[j, alpha, i]], {i, 1, dimG}, {j,1,dimG},
  {alpha,1,dimM}]
];


End[];


EndPackage[];
