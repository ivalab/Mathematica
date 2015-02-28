(*  configuration.m

  This package provides the necessary function and manipulations required
  to implement the dynamics of a object on some configuration manifold.  
  The configuration object stores the manifold, the equations of motion,
  and other elements of the system.

  Patricio Vela
*)
BeginPackage["Configuration`"];

Needs["Objects`","mathematica/libs/objects.m"];
Needs["Euclidean`","mathematica/libs/euclidean.m"];
Needs["Manifolds`","mathematica/libs/manifolds.m"];
Needs["LieGroups`","mathematica/libs/liegroups.m"];
Needs["LagrangianMechanics`","mathematica/libs/lagrangian.m"];
Needs["Equations`","mathematica/libs/equations.m"];
Needs["PrincipalBundle`","mathematica/libs/principal.m"];

Configuration::usage="Provides functionality of principal bundle Q = R^n x SE(2)";

dManifold::usage="Define the manifold object.";
gManifold::usage="Returns the manifold object.";

(*dEOM::usage="Define the equations of motion.";
gEOM::usage="Get the equations of motion.";*)

dImage::usage="dImage[Obj, Grafix]";
gImage::usage="gImage[Obj]";


Begin["`Private`"];

iNAME		= 1;    (* Name of the object for Output *)
iQ		= 2;	(* The manifold representation *)
iEOM		= 3;	(* The equations of motion of the system. *)
iGRAFIX		= 4	(* Graphical description of object. *)

lenConfiguration = 4;

(*--------------------oInit--------------------*)
(*--  Initialize the configuration object.
*)

oInit[Configuration] ^:= Configuration[Table[Null, {i, lenConfiguration}]];

oInit[Configuration, name_] ^:= dName[oInit[Configuration], name];

(*--------------------Naming--------------------*)

dName[Configuration[obj_], name_] ^:=
  Configuration[ReplacePart[obj,name, iNAME]];

gName[Configuration[obj_]] ^:= Part[obj, iNAME];

(*--------------------Manifold Representations--------------------*)

dManifold[Configuration[obj_], q_] := 
  Configuration[ReplacePart[obj, q, iQ]];

gManifold[Configuration[obj_]] := Part[obj,iQ];

(*--------------------Equations of Motion--------------------*)

dEOM[Configuration[obj_], eqns_] ^:=
  Configuration[ReplacePart[obj, eqns, iEOM]];

gEOM[Configuration[obj_]] := Part[obj,iEOM];

(*--------------------Graphical Representation--------------------*)

dImage[Configuration[obj_], grafx_] ^:= 
  Configuration[ReplacePart[obj, grafx, iGRAFIX]];

gImage[Configuration[obj_]] ^:= Part[obj, iGRAFIX];

(*--------------------Mathematica Representation--------------------*)

Format[Configuration[object_]] := Part[object, iNAME]; 

End[];

EndPackage[];
