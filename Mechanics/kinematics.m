(*--------------------Kinematics--------------------*)
(*

  Implements kinematic motions of manipulators and rigid
  objects.  Principally this has the objects needed to
  implement forces and torque as wrenches.

  Author: Patricio Vela
  Date: June 6, 2002
*)
BeginPackage["Kinematics`"];

Needs["Objects`","mathematica/libs/objects.m"];
Needs["Euclidean`","mathematica/libs/euclidean.m"];
Needs["Vectors`","mathematica/libs/vectors.m"];
Needs["Tangents`","mathematica/libs/tangents.m"];
Needs["Manifolds`","mathematica/libs/manifolds.m"];
Needs["LieGroups`", "mathematica/libs/groups.m"];

Force::usage="Defines a Force object for distinction. A force is a covector.";

gForce::usage="gTorque[wrench] \n Return the forces associated with the \
  wrench.";
gTorque::usage="gTorque[wrench] \n Return the torques associated with the \
  wrench.";

fStripTorque::usage="fStripTorque[wrench] \n Strips the torque from the \
  wrench and returns a wrench with only forces.";


gTorque[wrench_] := gRotation[wrench];
gForce[wrench_] := gPosition[wrench];

oStripTorque[wrench_] := gVector[wrench];

Begin["`Private`"];

oInit[Force, F_] ^:= Force[F];

gLocal[Force[F]] ^:= gLocal[CoVector[F]];

Format[Force[F]] ^:= F;

End[];

EndPackage[];
