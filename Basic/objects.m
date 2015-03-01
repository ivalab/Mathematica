(*=============================== Objects ==============================*)
(*

  This package provides the necessary functions and
  manipulations required to implement an Object Oriented functionality
  in Mathematica.  All basic object functions are described here.



  Author:       Patricio A. Vela
  Created:      2002/XX/XX
  Modified:     2015/03/01

*)
(*=============================== Objects ==============================*)
BeginPackage["Objects`"];

(*=============== Function help declarations: ==============*)
Object::usage="Object is the most basic object class.";

oInit::usage="oInit[ObjType, params] \n Initialize the oject.";
oType::usage="oType[Obj] \n Returns the object type.";

oSet::usage="oSet[Obj] \n Sets the raw contents of the object.";
oGet::usage="oGet[Obj] \n Gets the raw contents of the object.";

oContents::usage="oContents[Obj] \n Returns the raw contents of the object.";

gDim::usage="gDim[obj] \n Returns the length or dimension of the object.";

oCompare::usage="oCompare[Type1, Type2] \n Compares to see if same type";


(*========================== Package Functions =========================*)
Begin["`Private`"];

oInit[Object, val_] ^:= Object[val];

oType[obj_] := obj[[0]];

oContents[obj_] := obj[[1]];

oSet[Object[obj_], val_] ^:= Object[val];
oGet[Object[obj_]] ^:= obj;

oCompare[type1_, type2_] := ToString[type1] == ToString[type2];

End[];

EndPackage[];
(*=============================== Objects ==============================*)
