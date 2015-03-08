(*--------------------Objects--------------------*)
(*

  binary.m - This package provides the necessary functions and
  manipulations required to implement an Binary objectfunctionality
  in Mathematica.  All basic binary functions are described here.
  The role, for now, is to distinguish the different packages
  that exist in my library and to load the appropriate ones
  when needed with minimal effort.

  Patricio Vela
*)
BeginPackage["Binary`"];

Needs["Objects`","ivamativa/Basic/objects.m"];

Binary::usage="Binary is used for loading specific classes.";
oTest::usage="oTest[b] \n Test to see if b has any True values.";

Begin["`Private`"];

iLEN = 1;
iVAL = 2;

(* Will give error saying that Integer tag is protected. *)
(* Don't know why. This is new. Ignoring for now. *)
oInit[Binary, len_Integer] ^:= Binary[{len, 0}];
oInit[Binary, len_Integer, val_Integer] ^:= Binary[{len, val}];

gDim[Binary[b_]] ^:= Part[b, iLEN];

oSet[Binary[b_] , val_Integer ] ^:= Binary[ReplacePart[b,val,iVAL]];
oGet[Binary[b_]] ^:= Part[b, iVAL];

oTest[Binary[b_]] ^:= oGet[Binary[b]] > 0;

Binary[b1_] && Binary[b2_] ^:= Module[{tb},
  tb = Binary[b1] Binary[b2];
  oTest[tb]
 ];

Equal[Binary[b1_] , Binary[b2_]] ^:= 
  Equal[ oGet[Binary[b1]], oGet[Binary[b2]] ];

Plus[Binary[b1_] , Binary[b2_]] ^:= Module[
  {tb},
  
  tb = oInit[Binary, Max[gDim[Binary[b1]] , gDim[Binary[b2]] ] ];
  tb = oSet[tb, BitOr[oGet[Binary[b1]], oGet[Binary[b2]]]];
  tb
 ];

Times[Binary[b1_] , Binary[b2_]] ^:= Module[
  {tb},
  
  tb = oInit[Binary, Max[gDim[Binary[b1]] , gDim[Binary[b2]] ] ];
  tb = oSet[tb, BitAnd[oGet[Binary[b1]], oGet[Binary[b2]]]];
  tb
 ];

Format[Binary[b_]] := IntegerDigits[oGet[Binary[b]], 2,gDim[Binary[b]] ];

End[];

EndPackage[];
