BeginPackage["Extras`"];

doSimp::usage="doSimp[simpType], where simpType is one of the following \n 
  stNone (or None), stSimp, stFull, stTrigE, stTrigF, stTrigR. stChop \n";

fLinearize::usage="fLinearize[X, q, eval]";

fSign::usage="fSign[x] \n Returns +1 if non-negative and -1 if negative."

ArcTanQ::usage="ArcTanQ[x,y] \n Returns angle in  +/- \[Pi].";
ArcTanH::usage="ArcTanQ[x,y] \n Returns angle in  +/- \[Pi]/2.";

stNone = 0;
stSimp = 1;
stFull = 2;
stTrigE = 3;
stTrigF = 4;
stTrigR = 5;
stChop = 6;


Begin["`Private`"];

doSimp[simpType_] := Which[
  simpType == stFull, FullSimplify,
  simpType == stTrigE, TrigExpand,
  simpType == stTrigF, TrigFactor,
  simpType == stTrigR, TrigReduce,
  simpType == stChop, Chop,
  simpType == True || simpType == stSimp, Simplify,
  simpType == False || simpType == stNone, Identity,
  True, Context[stSimp]];

fLinearize[X_, q_, p_] := Module[{reprule, ldim},
  ldim = Length[q];
  reprule = Table[ q[[i]] -> p[[i]] , {i, ldim}];
  Table[ (D[X, q[[i]]]) /. reprule , {i, ldim}]
];

fSign[el_] := If[ el >= 0 , 1 , -1];

ArcTanH[x_, y_] := If[ x == 0 , Sign[y] Pi/2 , ArcTan[y/x] ]

ArcTanQ[x_, y_] := Module[ {angle},
  If[ x == 0 , angle = Sign[y] Pi/2 , angle = ArcTan[y/x] ];
  Which[ x < 0 && angle >= 0 , -Pi + angle, x < 0 && angle < 0, Pi +
  angle , True , angle ]
];

End[];

EndPackage[];
