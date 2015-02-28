(*--------------------Euclidean--------------------*)
(*

  euclidean.m - This package provides the necessary functions and
  manipulations required to implement Euclidean space as a manifold object 
  in Mathematica.  

  Author: Patricio Vela
  Date: 6/17/2002
*)
BeginPackage["Euclidean`"];

Needs["Objects`","mathematica/libs/objects.m"];
Needs["Vectors`","mathematica/libs/vectors.m"];

Euclidean::usage="Provides object functionality of Euclidean space."; 


gEOM::usage="Get equations of motion for system.";


Begin["`Private`"];

iBASE = 1;
iVECTOR = 2;

(*-----Euclidean object functions-----*)

oInit[Euclidean, q_] ^:= Euclidean[q];

gDim[Euclidean[q_]] ^:= Length[q];
gLocal[Euclidean[q_]] ^:= oGet[Object[q]];

Euclidean[q1_] == Euclidean[q2_] ^:=
    If[ Length[q1] == Length[q2] ,
      Module[{alist, aval},
        alist = Table[ q1[[i]] == q2[[i]], {i, Length[q1]} ];
	aval = True;
	For[ i=1, i <= Length[alist], i++,
	  aval = aval && alist[[i]] ];
	aval ]
    ,
     False]

(*----Equations of Motion-----*)

gEOM[Vector[v_], Euclidean[q_]] ^:= 
  Table[ q[[i]]' == v[[i]] , {i, gDim[Euclidean[q]]}];

gEOM[Vector[v_], Euclidean[q_], time_] ^:= 
  Table[ q[[i]]'[time] == v[[i]] , {i, gDim[Euclidean[q]]}];

gEOM[Vector[v_], Euclidean[q0_], Euclidean[q_], time_] ^:= 
  Join[Table[ q[[i]]'[time[[1]]] == v[[i]] , {i, gDim[Euclidean[q]]}],
    Table[ q[[i]][time[[2]]] == q0[[i]] , {i, gDim[Euclidean[q]]}]];

DSolve[Vector[v_], Euclidean[q_], time_] ^:=
  Module[{eom},

  eom = gEOM[Vector[v], Euclidean[q], time];
  DSolve[eom, q, time]
];

NDSolve[Vector[v_], Euclidean[q0_], Euclidean[q_], time_, params___] ^:=
  Module[{eom},

  eom = gEOM[Vector[v], Euclidean[q0], Euclidean[q], time];
  NDSolve[ eom, q, time, params]
];

(*-----pairing between points and vectors-----*)

Euclidean[q_] + Vector[v_] ^:=
  oInit[Euclidean, q + v];


(*-----Outputting the Euclidean object.-----*)

Format[Euclidean[q_]] := gLocal[Euclidean[q]];

End[];


EndPackage[];
