(*
  SE2.m - This package provides the necessary functions and
  manipulations required to implement the Lie groups SE(2) in 
  Mathematica.

  Patricio Vela
*)
(* BeginPackage["SE2`"];*)

(* Needs["Euclidean`","ivamatica/DiffGeometry/euclidean.m"];*)

SE2::usage="Provides the object functionality of the SE2 group.";
se2::usage="Provides the object functionality of the se2 Lie algebra.";
dse2::usage="Provides the object functionality of the dual to the Lie algebra se2.";

Begin["`Private`"];

(*--------------------SE(2)--------------------*)
(* SE(2) is given by g = (x,y,theta).  
   It's Lie algebra se(2) is xi = (vx, vy, vtheta)
   It's dual se(2)^* is alpha = (dx, dy, dtheta)
*)

(*--Constant definitions of the Lie group SE(2).--*)
gDim[SE2] ^:= 3;
gDim[SE2[g_]] ^:= 3;
gDim[se2] ^:= gDim[SE2];
gDim[se2[xi_]] ^:= gDim[SE2];
gDim[dse2] ^:= gDim[SE2];
gDim[dse2[a_]] ^:= gDim[SE2];

gGroup[SE2[obj_]] ^:= SE2;

gReduce[SE2] ^:= gGroupIdentity[SE2];

gGroupIdentity[SE2] ^:= SE2[{0, 0, 0}];
gGroupIdentity[SE2[g_]] ^:= gGroupIdentity[SE2];
gLieAlgebra[SE2] ^:= se2;
gLieAlgebra[SE2[g_]] ^:= se2;
gDualLieAlgebra[SE2] ^:= dse2;
gDualLieAlgebra[SE2[g_]] ^:= dse2;
gDual[se2] ^:= dse2;
gDual[se2[xi_]] ^:= dse2;

oInit[SE2, g1_ , g2_ , g3_] ^:= SE2[{g1,g2,g3}];
oInit[SE2, g_] ^:= SE2[g];
oInit[se2, xi1_ , xi2_ , xi3_] ^:= se2[{xi1,xi2,xi3}];
oInit[se2, xi_] ^:= se2[xi];
oInit[dse2, a1_ , a2_ , a3_] ^:= dse2[{a1,a2,a3}];
oInit[dse2, a_] ^:= dse2[a];

(*SE2[fInit, g1_ , g2_ , g3_] := {g1,g2,g3};*)

(*--Action of/on SE(2)--*)
Translate[SE2[g_], p_] ^:=
  SE2[g + Join[p, {0}]];

Left[SE2[g_]] ^:=
  { { Cos[g[[3]]] , -Sin[g[[3]]], g[[1]] },
    { Sin[g[[3]]] ,  Cos[g[[3]]], g[[2]] },
    { 0 , 0, 1 } };

gT[SE2[g_]] ^:= 
  { { Cos[g[[3]]] , -Sin[g[[3]]], 0 },
    { Sin[g[[3]]] ,  Cos[g[[3]]], 0 },
    { 0 , 0, 1 } };

Ad[SE2[g_]] ^:= 
  { { Cos[g[[3]]] , - Sin[g[[3]]] ,   g[[2]] } , 
    { Sin[g[[3]]] ,   Cos[g[[3]]] , - g[[1]] } ,
    { 0 , 0, 1} } ;

InfGen[se2[xi_], SE2[g_]] ^:=
  se2[{xi[[1]] - xi[[3]] g[[2]] , xi[[2]] + xi[[3]] g[[1]] , xi[[3]]}];

Left[SE2[g_], se2[xi_]] ^:= oInit[Vector, gT[SE2[g]].xi];

(*--------------------Inverse--------------------*)

Inverse[SE2[g_]] := Module[ {A} ,
  A = - Inverse[gT[SE2[g]]];
  SE2[A.g]
];


(*--------------------Reduce--------------------*)

cReduce[func_, SE2[g_]] ^:= Module[{e},
  e = gLocal[gGroupIdentity[SE2]];
  func /. Table[ g[[i]] -> e[[i]], {i, gDim[SE2]} ]
];


(*--------------------Other functions--------------------*)
gVector[SE2[g_]] ^:= { g[[1]], g[[2]], 0};
gVector[se2[xi_]] ^:= se2[gVector[SE2[xi]]];
gVector[dse2[alpha_]] ^:= dse2[gVector[SE2[alpha]]];

gRotation[SE2[g_]] ^:= g[[3]];
gRotation[se2[xi_]] ^:= gRotation[SE2[xi]];
gRotation[dse2[alpha_]] ^:= gRotation[SE2[alpha]];

gPosition[SE2[g_]] ^:= Table[g[[i]], {i, 2}];
gPosition[se2[xi_]] ^:= gPosition[SE2[xi]];
gPosition[dse2[alpha_]] ^:= gPosition[SE2[alpha]];

gLocal[SE2[g_]] ^:= g;
gLocal[se2[xi_]] ^:= xi;
gLocal[dse2[alpha_]] ^:= alpha;

gMag[SE2[g_]] ^:=  Sqrt[(g[[1]])^2 + (g[[2]])^2];
gAngle[SE2[g_]] ^:= If[ (g[[1]] == 0) && (g[[2]] == 0) , 0, 
  ArcTan[ g[[1]], g[[2]] ] ];
gBranch[SE2[g_], bc_] ^:= Module[{amod},
   amod = Mod[g[[3]], 2 Pi];
   Which [ amod>Max[bc] , amod-= 2 Pi , amod<Min[bc], amod+= 2 Pi];
   amod
  ];

(*--------------------se(2)--------------------*)
Wedge[ se2[xi_], se2[eta_] ] ^:= 
  se2[{  xi[[2]] eta[[3]] - xi[[3]] eta[[2]] ,
    - xi[[1]] eta[[3]] + xi[[3]] eta[[1]] ,
     0  }];
Ad[SE2[g_], se2[xi_]] := 
  se2[{  xi[[3]] g[[2]] + xi[[1]] Cos[g[[3]]] - xi[[2]] Sin[g[[3]]] ,
    -xi[[3]] g[[1]] + xi[[1]] Sin[g[[3]]] + xi[[2]] Cos[g[[3]]] ,
     xi[[3]] }];
ad[ se2[xi_], se2[eta_] ] ^:= Wedge[ se2[xi], se2[eta] ];

(*--------------------se(2)^*--------------------*)
Cross[dse2[alpha_], se2[xi_]] ^:=
  dse2[Cross[{alpha[[1]], alpha[[2]], 0} , {xi[[1]], xi[[2]], xi[[3]]}] ];

Ad[dse2[alpha_] , SE2[g_]] ^:= 
  dse2[ { alpha[[1]] Cos[g[[3]]] + alpha[[2]] Sin[g[[3]]] ,
    - alpha[[1]] Sin[g[[3]]] + alpha[[2]] Cos[g[[3]]] ,
    - alpha[[1]] g[[2]] + alpha[[2]] g[[1]] + alpha[[3]] } ];


(*==================================================================*)
(*------------------------Equations of Motion-----------------------*)
(*==================================================================*)

(*--------------------gEOM--------------------*)
(*
  Compute first order nonlinear matrix ODE of the form
    g' = F
  for use with NDSolve.
*)
gEOM[Vector[F_], SE2[g_]] ^:= Module [
  {gp} , 
  gp = Table[ gLocal[SE2[g]][[i]]' , {i, gDim[SE2]} ];
  Table[gp[[i]] == gLocal[Vector[F]][[i]] , {i,gDim[SE2]}]
]
gEOM[Vector[F_], SE2[g_], t_] ^:= Module [
  {gp, v} , 
  gp = Table[ gLocal[SE2[g]][[i]]'[t] , {i, gDim[SE2]} ];
  Table[gp[[i]] == gLocal[F][[i]] , {i,gDim[SE2]}]
]
gEOM[Vector[F_], SE2[g0_], SE2[g_], time_] ^:= Module [
  {gp, v} , 
  gp = Table[ gLocal[SE2[g]][[i]]'[time[[1]]] , {i, gDim[SE2]} ];
  Join[Table[gp[[i]] == gLocal[F][[i]] , {i,gDim[SE2]}], 
	        Table[gLocal[SE2[g]][[i]][time[[2]]] == g0[[i]] , {i, gDim[SE2]}] ]

]

(*--------------------gEOM--------------------*)
(*
  Compute first order nonlinear matrix ODE of the form
    g' = g F
  for use with NDSolve.
*)
gEOM[se2[F_], SE2[g_]] ^:= Module [
  {hp} , 
  hp = Table[ gLocal[SE2[g]][[i]]' , {i, gDim[SE2]} ];
  Table[hp[[i]] == (gT[SE2[g]].F)[[i]] , {i,gDim[SE2]}]
]

gEOM[se2[F_], SE2[g_], t_] ^:= Module [
  {h,hp} , 
  h = oInit[SE2, Table[ gLocal[SE2[g]][[i]][t], {i, gDim[SE2]}] ];
  hp = Table[ gLocal[SE2[g]][[i]]'[t] , {i, gDim[SE2]} ];
  Table[hp[[i]] == (gT[h].F)[[i]] , {i,gDim[SE2]}]
]

gEOM[se2[F_], SE2[g0_], SE2[g_], time_] ^:= Module [
  {h, hp, hh} , 
  h  = oInit[SE2,g[[1]][time[[1]]],g[[2]][time[[1]]],g[[3]][time[[1]]]];
  hp = {g[[1]]'[time[[1]]],g[[2]]'[time[[1]]],g[[3]]'[time[[1]]]};
  Join[ Table[hp[[i]] == (gT[h].F)[[i]] , {i,gDim[SE2]}], 
	        Table[g[[i]][time[[2]]] == g0[[i]] , {i, gDim[SE2]}] ]
]

(*
  Compute the first order nonlinear matrix ODE of the form
    xi' =  X(xi, t)
  for use with NDSolve
*)
gEOM[Vector[X_], se2[xi0_], se2[xi_], time_] ^:= Module [
  {xip, Xt} , 

  xip = {xi[[1]]'[time[[1]]],xi[[2]]'[time[[1]]],xi[[3]]'[time[[1]]]};
  Xt = X /. Table[ xi[[i]] -> xi[[i]][time[[1]]], {i, gDim[SE2]} ];

  Join[ Table[xip[[i]] == Xt[[i]] , {i,gDim[SE2]}], 
	Table[xi[[i]][time[[2]]] == xi0[[i]] , {i,gDim[SE2]}] ]
]


(*-------------------------NDSolve-------------------------*)

NDSolve[Vector[F_], SE2[g0_], SE2[g_], time_, params___] ^:=
  Module[{eom},

  eom = gEOM[Vector[F], SE2[g0], SE2[g], time];
  NDSolve[ eom, g, time, params]
];

NDSolve[se2[xi_], SE2[g0_], SE2[g_], time_, params___] ^:=
  Module[{eom},

  eom = gEOM[se2[xi], SE2[g0], SE2[g], time];
  NDSolve[ eom, g, time, params]
];

NDSolve[Vector[X_], se2[xi0_], se2[xi_], time_, params___] ^:=
  Module[ {eom},
  
  eom = gEOM[Vector[X], se2[xi0], se2[xi], time];
  NDSolve[eom, xi, time, params]
];


(*==================================================================*)
(*--------------------Exponential representation--------------------*)
(*==================================================================*)

Exp[se2[xi_]] ^:= Module[{w, sgn, v, m},
  w = Mod[xi[[3]], If[xi[[3]]> 0, 1, -1] 2 Pi];
  sgn = Sign[xi[[3]]];

  If[ xi[[3]] == 0 , oInit[SE2, {xi[[1]], xi[[2]], 0}],
    v = (1/xi[[3]])*
      {{1-Cos[w],Sin[w]}, {-Sin[w], 1-Cos[w]}}.{-xi[[2]], xi[[1]]};
    oInit[SE2,{v[[1]],v[[2]],xi[[3]]}] ]
];

Log[SE2[g_]] ^:= Module[{w,v,m},
  w = Mod[g[[3]], If[g[[3]]> 0, 1, -1] 2 Pi];
  
  If[ w == 0, oInit[se2, {g[[1]], g[[2]], 0}],
    v = g[[3]]*
      Inverse[{{Sin[w],-1+Cos[w]},{1-Cos[w],Sin[w]}}].{g[[1]],g[[2]]};
    oInit[se2,{v[[1]], v[[2]], g[[3]]}] ]
];


(*--Functions that support the configuration object.
    They are not bound to any instance in particular.  --*)

StructureConstants[LA_, T_] := Module[
  {eb, SC, dimG},

  dimG = gDim[LA];
  eb[i_] := Transpose[T][[i]];
  SC = Table[ Inverse[T].gLocal[ Wedge[LA[eb[j]], LA[eb[k]]] ],
      {j,1,dimG}, {k,1,dimG}];
  SC = Table[ -SC[[j,k,i]], {i, 1, dimG},{j,1,dimG},{k,1,dimG}]
];


Format[SE2[g_]] := MatrixForm[g];
Format[se2[xi_]] := xi;
Format[dse2[alpha_]] := alpha;

End[];

SE2[g_] == SE2[h_] ^:= Euclidean[g] == Euclidean[h];

SE2[g_] \[CenterDot] SE2[h_] ^:= 
  SE2[{  h[[1]] Cos[g[[3]]] - h[[2]] Sin[g[[3]]] + g[[1]],
     h[[1]] Sin[g[[3]]] + h[[2]] Cos[g[[3]]] + g[[2]],
     g[[3]] + h[[3]]  }];
SE2[g_] . SE2[h_] ^:= SE2[g] \[CenterDot] SE2[h];

se2[xi_] \[CenterDot] SE2[h_] ^:= InfGen[se2[xi],SE2[h]];

SE2[h_] \[CenterDot] se2[xi_]  ^:= Action[SE2[h],se2[xi]];

(*se2[xi_] \[Wedge] se2[eta_] ^:= Wedge[se2[xi], se2[eta]];*)

dse2[alpha_] \[Cross] se2[eta_] ^:= dse2[fCross, alpha, eta];

se2[xi_] + se2[eta_] ^:=  se2[xi + eta];
dse2[xi_] + dse2[eta_] ^:=  dse2[xi + eta];

a_ * se2[xi_] ^:= se2[a * xi];


(* EndPackage[]; *)
