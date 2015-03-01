(*--------------------TSE(2)--------------------

  TSO3.m - This package provides the necessary functions and
  manipulations required to implement the tangent and
  cotangent spaces to SE(2), including any reduced structures.

  Author: Patricio Vela
  Date: 07.12.2002
*)


(*--Define objects--*)
TSO3::uase="Provides the object functionality of the tangent to SO3.";
dSO3::uase="Provides the object functionality of the cotangent to SO3.";
SO3xso3::usage="Provides the object functionality of the SO3xso3 for
dynamics pruposes.";
SO3xdso3::usage="Provides the object functionality of the SO3xdso3 for
dynamics pruposes.";


Begin["`Private`"];

gTangent[SO3] ^:= TSO3;
gTangent[SO3[g_]] ^:= TSO3;

gCotangent[SO3] ^:= dSO3;
gCotangent[SO3[g_]] ^:= dSO3;

(*--------------------TSE(2) as a Tangent Space--------------------*)

oInit[TSO3, SO3[g_], Vector[gp_]] ^:= 
  TSO3[oContents[oInit[Bundle, SO3[g], Vector[gp]] ]];

Format[TSO3[gp_]] ^:= Format[TManifold[gp]];

gDim[TSO3] ^:= 2 gDim[SO3];

dBase[TSO3[qp_], SO3[q_]] ^:=
  TSO3[oContents[dBase[Bundle[qp], SO3[q]]]];
dVector[TSO3[qp_], Vector[v_]] ^:=
  TSO3[oContents[dVector[TManifold[qp], Vector[v]]]];

gBase[TSO3[qp_]] ^:= gBase[TManifold[qp]];
gVector[TSO3[qp_]] ^:= gVector[TManifold[qp]];

gLocal[TSO3[qp_]] ^:= gLocal[TManifold[qp]];

gDual[TSO3] ^:= dSO3;


(*--------------------TSE(2) Lie Group Structure--------------------*)

gReduce[TSO3] ^:= SO3xso3;

gGroupIdentity[TSO3] ^:= 
  oInit[TSO3, gGroupIdentity[SO3], 
     oInit[Vector, Table[ 0, {i, gDim[SO3]}]] ] ;
gGroupIdentity[TSO3[gp_]] ^:= gGroupIdentity[TSO3];


cReduce[func_, TSO3[gp_], so3[xi_]] ^:= Module[{e},
  e = gLocal[gGroupIdentity[SO3]];
  func /. Join[ Table[ gLocal[gVector[ TSO3[gp]] ][[i]] -> xi[[i]], 
    {i, gDim[so3]}] ,
    Table[ gLocal[gBase[ TSO3[gp] ]][[i]] -> e[[i]], {i, gDim[SO3]} ] ]
];


(*--------------------Differential Equations--------------------*)
(*
  Print the second order nonlinear matrix ODE of the form
    d  [g ]  = g'
   --- [  ]
   dt  [g'] =  X(g, g', t)
*)
gEOM[TSO3[gp_], X_, TSO3[gp0_], t_] ^:= Module[
  {v, gpt, gt, geqns, gpeqns, Xt},

  v = gVector[TSO3[gp]];
  gpt = oInit[Vector, Table[ gLocal[v][[i]][t], {i, gDim[v]}]];
  geqns = gEOM[SO3[g], gpt, gBase[TSO3[g0]], t];

  Xt = X /. Table[gLocal[TSO3[gp]][[i]] -> gLocal[TSO3[gp]][[i]][t], 
    {i, gDim[TSO3]} ];

  gpeqns = Join[ Table[ gLocal[v][[i]]'[t] == Xt[[i]], {i, gDim[v]}],
    Table[ gLocal[v][[i]][0] == gLocal[gVector[TSO3[g0]]][[i]],
      {i, gDim[v]} ] ];

  {geqns, gpeqns}
]



(*--------------------TSE(2) as a Cotangent Space--------------------*)

oInit[dSO3, SO3[g_], Covector[gp_]] ^:= 
  dSO3[oContents[oInit[Bundle, SO3[g], Covector[gp]] ]];

Format[dSO3[gp_]] ^:= Format[dManifold[gp]];

gDim[dSO3] ^:= 2 gDim[SO3];

dBase[dSO3[qp_], SO3[q_]] ^:=
  dSO3[oContents[dBase[Bundle[qp], SO3[q]]]];
dCovector[TSO3[qp_], Covector[v_]] ^:=
  dSO3[oContents[dVector[dManifold[qp], Covector[v]]]];

gBase[dSO3[qp_]] ^:= gBase[dManifold[qp]];
gCovector[dSO3[qp_]] ^:= gCovector[dManifold[qp]];

gLocal[dSO3[qp_]] ^:= gLocal[dManifold[qp]];

gDual[dSO3] ^:= TSO3;


(*--------------------dSE(2) Lie Group Structure--------------------*)

gReduce[dSO3] ^:= SO3xdso3;

gGroupIdentity[dSO3] ^:= 
  oInit[dSO3, gGroupIdentity[dE2], 
     oInit[Covector, Table[ 0, {i, gDim[SO3]}]] ] ;
gGroupIdentity[dSO3[gp_]] ^:= gGroupIdentity[dSO3];

cReduce[func_, dSO3[gp_], dso3[xi_]] ^:= Module[{e},
  e = gLocal[gGroupIdentity[SO3]];
  func /. Join[ Table[ gLocal[gCovector[ dSO3[gp]] ][[i]] -> xi[[i]], 
    {i, gDim[dso3]}] ,
    Table[ gLocal[gBase[ dSO3[gp] ]][[i]] -> e[[i]], {i, gDim[SO3]} ] ]
];





(*--------------------SO3xso3 as a Bundle Space--------------------*)

oInit[SO3xso3, SO3[g_], so3[xi_]] ^:=
  SO3xso3[oContents[oInit[Bundle, SO3[g], so3[xi]] ]];

Format[SO3xso3[gp_]] := Format[Bundle[gp]];

gDim[SO3xso3[gp_]] ^:= gDim[SO3] + gDim[so3];

dBase[SO3xso3[gp_], SO3[g_] ] ^:=
  SO3xso3[oContents[dBase[Bundle[gp], SO3[g]] ]];
dVector[SO3xso3[gp_], so3[xi_] ] ^:=
  SO3xso3[oContents[dVector[Bundle[gp], so3[xi]] ]];
dFiber[SO3xso3[gp_], so3[xi_] ] ^:= dVector[SO3xso3[gp], so3[xi]];

gBase[SO3xso3[q_]] ^:= gBase[Bundle[q]];
gFiber[SO3xso3[q_]] ^:= gFiber[Bundle[q]];
gVector[SO3xso3[q_]] ^:= gFiber[SO3xso3[q]];

gLocal[SO3xso3[q_]] ^:= gLocal[Bundle[q]];



(*--------------------SO3xdso3 as a Bundle Space--------------------*)

oInit[SO3xdso3, SO3[g_], dso3[a_]] ^:=
  SO3xdso3[oContents[oInit[Bundle, SO3[g], dso3[a]] ]];

Format[SO3xdso3[q_]] := Format[Bundle[q]];

gDim[SO3xdso3[q_]] ^:= gDim[SO3] + gDim[dso3];

dBase[SO3xdso3[gp_], SO3[g_] ] ^:=
  SO3xdso3[oContents[dBase[Bundle[gp], SO3[g]] ]];
dVector[SO3xso3[gp_], dso3[a_] ] ^:=
  SO3xdso3[oContents[dVector[Bundle[gp], so3[a]] ]];
dFiber[SO3xdso3[gp_], dso3[a_] ] ^:= dVector[SO3xdso3[gp], dso3[a]];

gBase[SO3xdso3[q_]] ^:= gBase[Bundle[q]];
gFiber[SO3xdso3[q_]] ^:= gFiber[Bundle[q]];
gVector[SO3xdso3[q_]] ^:= gFiber[SExdso3[q]];

gLocal[SO3xdso3[q_]] ^:= gLocal[Bundle[q]];



(*
  Print the second order nonlinear matrix ODE of the form
    g'  = g.xi
    xi' = F(g, xi, t)
*)
gEOM[SO3[g_], so3[xi_], X_, SO3[g0_], so3[xi0_], t_] ^:= Module[
  {xit, geqns, xieqns, Xt},

  xit = Table[ xi[[i]][t], {i, gDim[so3]}];
  geqns = gEOM[SO3[g], so3[xit], SO3[g0], t];

  Xt = X /. Table[g[[i]] -> g[[i]][t], {i, gDim[SO3]} ];

  xieqns = gEOM[so3[xi], Xt, so3[xi0], t];

  {geqns, xieqns}
]

gEOM[SO3xso3[q_], X_, SO3xso3[q0_], t_] ^:= Module[
  {g0, xi0, g, xi, Xt}, 

  g = Table[q[[i]], {i, gDim[SO3]}];
  xi = Table[q[[i]], {i, 1+gDim[SO3], gDim[SO3]+gDim[so3]}];
  g0 = Table[q0[[i]], {i, gDim[SO3]}];
  xi0 = Table[q0[[i]], {i, 1+gDim[SO3], gDim[SO3]+gDim[so3]}];

  gEOM[SO3[g], so3[xi], X, SO3[g0], so3[xi0], t]
]

(*-------------------------------------------------------------------*)
(*--------------------Operations into TSO3 or dSE--------------------*)
(*-------------------------------------------------------------------*)

(*--------------------D--------------------*)

D[f_, SO3[g_]] ^:= 
  dBase[dSO3[ oContents[ D[f, Euclidean[g]] ]], SO3[g]];

D[SO3[f_], SO3[g_]] ^:= 
  Table[ D[gLocal[SO3[f]][[i]], SO3[g]], {i, gDim[SO3]}];

(*--------------------------------------------------------*)
(*--------------------Other Operations--------------------*)
(*--------------------------------------------------------*)


dSO3[dg_] * TSO3[gv_] ^:=
  dManifold[dg] * TManifold[gv];

  (*If[ gBase[dSO3[dg]] == gBase[TSO3[qv]] ,
    Module[{cv, v},
    cv = gCovector[dSO3[dg]];
    v  = gVector[TSO3[gv]];
    cv.v ]
  ,
    "Not possible." , "Not possible." ];
  dEuclidean[alpha] * TEuclidean[v];*)


End[];
