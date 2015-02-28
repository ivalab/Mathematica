(*--------------------TSE(2)--------------------

  TSE2.m - This package provides the necessary functions and
  manipulations required to implement the tangent and
  cotangent spaces to SE(2), including any reduced structures.

  Author: Patricio Vela
  Date: 07.12.2002
*)


(*--Define objects--*)
TSE2::uase="Provides the object functionality of the tangent to SE2.";
dSE2::uase="Provides the object functionality of the cotangent to SE2.";
SE2xse2::usage="Provides the object functionality of the SE2xse2 for
dynamics pruposes.";
SE2xdse2::usage="Provides the object functionality of the SE2xdse2 for
dynamics pruposes.";


Begin["`Private`"];

gTangent[SE2] ^:= TSE2;
gTangent[SE2[g_]] ^:= TSE2;

gCotangent[SE2] ^:= dSE2;
gCotangent[SE2[g_]] ^:= dSE2;

(*--------------------TSE(2) as a Tangent Space--------------------*)

oInit[TSE2, SE2[g_], Vector[gp_]] ^:= 
  TSE2[oContents[oInit[Bundle, SE2[g], Vector[gp]] ]];

Format[TSE2[gp_]] ^:= Format[TManifold[gp]];

gDim[TSE2] ^:= 2 gDim[SE2];

dBase[TSE2[qp_], SE2[q_]] ^:=
  TSE2[oContents[dBase[Bundle[qp], SE2[q]]]];
dVector[TSE2[qp_], Vector[v_]] ^:=
  TSE2[oContents[dVector[TManifold[qp], Vector[v]]]];

gBase[TSE2[qp_]] ^:= gBase[TManifold[qp]];
gVector[TSE2[qp_]] ^:= gVector[TManifold[qp]];

gLocal[TSE2[qp_]] ^:= gLocal[TManifold[qp]];

gDual[TSE2] ^:= dSE2;


(*--------------------TSE(2) Lie Group Structure--------------------*)

gReduce[TSE2] ^:= SE2xse2;

gGroupIdentity[TSE2] ^:= 
  oInit[TSE2, gGroupIdentity[SE2], 
     oInit[Vector, Table[ 0, {i, gDim[SE2]}]] ] ;
gGroupIdentity[TSE2[gp_]] ^:= gGroupIdentity[TSE2];


cReduce[func_, TSE2[gp_], se2[xi_]] ^:= Module[{e},
  e = gLocal[gGroupIdentity[SE2]];
  func /. Join[ Table[ gLocal[gVector[ TSE2[gp]] ][[i]] -> xi[[i]], 
    {i, gDim[se2]}] ,
    Table[ gLocal[gBase[ TSE2[gp] ]][[i]] -> e[[i]], {i, gDim[SE2]} ] ]
];


(*--------------------Differential Equations--------------------*)
(*
  Print the second order nonlinear matrix ODE of the form
    d  [g ]  = g'
   --- [  ]
   dt  [g'] =  X(g, g', t)
*)
gEOM[TSE2[gp_], X_, TSE2[gp0_], t_] ^:= Module[
  {v, gpt, gt, geqns, gpeqns, Xt},

  v = gVector[TSE2[gp]];
  gpt = oInit[Vector, Table[ gLocal[v][[i]][t], {i, gDim[v]}]];
  geqns = gEOM[SE2[g], gpt, gBase[TSE2[g0]], t];

  Xt = X /. Table[gLocal[TSE2[gp]][[i]] -> gLocal[TSE2[gp]][[i]][t], 
    {i, gDim[TSE2]} ];

  gpeqns = Join[ Table[ gLocal[v][[i]]'[t] == Xt[[i]], {i, gDim[v]}],
    Table[ gLocal[v][[i]][0] == gLocal[gVector[TSE2[g0]]][[i]],
      {i, gDim[v]} ] ];

  {geqns, gpeqns}
]



(*--------------------TSE(2) as a Cotangent Space--------------------*)

oInit[dSE2, SE2[g_], Covector[gp_]] ^:= 
  dSE2[oContents[oInit[Bundle, SE2[g], Covector[gp]] ]];

Format[dSE2[gp_]] ^:= Format[dManifold[gp]];

gDim[dSE2] ^:= 2 gDim[SE2];

dBase[dSE2[qp_], SE2[q_]] ^:=
  dSE2[oContents[dBase[Bundle[qp], SE2[q]]]];
dCovector[TSE2[qp_], Covector[v_]] ^:=
  dSE2[oContents[dVector[dManifold[qp], Covector[v]]]];

gBase[dSE2[qp_]] ^:= gBase[dManifold[qp]];
gCovector[dSE2[qp_]] ^:= gCovector[dManifold[qp]];

gLocal[dSE2[qp_]] ^:= gLocal[dManifold[qp]];

gDual[dSE2] ^:= TSE2;


(*--------------------dSE(2) Lie Group Structure--------------------*)

gReduce[dSE2] ^:= SE2xdse2;

gGroupIdentity[dSE2] ^:= 
  oInit[dSE2, gGroupIdentity[dE2], 
     oInit[Covector, Table[ 0, {i, gDim[SE2]}]] ] ;
gGroupIdentity[dSE2[gp_]] ^:= gGroupIdentity[dSE2];

cReduce[func_, dSE2[gp_], dse2[xi_]] ^:= Module[{e},
  e = gLocal[gGroupIdentity[SE2]];
  func /. Join[ Table[ gLocal[gCovector[ dSE2[gp]] ][[i]] -> xi[[i]], 
    {i, gDim[dse2]}] ,
    Table[ gLocal[gBase[ dSE2[gp] ]][[i]] -> e[[i]], {i, gDim[SE2]} ] ]
];





(*----------SE2xse2 as a (Semidirect Product) Bundle Space----------*)

oInit[SE2xse2, SE2[g_], se2[xi_]] ^:=
  SE2xse2[oContents[oInit[Bundle, SE2[g], se2[xi]] ]];

Format[SE2xse2[gp_]] := Format[Bundle[gp]];

gDim[SE2xse2[gp_]] ^:= gDim[SE2] + gDim[se2];

dBase[SE2xse2[gp_], SE2[g_] ] ^:=
  SE2xse2[oContents[dBase[Bundle[gp], SE2[g]] ]];
dVector[SE2xse2[gp_], se2[xi_] ] ^:=
  SE2xse2[oContents[dVector[Bundle[gp], se2[xi]] ]];
dFiber[SE2xse2[gp_], se2[xi_] ] ^:= dVector[SE2xse2[gp], se2[xi]];

gBase[SE2xse2[q_]] ^:= gBase[Bundle[q]];
gFiber[SE2xse2[q_]] ^:= gFiber[Bundle[q]];
gVector[SE2xse2[q_]] ^:= gFiber[SE2xse2[q]];

gLocal[SE2xse2[q_]] ^:= gLocal[Bundle[q]];



(*--------------------SE2xdse2 as a Bundle Space--------------------*)

oInit[SE2xdse2, SE2[g_], dse2[a_]] ^:=
  SE2xdse2[oContents[oInit[Bundle, SE2[g], dse2[a]] ]];

Format[SE2xdse2[q_]] := Format[Bundle[q]];

gDim[SE2xdse2[q_]] ^:= gDim[SE2] + gDim[dse2];

dBase[SE2xdse2[gp_], SE2[g_] ] ^:=
  SE2xdse2[oContents[dBase[Bundle[gp], SE2[g]] ]];
dVector[SE2xse2[gp_], dse2[a_] ] ^:=
  SE2xdse2[oContents[dVector[Bundle[gp], se2[a]] ]];
dFiber[SE2xdse2[gp_], dse2[a_] ] ^:= dVector[SE2xdse2[gp], dse2[a]];

gBase[SE2xdse2[q_]] ^:= gBase[Bundle[q]];
gFiber[SE2xdse2[q_]] ^:= gFiber[Bundle[q]];
gVector[SE2xdse2[q_]] ^:= gFiber[SExdse2[q]];

gLocal[SE2xdse2[q_]] ^:= gLocal[Bundle[q]];



(*
  Print the second order nonlinear matrix ODE of the form
    g'  = g.xi
    xi' = F(g, xi, t)
*)
gEOM[SE2[g_], se2[xi_], X_, SE2[g0_], se2[xi0_], t_] ^:= Module[
  {xit, geqns, xieqns, Xt},

  xit = Table[ xi[[i]][t], {i, gDim[se2]}];
  geqns = gEOM[SE2[g], se2[xit], SE2[g0], t];

  Xt = X /. Table[g[[i]] -> g[[i]][t], {i, gDim[SE2]} ];

  xieqns = gEOM[se2[xi], Xt, se2[xi0], t];

  {geqns, xieqns}
]

gEOM[SE2xse2[q_], X_, SE2xse2[q0_], t_] ^:= Module[
  {g0, xi0, g, xi, Xt}, 

  g = Table[q[[i]], {i, gDim[SE2]}];
  xi = Table[q[[i]], {i, 1+gDim[SE2], gDim[SE2]+gDim[se2]}];
  g0 = Table[q0[[i]], {i, gDim[SE2]}];
  xi0 = Table[q0[[i]], {i, 1+gDim[SE2], gDim[SE2]+gDim[se2]}];

  gEOM[SE2[g], se2[xi], X, SE2[g0], se2[xi0], t]
]

(*-------------------------------------------------------------------*)
(*--------------------Operations into TSE2 or dSE--------------------*)
(*-------------------------------------------------------------------*)

(*--------------------D--------------------*)

D[f_, SE2[g_]] ^:= 
  dBase[dSE2[ oContents[ D[f, Euclidean[g]] ]], SE2[g]];

D[SE2[f_], SE2[g_]] ^:= 
  Table[ D[gLocal[SE2[f]][[i]], SE2[g]], {i, gDim[SE2]}];

(*--------------------------------------------------------*)
(*--------------------Other Operations--------------------*)
(*--------------------------------------------------------*)


dSE2[dg_] * TSE2[gv_] ^:=
  dManifold[dg] * TManifold[gv];

  (*If[ gBase[dSE2[dg]] == gBase[TSE2[qv]] ,
    Module[{cv, v},
    cv = gCovector[dSE2[dg]];
    v  = gVector[TSE2[gv]];
    cv.v ]
  ,
    "Not possible." , "Not possible." ];
  dEuclidean[alpha] * TEuclidean[v];*)


End[];
