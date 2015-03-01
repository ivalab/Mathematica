(*

  groups.m - This package provides the necessary functions and
  manipulations required to implement Lie groups in Mathematica.

  Patricio Vela
*)
BeginPackage["LieGroups`"];

Needs["Objects`","ivamatica/Basic/objects.m"];
Needs["Manifolds`","ivamatica/DiffGeometry/manifolds.m"];

SE2::usage="Provides the object functionality of the SE2 group.";
se2::usage="Provides the object functionality of the se2 Lie algebra.";
dse2::usage="Provides the object functionality of the dual to the Lie algebra se2.";
TSE2::uase="Provides the object functionality of the tangent to SE2.";
SE2xse2::usage="Provides the object functionality of the SE2xse2 for
dynamics pruposes.";

(*--SE(2)--*)
gGroup::usage"gGroup[Obj] \n Returns the Lie group associated to the object.";
gTangent::usage"gTangent[Obj] \n Returns the Tangent bundle of the object.";
gReduce::usage"gReduce[Obj] \n Gets the reduced structure of the object.";
gGroupIdentity::usage="Returns the group identity.";
gLieAlgebra::usage="Returns the Lie Algebra.";
gDualLieAlgebra::usage="Returns the dual to the Lie Algebra.";
gDual::usage="Returns the dual.";

oTranslate::usage="oTranslate[Obj, p]";
oLeft::usage="--";
gT::usage="--";
oInfGen::usage="--";
oInverse::usage="oInverse[g] \n Returns the inverse of g.";
gVector::usage="--";
gMag::usage="Gives the magnitude of the vector part.";
gAngle::usage="Gives the angle of the vector part.";
gRotation::usage="Gives the rotation angles for the group object";
gPosition::usage="Gives the positions for the group object";

gBase::usage="gBase[Obj] \n Gets the base point of the fibered object.";
gFiber::usage="gFiber[Obj] \n Gets the fiber of the fibered object.";

gBranch::usage="Gives the angle within the branch cut";
cReduce::usage="cReduce[g, func], cReduce[gp, xi, func] \n
Reduces the G-invariant function of G by the group G.\n Reduces the G-invariant function of TG by the group G.  Returns a function of the Lie algebra of G";

(*--se(2)--*)
oWedge::usage="--";
StructureConstants::usage="--";

(*--se(2)^*--*)
oCross::usage="--";

(*--SE(2), se(2), se(2)^*--*)
oAd::usage="oAd[element, g] \n Takes the Adjoint at g of the element. \ 
  May be in the group, lie algebra, or dual lie algebra. For now, \
  though, only valid for dual.";

(*--Solving dynamics--*)
gEOM::usage="Get equations of motion for system.";
oNDSolve::usage="oNDSolve[Obj, F, q0, int, params]";


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
gDim[SE2xse2] ^:= gDim[SE2] + gDim[se2];
gDim[SE2xse2[q_]] ^:= gDim[SE2] + gDim[se2];
gDim[TSE2] ^:=2 gDim[SE2];

gGroup[SE2[obj_]] ^:= SE2;

gTangent[SE2] ^:= TSE2;

gReduce[SE2] ^:= gGroupIdentity[SE2];
gReduce[TSE2] ^:= SE2xse2;

gGroupIdentity[SE2] ^:= SE2[{0, 0, 0}];
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
oInit[SE2xse2, SE2[g_], se2[xi_]] ^:= SE2xse2[Join[g, xi]];
oInit[TSE2, SE2[g_], gp_] ^:= TSE2[Join[g, gp]];

(*SE2[fInit, g1_ , g2_ , g3_] := {g1,g2,g3};*)

(*--Action of/on SE(2)--*)
oTranslate[SE2[g_], p_] ^:=
  SE2[g + Join[p, {0}]];

oLeft[SE2[g_]] ^:=
  { { Cos[g[[3]]] , -Sin[g[[3]]], g[[1]] },
    { Sin[g[[3]]] ,  Cos[g[[3]]], g[[2]] },
    { 0 , 0, 1 } };

gT[SE2[g_]] ^:= 
  { { Cos[g[[3]]] , -Sin[g[[3]]], 0 },
    { Sin[g[[3]]] ,  Cos[g[[3]]], 0 },
    { 0 , 0, 1 } };

oAd[SE2[g_]] ^:= 
  { { Cos[g[[3]]] , - Sin[g[[3]]] ,   g[[2]] } , 
    { Sin[g[[3]]] ,   Cos[g[[3]]] , - g[[1]] } ,
    { 0 , 0, 1} } ;

oInfGen[se2[xi_], SE2[g_]] ^:=
  se2[{xi[[1]] - xi[[3]] g[[2]] , xi[[2]] + xi[[3]] g[[1]] , xi[[3]]}];

(*--------------------Inverse--------------------*)

oInverse[SE2[g_]] := Module[ {A} ,
  A = - Inverse[SE2[fT,g]];
  SE2[A.g]
];


cReduce[SE2[g_], func_] ^:= Module[{e},
  e = gGroupIdentity[SE2];
  func /. Table[ g[[i]] -> e[[i]], {i, gDim[SE2]} ]
];

cReduce[TSE2[gp_], se2[xi_], func_] ^:= Module[{e},
  e = gGroupIdentity[SE2];
  func /. Join[ Table[ gp[[i+gDim[SE2]]] -> xi[[i]], {i, gDim[se2]}] ,
    Table[ gp[[i]] -> e[[i]], {i, gDim[SE2]} ] ]
];


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
gLocal[SE2xse2[q_]] ^:= q;
gLocal[TSE2[gp_]] ^:= gp;

gBase[SE2xse2[q_]] ^:= SE2[Table[q[[i]], {i, gDim[SE2]}]];
gBase[TSE2[gp_]] ^:= SE2[Table[gp[[i]], {i, gDim[SE2]}]];

gFiber[SE2xse2[q_]] ^:= 
  se2[Table[q[[i]], {i, 1+gDim[SE2], gDim[SE2]+gDim[se2]}]];
gFiber[TSE2[gp_]] ^:= 
  Table[gp[[i]], {i, 1+gDim[SE2], gDim[SE2]+gDim[se2]}];

gMag[SE2[g_]] ^:=  Sqrt[(g[[1]])^2 + (g[[2]])^2];
gAngle[SE2[g_]] ^:= If[ (el[[1]] == 0) && (el[[2]] == 0) , 0, 
  ArcTan[ el[[1]], el[[2]] ] ];
gBranch[SE2[g_], bc_] ^:= Module[{amod},
   amod = Mod[g[[3]], 2 Pi];
   Which [ amod>Max[bc] , amod-= 2 Pi , amod<Min[bc], amod+= 2 Pi];
   amod
  ];

(*--------------------se(2)--------------------*)
oWedge[ se2[xi_], se2[eta_] ] ^:= 
  se2[{  xi[[2]] eta[[3]] - xi[[3]] eta[[2]] ,
    - xi[[1]] eta[[3]] + xi[[3]] eta[[1]] ,
     0  }];
oAd[SE2[g_], se2[xi_]] := 
  se2[{ { xi[[3]] g[[2]] + xi[[1]] Cos[g[[3]]] - xi[[2]] Sin[g[[3]]] },
    {-xi[[3]] g[[1]] + xi[[1]] Sin[g[[3]]] + xi[[2]] Cos[g[[3]]] },
    { xi[[3]] } }];

(*--------------------se(2)^*--------------------*)
oCross[dse2[alpha_], se2[xi_]] ^:=
  dse2[Cross[{alpha[[1]], alpha[[2]], 0} , {xi[[1]], xi[[2]], xi[[3]]}] ];

oAd[dse2[alpha_] , SE2[g_]] ^:= 
  dse2[ { alpha[[1]] Cos[g[[3]]] + alpha[[2]] Sin[g[[3]]] ,
    - alpha[[1]] Sin[g[[3]]] + alpha[[2]] Cos[g[[3]]] ,
    - alpha[[1]] g[[2]] + alpha[[2]] g[[1]] + alpha[[3]] } ];


(*--------------------gEOM--------------------*)
(*
  Compute first order nonlinear matrix ODE of the form
    g' = g F
  for use with NDSolve.
*)
gEOM[SE2[g_], se2[F_]] ^:= Module [
  {hp} , 
  hp = Table[ gLocal[SE2[g]][[i]]' , {i, gDim[SE2]} ];
  Table[hp[[i]] == (gT[SE2[g]].F)[[i]] , {i,gDim[SE2]}]
]

gEOM[SE2[g_], se2[F_], t_] ^:= Module [
  {h,hp} , 
  h = oInit[SE2, Table[ gLocal[SE2[g]][[i]][t], {i, gDim[SE2]}] ];
  hp = Table[ gLocal[SE2[g]][[i]]'[t] , {i, gDim[SE2]} ];
  Table[hp[[i]] == (gT[h].F)[[i]] , {i,gDim[SE2]}]
]

gEOM[SE2[g_], se2[F_], SE2[g0_] , t_] ^:= Module [
  {h, hp, hh} , 
  h  = oInit[SE2,g[[1]][t],g[[2]][t],g[[3]][t]];
  hp = {g[[1]]'[t],g[[2]]'[t],g[[3]]'[t]};
  Join[ Table[hp[[i]] == (gT[h].F)[[i]] , {i,gDim[SE2]}], 
	        Table[g[[i]][0] == g0[[i]] , {i, gDim[SE2]}] ]
]

(*
  Compute the first order nonlinear matrix ODE of the form
    xi' =  X(xi, t)
  for use with NDSolve
*)
gEOM[se2[xi_], X_ , se2[xi0_] , t_] ^:= Module [
  {xip, Xt} , 

  xip = {xi[[1]]'[t],xi[[2]]'[t],xi[[3]]'[t]};
  Xt = X /. Table[ xi[[i]] -> xi[[i]][t], {i, gDim[SE2]} ];

  Join[ Table[xip[[i]] == Xt[[i]] , {i,gDim[SE2]}], 
	Table[xi[[i]][0] == xi0[[i]] , {i,gDim[SE2]}] ]
]

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

(*--------------------oNDSolve--------------------*)
(*
  Solve first order nonlinear matrix ODE for one of the
  following forms:
    1) g'  = g F
    2) xi' = X(xi, t)
    3)  d  [ g  ]    [ g xi ]
       --- [    ]  = [      ]
        dt [ xi ]    [  X   ]
  over the given interval.
*)
oNDSolve[Obj_, F_, g0_ , int_, params___] := Module [{Geqns},
  Geqns = gEOM[Obj, F, g0, int[[1]] ];
  NDSolve[Geqns, gLocal[Obj], int, params]
]

(*--Functions that support the configuration object.
    They are not bound to any instance in particular.  --*)

StructureConstants[LA_, T_] := Module[
  {eb, SC, dimG},

  dimG = gDim[LA];
  eb[i_] := Transpose[T][[i]];
  SC = Table[ Inverse[T].gLocal[ oWedge[LA[eb[j]], LA[eb[k]]] ],
      {j,1,dimG}, {k,1,dimG}];
  SC = Table[ -SC[[j,k,i]], {i, 1, dimG},{j,1,dimG},{k,1,dimG}]
];


Format[SE2[g_]] := g;
Format[se2[xi_]] := MatrixForm[xi];
Format[dse2[alpha_]] := MatrixForm[Identity[
  Table[{alpha[[i]]},{i,gDim[dse2]}] ]];
Format[SE2xse2[q_]] := q;

End[];

SE2[g_] \[CenterDot] SE2[h_] ^:= 
  {  h[[1]] Cos[g[[3]]] - h[[2]] Sin[g[[3]]] + g[[1]],
     h[[1]] Sin[g[[3]]] + h[[2]] Cos[g[[3]]] + g[[2]],
     g[[3]] + h[[3]]  };

se2[xi_] \[CenterDot] SE2[h_] ^:= oInfGen[se2[xi],SE2[h]];

se2[xi_] \[Wedge] se2[eta_] ^:= oWedge[se2[xi], se2[eta]];

dse2[alpha_] \[Cross] se2[eta_] ^:= dse2[fCross, alpha, eta];

se2[xi_] + se2[eta_] ^:=  se2[xi + eta];
dse2[xi_] + dse2[eta_] ^:=  dse2[xi + eta];

EndPackage[];
