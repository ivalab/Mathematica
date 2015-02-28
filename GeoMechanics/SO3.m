(*
  SO3.m - This package provides the necessary functions and
  manipulations required to implement the Lie group SO(3) in 
  Mathematica.

  Patricio Vela
  1/21/2003

  Note: Try to make this as general as possible so that it may
  be converted into a package for SO(n).
*)

SO3::usage="Provides the object functionality of the SO3 group.";
so3::usage="Provides the object functionality of the so3 Lie algebra.";
dso3::usage="Provides the object functionality of the dual to the Lie algebra so3.";

Begin["`Private`"];

(*--------------------SE(2)--------------------*)
(* S0(3) is given by g = (?,?,?,?,...).  
   It's Lie algebra s0(3) is xi = (w1, w2, w3)
   It's dual so(2)^* is alpha = (dw1, dw2, dw3)
*)

(*--Constant definitions of the Lie group SE(2).--*)
gDim[SO3] ^:= 3;
gDim[SO3[g_]] ^:= 3;
gDim[so3] ^:= gDim[SO3];
gDim[so3[xi_]] ^:= gDim[SO3];
gDim[dso3] ^:= gDim[SO3];
gDim[dso3[a_]] ^:= gDim[SO3];

gGroup[SO3[obj_]] ^:= SO3;

gReduce[SO3] ^:= gGroupIdentity[SO3];

gGroupIdentity[SO3] ^:= SO3[IdentityMatrix[gDim[SO3]]];
gGroupIdentity[SO3[g_]] ^:= gGroupIdentity[SO3];
gLieAlgebra[SO3] ^:= so3;
gLieAlgebra[SO3[g_]] ^:= so3;
gDualLieAlgebra[SO3] ^:= dso3;
gDualLieAlgebra[SO3[g_]] ^:= dso3;
gDual[so3] ^:= dso3;
gDual[so3[xi_]] ^:= dso3;

oInit[SO3, th1_ , th2_ , th3_] ^:= 
    SO3[{{Sec[th2] Cos[th3] , Sec[th2] Sin[th3], 0},
         {-Sin[th3], Cos[th3], 0},
	 {Cos[th3] Tan[th2], Sin[th3] Tan[th2], 1}}];

oInit[SO3, g_] ^:= SO3[g];
oInit[so3, xi1_ , xi2_ , xi3_] ^:= 
  so3[{{0, -xi3, xi2},{xi3, 0, -xi1},{-xi2, xi1, 0}}];
oInit[so3, xi_] ^:= so3[xi];
oInit[dso3, a1_ , a2_ , a3_] ^:= 
  dso3[{{0, xi3, -xi2},{-xi3, 0, xi1},{xi2, -xi1, 0}}];
oInit[dso3, a_] ^:= dso3[a];

oName[SO3, g_] ^:= Module[{dim},
  dim = gDim[SO3];
  SO3[Table[ g[i*dim+j], {i,0,dim-1},{j,dim}]]
];

(*--Group Actions--*)

Left[SO3[g_], SO3[h_]] ^:= SO3[g.h];

Right[SO3[g_], SO3[h_]] ^:= SO3[h.g];

Ad[SO3[g_], SO3[h_]] ^:= SO3[ Inverse[g].h.g ];

Left[so3[xi_], SO3[g_]] ^:= oInit[Vector, g.xi];
Right[so3[xi_], SO3[g_]] ^:= oInit[Vector, xi.g];

InfGen[so3[xi_], SO3[g_]] ^:= Action[so3[xi],SO3[g]];

(*--------------------Inverse--------------------*)

Inverse[SO3[g_]] := SO3[Inverse[g]];


(*--------------------Reduce--------------------*)

cReduce[func_, SO3[g_]] ^:= Module[{e},
  e = gGroupIdentity[SO3];
  func[e]
];

(*--------------------Other functions--------------------*)
gVector[SO3[g_]] ^:= Flatten[g];
gVector[so3[xi_]] ^:= {xi[[3,2]], xi[[1,3]], xi[[2,1]]};
gVector[dso3[alpha_]] ^:= {xi[[2,3]], xi[[1,2]], xi[[3,1]]};

gLocal[SO3[g_]] ^:= g;
gLocal[so3[xi_]] ^:= xi;
gLocal[dso3[alpha_]] ^:= alpha;

(*--------------------se(2)--------------------*)
Wedge[ so3[xi_], so3[eta_] ] ^:= so3[ xi.eta - eta.xi ];

Ad[SO3[g_], so3[xi_]] := so3[Inverse[g].xi.g];
ad[so3[xi_], so3[eta_] ] ^:= oWedge[ so3[xi], so3[eta] ];

(*--------------------se(2)^*--------------------*)
Cross[dso3[alpha_], so3[xi_]] ^:=
  dso3[{0,0,0}];

Ad[dso3[alpha_] , SO3[g_]] ^:= 
  dso3[ {0,0,0} ];


(*--------------------gEOM--------------------*)
(*
  Compute first order nonlinear matrix ODE of the form
    g' = F
  for use with NDSolve.
*)

gEOM[Vector[F_], SO3[g_]] ^:= Module [
  {dim} , 
  dim = gDim[SO3];
  Flatten[Table[g[[i,j]]' == gLocal[Vector[F]][[i,j]], 
    {i, dim}, {j, dim}]]
]

gEOM[Vector[F_], SO3[g_], t_] ^:= Module [
  {dim, Ft} , 
  dim = gDim[SO3];
  Ft = gLocal[Vector[F]] /. 
         Flatten[Table[g[[i,j]] -> g[[i,j]][t],{i, dim}, {j, dim}]];
  Flatten[Table[g[[i,j]]'[t] == Ft[[i,j]], {i, dim}, {j, dim}]]
]

gEOM[Vector[F_], SO3[g0_], SO3[g_], time_] ^:= Module [
  {eom, ic, dim} , 
  dim = gDim[SO3];
  eom = gEOM[Vector[F], SO3[g], time[[1]]];
  ic = Flatten[Table[g[[i,j]][time[[2]]] == g0[[i,j]], 
    {i, dim}, {j, dim}]];
  Join[eom, ic]
]

(*--------------------gEOM--------------------*)
(*
  Compute first order nonlinear matrix ODE of the form
    g' = F g
  for use with NDSolve.
*)
gEOM[so3[F_],SO3[g_]] ^:= Module [
  {hp} , 
  hp =  InfGen[so3[F],SO3[g]];
  gEOM[hp, SO3[g]]
]

gEOM[so3[F_], SO3[g_], t_] ^:= Module [
  {hp} , 
  hp =  InfGen[so3[F],SO3[g]];
  gEOM[hp,SO3[g], t]
]

gEOM[so3[F_], SO3[g0_] , SO3[g_], t_] ^:= Module [
  {hp} , 
  hp =  InfGen[so3[F],SO3[g]];
  gEOM[hp,SO3[g0],SO3[g], t]
]

(*
  Compute the first order nonlinear matrix ODE of the form
    xi' =  X(xi, t)
  for use with NDSolve
*)
(*--Not updated from SE2--*)
gEOM[so3[xi_], X_ , so3[xi0_] , t_] ^:= Module [
  {xip, Xt} , 

  xip = {xi[[1]]'[t],xi[[2]]'[t],xi[[3]]'[t]};
  Xt = X /. Table[ xi[[i]] -> xi[[i]][t], {i, gDim[SO3]} ];

  Join[ Table[xip[[i]] == Xt[[i]] , {i,gDim[SO3]}], 
	Table[xi[[i]][0] == xi0[[i]] , {i,gDim[SO3]}] ]
]

(*------------------------NDSolve-------------------------*)

NDSolve[so3[xi_], SO3[g0_], SO3[g_], time_,params___] ^:=
  Module[{eom},

  eom = gEOM[so3[xi], SO3[g0], SO3[g], time];
  NDSolve[eom, Flatten[g], time, params]
]



(*==================================================================*)
(*--------------------Exponential representation--------------------*)
(*==================================================================*)

(*Exp[so3[xi_]] ^:= oInit[SO3, MatrixExp[xi]];*)
Exp[so3[xi_]] ^:= Module[{mag},
  mag = PowerExpand[Sqrt[gVector[so3[xi]].gVector[so3[xi]] ]];
  oInit[SO3, IdentityMatrix[gDim[SO3]] + (Sin[mag]/mag) xi 
      + ((1-Cos[mag])/(mag^2)) xi.xi]
];

Log[SO3[IdentityMatrix[3]]] ^:= oInit[so3, 0, 0, 0];
Log[SO3[g_]] ^:= Module[{th = ArcCos[(Tr[g] - 1)/2]},
  oInit[so3, (th/(2 Sin[th])) (g - Transpose[g])]
];

(*--Functions that support the configuration object.
    They are not bound to any instance in particular.  --*)
(*--No need to have it here.--*)
(*
StructureConstants[LA_, T_] := Module[
  {eb, SC, dimG},

  dimG = gDim[LA];
  eb[i_] := Transpose[T][[i]];
  SC = Table[ Inverse[T].gLocal[ oWedge[LA[eb[j]], LA[eb[k]]] ],
      {j,1,dimG}, {k,1,dimG}];
  SC = Table[ -SC[[j,k,i]], {i, 1, dimG},{j,1,dimG},{k,1,dimG}]
];
*)

Format[SO3[g_]] := MatrixForm[g];
Format[so3[xi_]] := MatrixForm[xi];
Format[dso3[alpha_]] := MatrixForm[alpha];

End[];

SO3[g_] == SO3[h_] ^:=
  Euclidean[g] == Euclidean[h];

SO3[g_] \[CenterDot] SO3[h_] ^:= Action[SO3[g],SO3[h]];

SO3[g_] . SO3[h_] ^:= SO3[g] \[CenterDot] SO3[h];

so3[xi_] \[CenterDot] SO3[h_] ^:= InfGen[so3[xi],SO3[h]];

so3[xi_] \[Wedge] so3[eta_] ^:= Wedge[so3[xi], so3[eta]];

dso3[alpha_] \[Cross] so3[eta_] ^:= dso3[fCross, alpha, eta];

so3[xi_] + so3[eta_] ^:=  so3[xi + eta];
dso3[xi_] + dso3[eta_] ^:=  dso3[xi + eta];

a_ * so3[xi_] ^:= so3[a * xi];

