BeginPackage["LagrangianMechanics`"];

Needs["Objects`","mathematica/libs/objects.m"];
Needs["Euclidean`","mathematica/libs/euclidean.m"];
Needs["Bundles`","mathematica/libs/bundles.m"];
Needs["Vectors`","mathematica/libs/vectors.m"];
Needs["Tangents`","mathematica/libs/tangents.m"];
Needs["Manifolds`","mathematica/libs/manifolds.m"];
Needs["TangentManifolds`","mathematica/libs/tmanifolds.m"];
Needs["LieGroups`","mathematica/libs/liegroups.m"];


fLagrangeEqns::usage=
  "fLagrangeEqns[Lagrangian, Forces, q, q'] \n
   Solves for Lagrange's equations of motion with forces given a
Lagrangian system.";

fConsLagrangeEqns::usage=
  "fConsLagrangeEqns[Lagrangian, Forces,Constraints,q, q', Multiplier] \n
   Solves for Lagrange's equations of motion with forces given a
Lagrangian system with bilateral constraints.  ";

fConsLagrangeEqns::notes="This is a test.";
gLagrangianMetric::usage="--";
fGetLagrangianMetric::usage="--";
fReduceLagrangian::usage="fReduceLagrangian[L, g, e, xi]";
fReduceTransformation::usage="--";
fTransformLagrangian::usage="--";
fTransformLagrangianMatrix::usage="--";
fFDer::usage="Fiber derivative of the Lagrangian, L, at the point (q, q').";

fLift::usage="fLift[Y]\n Lift the vector Y to T2Q";


Begin["`Private`"];

(*--------------------fLagrangeEqns--------------------*)
(*Lagrange's equations of motion with forcing

  (d/dt) (d L / d q') - (d L / d q) = F 

  The form of the parameter group (F, q, qp) must agree in 
  dimension.
*)
fLagrangeEqns[L_ , F_ , q_, qp_] := Module[
  {Lnew, dqp, dq, lqp, qt, qpt} ,
  Lnew = L;
  If[ Length[q] > 0 , Module[ {} , 
    dqp = Table[ D[L,qp[[i]] ] , {i, Length[qp]}] ;
    For[i=1,i<=Length[qp],i++, Lnew = Lnew /. {qp[[i]] -> lqp[i]}];
    dq  = Table[ D[Lnew,q[[i]]] , {i, Length[q]} ]; 
    For[i=1,i<=Length[qp],i++, dq = dq /. {lqp[i] -> qp[[i]]}];
    qt = Table[ q[[i]][t] , {i, Length[q]} ];
    qpt = Table[ qp[[i]][t] , {i, Length[q]} ];
    For[i=1,i<=Length[q],i++, dqp = dqp /. {qp[[i]] -> lqp[i]}];
    For[i=1,i<=Length[q],i++, dqp = dqp /. {q[[i]] -> qt[[i]]}];
    For[i=1,i<=Length[q],i++, dqp = dqp /. {lqp[i] -> qpt[[i]]}];
    dqp = D[dqp,t];
    For[i=1,i<=Length[q],i++, dqp = dqp /. {qt[[i]] -> q[[i]]}];
    For[i=1,i<=Length[q],i++, dqp = dqp /. {qpt[[i]] -> qp[[i]]}];
    For[i=1,i<=Length[q],i++, dqp = dqp /. {q[[i]]''[t] -> q[[i]]''}];
    Table[ dqp[[i]] - dq[[i]] == F[[i]] , {i, Length[q]} ]
  ] , Module[ {} ,
    dqp = D[L,qp]; 
    Lnew = Lnew /. {qp -> lqp};
    dq  = (D[Lnew,q ]) /. {lqp -> qp};
    dqp = dqp /. {qp -> lqp};
    dqp = dqp /. {q -> q[t]};
    dqp = dqp /. {lqp -> qp[t]};
    dqp = D[dqp,t];
    dqp = dqp /. {q[t]-> q};
    dqp = dqp /. {qp[t] -> qp};
    dqp = dqp /. {q''[t] -> q''};
    {dqp - dq == F } 
  ] ]
];

(*--------------------fConsLagrangeEqns--------------------*)
(*Lagrange's equations of motion with forcing for bilaterally 
  constrained systems.

  (d/dt) (d L / d q') - (d L / d q) = F - lambda . (d phi / d q)

  The form of the parameter groups (F, q, qp) and (phi, lambda) 
  must agree in dimension.
*)
fConsLagrangeEqns[L_, F_, phi_, q_, qp_, lambda_] := 
  If[ Length[q] > 0 , 
    If[ Length[lambda] > 0 ,
      fLagrangeEqns[L, F - lambda . Transpose[Table[ D[phi, q[[i]]], 
	{i,Length[q]}]], q, qp]
    , (*else*)
      fLagrangeEqns[L, F - lambda Table[ D[phi, q[[i]] ], {i,Length[q]}], 
       q, qp] ]
  , (*else*)
    If[ Length[lambda] > 0 ,
      fLagrangeEqns[L, F - lambda . D[phi,q], q, qp]
    , (*else*)  
      fLagrangeEqns[L, F - lambda D[phi, q], q, qp] ]
  ]; 

(*--------------------gLagrangianMetric--------------------*)
(*Compute the Lagrangian kinetic energy metric for the Lagrangian
  given and the configuration tangent space variables passed along.
*)
gLagrangianMetric[L_, xp_] := 
  Table[ D[D[L,xp[[i]]],xp[[j]]] ,
  {i,Length[xp]},{j,Length[xp]}];

(*--------------------fReduceLagrangian--------------------*)
(*Reduce a Lagrangian under group symmetries.  This is done
  by solving for the Lagrangian at the group identity.
*)
fReduceLagrangian[L_, g_, e_, xi_] := Module[
  {l, lg, le, lxi} , 
   
  l = L;
  lg = gLocal[g];
  le = gLocal[e];
  lxi = gLocal[xi];
  l = l /. Table[ lg[[i]]' -> lgp[i] , {i, gDim[g]} ];
  l = l /. Table[ lg[[i]]' -> le[[i]] , {i, gDim[g]}];
  l = l /. Table[ lgp[i]   -> lxi[[i]] , {i, gDim[g]}];
  l
];

(*--------------------fReduceTransformation--------------------*)
(*Reduce a transformation of Lie algebra coordinates under group 
  symmetries.  This is done by solving for the transformation
  at the group identity.
*)
fReduceTransformation[T_, g_, e_] :=  Module[
  {t} , 
  t = T;
  For[i = gDim[g] , i > 0, i--, t = t /. {gLocal[g][[i]] -> gLocal[e][[i]]}]; 
  t
];

(*-------------------fTransformLagrangian---------------------*)
(*Transform the Lie algebra basis of the reduced Lagrangian.
*)
fTransformLagrangian [L_, T_, xi_, e_] := Module[
  {l , f , v} , 
   
  l = L;
  f = e.T;
  l = l /. Table[ xi[[i]] -> f[[i]] , {i, Length[xi]}];
  l
];

(*-------------------fTransformLagrangian---------------------*)
(*Transform the tangent space basis of the Lagrangian Matrix.
*)
fTransformLagrangiaMatrix [LM_, T_] := Module[
  {lm, f} , 

  lm = Transpose[T].LM.T
];

(*--------------------fFDer--------------------*)
(*Fiber derivative of the Lagrangian.

  |F L(q,q') = d L(q,q') / d q'
*)
fFDer[L_, q_, qdot_] := Module[
  {z, zdot, F, x, xdot} , 

  z = Table[ x[i], {i, Length[q]}];
  zdot = Table[ xdot[i] , {i, Length[q]}];
  F = Table[ D[L[z, zdot], zdot[[i]]] , {i, Length[qdot]}];
  For[i=Length[q], i>0, i--, F=F /. {z[[i]] -> q[[i]]}];
  For[i=Length[q], i>0, i--, F=F /. {zdot[[i]] -> qdot[[i]]}];
  F
];

fF[L_,  xi_] := Module[
  {F} , 

  F = Table[ D[L[xi], xi[[i]] ]  , {i, Length[qdot]}];
  F
];

(*---------------------fLift--------------------*)
(*Lift of Y to T2Q.
*)
fLift[Y_] := 
  Join[ Table[ 0, {i, Length[Y]}] , Y];


End[];

EndPackage[];
