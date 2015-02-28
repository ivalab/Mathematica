
(*

  sobstacles.m - This package provides the necessary functions and
  manipultions required to implement a sphere-world forumalation of
  a given robot and its configuration space.  The robot is a point mass
  and the obstacles are circles in the sphereworld.

*)

BeginPackage["SObstacles`"];

(* Constants *)
xi = 1;
eta = 1;
varphi = 100;
rho1 = 2;
rho2 = 0.01;
rho3 = 0.01;
delta = 0.02;
pow = 2;

(* 
   An s-obstacle is given by a position and radius.
   Obstacles have nonzero radius.
   The world has a nonzero radius.
   The robot and the goal have zero radius. 
*)


(*- Initialization of s-obstacles. -*)
SO[Init, x_, y_, r_] := {x,y,r};
SO[InitWorld, x_, y_, r_] := SO[Init, x, y, r];


(*- Information retrieval of the s-obstacle -*)
SO[Position, so_] := Take[so, 2];
SO[Radius, so_] :=  so[[3]] ;

(*- Set the data of the s-obstacle. -*)
SO[SetPosition, so_, pos_] := 
   SO[Init, pos[[1]], pos[[2]], SO[Radius, so] ];
SO[SetRadius, so_, radius_] := 
   SO[Init, so[[1]], so[[2]], radius ];

(*- Set the various parameters of the s-obstacle potential fields. -*)
SO[SetXi, newxi_] :=  xi = newxi;
SO[SetEta, neweta_] :=  eta = neweta;
SO[SetRho1, newrho1_] :=  rho1 = newrho1;
SO[SetRho2, newrho2_] :=  rho2 = newrho2;
SO[SetPow, newpow_] :=  pow = newpow;

(*- The distance functions and their gradients. -*)

(*- The distance of a robot to an s-obstacle. -*)
SO[d, ro_, so_] := Module[ {dist} , 
  dist = SO[Position, ro] - SO[Position, so];
  dist = Sqrt[dist.dist];
  If[ dist > SO[Radius, so] , 
    dist - SO[Radius, so]
  , (*else*)
    0
  ]
]

SO[Gradd, ro_, ob_] := Module[ {dist} , 
  dist = SO[d, ro, ob];
  If[ dist == 0 ,
    {0, 0} 
  ,  (*else*)
    (1/(dist+rho2)) * (SO[Position, ro] - SO[Position, ob])
  ]
];


(*- The distance of a robot to the edge of the world. -*)
SO[dedge, ro_, wo_] := Module[ {dist} , 
  dist = SO[Position, ro] - SO[Position, wo];
  dist = Sqrt[dist.dist];
  dist = SO[Radius, wo] - dist;
  If[ dist < 0 , -dist , dist]
]

SO[Graddedge, ro_, bo_] := Module[ {dist} , 
  dist = SO[dedge, ro, bo];
  If[ dist == 0 ,
    {0, 0} 
  ,  (*else*)
    - (1/(dist+rho2)) * (SO[Position, ro]-SO[Position, bo])
  ]
];

(*- The distance of a robot to the goal. -*)
SO[dgoal, ro_, go_] := Module[ {dist} , 
  dist = SO[Position, ro] - SO[Position, go];
  dist = Sqrt[dist.dist]
]

SO[Graddgoal, ro_, go_] := Module[ {dist} , 
  dist = SO[dedge, ro, go];
  If[ dist == 0 ,
    {0, 0} 
  ,  (*else*)
    (1/(dist+rho2)) * (SO[Position, ro] - SO[Position, go])
  ]
];


(*- The potentials and their gradients. -*)

(*- Attractive potential, robot -> goal. -*)
SO[Uattractive, ro_, go_] := (1/2)*xi*(SO[dgoal, ro, go])^2;

SO[GradUattractive, ro_, go_] := Module[ {dist}, 
  dist = SO[dgoal, ro, go];
  xi*dist*SO[Graddgoal, ro, go] 
];

(*- Repulsive potential, robot <- object. -*)
SO[Urepulsive, ro_, ob_] := Module[ {dist} , 
  dist = SO[d, ro, ob];
  If [ dist <= rho1^2 , 
    (eta/2) * (1/(dist+rho2) - 1/rho1)^2
  , (*else*)
    0
  ]
]

SO[GradUrepulsive, ro_, ob_] := Module[ {dist, pot}, 
  dist = SO[d, ro, ob];
  pot = SO[Urepulsive, ro, ob];
  If[ dist == 0 || pot == 0,
    {0, 0} 
  ,  (*else*)
    - eta * (1/(dist+rho2)^2) * (1/(dist+rho2) - 1/rho1) * SO[Gradd, ro, ob]
  ]
];

(*- Repulsive potential, robot <- world. -*)
SO[Uboundary, ro_, bo_] := Module[ {dist} , 
  dist = SO[dedge, ro, bo];
  If [ dist <= rho1^2 , 
    (eta/2) * (1/(dist+rho2) - 1/rho1)^2
  , (*else*)
    0
  ]
]

SO[GradUboundary, ro_, bo_] := Module[ {dist, pot}, 
  dist = SO[dedge, ro, bo];
  pot = SO[Uboundary, ro, bo];
  If[ dist == 0 || pot == 0,
    {0, 0} 
  ,  (*else*)
    - eta * (1/(dist+rho2)^2) * (1/(dist+rho2) - 1/rho1) 
      * (SO[Graddedge, ro, bo])
  ]
];


(*- The potentials that we use. -*)

SO[Urk , ro_, go_, bo_, objects_] := Module[  
  {dist, Uden, Unum} 
  , 
     (* Compute the distances. *)
  dist = Join[ {SO[dedge, ro, bo]} , 
	Table[ SO[d, ro, objects[[i]] ] , {i, 1, Length[objects]} ] ];
     (* Take the product of the distances.  Our denominator. *)
  Uden = Product [ dist[[i]] , {i, 1, Length[dist]} ] + rho2;

     (* The is the numerator * gradient of denominator *)
  Unum = (1/2)*xi*(SO[dgoal, ro, go])^pow;
  Unum/Uden
];

SO[GradUrk, ro_ , go_, bo_, objects_ ] := Module[ 
  {dist, graddist, Uden, U1, U2} 
  , 
     (* Compute the distances. *)
  dist = Join[ {SO[dedge, ro, bo]} , 
	Table[ SO[d, ro, objects[[i]] ] , {i, 1, Length[objects]} ] ];
     (* Compute the gradients of the distances. *)
  graddist = Join[ {SO[Graddedge, ro, bo]} , 
	Table[ SO[Gradd, ro, objects[[i]] ] , {i, 1, Length[objects]} ] ];
     (* Take the product of the distances.  Our denominator. *)
  Uden = Product [ dist[[i]] , {i, 1, Length[dist]} ] + rho2;

     (* This is the gradient of the numerator / denominator *)
  U1 = pow*(SO[dgoal, ro, go])^(pow-1) * SO[Graddgoal, ro, go] / Uden;

     (* The is the numerator * gradient of denominator *)
  U2 = Sum[ ( Product[ dist[[i]] , {i, 1, j-1} ] 
	 * graddist[[j]] 
	 * Product[ dist[[i]], {i, j+1, Length[dist]} ] ) 
       , {j, 1, Length[dist]} ];
  U2 *= ((SO[dgoal, ro, go])^pow)/(Uden^2);
  U1-U2
]

SO[Upv, ro_, go_, wo_ , objects_] := Module[
  {dist, U1, U2}
  ,
  dist = SO[dgoal, ro, go];
  U1 = SO[Uattractive, ro, go]  - varphi/(SO[dgoal,ro,go]+rho3)^2;
     (* Compute the distances. *)
  U2 = Join[ {SO[Uboundary, ro, wo]} , 
	  Table[ SO[Urepulsive, ro, objects[[i]] ] 
	    , {i, 1, Length[objects]} ] ];
     (* Take the product of the distances.  Our denominator. *)
  U1+Sum[ U2[[i]], {i, 1, Length[U2]} ]
]

SO[GradUpv, ro_, go_, wo_, objects_] := Module[
  {gradU1, gradU2, gradtemp}
  ,
  dist = SO[dgoal, ro, go];
  If[ dist == 0  ,
    gradtemp = {0, 0}
  ,  (*else*) 
    gradtemp = varphi*(1/(dist+rho3)^2)*SO[Graddgoal, ro, go]
  ];
  gradU1  = SO[GradUattractive, ro, go] + gradtemp;
  gradU2 = Join[ { SO[GradUboundary, ro, wo]} , 
	      Table[ SO[GradUrepulsive, ro, objects[[i]] ]
		, {i, 1, Length[objects]} ] ];
  gradU1 + Sum[ gradU2[[i]], {i, 1, Length[gradU2]} ]
]

(*- Draw the list of objects passed in. -*)
SO[Draw, objects_] := Module[ {gfx, object, i} , 
  
  gfx = {};
  Do[ Module[ {} , 
    object = objects[[i]];
    If[ SO[Radius, object] > 0 , 
      gfx = Append[gfx, Circle[ SO[Position, object] , SO[Radius, object] ] ] 
    , (*else*)
      gfx = Append[gfx, Point[ SO[Position, object] ] ] 
    ]; ]
, {i, Length[objects]}];

  Show[Graphics[gfx] ,{ AspectRatio->Automatic} ];

]

(*- Draw the lists of objects passed in plus the path given. -*)
SO[DrawPath, path_, objects_] := Module[ {gfx, object, i} ,

  gfx = {};
  Do[ Module[ {} , 
    object = objects[[i]];
    If[ SO[Radius, object] > 0 , 
      gfx = Append[gfx, Circle[ SO[Position, object] , SO[Radius, object] ] ] 
    , (*else*)
      gfx = Append[gfx, Point[ SO[Position, object] ] ] 
    ]; ]
, {i, Length[objects]}];

  gfx = Append[ gfx, Line[path] ];
  Show[Graphics[gfx] , {AspectRatio->Automatic, PlotRange->All} ];

]

EndPackage[];
