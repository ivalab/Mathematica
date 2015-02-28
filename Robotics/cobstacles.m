(*

  cobstacles.m - This package provides the necessary functions and
  manipulations required to implement a c-obstacle formulation of a
  given robot and the calculations to determine boundaries for two
  c-obstacles.

  Patricio Vela
*)

BeginPackage["CObstacles`"];

(* A c-obstacle is given by a list of vertices for the polygonal oject
   given with respect to the frame of reference of the object, and
   the location of the body reference frame in the global reference frame.

   The reference frame should be a Lie group, ie. SE(2) or SE(3), etc.
   and have its own set of definitions.  The points should also be given
   with respect to the definition of the Lie group.
*)

CO[Init, G_, Rf_, vert_] := {G, Rf,  vert};
CO[Length, co_] := Length[ co[[3]] ];
CO[Group, co_] := co[[1]];
CO[Frame, co_] := co[[2]];
CO[Vertices, co_] := co[[3]];

CO[Polygon, co_ ] := Module[
  {h,p} , 
  h = CO[Group, co][LieGroups`L, CO[Frame, co]];
  p = Transpose[h.Transpose[CO[Vertices, co]]];
  p = Table[ Take[p[[i]], 2] , {i , 1, Length[p]}];
  Polygon[p] 
]

CO[SetOrientation, co_ , theta_] := 
  ReplacePart[co, Append[ Take[CO[Frame, co], 2] , theta] , 2];

CO[Translate, co_ , g_] := Module[
  CO[Frame, co] = CO[Group, co][g] * CO[Group, co][CO[Frame, co]];
]

CO[Place, co_, h_] := ReplacePart[co, h, 2];

(*- Moves vertex number vn to the point given by p, in the reference
    frame.
*)
CO[Moveto, co_ , p_ , vn_] := Module[
  {h, v} , 
  h = CO[Group, co][LieGroups`L, CO[Frame, co]];
  v = p - h.CO[Vertices, co][[vn]];
  h = CO[Group, co][LieGroups`Translate, CO[Frame, co] , Take[v,2]];
  ReplacePart[co, h, 2 ] 
]

CO[Obstacle , co1_ , co2_] := Module[
  {i, in, j, b, v1,v2, h1, h2,l1,l2,f1,f2,V2,J, co} ,

    (*- The vertices of object #1 in the reference plane -*)
  h1 =  CO[Group, co1][LieGroups`L, CO[Frame, co1]];
  v1 = Transpose[h1.Transpose[CO[Vertices, co1]]];
  v1 = Table[ Take[v1[[i]], 2] , {i , 1, Length[v1]}]; 
    (*- The vertices of object #2 in the reference plane -*)
  h2 =  CO[Group, co2][LieGroups`L, CO[Frame, co2]];
  v2 = Transpose[h2.Transpose[CO[Vertices, co2]]];
  v2 = Table[ Take[v2[[i]], 2] , {i , 1, Length[v2]}]; 
    (*- Determine the tanget vectors along the edges of object #1 -*)
  l1 = CO[Length, co1];
  f1[t_] := v1[[t+1]] - v1[[t]];
  V1 = Table[ f1[i] , {i , 1 , l1-1} ];
  V1 = Append[ V1, v1[[1]] - v1[[l1]] ]; 
    (*- Determine the vector normals to the edges of object #2 -*)
  l2 = CO[Length, co2];
  f2[t_] := v2[[t+1]] - v2[[t]];
  V2 = Table[ f2[i] , {i , 1 , l2-1} ];
  V2 = Append[ V2, v2[[1]] - v2[[l2]] ]; 
  J = { {0 , 1} , {-1 , 0} };
  V2 = Transpose[J.Transpose[V2]];
    (*- Go through the vertices of object #1 to figure out where
	to begin the c-obstacle outline algorithm.  We want a vertice
	that is to one side of an edge of object #2.  In particular
	we will begin with edge #1 of object #2, and vertice #2 of
	object #1. -*)
  f1[t_] := Mod[t,l1+1] + Floor[t/(l1+1)];
  f2[t_] := Mod[t,l2+1] + Floor[t/(l2+1)];
  i = 2;  
  j = 1;
  While[ V1[[f1[i]]].V2[[f2[j]]] < 0 
    || -V1[[f1[i-1]]].V2[[f2[j]]] < 0, i++ ]; 
  co = CO[Moveto, co1 , Append[ v2[[f2[j]]] ,1] , f1[i] ];
  b = { Take[CO[Frame, co], 2] };
  i = f1[i];

  Do[ Module[ {},  
    (* Print[{i,j, CO[Frame, co]}]; *)
    co = CO[Moveto, co, Append[ v2[[f2[j]]] , 1] , f1[i]];
    b = Append[b, Take[CO[Frame, co], 2] ];
      (*- find vertice to one side of edge. -*)
    in = i;
    While[ V1[[f1[in]]].V2[[f2[j]]] < 0 
      || -V1[[f1[in-1]]].V2[[f2[j]]] < 0, in++ ]; 
    (* Print[ { V1[[f1[in]]].V2[[f2[j]]] , -V1[[f1[in]]].V2[[f2[j]]] }]; *)
      (*- Check to see how it compares to current vertice -*)
    s = Sign[in-i];
    (* Print[{in,i,s}]; *)
    While[ in != i , Module[ {},
      i+=s; 
     (*  Print[i]; *)
      co = CO[Moveto, co, Append[ v2[[f2[j]]] , 1] , f1[i]];
      b = Append[b, Take[CO[Frame, co], 2] ];
    ] ] ; ], 
  {j , 2, l2+1}]; 
  (* {v1 , v2, h1, h2,l1, l2,V1,V2,i,j, co, *)
  b
]

CO[Boundary, co1_, co2_, min_, max_, step_] := Module[
  {theta, p, b, co} , 

  b = {};
  For[ theta = min , theta <= max , theta+=step ; 
    Module[ {},
    co = CO[SetOrientation, co1, Pi*theta/180];
    p = CO[Obstacle, co, co2];
    p = Table[ Append[p[[i]] , theta ] , {i, 1, Length[p]} ];
    b = Append[b, p]; ] ];
  b
]
    
CO[Slice, co1_, co2_, vertex_] := Module[
  {i, in, j, b, f, com, last, 
   tv, nv, h, tl, nl,		(*- tangent and normal vertices -*)
   ti, ni, T, N, J} ,		(*- tangent and normal vectors  -*)

    (*- Functions to ensure that we don't index incorrectly. -*)
  nl = CO[Length, co1];				(*- #vertices. -*)
  tl = CO[Length, co2];				(*- #vertices. -*)
  ni[t_] := Mod[t,nl+1] + Floor[t/(nl+1)];	(*- index function -*)
  ti[t_] := Mod[t,tl+1] + Floor[t/(tl+1)];	(*- index function -*)

    (*- The vertices of object #1 in the reference plane -*)
  h  = CO[Group, co1][LieGroups`L, CO[Frame, co1]];
  nv = Transpose[h.Transpose[CO[Vertices, co1]]];
  nv = Table[ Take[nv[[i]], 2] , {i , 1, Length[nv]}]; 

    (*- The vertices of object #2 in the reference plane -*)
  h  = CO[Group, co2][LieGroups`L, CO[Frame, co2]];
  tv = Transpose[h.Transpose[CO[Vertices, co2]]];
  tv = Table[ Take[tv[[i]], 2] , {i , 1, Length[tv]}]; 

    (*- Determine the normal vectors along the edges of object #1 -*)
  f[t_] := nv[[ni[t+1]]] - nv[[ni[t]]];
  N = Table[ f[i] , {i , 1 , nl} ];
  J = { {0 , 1} , {-1 , 0} };
  N = Transpose[J.Transpose[N]];

    (*- Determine the tangent vectors to the vertices of object #2 -*)
  f[t_] := tv[[ti[t+1]]] - tv[[ti[t]]];
  T = Table[ f[i] , {i , 1 , tl} ];

    (*- Go through the edges of object #1 to figure out which
	vertices are allowed to make contact with the given vertex
	of object #2. -*)

  last = 0;
  b = {};
  Do[ Module[ {},  
    If[ N[[ni[i]]].T[[ti[vertex]]] >= 0 
	&& -N[[ni[i]]].T[[ti[vertex+tl-1]]] >= 0, Module[ {},
      co = CO[Moveto, co1, Append[ tv[[ti[vertex]]], 1] , ni[i] ];
      b = Append[b, CO[Frame, co]]; 
      last = 1; ] 
    , (*else*) If[ last == 1 , Module[ {}, 
      co = CO[Moveto, co1, Append[ tv[[ti[vertex]]], 1] , ni[i] ];
      b = Append[b, CO[Frame, co]]; 
      last = 0; ] ]
    , (*else*)
      last=0 ];  
    ] ,
  {i, 1, nl}];
  i = nl+1;
  If[ N[[ni[i]]].T[[ti[vertex]]] >= 0 
      && -N[[ni[i]]].T[[ti[vertex+tl-1]]] >= 0, Module[ {},
    co = CO[Moveto, co1, Append[ tv[[ti[vertex]]], 1] , ni[i] ];
    b = Prepend[b, CO[Frame, co]]; 
    last = 1; ] 
  , (*else*) If[ last == 1 , Module[ {}, 
    co = CO[Moveto, co1, Append[ tv[[ti[vertex]]], 1] , ni[i] ];
    b = Prepend[b, CO[Frame, co]]; 
    last = 0; ] ]
  , (*else*)
    last=0 ];  

  b
]

CO[Slice, co1_, co2_, edge_, vertex_] := Module[
  {i, in, j, b, f, com, last, 
   tv, nv, h, tl, nl,		(*- tangent and normal vertices -*)
   ti, ni, T, N, J} ,		(*- tangent and normal vectors  -*)

    (*- Functions to ensure that we don't index incorrectly. -*)
  nl = CO[Length, co1];				(*- #vertices. -*)
  tl = CO[Length, co2];				(*- #vertices. -*)
  ni[t_] := Mod[t,nl+1] + Floor[t/(nl+1)];	(*- index function -*)
  ti[t_] := Mod[t,tl+1] + Floor[t/(tl+1)];	(*- index function -*)

    (*- The vertices of object #1 in the reference plane -*)
  h  = CO[Group, co1][LieGroups`L, CO[Frame, co1]];
  nv = Transpose[h.Transpose[CO[Vertices, co1]]];
  nv = Table[ Take[nv[[i]], 2] , {i , 1, Length[nv]}]; 

    (*- The vertices of object #2 in the reference plane -*)
  h  = CO[Group, co2][LieGroups`L, CO[Frame, co2]];
  tv = Transpose[h.Transpose[CO[Vertices, co2]]];
  tv = Table[ Take[tv[[i]], 2] , {i , 1, Length[tv]}]; 

    (*- Determine the normal vectors along the edges of object #1 -*)
  f[t_] := nv[[ni[t+1]]] - nv[[ni[t]]];
  N = Table[ f[i] , {i , 1 , nl} ];
  J = { {0 , 1} , {-1 , 0} };
  N = Transpose[J.Transpose[N]];

    (*- Determine the tangent vectors to the vertices of object #2 -*)
  f[t_] := tv[[ti[t+1]]] - tv[[ti[t]]];
  T = Table[ f[i] , {i , 1 , tl} ];

    (*- Go through the edges of object #1 to figure out which
	vertices are allowed to make contact with the given vertex
	of object #2. -*)

  last = 0;
  b = {};

  If[ N[[ni[edge]]].T[[ti[vertex]]] >= 0 
      && -N[[ni[edge]]].T[[ti[vertex+tl-1]]] >= 0, Module[ {},
    co = CO[Moveto, co1, Append[ tv[[ti[vertex]]], 1] , ni[edge] ];
    b = Append[b, CO[Frame, co]]; 
    co = CO[Moveto, co1, Append[ tv[[ti[vertex]]], 1] , ni[edge+1] ];
    b = Append[b, CO[Frame, co]]; 
  ] ]; 
  b
]

CO[Patch, co1_, co2_, vertex_, min_, max_, step_] := Module[
  {theta, p, b, co} , 

  b = {};
  For[ theta = min , theta <= max , theta+=step , 
    Module[ {},
    co = CO[SetOrientation, co1, Pi*theta/180];
    p = CO[Slice, co, co2, vertex];
    b = Append[b, p]; 
  ] ];
  b
]

CO[Patch, co1_, co2_, edge_, vertex_, min_, max_, step_] := Module[
  {theta, p, b, co} , 

  b = {};
  For[ theta = min , theta <= max , theta+=step , 
    Module[ {},
    co = CO[SetOrientation, co1, Pi*theta/180];
    p = CO[Slice, co, co2, edge, vertex];
    b = Append[b, p]; 
  ] ];
  b
]

CO[DrawPatch, co1_, co2_, vertex_, min_, max_, step_] := Module[
  {b, gfx, currlen, lastlen, curr, last} ,

  b = CO[Patch, co1, co2, vertex, min, max, step];
  plast = 0;
  gfx = {};
  curr = b[[1]];
  currlen = Length[curr];

  Do[ Module[ {},
    lastlen = currlen;
    last = curr;
    curr = b[[i]];
    currlen = Length[curr];

    If[lastlen == 1 , 
      If[currlen == 1, 
	gfx = Append[ gfx, Line[ Join[last, curr] ] ];
      , (*else*) If[ currlen > 1 ,
	gfx = Append[ gfx, Polygon[ Join[last, Reverse[curr]] ] ];
      ] ]
    , (* else *) If[lastlen > 1, 
      gfx = Append[ gfx, Polygon[ Join[last, Reverse[curr]] ] ];
    ] ]; ]
  , {i, 2, Length[b]} ];

  Show[Graphics3D[gfx], 
    {BoxRatios -> {2 , 2 , 5} , Axes -> True, 
    PlotRange -> {{-15, 15}, {-15, 15}, {min*Pi/180, max*Pi/180}}, 
    ViewPoint -> {2, -2, 3}, Ticks->Automatic,
    FaceGrids->{{0,0,-1}, {-1,0,0}, {0,1,0}},
    AxesLabel -> {"x", "y", "\[Theta]"}}];
]

CO[DrawPatch, co1_, co2_, edge_, vertex_, min_, max_, step_] := Module[
  {b, gfx, currlen, lastlen, curr, last} ,

  b = CO[Patch, co1, co2, edge, vertex, min, max, step];
  last = 0;
  gfx = {};
  curr = b[[1]];
  currlen = Length[curr];

  Do[ Module[ {},
    lastlen = currlen;
    last = curr;
    curr = b[[i]];
    currlen = Length[curr];

    If[lastlen == 1 , 
      If[currlen == 1, 
	gfx = Append[ gfx, Line[ Join[last, curr] ] ];
      , (*else*) If[ currlen > 1 ,
	gfx = Append[ gfx, Polygon[ Join[last, Reverse[curr]] ] ];
      ] ]
    , (* else *) If[lastlen > 1, 
      gfx = Append[ gfx, Polygon[ Join[last, Reverse[curr]] ] ];
    ] ]; ]
  , {i, 2, Length[b]} ];

  Show[Graphics3D[gfx], 
    {BoxRatios -> {2 , 2 , 5} , Axes -> True, 
    PlotRange -> {{-15, 15}, {-15, 15}, {(1-.1)min*Pi/180, (1+.1)max*Pi/180}}, 
    ViewPoint -> {2, -2, 3}, Ticks->Automatic,
    FaceGrids->{{0,0,-1}, {-1,0,0}, {0,1,0}},
    AxesLabel -> {"x", "y", "\[Theta]"}}];
]


EndPackage[];
