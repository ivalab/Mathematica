(*============================== Animation =============================*)
(*

  This package implements an object that keeps track of the
  equations of motion for a dynamical system.  These are typically
  tied in to a particular manifold object and do not necessarily
  function alone.


  Author:       Patricio Vela
  Created:      2002/06/07
  Modified:     2015/03/01

  Notes: 
    - At some point, the different objects need to be separated.

*)
(*============================== Animation =============================*)
BeginPackage["Anim`"];

(*================================ Usage ===============================*)
fMakeGrid::usage="";
fMakeGridWindow::usage="";
fGetTicks::usage="fGetTicks[point_, bounds_, prange_, step_, tolerance_]";
fMakeMovie::usage="fMakeMove[basedir, basename, thegrafx]";
fMakeMovie2::usage="";

(*========================== Package Functions =========================*)
Begin["`Private`"];

(*----------------------------- fMakeGrid ----------------------------*)
(*
    Creates a grid for insertion into a graphic object.
*)
fMakeGrid[xbounds_, ybounds_, step_] :=  Module[{i, cx, cy, gfx},
  gfx = {};
  For[i = xbounds[[1]] , i <= xbounds[[2]], i+=step,
    gfx = Append[gfx, Line[{{i, ybounds[[1]]},{i,ybounds[[2]]}}]] ];
  For[i = ybounds[[1]] , i <= ybounds[[2]], i+=step,
    gfx = Append[gfx, Line[{{xbounds[[1]], i},{xbounds[[2]],i}}]] ];
  gfx
];

(*-------------------------- fMakeGridWindow -------------------------*)
(*
    Creates a grid for insertion into a graphic object, but about a
    point.  The grid is decimated about that point.  Basically creates
    a roving window around a moving object that keeps a consistent
    grid for the object as it moves.  For a moving object centered at 
    the point, the grid will displace therefore providing the real
    illusion of motion.
*)
fMakeGridWindow[point_, xbounds_, ybounds_, step_] := Module[{i, cx, cy, gfx},
  gfx = {};
  For[i = point[[1]] + xbounds[[1]] , i <= point[[1]] + xbounds[[2]], i+=step,
    gfx = Append[gfx, Line[{{i, point[[2]] + ybounds[[1]]},
      {i,point[[2]] + ybounds[[2]]}}]] ];
  For[i = point[[2]] + ybounds[[1]] , i <= point[[2]] +ybounds[[2]], i+=step,
    gfx = Append[gfx, Line[{{point[[1]] + xbounds[[1]], i},
      {point[[1]] + xbounds[[2]],i}}]] ];
  gfx
];

(*----------------------------- fGetTicks ----------------------------*)
fGetTicks[point_, bounds_, prange_, step_, tolerance_] := Module[
  {i, xticks, yticks, xmin, xmax, ymin, ymax},
  
  xticks = {};
  yticks = {};
  xmin = point[[1]] + prange[[1,1]] + tolerance;
  xmax = point[[1]] + prange[[1,2]] - tolerance;
  ymin = point[[2]] + prange[[2,1]] + tolerance;
  ymax = point[[2]] + prange[[2,2]] - tolerance;

  For[i = bounds[[1,1]] , i <= bounds[[1,2]], 
    i+=step,
    If[ (i >= xmin) && (i <= xmax) ,  xticks = Append[xticks, i] ] ];
  For[i = bounds[[2,1]] , i <= bounds[[2,2]], 
    i+=step,
    If[ (i >= ymin) && (i <= ymax) ,  yticks = Append[yticks, i] ] ];
 Chop[{xticks, yticks, {}, {}}]
];

(*---------------------------- fMakeMovie ----------------------------*)
(*
    Long time ago, Mathematica did not have the ability to export a
    movie.  This output images, then refactored them into a movie
    using ffmpeg (or maybe mplayer) within a system executable script.

    If Mathematica can now output video, then this is unnecessary.
*)
fMakeMovie[basedir_, basename_, thegrafx_] := Module[
  {fnames, jnames, mname, mlen},
  mlen = Length[thegrafx];
  fnames = Table[ basedir <> "/" <> basename <> "." <> 
    ToString[PaddedForm[i, 4, NumberPadding -> {"0","0"},
    NumberSigns -> {"", ""}]] <> ".eps",  
    {i,mlen}];
  jnames = Table[ basedir <> "/" <> basename <> "." <> 
    ToString[PaddedForm[i, 4, NumberPadding -> {"0","0"},
    NumberSigns -> {"", ""}]] <> ".jpg",  
    {i,mlen}];
  For[i=1, i <= mlen, i++, 
    Export[ fnames[[i]], thegrafx[[i]] ];
    Run["convert " <> fnames[[i]] <> " " <> jnames[[i]] ] ];
  Run["cd " <> basedir <> "; ~/mathematica/libs/makempeg -fs 1 -fe " <> ToString[mlen] <> 
    " -fi 1 -base " <> basename <> " -ext jpg"];
];

End[];

EndPackage[];
