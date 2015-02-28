BeginPackage["Anim`"];

fMakeGrid::usage="";
fMakeGridWindow::usage="";
fGetTicks::usage="fGetTicks[point_, bounds_, prange_, step_, tolerance_]";
fMakeMovie::usage="fMakeMove[basedir, basename, thegrafx]";
fMakeMovie2::usage="";

Begin["`Private`"];

fMakeGrid[xbounds_, ybounds_, step_] :=  Module[{i, cx, cy, gfx},
  gfx = {};
  For[i = xbounds[[1]] , i <= xbounds[[2]], i+=step,
    gfx = Append[gfx, Line[{{i, ybounds[[1]]},{i,ybounds[[2]]}}]] ];
  For[i = ybounds[[1]] , i <= ybounds[[2]], i+=step,
    gfx = Append[gfx, Line[{{xbounds[[1]], i},{xbounds[[2]],i}}]] ];
  gfx
];

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
