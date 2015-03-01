(*========================== Averaging Theory ==========================*)
(*

  This package implements some of the basic equations associated to
  averaging theory.  At least first order is done, but possible second
  or higher.  


  Author:       Patricio A. Vela
  Created:      2002/06/07 
  Modified:     2015/03/01

  NOTES: 
    - At some point, the different objects need to be separated.
    - Definitely the code needs some cleaning up and documentation.

*)
(*========================== Averaging Theory ==========================*)
BeginPackage["Averaging`"];

fAverage::usage=
  "fAverage[f,{t, T}] \n\nFinds the average of the function over the time interval (0,T) with respect to the variable t. \n Of course f should be a T-periodic function";

fNAverage::usage=
  "fNAverage[f,{t, T}] \n\nFinds the average of the function over the
  time interval (0,T) with respect to the variable t. \n Of course f
  should be a T-periodic function. Uses numerical integration";

fAvgPotential::usage=
  "fAvgProduct[U, {t,T}] \n\nFind the averaged potential matrix corresponding to the oscillatory input vector U with period T and with respect to variable t.";

fAvgProduct::usage=
  "fAvgProduct[U, indices, powers, {t, T}] \nfAvgProduct[U, pow_, {t, T}] \n\nFinds the...";

Begin["`Private`"];

fAverage[ f_ , spec_ ] := (1/spec[[2]]) Integrate[ f , {spec[[1]], 0, spec[[2]]} ];

fNAverage[ f_ , t_ ,lim_ ] := 
  (1/lim) NIntegrate[ f , {t, -Infinity, lim}, Method->Oscillatory ];

fAvgPotential[U_ , spec_ ] := Module[
  {UU, V, VM, t, T} , 
  t = spec[[1]];
  T = spec[[2]];
  UU[tt_] := Integrate[ U /. {t -> tt}, {tt, 0, t} ];
  V[tt_] := Outer[Times, UU[tt],  UU[tt]];
  VM = (1/2) fAverage[ V[t] , {t,T} ] ;
  VM
];

IteratedIntegral[f_, int_, n_] := Module[
  {x} , 
  If[ n < 1 , 
    f 
  ,
    If[ n == 1 ,
      Integrate[f, int]
    ,
      If[Length[int] > 0,
        Integrate[IteratedIntegral[ f /. {int[[1]] -> x[n]} , {x[n],
          int[[2]] , int[[1]]} , n-1] , int],
       Integrate[IteratedIntegral[f /. {int -> x[n]} , x[n], n-1] /.
       {x[n] -> int} , int ] 
      ]
    ]
  ]
];

fAvgProduct[U_, ind_, pow_, int_] := 
  Product[ IteratedIntegral[ U[[ind[[i]]]] , int , pow[[i]] ] , 
    {i, Length[ind]}];


fAvgProduct[U_, pow_, int_] := Module[
  {i, j, str, istr, Uc}, 
  
  j = Length[pow];
  str =  ToString[ Uc ] <> "," <> ToString[ Array[x,j] ];
  str = str <> ", " <> ToString[ Unevaluated[pow]]; 
  str = str <> ", " <> ToString[ Unevaluated[int]] <> "]";
  istr = "Table[ fAvgProduct[" <> str;
  str = Table[ ToString[{x[i] , 1, Length[U]}] , {i, j} ];
  For[i = 1, i <= Length[str] , i++ ,
    istr = istr <> "," <> str[[i]] ];
  istr = istr <> "]";
  Uc = U;
  ToExpression[istr]
];


End[];

EndPackage[];

