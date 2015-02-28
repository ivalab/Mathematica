(*

  matrices.m - allows one to specify certain function with matrices,
  such as equation solving.

*)

BeginPackage["matrices`"];
<<LinearAlgebra`MatrixManipulation`;

(* 
  Solve second order linear matrix ODE of the form
  My'' + Ny' + Oy = F
*)
MLODE2 [M_ , N_ , O_ , F_ , int_] := Module [
  {ANS, Y,Yp,Yeqns,Minv,MN,MO,x,ROWS,COLS} , 
  Minv = Inverse[M];
  Yp = Table[{y[i]'[int[[1]]]} , {i,Length[M]}];  
  Y  = Table[{y[i][int[[1]]]}  , {i,Length[M]}];
  MN = Minv N;
  MO = Minv O;
  Yeqns = Join [ Table [ y[i]''[int[[1]]] + (MN[[i]] . Yp )[[1]]
	 + (MO[[i]] . Y)[[1]] == (Minv[[i]] . F)[[1]]  
	 , {i,Length[M]}],
         Table [ y[i][0] == 0 , {i,Length[M]} ],
         Table [ y[i]'[0] == 0 , {i,Length[M]} ]
       ];
  ANS = NDSolve [ Yeqns , Table[y[i] , {i, Length[M]}],int]
] 

(* 
  Solve first order linear matrix ODE of the form
  My' + Ny = F
*)
MLODE1 [M_ , N_ , F_ , int_] := Module [
  {ANS, Y,Yeqns,Minv,MN,MO,x,ROWS,COLS} , 
  ANS = 0;
  If [ SquareMatrixQ[M] && Det[M] != 0 , 
    Y  = Table[{y[i][int[[1]]]}  , {i,Length[M]}];
    Minv = Inverse[N];
    MO = Minv O;
    Yeqns = Join [ Table [ y[i]'[int[[1]]] 
            + (MO[[i]] . Y)[[1]] == (Minv[[i]] . F)[[1]]  
	    , {i,Length[M]}],
            Table [ y[i][0] == 0 , {i,Length[M]} ]
           ];
    ANS = NDSolve [ Yeqns , Table[y[i] , {i, Length[M]}],int]
  ];
  ANS
]


EndPackage[];
