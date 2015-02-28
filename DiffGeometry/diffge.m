BeginPackage["DifferentialGeometry`"];

Lie::usage=
  "Lie[X1, X2, q] \n
   Compute the Lie bracket of X1 and X2 with respect to variables q in R^n";

SymProd::usage=
  "SymProd[X1, X2, S, q] \n
   Compute the Symmetric product of X1 and X2 with respect to variables q in R^n";

Begin["`Private`"];

Lie[X_, Y_, q_] := Outer[D, Y, q].X - Outer[D, X, q].Y;

SymProd[Y1_, Y2_, S_, q_] := Lie[Y2, Lie[S, Y1, q], q];

End[];

EndPackage[];
