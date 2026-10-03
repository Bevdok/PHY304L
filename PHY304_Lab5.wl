(* ::Package:: *)

Plot[y=x , {x, 0, 5}]


w[t] = F[t] = \[Alpha]*Log[\[Alpha]*t]*Cos[\[Omega]*t]*Sin[\[Omega]*t]


D[w[t], t]


w[x] = f[x]*G[x]*h[x]


D[w[x], x]



g = 9.81;
l = 0.1;
Omega = Sqrt[g/l];

sol = NDSolve[{Phi''[t] == -Omega^2*Sin[Phi[t]], Phi[0] == Pi/6, Phi'[0] == 0}, Phi, {t, 0, 2}];

Plot[Evaluate[Phi[t] /. sol], {t, 0, 2}, AxesLabel -> {"t (s)", "\[Phi](t)"}, PlotLabel -> "Pendulum Motion"]



Graphics3D[{Cylinder[{{0, 0, 0}, {0, 0, 5}}, .1], 
Cylinder[{{0, 0, 0}, {0, 0, 0.2}}, 1]},
 Axes -> True,
 AxesLabel -> {"x", "y", "z"},
 PlotRange -> {{-2, 2}, {-2, 2}, {0, 5}}]
