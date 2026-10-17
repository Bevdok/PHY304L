m = 1;
k = 4;
x0 = 1;
v0 = 0;
\[Omega] = Sqrt[k/m];
\[Gamma] = b/m;

#Underdampedcasebequals2

SolUnder = NDSolve[{x''[t] +  2*x'[t] + \[Omega]^2 x[t] == 0, x[0] == x0, x'[0] == v0}, x, {t, 0, 10}]

#Criticallydampedcasebquals4

SolCrit = NDSolve[{x''[t] + 4*x'[t] + \[Omega]^2 x[t] == 0, x[0] == x0, x'[0] == v0}, x, {t, 0, 10}]

#Overdampedcasebequals6

SolOver = NDSolve[{x''[t] + 6*x'[t] + \[Omega]^2 x[t] == 0, x[0] == x0, x'[0] == v0}, x, {t, 0, 10}]


Plot[{Evaluate[x[t] /. SolUnder], Evaluate[x[t] /. SolCrit], Evaluate[x[t] /. SolOver]}, {t, 0, 10}, PlotLegends -> {"Underdamped", "Critically Damped", "Overdamped"}, PlotStyle -> {Red, Blue, Green}, AxesLabel -> {"Time (s)", "Displacement (m)"}, PlotLabel -> "Damped Harmonic Oscillator Cases"]

m = 1;
k = 4;
\[Gamma] = 0.5;
F0 = 1;
\[Omega] = Sqrt[k/m];
\[Omega]d = 1.5;

x0 = 2;
v0 = 0;

SolA = NDSolve[{x''[t] + \[Gamma] x'[t] + \[Omega]^2 x[t] == F0 Cos[\[Omega]d t], x[0] == x0, x'[0] == v0}, x, {t, 0, 30}];

Plot[Evaluate[x[t] /. SolA], {t, 0, 30}, PlotStyle -> {Red}, AxesLabel -> {"Time (s)", "Displacement (m)"}, PlotLabel -> "Driven Damped Harmonic Oscillator", Epilog -> {Red, Dashed, Thick, Line[{{15, -5}, {15, 5}}]}]

m = 1;
k = 4;
\[Gamma] = 0.5;
F0 = 1;
\[Omega] = Sqrt[k/m];

\[Omega]0 = 1;
\[Alpha] = 0.05;
Tmax = 50;

x0 = 0;
v0 = 0;

SolB = NDSolve[{x''[t] + \[Gamma] x'[t] + \[Omega]^2 x[t] == F0 Cos[\[Omega]0*t + 0.5*\[Alpha]*t^2], x[0] == x0, x'[0] == v0}, x, {t, 0, Tmax}];

Plot[Evaluate[x[t] /. SolB], {t, 0, Tmax}, PlotStyle -> {Red}, AxesLabel -> {"Time (s)", "Displacement (m)"}, PlotLabel -> "Driven Damped Harmonic Oscillator with Chirped Driving Force"]

m = 1;
F0 = 1;
x0 = 0;
v0 = 1;

\[Tau] = (2*m*v0)/F0;

x[t] := Module[{t = Mod[t, 2*\[Tau]]}, Piecewise[{{v0*t -(F0/(2*m))*t^2, 0 <= t <= \[Tau]}, {-v0*(t - \[Tau]) + (F0/(2*m))*(t - \[Tau])^2, \[Tau] <= t <= 2*\[Tau]}}]];


Plot[Evaluate[x[t]], {t, 0, 2*\[Tau]}, PlotStyle -> {Red}, AxesLabel -> {"Time (s)", "Displacement (m)"}, PlotLabel -> "Piecewise Motion of a Particle under Constant Force"]

As the negative force is applied, the displacement reaches its peak positive displacement, before a restorative force, the spring, brings it back to its original displacement, where the positive force is applied, and the displacement reaches its peak negative displacement and reverses back to its original displacement. A scenario where we see similair motion is when a ball is rolling up the two sides of a V-shaped well, it reaches its max displacement on both sides, slows down and "restores" its position back to its starting displacment, the bottom of the well.