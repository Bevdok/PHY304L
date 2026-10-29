One way in which we can compare the data between our pendulum simulation and an analytical model is by using the small angle appoximation on 
our analytical model and compare it to our built computational pendulum solver. We could then do a simple percent error caluclation between 
the compuational and analytical results to compare, with that percent error being our figure of merit.
Percent error can be calculated using the formula: Percent Error = |(Analytical Value - Computational Value) / Analytical Value| * 100

As the initial angle of the pendulum increases, we expect the percent error to also increase due to the limitations of the small angle approximation. 
The small angle approximation assumes that sin(\[Phi]) ≈ \[Phi] when \[Phi] is measured in radians and is small.
So in general we can expect our analytical results to deviate more from the computational results as the initial angle increases, 
leading to a higher percent error.

g = 9.81;
l = 0.1;
Omega = Sqrt[g/l];
Tanalytical = 2*Pi*Sqrt[l/g];
angles = {5, 10, 15, 30, 60};
Tcomputation = Table[
  sol = NDSolve[
    {Phi''[t] == -Omega^2*Sin[Phi[t]],
     Phi[0] == angles[[i]]*Pi/180,
     Phi'[0] == 0},
    Phi, {t, 0, 2}
  ];
  FindRoot[
    (Phi[t] /. sol[[1]]) == angles[[i]]*Pi/180,
    {t, 0.6}
  ][[1, 2]],
  {i, 1, Length[angles]}
];
PercentError = 100*(Tcomputation - Tanalytical)/Tanalytical;
ListPlot[
  Transpose[{angles, PercentError}],
  Joined -> True,
  PlotMarkers -> Automatic,
  PlotRange -> All,
  AxesLabel -> {"Initial angle (degrees)", "Period Error (%)"},
  PlotLabel -> "Small-Angle Approximation Error",
  GridLines -> Automatic,
  PlotStyle -> Red
]
We can see that at 60 degrees, our percent error for the small angle approximation is slightly above 7%, so we can mark this angle as the point of significant deviation.

An increase in the mass of the pendulum bob will not change our systems motion as we can see through a Newtons 2nd law equation, that the mass is canceled out. If we set the Force = -mgsin(Phi), which we can turn into ma = -mgsin(Phi). with mass on both sides, we can cancel it out, proving that our systems motion does not depend on the mass on the pendulum bob.