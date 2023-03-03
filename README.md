# boussinesq-waves
This simulation solves the Boussinesq equation for beach waves using the Potential Flow model described by Wei &amp; Kirby 1995.

# Running the simulation
`A = Boussinesq([200, 0.01, 0.045, 0.45, 0.05, 0.05, 0.05, 10, 10, FloorProfile.FLAT, InitialCondition.EXPONENTIAL]);`
`A = A.solve();`
`A = A.displayMeshes();`
