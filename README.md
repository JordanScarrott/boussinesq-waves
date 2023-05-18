# boussinesq-waves
This simulation solves the Boussinesq equation for beach waves using the Potential Flow model described by Wei &amp; Kirby 1995.

# Running the simulation
`A = Boussinesq([200, 0.01, 0.045, 0.45, 0.05, 0.05, 0.05, 10, 10, FloorProfile.FLAT, InitialCondition.EXPONENTIAL]);`
`A = A.solve();`
`A = A.displayMeshes();`

# Other Interesting parameter sets:
Here are a few interesting sets of parameters I've found when playing around with the program.

## Water down a channel:
I found that you can actually increase dx and dy quite drastically while keeping the time set the same as above.
This way I was able to adjust the length of the channel to allow waves to propagate up a slope for longer.

`A = Boussinesq([2000, 0.01, 0.045, 0.45, 0.5, 0.5, 0.05, 300, 10, FloorProfile.FLAT, InitialCondition.SECH]);`
InitialCondition.SECH sets up a plane wave initial distribution that falls and creates two oppositely directed child waves.
I found that the SINGLE_BAR floor profile wasnt working so I adapted the FLAT one to actually be a slope.
I can get easily get 2000 iterations in a short amount of time.
It seems like the boundary errors are greatly reduced when dx and dy are larger than the original 0.05.
