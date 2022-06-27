### Simple crystalline silicon example.

1. Run the lammps input script `in.run`.

2. This will dump velocities into a file `dump.velocities`, which is read by the `postproc` program.

3. Do `mpirun -np 4 postproc powerspec 2000 640 0.0025` to calculate the density of states. 

The arguments of `postproc` are explained below.

`powerspec` - this option tells us to calculate the power spectrum (i.e. density of states).

`2000` - this is the number of timesteps (200000) divided by the frequency of dumping (100).

`640` - number of atoms

`0.0025` - sampling interval in ps, which is determined by timestep (0.5 fs = 0.0005 ps) times 5 because we dumped every 5 timesteps.


