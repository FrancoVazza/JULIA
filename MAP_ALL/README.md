Collection of functions for efficient map making of large (up to 2400^3) datasets from ENZO cosmological simulations.

This is a (work in progress) sort of universal macro that does multiple map making in Julia of ENZO simulations, for the moment:

- detect shocks on the grid;
- projected field quantities;
- synchrotron radio emission from shocks;
- HI emission (it can be still speeded up, probably);
- X-ray emission from an bAPEC model (externally generated with IDL); 
(-additionally, also the filtering of local velocity field to get a turbulence estimate is possible)

It can be executed 
- in serial, by cutting the volume into multiple HDF5 blocks;
- in parallel, by distributing the map making to different cores;
- in mixed parallel way (suitable for very large runs), by distributing tasks to different cores and allowing each of them to cut the volume in multiple HDF5 blocks.

PABLO_opt.jl is the workhorse, and it assumes a parallel distribution of data via DistributedArrays (require Julia v0.6), from 1 to N procs. 
