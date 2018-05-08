# JULIA
(Growing) collection of programs coded in Julia, for postprocessing of grid (mostly Enzo) data

Updated list

1. map_par.jl   =  map making of cubic grid data with uniform resolution, written in hdf5. Works in parallel using Distributed Arrays. Tested up to 64 cores on Piz Daint, and with grids up to 1024^3. Necessary packages: Distributed Arrays, HDF5, FITSIO 
 
2. map_serial.jl = simple map making tool for ENZO simulations (using hdf5 datasets) and producing maps with Winston. 

3. miro_par4.jl   = parallel version of the MIRO code (here the IDL version https://github.com/FrancoVazza/MIRO). It feature the generation of random 3D magnetic field with components obeying an input power-law spectrum and keeping div(B)=0 through the use of vector potential. Can work in parallel using DistributedArrays. Necessary packages: FFTViews and DistributedArrays
