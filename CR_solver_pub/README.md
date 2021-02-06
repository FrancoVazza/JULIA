
Public serial version of the cosmic ray solver (e.g. Fokker Planck without diffusion terms) used in Vazza, Wittor, Brunetti & Bruggen 2021.

The code requires as input a sequence of tracers data, written in HDF5 format (see /tracers folder), and computes particle spectra under radiative and coulomb losses, adiabatic compression/rarefaction, and injection+re-acceleration by shock waves. 

- The input tracer files contain particle ID information and gas physical quantities like gas density [g/cm3], gas temeperature[K], magnetic fields strength [microGauss], 
3D vorticity and divergence [1/s] and redshift). A file sample of 25 tracers is given in this the /tracers folder

- The main code is FP_public.jl

- Additional necessary functions are: 

      param_spectra.jl  > containing  parameters of input spectra and momentum binning, which can be changed here (pmin, pmax, dp).

      loss_gain1.jl   > containing all acceleration and loss terms, necessary for the Fokker-Planck evolution (Chang & Cooper 1970 - like). 

- FP_public.jl generates the simultaneous evolution of the same set of electron spectra, over 50 time snapshots, for three different scenarios of cosmic ray acceleration and losses.

- Details on the method and on its application to cosmological simulations are described in Vazza, Wittor, Brunetti & Bruggen, 2021 A&A (submitted, as of this writing). 
