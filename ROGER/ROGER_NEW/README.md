Updated version of the ROGER code to simulate the evolution of relativistic electron spectra carried by tracers, solving the time-dependent diffusion-loss equation of relativistic electrons with a parallel Fokker-Planck solver. 

The solver includes:

* adiabatic changes;
* collisional and ionisation losses;
* synchrotron and inverse Compton losses;
* shock acceleration assuming diffusive shock acceleration from the thermal pool;
* shock re-acceleration from diffusive shock acceleration;
* Fermi II turbulent re-acceleration based on the adiabatic stochastic acceleration model (ASA).
    

Used already for production in Beduzzi et al. 2024 (https://arxiv.org/pdf/2406.09859); Vazza et al. 2025 (https://arxiv.org/abs/2501.19041)

Example of the evolving map of gas density and of radio emission at 150 MHz produced by ROGER:

<img src="E1_dens_radio1.gif" alt="alt text" width="whatever" height="whatever">

Example of the evolving momentun spectra of relativistic electrons (for three different models simulated at the same time by ROGER) for a full tracer simulation.

<img src="E62_spectra.gif" alt="alt text" width="whatever" height="whatever">
