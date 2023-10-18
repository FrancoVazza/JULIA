Collection of Julia routines to produce final datasets for the SKA's DC.
This routines do not allow users to actually see the underlying physical models and correct answers of the challange. 

* [RM_peaks.jl](https://github.com/FrancoVazza/JULIA/SKA_DC/RM_peaks.jl) - locates the maximum of the RM along the LOS in a simulated box (100^3 Mpc^3 volume).

* [interpolate_frequencies.jl](https://github.com/FrancoVazza/JULIA/SKA_DC/interpolate_frequencies.jl) - produces N additional IQU maps by interpolating 5 fundamental frequencies.

* [RM_redshift.jl](https://github.com/FrancoVazza/JULIA/SKA_DC/RM_redshift.jl) - produces RM in a long redshift interval, by replicating and redshift-correcting a set of input boxes (85^3 Mpc^3 volume).
