#...parameters of the simulation

@everywhere run_model = ["E15B", "E62", "E21", "E1", "E11", "E7", "E14", "E16A", "E18A", "E18B"]  #....list of runs to analyse (just one for the moment )

@everywhere const dx = 40.69    #resolution of the eulerian simulation at this level 
@everywhere xc = 1 / (dx)
@everywhere snap_in0 = [35, 35, 34, 25, 32, 25, 35, 51, 47, 40]  #..initial snapshot
@everywhere ntr00 = [262440, 249028, 155261, 158595, 277240, 199856, 255474, 350583, 312260, 206980] #...maximum possible number of tracers in each run (can be reduced by selecting based on density)
@everywhere snap_fin0 = [199, 214, 206, 198, 237, 197, 202, 283, 274, 236]#  
@everywhere dmin = [1e-29, 1e-29, 1e-29, 1e-29, 1e-29, 1e-29, 1e-29, 1e-29, 3e-29, 1e-29]   #...low density threshold 

#...cosmological parameters
@everywhere using Cosmology
@everywhere cOmegaM = 0.258
@everywhere ch = 0.72
@everywhere cosmo = cosmology(OmegaM=cOmegaM, h=ch)


#...parameters for spectra
@everywhere using SpecialFunctions
@everywhere const p_max = log10(1.1e6)
@everywhere const p_min = log10(0.1)
@everywhere const dp = 0.1  # log spacing of momenta in the spectra
@everywhere const part = 2  #...1=proton  2=electron
@everywhere const np = 1 + (floor(Int64, (p_max - p_min) / dp))
@everywhere const pend = np
@everywhere pval = collect(p_min:dp:p_max)
# parameters for shock modelling
@everywhere const norm = 30 #to have all in 1e30 erg or 1e30 erg/s
@everywhere const lcorr = 20.0 #...in kpc
@everywhere const mthr = 2.0  #...minimum Mach number for shock re-acceleration (i.e. no injection of fresh new electrons)
@everywhere const minj = 2.5  #...minimum Mach number for the shock injection (new fresh electrons are added to spectra)
#..other parameters
@everywhere const ηB = 0.05 #...was 0.035   #...conversion efficiency of turbulent energy into magnetic energy 
@everywhere const scale_turbo = 2 * dx   #....stencil used to compute vorticity
@everywhere const normcr = 0.000696

#....radio frequencies used to compute synchrotron emission
@everywhere freqf = [5e7, 1.5e8, 6.5e8, 1.4e9, 5e9]
@everywhere nfreq = size(freqf)
