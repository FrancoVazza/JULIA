#....packages
@everywhere using LaTeXStrings
@everywhere using Plots
@everywhere Plots.PyPlotBackend()
@everywhere using FITSIO
@everywhere using Statistics
@everywhere using Distributed
@everywhere using SharedArrays
@everywhere using SpecialFunctions
@everywhere using DelimitedFiles
@everywhere using HDF5
@everywhere using Unitful

@everywhere main = "./"  #..main folder containing ROGER functions
@everywhere include(string(main, "/constants.jl"))   #...constants 
@everywhere include(string(main, "/param.jl"))       #...input simulation parameters & cosmology
@everywhere include(string(main, "/loss_gain_par_log.jl"))    #...routines for loss and gain processes
@everywhere include(string(main, "/synchrotron_par_log3.jl"))   #...synchrotron modules 
@everywhere include(string(main, "/plot_maps.jl"))              #...routines for map making 

#...custom label to add to outputs of ROGER
@everywhere ftag = "_newTz1.5"
@everywhere test = "_newS"
@everywhere rrr = 0

for ro in 1:1        #....loops over different simulations (list of runs to be set in param.jl) 
   @everywhere rrr += 1
   @everywhere ntr0 = ntr00[rrr]
   @everywhere snap_in = snap_in0[rrr]
   @everywhere snap_fin = snap_fin0[rrr]
   @everywhere initial_snap = snap_in
   @everywhere runm = run_model[rrr]
   @everywhere rootf = string("/Users/francovazza/Desktop/data/DATA/MH_LEO/", runm, "/trac/")           #..folder containing the hdf5 files of tracers
   @everywhere root_map = string("/Users/francovazza/Desktop/data/DATA/MH_LEO/", runm)
   @everywhere root_out = string("/Users/francovazza/Desktop/data/DATA/MH_LEO/", runm, "/trac/roger/")  #..folder collecting ROGER outputs 

   @everywhere n = 384            #grid size of input simulation 
   println("doing model ", runm)
   @everywhere ntr0 = ntr00[rrr]
   @everywhere snap_in = snap_in0[rrr]
   @everywhere snap_fin = snap_fin0[rrr]

   @everywhere Dmin = dmin[rrr]  #...minimum physical density to start the injection of tracers (set in param.jl)


   println("° 1st SERVE : selecting tracers°°°°°°°°°°°°°°°°°°°")
   #.....this first part is to select which tracers, out of the total contained in HDF5 files, are going to be
   #.....initialised, and with which initial spectrum
   @everywhere lag = 1
   @everywhere snapn = string("0", snap_in)
   @everywhere snapn2 = string("0", snap_in + lag)
   if snap_in < 100
      @everywhere snapn = string("0", snap_in)
   end
   if snap_in < 10
      @everywhere snapn = string("00", snap_in)
   end
   @everywhere snapn2 = string(snap_in + lag)
   if snap_in + lag < 100
      @everywhere snapn2 = string("0", snap_in + lag)
   end
   if snap_in + lag < 10
      @everywhere snapn2 = string("00", snap_in + lag)
   end

   @everywhere file_trac = string(rootf, "_tracers_", snapn, ftag, ".hdf5")   #file with Lagrangian tracers information


   #.....read conversion factors from numerical to cgs and redshift - otherwise the latter should be provided in a different way 
   @everywhere cdd = h5read(file_trac, "DensityConversionFactor")
   @everywhere z = h5read(file_trac, "Redshift")
   @everywhere cv = h5read(file_trac, "VelocityConversionFactor")
   cb = sqrt(cdd * 4.0 * pi) * cv * (1 + z)^2  #b-field in proper Gauss   
   println("redshift=", z)
   @everywhere zz = z

   px = h5read(file_trac, "xcoord")
   pn = size(px)
   np1 = pn[1]
   ntr0 = np1

   tr = Array{Float64}(undef, 16, ntr0)   #...main tracer array with all fields
   tr .= 0.0

   @views tr[1, 1:ntr0] = px[1:ntr0]
   py = h5read(file_trac, "ycoord")
   @views tr[2, 1:ntr0] = py[1:ntr0]
   pz = h5read(file_trac, "zcoord")
   @views tr[3, 1:ntr0] = pz[1:ntr0]
   pd = h5read(file_trac, "density")
   @views tr[7, 1:ntr0] = pd[1:ntr0]


   normcr = 1.0

   pv = h5read(file_trac, "volume")
   @views tr[15, 1:ntr0] = pv[1:ntr0]
   #....selection based on density of tracers, only tracers with rho>Dmin are selected
   sel = findall(x -> (x >= Dmin), pd[1:ntr0])
   ns = size(sel)
   ntr = ns[1]
   @everywhere lsize = 1.0
   println("going to evolve ", ntr, " tracers")

   println("° °2nd SERVE : generating starting spectra°°°°°°°°°°°°°°")
   #.....array for momentum and radio spectra are created only for the selected particles
   #.....ROGER evolves 3 model of spectra at the same time
   #.....each model evoles the exact same tracers, under different scenarios for losses and gains (to be specified by the user)
   global pe = fill(1e-40, np, ntr)
   global pradio1 = fill(1e0, ntr, nfreq[1])
   global pe2 = fill(1e-40, np, ntr)
   global pradio2 = fill(1e0, ntr, nfreq[1])
   global pe3 = fill(1e-40, np, ntr)
   global pradio3 = fill(1e0, ntr, nfreq[1])



   #...CIOFF model
   #...initialisation using the CIOFF model - assuming particles are initially injected by radiogalaxies (-> Vazza et al. 2022 A&A for details)
   #...replace with different spectral initialisations for different problems
   #...avoid initialising particles with 0 values as this produces NaN.
   @everywhere φ = 1e-4   #.....initial density of relativistic electrons >pmin wtr thermal electron density. Be careful - have a physical model to decide what this should be, depending on pmin 
   @everywhere mass_gas = msol * 5e8 #....gas mass sampled by each tracer. the CR electron mass is a tiny fraction, φ, of this!
   @everywhere mass_trac = mass_gas * φ
   @everywhere t_on = 1e6 * yr  #...t_on [s]   parameters for the CIOFF injection model - see Harwood et al. 
   @everywhere t_off = 1e4 * yr #...t_off [s]
   @everywhere t_tot = t_on + t_off
   @everywhere bjet = 1e-6    #..B-field inside lobes at injection [in G]
   @everywhere β = (b_syn * bjet^2.0 + (1.0 + zz)^4 * b_IC)
   @everywhere γc1 = 1 / (β * (t_tot))
   @everywhere γc2 = 1 / (β * (t_off))
   global in_spec = 2.2      #input spectrum 
   global pmin_jet = 10^pval[1]   #...minimum electron momentum in the initial spectra 
   @inbounds for g in 1:np
      ppv = 10^pval[g]
      ppmax = 10^pval[np-1]
      ga = sqrt(1 + (ppv)^2)
      gamma_max = sqrt(1 + (ppmax)^2)
      @inbounds @simd for i in 1:ntr
         mass0 = φ * pd[i] * pv[i] * (dx * kpctocm)^3  #mass using comoving volume and comoving density 
         massCRe = mass0

         peg = 0
         if ga <= γc1
            peg = (1 - β * ga * t_off)^(in_spec - 1) - (1 - β * ga * (t_tot))^(in_spec - 1)
         end
         if ga >= γc1 && ga <= γc2
            peg = (1 - β * ga * t_off)^(in_spec - 1)
         end
         if ga > gamma_max
            peg = 0.0
         end

         #...possible way of initialising three different populations , which will be evolved in parallel 
         pe[g, i] = (mass0 * pmin_jet^(in_spec - 1) / (0.6 * mp)) * (ppv)^(-in_spec - 1) * peg
         pe2[g, i] = (massCRe * pmin_jet^(in_spec - 1) / (0.6 * mp)) * (ppv)^(-in_spec - 1) * peg
         pe3[g, i] = pe[g, i]
         pe[g, i] *= 1e-20
      end
   end

   #....plotting of the initial input spectrum (before losses)
   pe_tot = fill(0.0, np)
   pe_tot2 = fill(0.0, np)
   pe_tot3 = fill(0.0, np)
   @inbounds @simd for i in 1:ntr
      @inbounds @simd for g in 1:np
         pe_tot[g] += (pe[g, i] * 10^pval[g])
         pe_tot2[g] += (pe2[g, i] * 10^pval[g])
         pe_tot3[g] += (pe3[g, i] * 10^pval[g])
      end
   end

   println("plotting injection spectra")
   a = plot_spectra_ini(pe_tot, pe_tot2, pe_tot3, pval)

   @everywhere maxp = np - 1
   @everywhere so = 0


   println("° ° °THE EXCHANGE BEGINS : evolving the spectra°°°°°°°°°")
   @everywhere sss = snap_in - lag
   @inbounds for snap in snap_in:lag:snap_fin   #..."lag" allows to use tracers file data every lag-steps (for too refined time sequences)
      @everywhere sss += lag
      @everywhere so += 1

      write_spectra_hdf5 = 0 #...1 if we want to write spectra on disk as hdf5 files (can use a lot of memory )
      t00 = time()

      @everywhere snapn = string(sss)

      #...the redshift of two snapshots are read to compute the dt time in between
      if snap < 100
         @everywhere snapn = string("0", sss)
      end
      if snap < 10
         @everywhere snapn = string("00", sss)
      end
      @everywhere snapn2 = string(sss + lag)
      if sss + lag < 100
         @everywhere snapn2 = string("0", sss + lag)
      end
      if snap + lag < 10
         @everywhere snapn2 = string("00", sss + lag)
      end

      @everywhere file_trac = string(rootf, "_tracers_", snapn, ftag, ".hdf5")
      @everywhere z = h5read(file_trac, "Redshift")
      @everywhere file_trac2 = string(rootf, "_tracers_", snapn2, ftag, ".hdf5")
      @everywhere z2 = h5read(file_trac2, "Redshift")
      @everywhere time1 = ustrip(lookback_time(cosmo, z)::Number)   #...use cosmology to compute time elapsed between snapshots 
      @everywhere time2 = ustrip(lookback_time(cosmo, z2)::Number)
      @everywhere dt = abs(1e9 * (time2 - time1))
      println("evolving snapshots", snapn, " with timestep dt=", dt, " yr and redshift z=", z)
      @everywhere zz = convert(Float64, z)

      cdd = h5read(file_trac, "DensityConversionFactor")
      cv = h5read(file_trac, "VelocityConversionFactor")
      cb = sqrt(cdd * 4.0 * pi) * cv * (1 + z)^2  #b-field in proper Gauss   \         \
      p1 = h5read(file_trac, "xcoord")
      @views tr[1, 1:ntr0] .= p1[1:ntr0]
      p1 = h5read(file_trac, "ycoord")
      @views tr[2, 1:ntr0] .= p1[1:ntr0]
      p1 = h5read(file_trac, "zcoord")
      @views tr[3, 1:ntr0] .= p1[1:ntr0]
      p1 = h5read(file_trac, "density")
      @views tr[7, 1:ntr0] .= p1[1:ntr0]
      p1 = h5read(file_trac, "temperature")
      @views tr[8, 1:ntr0] .= p1[1:ntr0]
      p1 = h5read(file_trac, "mach")
      @views tr[9, 1:ntr0] .= p1[1:ntr0]
      p1 = h5read(file_trac, "B")
      @views tr[10, 1:ntr0] .= (p1[1:ntr0]) #....B is in Gauss
      p1 = h5read(file_trac, "div_v")
      @views tr[11, 1:ntr0] .= p1[1:ntr0]
      p1 = h5read(file_trac, "curl_v")
      @views tr[12, 1:ntr0] .= p1[1:ntr0]
      p1 = h5read(file_trac, "vturb")
      @views tr[14, 1:ntr0] .= p1[1:ntr0]

      #....checking whether amplified field, bturb, should be larger than the MHD one (see Beduzzi et al. 2024 for details)
      @inbounds for i in 1:ntr0
         @fastmath turb = tr[12, i] * (scale_turbo * kpc) #...physical units  
         @fastmath bturb = (4 * pi * ηB * tr[7, i] * turb^2)^0.5    #...physical B
         if bturb >= tr[10, i]
            tr[10, i] = bturb
         end
         end

          @inbounds for i in 1:ntr0   #...generated particle ID 
         tr[13, i] = i * 1.0
          end

      #....PARALLEL CALL TO THE FP SOLVER
      #...the user has to select a number of tracers before hand, with addprocs(N)
      nw = nworkers()
      dax = convert(Int64, trunc(ntr / nw))
      proc_0 = 2 #...
      pid = collect(proc_0:proc_0+nw-1)
      intern_id = collect(1:ntr0)
      use_particle = intern_id[sel]
      use_spectra = intern_id[1:ntr]
      chunk_particle = chunk(use_particle, nw)
      chunk_spectra = chunk(use_spectra, nw)
      nw_chunk = nw

      println("using ", nw_chunk, " cores")

      if length(chunk_particle) < nw
         nw_chunk = length(chunk_particle)
      end

      @sync for i in 1:nw_chunk
         @async pe[:, chunk_spectra[i]] = remotecall_fetch(spectra, pid[i], tr[:, chunk_particle[i]], pe[:, chunk_spectra[i]], dt * 1e-9, zz, lsize, snap, pradio1[chunk_spectra[i], :], "model1", maxp)
      end

      @sync for i in 1:nw_chunk
         @async pradio1[chunk_spectra[i], :] = remotecall_fetch(sync_wrapper, pid[i], pval, tr[10, chunk_particle[i]], pe[:, chunk_spectra[i]], freqf, pradio1[chunk_spectra[i], :])
      end
      @sync for i in 1:nw_chunk
         @async pe2[:, chunk_spectra[i]] = remotecall_fetch(spectra, pid[i], tr[:, chunk_particle[i]], pe2[:, chunk_spectra[i]], dt * 1e-9, zz, lsize, snap, pradio2[chunk_spectra[i], :], "model2", maxp)
      end
      @sync for i in 1:nw_chunk
         @async pradio2[chunk_spectra[i], :] = remotecall_fetch(sync_wrapper, pid[i], pval, tr[10, chunk_particle[i]], pe2[:, chunk_spectra[i]], freqf, pradio2[chunk_spectra[i], :])
      end
      @sync for i in 1:nw_chunk
         @async pe3[:, chunk_spectra[i]] = remotecall_fetch(spectra, pid[i], tr[:, chunk_particle[i]], pe3[:, chunk_spectra[i]], dt * 1e-9, zz, lsize, snap, pradio3[chunk_spectra[i], :], "model3", maxp)
      end
      @sync for i in 1:nw_chunk
         @async pradio3[chunk_spectra[i], :] = remotecall_fetch(sync_wrapper, pid[i], pval, tr[10, chunk_particle[i]], pe3[:, chunk_spectra[i]], freqf, pradio3[chunk_spectra[i], :])
      end

      #....OUTPUTS : MAPS, PARTICLE AND RADIO SPECTRA
      dezoom = 1.0 #....this makes the side maps a factor "dezoom" smaller than the original n,n resolution - to save space, since most of the maps are empty
      nima = convert(Int64, trunc(n * dezoom))
      map_radioy = Array{Float64}(undef, nima, nima, nfreq[1] + 1)
      map_radioy[:] .= 0.0
      map_radioy2 = Array{Float64}(undef, nima, nima, nfreq[1] + 1)
      map_radioy2[:] .= 0.0
      map_radioy3 = Array{Float64}(undef, nima, nima, nfreq[1] + 1)
      map_radioy3[:] .= 0.0

      pe_tot = fill(0.0, np)
      pe_tot2 = fill(0.0, np)
      pe_tot3 = fill(0.0, np)
      nn1 = 1
      nn2 = ntr - 1

      rspec = 3000.0  #.....useful to select tracers using a radial range - presently not used

      total_radio1 = Array{Float64}(undef, nfreq[1])
      total_radio2 = Array{Float64}(undef, nfreq[1])
      total_radio3 = Array{Float64}(undef, nfreq[1])

      #..assign tasks to processors
      ip1 = [1, 2, 3, 4, 5, 6]

      if nw < 6 #...we have more tasks than processors
         for i in nw+1:6
            ip1[i] = i - nw
         end
      end

      #.....parallel calls for different analysis tasks
      @sync @async pe_tot[:] = remotecall_fetch(total_spectra, pid[ip1[1]], pval, pe[:, :], pe_tot[:], nn1, nn2, tr[1, sel[nn1:nn2]], tr[2, sel[nn1:nn2]], tr[3, sel[nn1:nn2]], rspec, xc)
      @sync @async pe_tot2[:] = remotecall_fetch(total_spectra, pid[ip1[2]], pval, pe2[:, :], pe_tot2[:], nn1, nn2, tr[1, sel[nn1:nn2]], tr[2, sel[nn1:nn2]], tr[3, sel[nn1:nn2]], rspec, xc)
      @sync @async pe_tot3[:] = remotecall_fetch(total_spectra, pid[ip1[3]], pval, pe3[:, :], pe_tot3[:], nn1, nn2, tr[1, sel[nn1:nn2]], tr[2, sel[nn1:nn2]], tr[3, sel[nn1:nn2]], rspec, xc)
      #.....maps along Y - different projections are generated just by selecting different coordinates of tracers here below 
      @sync @async map_radioy = remotecall_fetch(map2d, pid[ip1[4]], pradio1[:, :], tr[1, sel[nn1:nn2]], tr[3, sel[nn1:nn2]], tr[2, sel[nn1:nn2]], nn1, nn2, freqf, map_radioy, xc, pe[:, :], dezoom)
      @sync @async map_radioy2 = remotecall_fetch(map2d, pid[ip1[5]], pradio2[:, :], tr[1, sel[nn1:nn2]], tr[3, sel[nn1:nn2]], tr[2, sel[nn1:nn2]], nn1, nn2, freqf, map_radioy2, xc, pe2[:, :], dezoom)
      @sync @async map_radioy3 = remotecall_fetch(map2d, pid[ip1[6]], pradio3[:, :], tr[1, sel[nn1:nn2]], tr[3, sel[nn1:nn2]], tr[2, sel[nn1:nn2]], nn1, nn2, freqf, map_radioy3, xc, pe3[:, :], dezoom)
    
         @everywhere maxp = 1
         @inbounds for j in 1:np
         if pe_tot3[j] > 1e50
            @everywhere maxp += 1
         end
         end
      #...avoiding the very outer frame where some tracer can be artificially stuck
      nima = size(map_radioy)
      d1 = 1 #5
      d2 = nima[2] #nima[2]-5
      maxm = maximum(map_radioy3[d1:d2, d1:d2, nfreq[1]])
      @inbounds for j in 1:nima[2]
         @inbounds for i in 1:nima[2]
            if map_radioy3[i, j, nfreq[1]] == maxm
               imax = i
               jmax = j
               global imax, jmax
            end
         end
      end
      dr = convert(Int64, trunc(dezoom * 500 / dx))
      i1 = imax - dr
      i2 = imax + dr
      j1 = jmax - dr
      j2 = jmax + dr
      if i1 <= 1
         i1 = 1
      end
      if j1 <= 1
         j1 = 1
      end
      if i2 >= nima[2]
         i2 = nima[2]
      end
      if j2 >= nima[2]
         j2 = nima[2]
      end
      @inbounds for f in 1:nfreq[1]
         total_radio1[f] = sum(map_radioy[i1:i2, j1:j2, f])
         total_radio2[f] = sum(map_radioy2[i1:i2, j1:j2, f])
         total_radio3[f] = sum(map_radioy3[i1:i2, j1:j2, f])
      end

      a = plot_spectra(pe_tot, pe_tot2, pe_tot3, pval)
      a = plot_radio_spectra(total_radio1, total_radio2, total_radio3, freqf)
      model1 = "1"
      model2 = "2"
      model3 = "3"

      #...fits maps of radio emission 
      a = write_maps_fits(map_radioy, model1)
      a = write_maps_fits(map_radioy2, model2)
      a = write_maps_fits(map_radioy3, model3)

      #....png maps with projected density and contours of radio emission 
      a = do_contours(map_radioy, map_radioy2, map_radioy3, so)

      #.....writing all HDF5 data is memory expensive so we activate write_data_sync only sometime, when ss=si
      write_spectra_hdf5 = 0
      ss = 20 * (snap - snap_in) / (90)      #....this activates it only every few timesteps
      si = convert(Int64, trunc(ss))
      println(ss, " ", si)
      if ss == si
         write_spectra_hdf5 = 1
      end

      if write_spectra_hdf5 == 10 #...n
         a = write_hdf5(pe, pe2, pe3, pval, tr, sel, pradio1, pradio2, pradio3, z, freqf, nn1, nn2)
      end

   end
   println("° ° ° ° GAME, SET AND MATCH : ROGER WON°°°°°°°°°")
end
