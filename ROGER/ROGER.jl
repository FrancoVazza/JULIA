   @everywhere using LaTeXStrings
   @everywhere using Plots
   @everywhere pyplot()
   @everywhere Plots.PyPlotBackend()

   @everywhere  using FITSIO
   @everywhere using Statistics
   @everywhere using Distributed
   @everywhere using SharedArrays
   @everywhere using SpecialFunctions
   @everywhere using DelimitedFiles
   @everywhere using HDF5
   # SIMULATION PARAMETERS

   @everywhere const n=640        #grid size of parent simulation
   @everywhere const dx=0.001     #cell resolution in kpc
   @everywhere xc=1/(dx)
   @everywhere const lsize=dx   #...lsize is also assumed to be the volume associated with each tracer - users must adjust this to their specific needs
   @everywhere ntr0=80        #...number of tracer particles
   @everywhere z_in=0.5        #...initial snapshot
   @everywhere snap_in=103       #...initial snapshot
   @everywhere snap_fin=199       #...final snapshot
   @everywhere dt=5.0e7          #..timestep [yr] - simple constant case. 
   @everywhere scale = 15.0         #...[kpc] - spatial scale used to measure the vorticity - goes into the turbulent energy flux and Fermi II reacceleration. 
   @everywhere test=""            #....name of this test run
   @everywhere run="test"         #....run name

   @everywhere main=string("/Users/francovazza/Dropbox/Julia/FOKKER_DEBUG/public/")
   @everywhere root_out=string(main)
   @everywhere root_tracers=string(main,"tracers/")

   @everywhere include(string(main,"/param_spectra_par_log.jl"))   #...parameters of spectra
   @everywhere include(string(main,"/loss_gain_par_log.jl"))       #...loss and gain terms
   @everywhere include(string(main,"/analysis_log.jl"))

   
   #....TRACER INITIALISATION
   @everywhere  snapn=string("",snap_in)

   @everywhere  file_trac=string(root_tracers,"_tracers_",snapn,"_public.hdf5")
   tr=Array{Float64}(undef,13,ntr0)
   tr.=0.
   @everywhere fi_jet=1e-4   #.....initial density of relativistic electrons >pmin wtr thermal electron density
   @everywhere mass_gas=msol*1e8
   @everywhere mass_trac=mass_gas*fi_jet



   global   pe=fill(1e-30,np,ntr0)
   global   pe2=fill(1e-30,np,ntr0)
   global   pe3=fill(1e-30,np,ntr0)


   #...initialisation using the CIOFF model - assuming particles are initially injected by radiogalaxies (-> Vazza et al. 2022 A&A for details)
   #...replace with different spectral initialisations for different problems
   #...avoid initialising particles with 0 values as this produces NaN.

   @everywhere t_on=4.e7*3.154e7    #...t_on [s]
   @everywhere t_off=1.1e5*3.154e7  #...t_off [s]
   @everywhere t_tot=t_on+t_off
   @everywhere bjet =1e-5    #..B-field inside lobes at injection [in G]
   @everywhere β=b2*(bjet^2. +(1. +z_in)^4*b3)
   @everywhere γc1=1/(β*(t_tot))
   @everywhere γc2=1/(β*(t_off))
   global in_spec=2.2
   global pmin_jet=10^pval[1]   #...minimum electron momentum in jets = minimum momentum set in parameter file.

   @inbounds @simd    for g in 1:np
    ppv=10^pval[g]
    ppmax=10^pval[np-1]
    ga=sqrt(1+(ppv)^2)
    gamma_max=sqrt(1+(ppmax)^2)
    if ga <= γc1
       pe[g,:].=(1-β*ga*t_off)^(in_spec-1)-(1-β*ga*(t_tot))^(in_spec-1)
    end
    if ga >= γc1 && ga <=γc2
       pe[g,:].=(1-β*ga*t_off)^(in_spec-1)
    end
    if ga > gamma_max
       pe[g,:].=0.0
    end

    pe[g,:].=(mass_trac*pmin_jet^(in_spec-1)/(0.6*mp))*(ppv)^(-in_spec-1)*pe[g,:]
    pe2[g,:].=pe[g,:]
    pe3[g,:].=pe[g,:]
   end


#....plotting of the initial input spectrum (before losses)
      pe_tot=fill(0.0,np)
      pe_tot2=fill(0.0,np)
      pe_tot3=fill(0.0,np)
   @inbounds @simd   for i in 1:ntr0
     @inbounds @simd for g in 1:np
   pe_tot[g]+=(pe[g,i]*10^pval[g])
   pe_tot2[g]+=(pe2[g,i]*10^pval[g])
   pe_tot3[g]+=(pe3[g,i]*10^pval[g])
   end
   end

   lfs=14
   xfs=12
   xtfs=12
           plot(pval,pe_tot[:],label="losses",color="black",line=:dash,linewidth=1,dpi=500,alpha=0.7,grid=false,legendfontsize=lfs,yguidefontsize=xfs,xguidefontsize=xfs,xtickfontsize=xtfs,ytickfontsize=xtfs)
           plot!(pval,pe_tot2[:],color="green",label="losses+shocks",line=:dash,linewidth=2,alpha=0.6)
           plot!(pval,pe_tot3[:],color="red",label="losses+shocks+turb",line=:dash,linewidth=3,alpha=0.5)
           yaxis!(L"pN(p)",:log10,(1e40,1e65),fonts=20)

       xaxis!(L"\log(P/(m_ec)",(p_min,p_max),fonts=20)


       filep1=string(root_out,run,"initial.png")
       savefig(filep1)


   println("injection of electrons done, starting simulation")

   @inbounds for snap in snap_in:snap_fin

   t00=time()
   snapn=string(snap)

   file_trac=string(root_tracers,"_tracers_",snapn,"_public.hdf5")
   z=h5read(file_trac,"Redshift")
   #for cosmological applications, the link between time and redshift might be included here.

    p1=h5read(file_trac, "xcoord")
    @views     tr[1,1:ntr0].=p1[1:ntr0]
    p1=h5read(file_trac, "ycoord")
    @views     tr[2,1:ntr0].=p1[1:ntr0]
    p1=h5read(file_trac, "zcoord")
    @views     tr[3,1:ntr0].=p1[1:ntr0]
    p1=h5read(file_trac, "density")
    @views     tr[7,1:ntr0].=p1[1:ntr0]
    p1=h5read(file_trac, "temperature")
    @views     tr[8,1:ntr0].=p1[1:ntr0]
    p1=h5read(file_trac, "mach")
    @views     tr[9,1:ntr0].=p1[1:ntr0]

   m_max=maximum(p1[1:ntr0])
   if m_max > mthr
   println("Mmax=",m_max)
    end
    p1=h5read(file_trac, "B")

    @views      tr[10,1:ntr0].=(1e-6*p1[1:ntr0])   #..B-field converted into G
    p1=h5read(file_trac, "div_v")
    @views      tr[11,1:ntr0].=p1[1:ntr0]
    p1=h5read(file_trac, "curl_v")
    @views      tr[12,1:ntr0].=p1[1:ntr0]
    p1=h5read(file_trac, "particleID")
    @views      tr[13,1:ntr0].=p1[1:ntr0]
    p1=nothing     #...only first injection
     tr2=similar(tr)  #...shock (re)acceleration only
     tr3=similar(tr)  #...also ASA
     tr2[:].=tr[:]
     tr3[:].=tr[:]
   @views   tr2[12,:].=0.   #model CS has no turbulent re-acceleration
   @views   tr[12,:].=0.    #model C has no turbulent and shock re-acceleration
   @views   tr[9,:].=0.

   if snap <= snap_in+3   #in case particle are initially in fast jets (as in its first application), we prefer to manually switch off turbulent re-acceleration from spurious vorticity
   @views   tr3[12,:].=0.
   end

    println(snapn," z=",z," dt=",dt)
    zz=convert(Float64,z)
    #....PARALLEL CALL TO THE FP SOLVER
    nw=nworkers()
    dax=convert(Int64,trunc(ntr0/nw))


  p_0=2
   pid=collect(p_0:p_0+nw-1)
        i1=1
        i2=dax
        i3=2*dax
        i4=3*dax
        i5=4*dax

        r1=remotecall(spectra,pid[1],tr[:,i1:i2-1],pe[:,i1:i2-1],dt*1e-9,zz,lsize,snap)
        r2=remotecall(spectra,pid[2],tr[:,i2:i3-1],pe[:,i2:i3-1],dt*1e-9,zz,lsize,snap)
        r3=remotecall(spectra,pid[3],tr[:,i3:i4-1],pe[:,i3:i4-1],dt*1e-9,zz,lsize,snap)
        r4=remotecall(spectra,pid[4],tr[:,i4:i5],pe[:,i4:i5],dt*1e-9,zz,lsize,snap)

        @views    pe[:,i1:i2-1]=fetch(r1)
        @views    pe[:,i2:i3-1]=fetch(r2)
        @views    pe[:,i3:i4-1]=fetch(r3)
        @views    pe[:,i4:i5]=fetch(r4)

       r1=remotecall(spectra,pid[1],tr2[:,i1:i2-1],pe2[:,i1:i2-1],dt*1e-9,zz,lsize,snap)
       r2=remotecall(spectra,pid[2],tr2[:,i2:i3-1],pe2[:,i2:i3-1],dt*1e-9,zz,lsize,snap)
       r3=remotecall(spectra,pid[3],tr2[:,i3:i4-1],pe2[:,i3:i4-1],dt*1e-9,zz,lsize,snap)
       r4=remotecall(spectra,pid[4],tr2[:,i4:i5],pe2[:,i4:i5],dt*1e-9,zz,lsize,snap)

       @views    pe2[:,i1:i2-1]=fetch(r1)
       @views    pe2[:,i2:i3-1]=fetch(r2)
       @views    pe2[:,i3:i4-1]=fetch(r3)
       @views    pe2[:,i4:i5]=fetch(r4)

      r1=remotecall(spectra,pid[1],tr3[:,i1:i2-1],pe3[:,i1:i2-1],dt*1e-9,zz,lsize,snap)
      r2=remotecall(spectra,pid[2],tr3[:,i2:i3-1],pe3[:,i2:i3-1],dt*1e-9,zz,lsize,snap)
      r3=remotecall(spectra,pid[3],tr3[:,i3:i4-1],pe3[:,i3:i4-1],dt*1e-9,zz,lsize,snap)
      r4=remotecall(spectra,pid[4],tr3[:,i4:i5],pe3[:,i4:i5],dt*1e-9,zz,lsize,snap)

      @views    pe3[:,i1:i2-1]=fetch(r1)
      @views    pe3[:,i2:i3-1]=fetch(r2)
      @views    pe3[:,i3:i4-1]=fetch(r3)
      @views    pe3[:,i4:i5]=fetch(r4)

     pe_tot=fill(0.0,np)
     pe_tot2=fill(0.0,np)
     pe_tot3=fill(0.0,np)

        nn1=1
        nn2=ntr0

   r1=remotecall(total_spectra,pid[1],pval,pe[:,:],pe_tot,nn1,nn2,tr[1,nn1:nn2],tr[2,nn1:nn2],tr[3,nn1:nn2])
   r2=remotecall(total_spectra,pid[2],pval,pe2[:,:],pe_tot2,nn1,nn2,tr2[1,nn1:nn2],tr2[2,nn1:nn2],tr2[3,nn1:nn2])
   r3=remotecall(total_spectra,pid[3],pval,pe3[:,:],pe_tot3,nn1,nn2,tr3[1,nn1:nn2],tr3[2,nn1:nn2],tr3[3,nn1:nn2])

   pe_tot[:]=fetch(r1)#
   pe_tot2[:]=fetch(r2)
   pe_tot3[:]=fetch(r3)


lfs=14
xfs=12
xtfs=12

        plot(pval,pe_tot[:],label="losses",color="black",line=:dash,dpi=500,linewidth=1,alpha=0.7,grid=false,legendfontsize=lfs,yguidefontsize=xfs,xguidefontsize=xfs,xtickfontsize=xtfs,ytickfontsize=xtfs)
        plot!(pval,pe_tot2[:],color="green",label="losses+shocks",line=:dash,linewidth=2,alpha=0.6)
        plot!(pval,pe_tot3[:],color="red",label="losses+shocks+turb",line=:dash,linewidth=3,alpha=0.5)

println("time per timestep of parallel run=",time()-t00)

    title!(string("z=",zz),fonts=20)
    yaxis!(L"pN(p)",:log10,(1e40,1e65),fonts=20)

    xaxis!(L"\log(P/(m_ec)",(p_min,p_max),fonts=20)


    filep1=string(root_out,run,"spectra_new_momenta_radial",snapn,"_",test,".png")
    savefig(filep1)

   #   uncomment to write particle spectra on disk for each snapshot
   #   warning: quite memory consuming since all particles and all momentum bins will be saved
   #   notice that it won't allow overwriting previously generated HDF5 files with the same names - they must be manually removed
   # filep1=string(root_out,"_radio_tracer_",run,"_",snapn,"_pradio_",test,".hdf5")
   #   h5write(filep1,"N(P)_A",pe[:,:])
   #   h5write(filep1,"N(P)_B",pe2[:,:])
   #   h5write(filep1,"N(P)_C",pe3[:,:])
   #   h5write(filep1,"P",pval[:])
   #   h5write(filep1,"Redshift",z)

end
