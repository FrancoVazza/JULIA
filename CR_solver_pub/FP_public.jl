using HDF5
using FITSIO
using Statistics
using Plots
using SpecialFunctions
pyplot()
Plots.PyPlotBackend()
using LaTeXStrings

        #folders and includes
        main=""
        include(string(main,"/param_spectra.jl"))
        include(string(main,"/loss_gain1.jl"))
        const root=string("/tracers/")
        const root_map=root
        const root_out=root
        const root1=root


        # Run paramters
        const dx=8.860   #constant resolution of the simulations [ckpc]
        snap_in=103   #initial timestep
        snap_fin=153  #final timestep
        ntr=25        #total number of tracers
        dt=3.98e7     #time difference between snapshots [yr] - assumed to be constant here
        mass_trac=msol*1e2 #..mass of relativistic electrons injected by jets for each tracer

        tr=Array{Float64}(undef,13,ntr) #..tracer fields
        tr.=0.

        #.....particle spectra, for 2 simultanously tested scenarios
        #....parameters for the minimum/maximum momenta and momentum bins have to be set in param_spectra.jl function

        pe=fill(1e0,np,ntr)    #...only cooling
        pe2=fill(1e0,np,ntr)   #...cooling and shocks injection
        pe3=fill(1e0,np,ntr)   #...cooling and shocks injection & reacceleration
        pe_total=fill(1e0,np)  #...total particle spectra
        pe_total2=fill(1e0,np)
        pe_total3=fill(1e0,np)
      #...initialisation of particle spectra from radiogalaxies (assuming an in_spec energy spectrum)
      in_spec=2.00
      @inbounds @simd   for g in 1:np-1
            @views pe[g,:].=(mass_trac*p_min^(in_spec-1)/(0.6*mp))*(pval[g])^(-in_spec)   #...number of relativistic electrons injected by radiogalaxy
            @views pe2[g,:].=(mass_trac*p_min^(in_spec-1)/(0.6*mp))*(pval[g])^(-in_spec)
            @views pe3[g,:].=(mass_trac*p_min^(in_spec-1)/(0.6*mp))*(pval[g])^(-in_spec)
      end


      #....main simulation loop
@inbounds for snap in snap_in:snap_fin
        println("doing snapshot",snap)
        snapn=string(snap)

        file_trac=string(root,"_tracers_",snapn,"_test.hdf5")
        z=h5read(file_trac,"Redshift")
        p1=h5read(file_trac, "xcoord")
        @views     tr[1,1:ntr].=p1[1:ntr]
        p1=h5read(file_trac, "ycoord")
        @views     tr[2,1:ntr].=p1[1:ntr]
        p1=h5read(file_trac, "zcoord")
        @views     tr[3,1:ntr].=p1[1:ntr]
        p1=h5read(file_trac, "density")
        @views     tr[7,1:ntr].=p1[1:ntr]
        p1=h5read(file_trac, "temperature")
        @views     tr[8,1:ntr].=p1[1:ntr]
        p1=h5read(file_trac, "mach")
        @views     tr[9,1:ntr].=p1[1:ntr]
        p1=h5read(file_trac, "B")
        @views      tr[10,1:ntr].=p1[1:ntr]
        p1=h5read(file_trac, "div_v")
        @views      tr[11,1:ntr].=p1[1:ntr]
        p1=h5read(file_trac, "curl_v")
        @views      tr[12,1:ntr].=p1[1:ntr]
        p1=h5read(file_trac, "particleID")
        @views      tr[13,1:ntr].=p1[1:ntr]
        p1=nothing


        tr2=similar(tr)
        tr3=similar(tr2)
        tr2[:].=tr[:]
        tr3[:].=tr[:]

         sel=findall( x->(x >= 2.0),tr[9,:])
        @views   tr[9,:].=0.   #all shocks in the first population are removed - to compare the effect of shocks in the tr2 population
        @views   tr2[9,sel].=0. #M>2.0 shocks are removed - only shock-recceleration is relevant in this regime

        lsize=dx*0.1 #linear size of the equivalent volume associated with each tracer. Can be fixed of set variable with the changing density, depending on the setup.                                          \

        spectra(tr,pe,dt*1e-9,z,lsize,snap)
        spectra(tr2,pe2,dt*1e-9,z,lsize,snap)
        spectra(tr3,pe3,dt*1e-9,z,lsize,snap)

                @inbounds  @simd   for g in 1:np
                pe_total[g]=sum(pe[g,1:ntr])*pval[g]
                pe_total2[g]=sum(pe2[g,1:ntr])*pval[g]
                pe_total3[g]=sum(pe3[g,1:ntr])*pval[g]
                end

     plot(pval,pe_total[:],label="losses",color="black",alpha=1.0,line=:solid)
     plot!(pval,pe_total2[:],color="blue",label="losses+shocks inj.",line=:dot,alpha=0.5)
     plot!(pval,pe_total3[:],color="red",label="losses+shocks (inj.+re.)",line=:dash,alpha=0.5)
     title!(string("z=",z),fonts=20)
     yaxis!(L"p \cdot N(p)",:log10,(1e40,1e70),fonts=20)
     xaxis!(L"P/(m_ec)",:log10,(p_min,p_max*0.99),fonts=20)

    filep1=string(root_out,"spectra_out",snapn,".png")
    savefig(filep1)

end

end#
