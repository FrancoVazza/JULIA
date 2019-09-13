#.....JULIA SHOCK FINDER FOR UNIGRID DATA
#.....INPUT DATA ARE SUPPOSED TO BE REGULAR GRID, WRITTEN IN HDF5 FORMAT
#.....THE CODE CAN ALSO RUN IN PARALLEL JUST BY ADDING PROCS


@everywhere main="/home/PERSONALE/franco.vazza2/Dropbox/Julia/"   #Julia folder
#@everywhere root="/home/PERSONALE/franco.vazza2/Desktop/DATA/CHRONOS++/85Mpc/1024/"   #input files folder 
@everywhere root="/home/PERSONALE/franco.vazza2/Desktop/DATA/CHRONOS3/800/"

#@everywhere mod="newBH12b_csf3dR"  #filename
@everywhere mod="csfR"  #filename
@everywhere snap="004"             #snapshot


@everywhere file1=string(root,"s",snap,"/full_dtb_",mod,"_",snap)    #input files    *dt* contains gas density, dark matter density and gas temperature
@everywhere file2=string(root,"s",snap,"/full_v_",mod,"_",snap)   # *b* contains the 3-d components of the v-field,
@everywhere file3=string(root,"s",snap,"/full_dtb_",mod,"_",snap)   # *b* contains the 3-d components of the B-field,
@everywhere file4=string(root,"s",snap,"/full_chem_",mod,"_",snap)   #chemistry files

 #...conversion factors
#@everywhere filec=string("/home/PERSONALE/franco.vazza2/Desktop/DATA/CHRONOS++/85Mpc/1024/conv/RD0",snap,".conv2")
@everywhere filec=string("/home/PERSONALE/franco.vazza2/Desktop/DATA/CHRONOS3/800/conv/RD0",snap,".conv2")
@everywhere a=readdlm(filec)
@everywhere const z=a[3]     #redshift                                         
@everywhere const cd=a[4] #conversion factor for density                    \  
@everywhere const cv=a[5] #conversion factor for velocity                   \  
@everywhere const cb=sqrt(cd*4.0*pi)*cv*(1+z)^2.  #b-field in Gauss   \                          
 

  #...modules related to cosmology and contants
@everywhere include(string(main,"constants.jl"))
@everywhere include(string(main,"cosmological_param.jl"))

  #...modules for different physical effects
@everywhere include(string(main,"shocks.jl")) #shock finder 
@everywhere include(string(main,"HI.jl"))  #neutral Hydrogen
@everywhere include(string(main,"SZ.jl"))  #SZ decrement
@everywhere include(string(main,"radio_power.jl")) #synchrotron radio power
@everywhere include(string(main,"Xray_em.jl"))  #Xray emission in band
@everywhere include(string(main,"turb.jl")) #turbulence estimate

@everywhere using HDF5    #necessary to open hdf5 files	
@everywhere using LaTeXStrings  # to include Latex-style writing in captions 



###############################################
#.......Reads in array fields and calls shock finder
@everywhere function make_map_parallel2(dat1)

n1=(size(dat1[1], 1), size(dat1[2], 1)) #size(dat1[3],1))   #this builds up the 3D structure from the dat1 array which is passed to each processor

#...here below we compute the min and max in each direction of the array.
#...notice in this implementation the domain is decomposed in rectangles which are long as the entire domain in the z direction (best for map-making)

	if los == "Z"
	i1=dat1[1][1]
	j1=dat1[2][1]
	i2=dat1[1][n1[1]]
	j2=dat1[2][n1[2]]
       
	end
	
	if los == "Y"
	i1=dat1[1][1]
	l1=dat1[2][1]
	i2=dat1[1][n1[1]]
	l2=dat1[2][n1[2]]
        
	end

	    
	if los == "X"
	j1=dat1[1][1]
	l1=dat1[2][1]
	j2=dat1[1][n1[1]]
	l2=dat1[2][n1[2]]
       
	end

	    map=Array{Float64}(n1[1],n1[2],nfields)
	    map[:]=0.

	const dl=20   # thickness of each tile (must be calibrated depending on available memory/node)
	const nstepsl=convert(Int64,trunc((n)/dl))  #number of tiles

	@inbounds @simd for kz in 1:nstepsl

    if los == "X"
    i1=dl*(kz-1)+1
    i2=dl*(kz)-1
    s2=i2
    s1=i1
    end

      if los == "Y"
    j1=dl*(kz-1)+1
    j2=dl*(kz)-1
    s2=j2
    s1=j1
    end

      if los == "Z"
    l1=dl*(kz-1)+1
    l2=dl*(kz)-1
    s2=l2
    s1=l1
      end

#...each processor reads the input files file1,file2 in parallel, i.e. each processor only reads a part of the entire domain, depending on the bounds defined above.

    #...reading raw HDF5 ENZO files.
    
  d=cd*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
 
  t=h5read(file1,"Temperature",(i1:i2,j1:j2,l1:l2))
  
  t=convert(Array{Float64,3},t)
  
  
  lx1=similar(t)
  lx2=similar(t)
   
  Xray_em(d,t,lx1,lx2) 
  
  sz_dec=similar(t)

  freqSZ=220e9 #...GHz (SPT-like)
  SZ(d,t,freqSZ,sz_dec)

  dm=cd*h5read(file1,"Dark_Matter_Density",(i1:i2,j1:j2,l1:l2))
                                      
   
    #...first round of map making (to spare memory)
@inbounds        @simd    for i in 1:s2-s1+1 

   if los == "Z"
  @views       map[:,:,1]+=d[:,:,i] #av.density map                                                                 
  @views       map[:,:,9]+=dm[:,:,i] 
  @views       map[:,:,10]+=(lx1[:,:,i]+lx2[:,:,i])
  @views       map[:,:,11]+=lx1[:,:,i]
  @views       map[:,:,12]+=lx2[:,:,i]
  @views       map[:,:,13]+=sz_dec[:,:,i]
   end

      if los == "Y"
   @views      map[:,:,1]+=d[:,i,:] #av.density map                                       
   @views      map[:,:,9]+=dm[:,i,:] 
   @views      map[:,:,10]+=(lx1[:,i,:]+lx2[:,i,:])
   @views      map[:,:,11]+=lx1[:,i,:]
   @views      map[:,:,12]+=lx2[:,i,:]
   @views      map[:,:,13]+=sz_dec[:,i,:]
   end

    if los == "X"
    @views     map[:,:,1]+=d[i,:,:] #av.density map                                       
    @views     map[:,:,9]+=dm[i,:,:]  #av. temperature map                                 
    @views     map[:,:,10]+=(lx1[i,:,:]+lx2[i,:,:])
    @views     map[:,:,11]+=lx1[i,:,:]
    @views     map[:,:,12]+=lx2[i,:,:]
    @views     map[:,:,13]+=sz_dec[i,:,:]
    
   end
   end
println("proc ", myid(), "has done first maps within ",i1," ",i2," ",j1," ",j2," ",l1," ",l2)
 
   lx1=nothing
   lx2=nothing
   sz_dec=nothing
   dm=nothing

  vx=cv*h5read(file2,"x-velocity",(i1:i2,j1:j2,l1:l2))
  vy=cv*h5read(file2,"y-velocity",(i1:i2,j1:j2,l1:l2))
  vz=cv*h5read(file2,"z-velocity",(i1:i2,j1:j2,l1:l2))

  n3=size(d)
  
    macx=similar(t)
    macy=similar(t)
    macz=similar(t)
    mach=similar(t)
    mach[:]=0.
    macx[:]=0
    macz[:]=0
    macy[:]=0
    
  mach=shocks(d,t,vx,vy,vz,macx,macy,macz,mach)   #...SHOCK FINDER

  vx=nothing
  vy=nothing
  vz=nothing

    
    tagnan=find( x->(x == NaN), mach)
    mach[tagnan]=0.
    tagnan=find( x->(x == Inf), mach)
    mach[tagnan]=0.
    
    tag0=find( x->(x >= mthr), mach)
    nshock=size(tag0)

    
    bz=cb*h5read(file3,"Bz",(i1:i2,j1:j2,l1:l2))
    by=cb*h5read(file3,"By",(i1:i2,j1:j2,l1:l2))
    bx=cb*h5read(file3,"Bx",(i1:i2,j1:j2,l1:l2))
    
    #...synchrotron map making
    #....pradio is a vector generated only using shocked cells
    #....pradio3D is a 3D array used to generate maps (this can be done more efficiently probably)
    nu=1.4 #...GHz  observing frequency
    pradio=0.
    
    pradio3D=similar(d)
    pradio3D[:]=0.

   @inbounds @simd for i in eachindex(tag0) #loop only over shocked cells
        ta=tag0[i]
        d1=d[ta]
        t1=t[ta]
        m1=mach[ta]
        bx1=bx[ta]
        by1=by[ta]
        bz1=bz[ta]
        prad=radio(d1,t1,m1,bx1,by1,bz1,nu,pradio)
       
        pradio3D[ta]=prad

     end

   rm=similar(d)
   
    for ir in eachindex(d)

        if los == "X"
            rm[ir]=dot(bx[ir],d[ir])
        end

        if los == "Y"
            rm[ir]=dot(by[ir],d[ir])
           end

        if los == "Z"
            rm[ir]=dot(bz[ir],d[ir])
        end
  
   # lx[ir]=dot(d[ir]^2.,sqrt(t[ir]))
    end
#second round of map making
@inbounds        @simd    for i in 1:s2-s1+1 

   if los == "Z"
  @views       map[:,:,3]+=abs(bz[:,:,i])  #projected Bz map
  @views       map[:,:,4]+=rm[:,:,i]  #RM map  
  @views       map[:,:,5]+=mach[:,:,i]    #av. Mach number 
  @views       map[:,:,8]+=pradio3D[:,:,i]
 
   end

      if los == "Y"  
   @views      map[:,:,3]+=abs(by[:,i,:])  #projected Bz map
   @views      map[:,:,4]+=rm[:,i,:]  #RM map  
   @views      map[:,:,5]+=mach[:,i,:]    #av. Mach number 
   @views      map[:,:,8]+=pradio3D[:,i,:]
   end

    if los == "X"
    @views     map[:,:,3]+=abs(bx[i,:,:])  #projected Bz map
    @views     map[:,:,4]+=rm[i,:,:]  #RM map  
    @views     map[:,:,5]+=mach[i,:,:]    #av. Mach number 
    @views     map[:,:,8]+=pradio3D[i,:,:]
    end
    end

println("proc ", myid(), "has done second maps within ",i1," ",i2," ",j1," ",j2," ",l1," ",l2)
 
    bx=nothing
    by=nothing
    bz=nothing
    mach=nothing
    pradio3D=nothing
    rm=nothing
   
    #...HI map making
#     dHI=cd*h5read(file4,"HI_Density",(i1:i2,j1:j2,l1:l2))
if isfile(file4) == true 
  dHI=cd*h5read(file4,"HI_Density",(i1:i2,j1:j2,l1:l2))
end

if isfile(file4) != true
     dHI=d*1e-6
end
     Ts=HI_compute(d,t,dHI)
    ####################################
    
          

 @inbounds @simd    for i in eachindex(d)
   @fastmath t[i]=dot(d[i],t[i])
    end
    
@inbounds        @simd    for i in 1:s2-s1+1 

   if los == "Z"                                  
  @views       map[:,:,2]+=t[:,:,i]  #mass weighted temperature map                                 
  @views       map[:,:,6]+=Ts[:,:,i]                                
  @views       map[:,:,7]+=dHI[:,:,i]
  end

      if los == "Y"
                                   
   @views      map[:,:,2]+=t[:,i,:]  #av. temperature map                                 
   @views      map[:,:,6]+=Ts[:,i,:]                                
   @views      map[:,:,7]+=dHI[:,i,:]
  
   end

    if los == "X"
                                      
    @views    map[:,:,2]+=t[i,:,:]  #av. temperature map                                 
    @views    map[:,:,6]+=Ts[i,:,:]                                
    @views     map[:,:,7]+=dHI[i,:,:]
   
   end
   end
    
println("proc ", myid(), "has done third maps within ",i1," ",i2," ",j1," ",j2," ",l1," ",l2)
 
    d=nothing
    t=nothing
    Ts=nothing
    dHI=nothing
   
end

    return map
    
end









######################################################
##....................Pablo MAIN PROGRAM..................#
#...output: a single FITS fle with many maps
######################################################


	for lll in 2:3   #selects a line of sight (Z is faster, X is slower
    	if lll == 1
        @everywhere los="Z"
	    end
	    if lll == 2
        @everywhere los="Y"
	    end
	     if lll == 3
        @everywhere los="X"
	    end

@everywhere const nfields=13
    
@everywhere using DistributedArrays
 map=DArray(make_map_parallel2,(n,n,nfields))  #...this builds up a Distributed Array (see Julia's documentation), i.e. an array made of chunks, each passed to a different processor. The distributed array is directly passed to the make_map function defined above

 map_total=convert(Array,map)  #....this step assembles a unique array by combining the different pieces of the distributed array 
println("MAPS DONE")
println("The minimum and maximum density in the map are ",minimum(map_total[:,:,1])," ",maximum(map_total[:,:,1]))
println("The minimum and maximum Mach  in the map are ",minimum(map_total[:,:,4])," ",maximum(map_total[:,:,4]))



#.....total projected values along z-axis
id=1
itmw=2
ib=3
irm=4
imach=5
iTH=6
idH=7
ipradio=8
idm=9
ilx=10
ilx1=11
ilx2=12
isz=13

@views    map_total[:,:,ib]/=n #...Bz
@views    map_total[:,:,irm]*=(812.*res/(1e-6*1e-3*mu*mp))
@views    map_total[:,:,imach]=(map_total[:,:,imach]+e-10)/n #...Bz
@views    map_total[:,:,iTH]/=n 
@views    map_total[:,:,idH]/=n
@views    map_total[:,:,ipradio]+=1e-30
@views    map_total[:,:,idm]/=n

@inbounds for j in 1:n
@simd for i in 1:n
    map_total[i,j,itmw]=map_total[i,j,itmw]/(map_total[i,j,id]) 
 end
 end

@views    map_total[:,:,id]/=n  #...density 


#....writing FITS files with images
@everywhere using FITSIO

filep1=string(root,"map_all",mod,"_",snap,"_800_",los,"new.fits")
 f = FITS(filep1, "w");
write(f,map_total)
close(f)

end


