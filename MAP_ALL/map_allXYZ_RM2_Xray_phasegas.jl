
@everywhere main="/home/PERSONALE/franco.vazza2/Dropbox/Julia/"    #main folder for Julia programs
#....files to read in
#@everywhere root="/Users/francovazza/Deskto   /data/DATA/CHRONOS/test_dataset/"
@everywhere root="/home/PERSONALE/franco.vazza2/Desktop/DATA/CHRONOS++/85Mpc/1024/s007/"
#root="/p/scratch/stressicm/highZ/A/7DD0051/"

#@everywhere root="/home/PERSONALE/franco.vazza2/Desktop/DATA/CHRONOS3/800/"

@everywhere mod="newBH12b_csf3dR"

@everywhere snap="007"
@everywhere sna=snap

@everywhere file1=string(root,"full_dtb_",mod,"_",snap)    #input files    *dt* contains gas density, dark matter density and gas temperature#
@everywhere file2=string(root,"full_v_",mod,"_",snap)   # *b* contains the 3-d components of the v-field,
@everywhere file3=string(root,"full_dtb_",mod,"_",snap)   # *b* contains the 3-d components of the B-field,

#@everywhere file1=string(root,"full_",mod,"_dt_",snap)    #input files    *dt* contains gas density, dark matter density and gas temperature
#@everywhere file2=string(root,"full_",mod,"_v_",snap)   # *b* contains the 3-d components of the v-field,
#@everywhere file3=string(root,"full_",mod,"_b_",snap)

@everywhere file4=string(root,"s",snap,"/full_chem_",mod,"_",snap)

#@everywhere file4=string(root,"full_chem_",mod,"_",snap)

@everywhere filec=string("/home/PERSONALE/franco.vazza2/Desktop/DATA/CHRONOS++/85Mpc/1024/conv/RD0",sna,".conv2")
@everywhere a=readdlm(filec)
@everywhere const z=a[3]     #redshift                                         \
                                                                                
@everywhere const cd=a[4] #conversion factor for density                    \
                                                                                
@everywhere const cv=a[5] #conversion factor for velocity                   \
                                                                                
@everywhere const cb=sqrt(cd*4.0*pi)*cv*(1+z)^2.  #b-field in Gauss   \
                                                                                
#@everywhere const dx=14 #kpc                                                  
@everywhere const mp=1.67e-24 #proton mass                                      
#@everywhere const n=1536  #...the 1D size of input grid in the HDF5 file.      

@everywhere include(string(main,"constants.jl"))
@everywhere include(string(main,"cosmological_param.jl"))
@everywhere dx=res
@everywhere include(string(main,"shocks.jl"))
@everywhere include(string(main,"HI.jl"))
@everywhere include(string(main,"radio_power.jl"))
@everywhere include(string(main,"Xray_em.jl"))
@everywhere using Interpolations
@everywhere include(string(main,"turb.jl"))

#@everywhere using DistributedArrays    #necessary for parallel map making
@everywhere using HDF5    #necessary to open hdf5 files	
#@everywhere using LaTeXStrings  # to include Latex-style writing in captions 
#@everywhere using Plot   #...to use PyPlot routines

@everywhere const mthr=1.5  #...threshold for Mach number estimation




###############################################
#.......Reads in array fields and calls shock finder
@everywhere function make_map_parallel2(dat1)

n1=(size(dat1[1], 1), size(dat1[2], 1)) #size(dat1[3],1))   #this builds up the 3D structure from the dat1 array which is passed to each processor

#...here below we compute the min and max in each direction of the array.
#...notice in this implementation the domain is decomposed in rectangles which are long as the entire domain in the z direction (best for map-making)

dz=0
if los == "Z"
i1=dz+dat1[1][1]
j1=dz+dat1[2][1]
i2=dz+dat1[1][n1[1]]
j2=dz+dat1[2][n1[2]]
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
    map=Array{Float64}(n1[1],n1[2],11)
    map[:]=0.

const dl=16
const nstepsl=convert(Int64,trunc((n)/dl))
@inbounds @simd for kz in 1:nstepsl

    if los == "X"
    i1=dl*(kz-1)+1
    i2=dl*(kz)-1
    end

      if los == "Y"
    j1=dl*(kz-1)+1
    j2=dl*(kz)-1
      end

      if los == "Z"
    l1=dz+dl*(kz-1)+1
    l2=dz+dl*(kz)-1
      end
    
println("processor number ", myid(), " is reading the following dataset:")
println(i1," ",i2," ",j1," ",j2," ",l1," ",l2)
#...each processor reads the input files file1,file2 in parallel, i.e. each processor only reads a part of the entire domain, depending on the bounds defined above.

    #...reading raw HDF5 ENZO files.
    #...the cd,cv,cb in front are ENZO conversion units (from code to cgs)
    
  d=cd*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
 
  t=h5read(file1,"Temperature",(i1:i2,j1:j2,l1:l2))
#  vx=cv*h5read(file2,"x-velocity",(i1:i2,j1:j2,l1:l2))
#  vy=cv*h5read(file2,"y-velocity",(i1:i2,j1:j2,l1:l2))
#  vz=cv*h5read(file2,"z-velocity",(i1:i2,j1:j2,l1:l2))

  n3=size(d)
  
  t=convert(Array{Float64,3},t)
  lx1=similar(t)
  lx2=similar(t)
  Xray_em(d,t,lx1,lx2)  
 #   macx=similar(t)
 #   macy=similar(t)
 #   macz=similar(t)
 #   mach=similar(t)
 #   mach[:]=0.
 #   macx[:]=0
 #   macz[:]=0
 #   macy[:]=0
    
#  mach=shocks(d,t,vx,vy,vz,macx,macy,macz,mach)   #...SHOCK FINDER
  
#  vx=nothing
#  vy=nothing
#  vz=nothing

    
#    tagnan=find( x->(x == NaN), mach)
#    mach[tagnan]=0.
#    tagnan=find( x->(x == Inf), mach)
#    mach[tagnan]=0.
    
#    tag0=find( x->(x >= mthr), mach)
#    nshock=size(tag0)

    
    bz=cb*h5read(file3,"Bz",(i1:i2,j1:j2,l1:l2))
    by=cb*h5read(file3,"By",(i1:i2,j1:j2,l1:l2))
    bx=cb*h5read(file3,"Bx",(i1:i2,j1:j2,l1:l2))
    
    #...synchrotron map making
    #....pradio is a vector generated only using shocked cells
    #....pradio3D is a 3D array used to generate maps (this can be done more efficiently probably)
#    nu=1.4 #...GHz  observing frequency
#    pradio=0.
#    println("1")
#    pradio3D=similar(d)
#    pradio3D[:]=0.
#   @inbounds @simd for i in eachindex(tag0)
#        ta=tag0[i]
#        d1=d[ta]
#        t1=t[ta]
#        m1=mach[ta]
#        bx1=bx[ta]
#        by1=by[ta]
#        bz1=bz[ta]
#        prad=radio(d1,t1,m1,bx1,by1,bz1,nu,pradio)
#       
#        pradio3D[ta]=prad

#     end

    
   #  bx=nothing
   #  by=nothing
   
    #...HI map making
 #    dHI=cd*h5read(file4,"HI_Density",(i1:i2,j1:j2,l1:l2))

#     dHI=d*1e-6
 #    Ts=HI_compute(d,t,dHI)
    ####################################
    
#     dm=cd*h5read(file1,"Dark_Matter_Density",(i1:i2,j1:j2,l1:l2))
#...main loop that produces the map                                             

    rm=similar(d)
    sz=similar(d)
   
 
#....compute turbulence
#  vx=cv*h5read(file2,"x-velocity",(i1:i2,j1:j2,l1:l2))                         
#  vy=cv*h5read(file2,"y-velocity",(i1:i2,j1:j2,l1:l2))                         
#  vz=cv*h5read(file2,"z-velocity",(i1:i2,j1:j2,l1:l2))    
#  sigma=similar(vx) 
#  turb2(vx,vy,vz,sigma)
   

#...delete stuff and select indices based on the LOS we are processing
if los=="Z" 
 s2=l2
 s1=l1
 b=bz
# v=vz
end

if los=="Y"
 s2=j2
 s1=j1
 b=by
# v=vy
end
if los=="X"
 s2=i2
 s1=i1
 b=bx
# v=vx
end
# vx=nothing
# vy=nothing
# vz=nothing
 bx=nothing
 by=nothing
 bz=nothing


    rm=similar(d)  #...Faraday Rotation
    sz=similar(d)  #...SZ proxy
    tlx1=similar(d) #....X-ray weighted temperature (free-free)
    tlx2=similar(d)  #.....X-ray weighted temperature (line emission)
    dlx1=similar(d)   #.....X-ray weighted density (free-free)
    dlx2=similar(d)   #....X-ray weighted density (line emission)


 for ir in eachindex(d)
     rm[ir]=dot(b[ir],d[ir])
     sz[ir]=dot(d[ir],t[ir])
#    lx[ir]=dot(d[ir]^2.,sqrt(t[ir]))   #...simple bolometric, not used       
     tlx1[ir]=dot(t[ir],lx1[ir])
     tlx2[ir]=dot(t[ir],lx2[ir])
      dlx1[ir]=dot(d[ir],lx1[ir])
      dlx2[ir]=dot(d[ir],lx2[ir])
     

 #    v[ir]=dot(v[ir],d[ir])
 #    sigma[ir]=dot(sigma[ir],d[ir])
    end


@inbounds        @simd    for i in 1:s2-s1+1 

   if los == "Z"
@views  @.       map[:,:,1]+=d[:,:,i] #av.density map                                       
@views  @.       map[:,:,2]+=t[:,:,i]  #av. temperature map                                 
@views  @.       map[:,:,3]+=abs(b[:,:,i])  #projected Bz map
@views  @.       map[:,:,4]+=rm[:,:,i]  #RM map  
#  @.       map[:,:,5]+=mach[:,:,i]    #av. Mach number 
#  @.       map[:,:,6]+=Ts[:,:,i]                                
#  @.       map[:,:,7]+=dHI[:,:,i]
#  @.       map[:,:,8]+=pradio3D[:,:,i]
#  @.       map[:,:,5]+=dm[:,:,i] 
@views  @.       map[:,:,5]+=sz[:,:,i]
@views  @.       map[:,:,6]+=lx1[:,:,i]
@views  @.       map[:,:,7]+=lx2[:,:,i]
@views  @.       map[:,:,8]+=tlx1[:,:,i]
@views  @.       map[:,:,9]+=tlx2[:,:,i]
@views  @.       map[:,:,10]+=dlx1[:,:,i]
@views  @.       map[:,:,11]+=dlx2[:,:,i]

  end

      if los == "Y"
@views   @.      map[:,:,1]+=d[:,i,:] #av.density map                                       
@views   @.      map[:,:,2]+=t[:,i,:]  #av. temperature map                                 
@views   @.      map[:,:,3]+=abs(b[:,i,:])  #projected Bz map
@views   @.      map[:,:,4]+=rm[:,i,:]  #RM map  
#   @.      map[:,:,5]+=mach[:,i,:]    #av. Mach number 
#   @.      map[:,:,6]+=Ts[:,i,:]                                
#   @.      map[:,:,7]+=dHI[:,i,:]
#   @.      map[:,:,8]+=pradio3D[:,i,:]
#   @.      map[:,:,5]+=dm[:,i,:] 
@views   @.      map[:,:,5]+=sz[:,i,:]
@views   @.      map[:,:,6]+=lx1[:,i,:]
@views   @.      map[:,:,7]+=lx2[:,i,:]
@views  @.       map[:,:,8]+=tlx1[:,i,:]
@views  @.       map[:,:,9]+=tlx2[:,i,:]
@views  @.       map[:,:,10]+=dlx1[:,i,:]
@views  @.       map[:,:,11]+=dlx2[:,i,:]
  

   end

    if los == "X"
@views    @.     map[:,:,1]+=d[i,:,:] #av.density map                                       
@views    @.     map[:,:,2]+=t[i,:,:]  #av. temperature map                                 
@views    @.     map[:,:,3]+=abs(b[i,:,:])  #projected Bz map
@views    @.     map[:,:,4]+=rm[i,:,:]  #RM map  
#    @.     map[:,:,5]+=mach[i,:,:]    #av. Mach number 
#    @.     map[:,:,6]+=Ts[i,:,:]                                
#    @.     map[:,:,7]+=dHI[i,:,:]
#    @.     map[:,:,8]+=pradio3D[i,:,:]
#    @.     map[:,:,5]+=dm[i,:,:] 
@views    @.     map[:,:,5]+=sz[i,:,:]
@views    @.     map[:,:,6]+=lx1[i,:,:]
@views    @.     map[:,:,7]+=lx2[i,:,:]
@views  @.       map[:,:,8]+=tlx1[i,:,:]
@views  @.       map[:,:,9]+=tlx2[i,:,:]
@views  @.       map[:,:,10]+=dlx1[i,:,:]
@views  @.       map[:,:,11]+=dlx2[i,:,:]
  
   end

    
end
    d=nothing
    dm=nothing
    t=nothing
    b=nothing
    Ts=nothing
    dHI=nothing
    mach=nothing
    pradio3D=nothing
    sz=nothing
    lx1=nothing
    v=nothing
    lx2=nothing 
    tlx1=nothing
    tlx2=nothing
    dlx1=nothing
    dlx2=nothing
end

    return map
    
end









######################################################
##....................MAIN PROGRAM..................#
#...basic output: a map with projected density, temperature, RM and mean projected Mach number
#...it assumes reading HDF5 files (ENZO-style), also including B-fields.
######################################################



for lll in 2:3
    if lll == 1
        @everywhere los="Z"
    end
    if lll == 2
        @everywhere los="Y"
    end
     if lll == 3
        @everywhere los="X"
    end


@everywhere using DistributedArrays
 map=DArray(make_map_parallel2,(n,n,11))  #...this builds up a Distributed Array (see Julia's documentation), i.e. an array made of chunks, each passed to a different processor. The distributed array is directly passed to the make_map function defined above

 map_total=convert(Array,map)  #....this step assembles a unique array by combining the different pieces of the distributed array 

println("The minimum and maximum density in the map are ",minimum(map_total[:,:,1])," ",maximum(map_total[:,:,1]))
println("The minimum and maximum Mach  in the map are ",minimum(map_total[:,:,4])," ",maximum(map_total[:,:,4]))

#.....total projected values along z-axis
irm=4
ib=3
itvw=2
id=1
itlw1=8
itlw2=9
idlw1=10
idlw2=11
itmw=5
ilx1=6
ilx2=7
@views    map_total[:,:,itvw]/=n  #...temperature
@views    map_total[:,:,ib]/=n #...Bz
@views    map_total[:,:,irm]*=(812.*res/(1e-6*1e-3*mu*mp))

@inbounds for j in 1:n
@simd for i in 1:n
    map_total[i,j,itlw1]=map_total[i,j,itlw1]/(map_total[i,j,ilx1]) #+map[i,j,ilx2])
    map_total[i,j,itlw2]=map_total[i,j,itlw2]/(map_total[i,j,ilx2])
    map_total[i,j,idlw1]=map_total[i,j,idlw1]/map_total[i,j,ilx1]
    map_total[i,j,idlw2]=map_total[i,j,idlw2]/map_total[i,j,ilx2]
    map_total[i,j,itmw]=map_total[i,j,itmw]/map_total[i,j,id]
#    map_total[i,j,iv]/=map_total[i,j,id]
#    map_total[i,j,isv]/=map_total[i,j,id]
 end
 end
@views    map_total[:,:,id]/=n  #...density 

#    map_total[:,:,5]/=n
 #   map_total[:,:,5]+=1e-10  #...Mach number
 #   map_total[:,:,6]/=n  #...HI  brightness temperature
 #   map_total[:,:,7]/=n    #..HI density
 #   map_total[:,:,8]+=1e-30   #...radio map
 #   map_total[:,:,9]/=n

#....writing FITS files with images
@everywhere using FITSIO

filep1=string(root,"map_all",mod,"_",snap,"_RM_X_weight1024",los,".fits")
 f = FITS(filep1, "w");
write(f,map_total)
close(f)

end


