#.....JULIA SHOCK FINDER FOR UNIGRID DATA
#.....INPUT DATA ARE SUPPOSED TO BE REGULAR GRID, WRITTEN IN HDF5 FORMAT
#.....THE CODE CAN ALSO RUN IN PARALLEL JUST BY ADDING PROCS


@everywhere main="/home/PERSONALE/franco.vazza2/Dropbox/Julia/"
#@everywhere main="/Users/francovazza/Dropbox/Julia/"
#....files to read in
#@everywhere root="/Users/francovazza/Desktop/data/DATA/CHRONOS/test_dataset/"
@everywhere root="/home/PERSONALE/franco.vazza2/Desktop/DATA/CHRONOS++/85Mpc/1024/"

@everywhere mod="newBH12b_csf3dR"
@everywhere snap="007"


@everywhere file1=string(root,"s",snap,"/full_dtb_",mod,"_",snap)    #input files    *dt* contains gas density, dark matter density and gas temperature
@everywhere file2=string(root,"s",snap,"/full_v_",mod,"_",snap)   # *b* contains the 3-d components of the v-field,
@everywhere file3=string(root,"s",snap,"/full_dtb_",mod,"_",snap)   # *b* contains the 3-d components of the B-field,


@everywhere file4=string(root,"s",snap,"/full_chem_",mod,"_",snap)
#@everywhere file_conv=string(root,"conv/",mod,snap,".conv2")  #file with conversion factors from Enzo to cgs


#@everywhere conv=readdlm(file_conv,'\t',Float64,'\n')
#@everywhere const cd=conv[4]  #...from code density to g/cm^3                                                                              
#@everywhere const cv=conv[5]  #...from code velocity to cm/s                                                                               
#@everywhere const cb=sqrt(cd*4*π)*cv  #...from code b-field to Gauss                                                                      
#@everywhere const z=conv[3]  #...redshift                                                                                                  

@everywhere cd=2.822e-30
@everywhere cv=2.5167e9
@everywhere cb=sqrt(cd*4*π)*cv

@everywhere include(string(main,"constants.jl"))
@everywhere include(string(main,"cosmological_param.jl"))
@everywhere include(string(main,"shocks.jl"))
@everywhere include(string(main,"HI.jl"))
@everywhere include(string(main,"radio_power.jl"))


#@everywhere using DistributedArrays    #necessary for parallel map making
@everywhere using HDF5    #necessary to open hdf5 files	


#.....HAVING THE APPROPRIATE CONVERSIONS IS CRUCIAL HERE, as Mach=v/sound_speed

#@everywhere const cd= 8.48e-30  #physical density 
#@everywhere const cv= 1.73e+09    #physical velocity
#@everywhere const cb=sqrt(cd*4*pi)*cv   #physical B-fild Gauss

@everywhere const mthr=3.  #...threshold for Mach number estimation

@everywhere function make_map_serial(n,kz,dz)

i1=1
j1=1
i2=n
j2=n
l1=dz*(kz-1)+1
l2=dz*(kz)-1

println(i1," ",i2," ",j1," ",j2," ",l1," ",l2)

    #...reading raw HDF5 ENZO files.                                            
    #...the cd,cv,cb in front are ENZO conversion units (from code to cgs)      

  d=cd*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
#  dm=cd*h5read(file1,"Dark_Matter_Density",(i1:i2,j1:j2,l1:l2))
  t=h5read(file1,"Temperature",(i1:i2,j1:j2,l1:l2))
  vx=cv*h5read(file2,"x-velocity",(i1:i2,j1:j2,l1:l2))
  vy=cv*h5read(file2,"y-velocity",(i1:i2,j1:j2,l1:l2))
  vz=cv*h5read(file2,"z-velocity",(i1:i2,j1:j2,l1:l2))
  bz=cb*h5read(file3,"Bz",(i1:i2,j1:j2,l1:l2))
#  by=cb*h5read(file3,"By",(i1:i2,j1:j2,l1:l2))
#  bx=cb*h5read(file3,"Bx",(i1:i2,j1:j2,l1:l2))
#  dHI=cd*h5read(file4,"HI_Density",(i1:i2,j1:j2,l1:l2))
  t=convert(Array{Float64,3},t)
#  Ts=HI_compute(d,t,dHI)
    macx=similar(t)                                                                                                        
    macy=similar(t)                                                                                                        
    macz=similar(t)                                                                                                        
    mach=similar(t)                                                                                                        
    mach[:]=0.                                                                                                             
    macx[:]=0                                                                                                              
    macz[:]=0                                                                                                              
    macy[:]=0                                                                                                              


 mach=shocks(d,t,vx,vy,vz,macx,macy,macz,mach)

# mach=shocks(d,t,vx,vy,vz)#,bx,by,bz)   #...SHOCK FINDER   
 vx=nothing
 vy=nothing 
 vz=nothing 

# dHI=cd*h5read(file4,"HI_Density",(i1:i2,j1:j2,l1:l2))
 dHI=d*1e-6
 Ts=HI_compute(d,t,dHI)
#...this is going to be the final map                                           
 map=Array{Float64}(i2-i1+1,j2-j1+1,6)
 map[:]=0.
#...main loop that produces the map                                             

@inbounds        @simd    for i in 1:l2-l1+1 
  @views          map[:,:,1]+=d[:,:,i] #av.density map                                       
  @views          map[:,:,2]+=t[:,:,i]  #av. temperature map                                 
  @views          map[:,:,3]+=abs(bz[:,:,i])  #RM map                                             
  @views          map[:,:,4]+=mach[:,:,i]    #av. Mach number 
  @views          map[:,:,5]+=Ts[:,:,i]                                
  @views          map[:,:,6]+=dHI[:,:,i]
end
    d=nothing
    t=nothing
    bz=nothing
    Ts=nothing
    dHI=nothing
    mach=nothing
 return map
end


###############################################
#.......Reads in array fields and calls shock finder
@everywhere function make_map_parallel(dat1)

n1=(size(dat1[1], 1), size(dat1[2], 1)) #size(dat1[3],1))   #this builds up the 3D structure from the dat1 array which is passed to each processor

#...here below we compute the min and max in each direction of the array.
#...notice in this implementation the domain is decomposed in rectangles which are long as the entire domain in the z direction (best for map-making)

i1=dat1[1][1]
j1=dat1[2][1]
i2=dat1[1][n1[1]]
j2=dat1[2][n1[2]]
l1=1
l2=n

println("processor number ", myid(), " is reading the following dataset:")
println(i1," ",i2," ",j1," ",j2," ",l1," ",l2)
#...each processor reads the input files file1,file2 in parallel, i.e. each processor only reads a part of the entire domain, depending on the bounds defined above.

    #...reading raw HDF5 ENZO files.
    #...the cd,cv,cb in front are ENZO conversion units (from code to cgs)
    
  d=cd*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
#  dm=cd*h5read(file1,"Dark_Matter_Density",(i1:i2,j1:j2,l1:l2))
  t=h5read(file1,"Temperature",(i1:i2,j1:j2,l1:l2))
#  vx=cv*h5read(file2,"x-velocity",(i1:i2,j1:j2,l1:l2))
#  vy=cv*h5read(file2,"y-velocity",(i1:i2,j1:j2,l1:l2))
#  vz=cv*h5read(file2,"z-velocity",(i1:i2,j1:j2,l1:l2))
#  bz=cb*h5read(file3,"Bz",(i1:i2,j1:j2,l1:l2))
#  bx=bz
#  by=bz
#  by=cb*h5read(file3,"By",(i1:i2,j1:j2,l1:l2))
#  bx=cb*h5read(file3,"Bx",(i1:i2,j1:j2,l1:l2))
#  dHI=cd*h5read(file4,"HI_Density",(i1:i2,j1:j2,l1:l2))
  t=convert(Array{Float64,3},t)

#    macx=similar(t)
#    macy=similar(t)
#    macz=similar(t)
#    mach=similar(t)
#    mach[:]=0.
#    macx[:]=0
#    macz[:]=0
#    macy[:]=0
    
#  mach=shocks(d,t,vx,vy,vz,macx,macy,macz,mach)   #...SHOCK FINDER
#  vx=nothing
#  vy=nothing
#  vz=nothing
  
     dHI=cd*h5read(file4,"HI_Density",(i1:i2,j1:j2,l1:l2))

     Ts=HI_compute(d,t,dHI)
 

     

    
    
#...this is going to be the final map
 map=Array{Float64}(n1[1],n1[2],6)
 map[:]=0.
#...main loop that produces the map 
 bz=cb*h5read(file3,"Bz",(i1:i2,j1:j2,l1:l2))

    
#...main loop that produces the map                                             

@inbounds        @simd    for i in 1:l2-l1+1 
            map[:,:,1]+=d[:,:,i] #av.density map                                       
            map[:,:,2]+=t[:,:,i]  #av. temperature map                                 
            map[:,:,3]+=abs(bz[:,:,i])  #RM map                                             
#            map[:,:,4]+=mach[:,:,i]    #av. Mach number 
            map[:,:,5]+=Ts[:,:,i]                                
            map[:,:,6]+=dHI[:,:,i]
end
    d=nothing
    t=nothing
    bz=nothing
    Ts=nothing
    dHI=nothing
    mach=nothing
    return map
    
end


###############################################
#.......Reads in array fields and calls shock finder
@everywhere function make_map_parallel2(dat1)

n1=(size(dat1[1], 1), size(dat1[2], 1)) #size(dat1[3],1))   #this builds up the 3D structure from the dat1 array which is passed to each processor

#...here below we compute the min and max in each direction of the array.
#...notice in this implementation the domain is decomposed in rectangles which are long as the entire domain in the z direction (best for map-making)

i1=dat1[1][1]
j1=dat1[2][1]
i2=dat1[1][n1[1]]
j2=dat1[2][n1[2]]
    map=Array{Float64}(n1[1],n1[2],7)
    map[:]=0.

const dl=64
const nstepsl=convert(Int64,trunc((n)/dl))
@inbounds @simd for kz in 1:nstepsl
l1=dl*(kz-1)+1
l2=dl*(kz)-1

println("processor number ", myid(), " is reading the following dataset:")
println(i1," ",i2," ",j1," ",j2," ",l1," ",l2)
#...each processor reads the input files file1,file2 in parallel, i.e. each processor only reads a part of the entire domain, depending on the bounds defined above.

    #...reading raw HDF5 ENZO files.
    #...the cd,cv,cb in front are ENZO conversion units (from code to cgs)
    
  d=cd*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
#  dm=cd*h5read(file1,"Dark_Matter_Density",(i1:i2,j1:j2,l1:l2))
  t=h5read(file1,"Temperature",(i1:i2,j1:j2,l1:l2))
  vx=cv*h5read(file2,"x-velocity",(i1:i2,j1:j2,l1:l2))
  vy=cv*h5read(file2,"y-velocity",(i1:i2,j1:j2,l1:l2))
  vz=cv*h5read(file2,"z-velocity",(i1:i2,j1:j2,l1:l2))
#  bz=cb*h5read(file3,"Bz",(i1:i2,j1:j2,l1:l2))

  n3=size(d)
    
#  dHI=cd*h5read(file4,"HI_Density",(i1:i2,j1:j2,l1:l2))
  t=convert(Array{Float64,3},t)

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
    println("1")
    pradio3D=similar(d)
    pradio3D[:]=0.
   @inbounds @simd for i in eachindex(tag0)
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

    
     bx=nothing
     by=nothing

    #...HI map making
     dHI=cd*h5read(file4,"HI_Density",(i1:i2,j1:j2,l1:l2))
    
     Ts=HI_compute(d,t,dHI)
    ####################################
    
    
#...main loop that produces the map                                             

@inbounds        @simd    for i in 1:l2-l1+1 
            map[:,:,1]+=d[:,:,i] #av.density map                                       
            map[:,:,2]+=t[:,:,i]  #av. temperature map                                 
            map[:,:,3]+=abs(bz[:,:,i])  #RM map                                             
            map[:,:,4]+=mach[:,:,i]    #av. Mach number 
            map[:,:,5]+=Ts[:,:,i]                                
            map[:,:,6]+=dHI[:,:,i]
            map[:,:,7]+=pradio3D[:,:,i]
end
    d=nothing
    t=nothing
    bz=nothing
    Ts=nothing
    dHI=nothing
    mach=nothing
    pradio3D=nothing

end

    return map
    
end









######################################################
##....................MAIN PROGRAM..................#
#...basic output: a map with projected density, temperature, RM and mean projected Mach number
#...it assumes reading HDF5 files (ENZO-style), also including B-fields.
######################################################

tic()   
what="parallel"

if what == "serial"
const dz=32
const nsteps=convert(Int64,trunc(n/dz))
map_total=Array{Float64}(n,n,6)
map_total[:]=0.
for kz in 1:nsteps
tic()
map_total0=make_map_serial(n,kz,dz)
print(size(map_total0))
map_total+=map_total0
toc()
end
end

if what == "parallel"
@everywhere using DistributedArrays
 map=DArray(make_map_parallel2,(n,n,7))  #...this builds up a Distributed Array (see Julia's documentation), i.e. an array made of chunks, each passed to a different processor. The distributed array is directly passed to the make_map function defined above

 map_total=convert(Array,map)  #....this step assembles a unique array by combining the different pieces of the distributed array 

end
println("The minimum and maximum density in the map are ",minimum(map_total[:,:,1])," ",maximum(map_total[:,:,1]))
println("The minimum and maximum Mach  in the map are ",minimum(map_total[:,:,4])," ",maximum(map_total[:,:,4]))

#.....total projected values along z-axis
    ma1=(map_total[:,:,1])/n  #...density
    ma2=(map_total[:,:,2])/n  #...temperature
    ma3=(map_total[:,:,3])    #...RM
    ma4=(map_total[:,:,4])/n+1e-10  #...Mach number
    ma5=(map_total[:,:,5])/n  #...HI  brightness temperature
    ma6=map_total[:,:,6]/n    #..HI density
    ma7=map_total[:,:,7]+1e-30   #...radio map
toc()

#....writing FITS files with images
@everywhere using FITSIO

filep1=string(root,"map_dens",mod,"_",snap,"_1024.fits")
 f = FITS(filep1, "w");
 write(f,ma1)
 close(f)

 filep2=string(root,"map_temp",mod,"_",snap,"_1024.fits")
 f = FITS(filep2, "w");
 write(f,ma2)
 close(f)

filep5=string(root,"map_dT2HI",mod,"_",snap,"_1024.fits")
 f = FITS(filep5, "w");
 write(f,ma5)
 close(f)

 filep3=string(root,"map_Bz",mod,"_",snap,"_1024.fits")
 f = FITS(filep3, "w");
 write(f,ma3)
 close(f)

filep4=string(root,"map_Mach",mod,"_",snap,"_1024_clean.fits")
 f = FITS(filep4, "w");
 write(f,ma4)
 close(f)

filep6=string(root,"map_densHI",mod,"_",snap,"_1024.fits")
 f = FITS(filep6, "w");
 write(f,ma6)
 close(f)


filep7=string(root,"radio_power",mod,"_",snap,"_1024.fits")
 f = FITS(filep7, "w");
 write(f,log10(ma7))
 close(f)





