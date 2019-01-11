# SIMULATION PARAMETERS
@everywhere const h0=1.0
@everywhere const Lbox=0.1 #Mpc
@everywhere const n=120 #grid size                                              
@everywhere const dx=0.195   #kpc
@everywhere const root="/Users/francovazza/Desktop/data/DATA/BOX2018/512/weakC/"

#CONSTANTS
@everywhere const kpc=3.086e21 #cm in one kpc
@everywhere const yr=3.154e7 #s in one yr
                                                                                                                              
@everywhere const  mp=1.67e-24
@everywhere const  me=0.1093897e-28   #g                                                                                          
@everywhere const  vc=2.99792458e10   #cm s-                                  
@everywhere const  kpctocm=3.086e21
@everywhere const  cmtoMpc=3.086e24
@everywhere const  kb=1.380658e-16    #erg k-1                                                
#COSMOLOGY PARAMETERS
@everywhere const ΩΛ=0.742
@everywhere const ΩM=0.258
@everywhere const Ωb=0.0441



#MODULES
#@everywhere using HDF5
@everywhere using FITSIO
#@everywhere using PyPlot
@everywhere using Optim
@everywhere using Base
#@everywhere using Devectorize

#first we define all functions that are called by the main algorithm (down below). Making them @everywhere defines them on all processors.

#...DOMAIN EDGES and TRACERS NUMBER
@everywhere i1=1
@everywhere j1=1
@everywhere l1=1
@everywhere i2=n
@everywhere j2=n
@everywhere l2=n
@everywhere np=40000
@everywhere model="test40000"
@everywhere ng1=i2-i1+1
@everywhere ng2=j2-j1+1
@everywhere ng3=l2-l1+1
#...SPECTRAL PARAMETERS

@everywhere const  part=2  #...1=proton  2=electron                                                                                             
#@everywhere const  volume=3*log10(vol_trac*kpctocm) #vol im cm^3case of                                                                        
@everywhere const norm=30 #to have all in 1e30 erg or 1e30 erg/s                                                                               
@everywhere const erest=0.510988946 #electron rest mass energy                                                                                 
@everywhere const evtoerg=1.60218e-12 #                                                                                                        
@everywhere const  b1=1.37e-20 #...1/s                                                                                                         
@everywhere const  b2=1.9e-9  #...1/s                                                                                             
@everywhere const  b3=7.2e-12 #...1/s                                                                                             
@everywhere const xi=0.5 # acceleration efficiency               
@everywhere const lcorr=20.0 #...in[100]kpc                                                                                                   
@everywhere const mthr=1.3
@everywhere const gyrtosec=1.0e9*3.0e7
@everywhere const fi=1e-3 #....


@everywhere include("/Users/francovazza/Desktop/Julia_prog/FOKKER/bo/spectra_lab/loss_gain.jl")
@everywhere include("/Users/francovazza/Desktop/Julia_prog/FOKKER/bo/spectra_lab/advect_assign.jl")
@everywhere include("/Users/francovazza/Desktop/Julia_prog/FOKKER/bo/spectra_lab/map_making.jl")


@everywhere const ngamma=(floor(Int64,(g_max-g_min)/dg))

@everywhere const gend=ngamma
@everywhere const gammaval=collect(g_min:dg:g_max)


##############
#MAIN PROGRAM#
##############
tic()


dt=0.0027374075593099984*1.543e+15/3e7/10  #yr .....fixed timestep. This is changed during the simulation after reading the redshift value for each snap-

sna=7
snap_name=string(sna)
snap_name="0224"
mod="214A"
fil="fits"
#...INFO dataset
#file1=string(root,"DD0A0_dt_",snap_name)    #input files    *dt* contains gas density, dark matter \density and gas temperature + 3 mag.field components                                                
#file2=string(root,"DD0A0_v_",snap_name)   # *v* contains the 3-d components of the v-field          
 #...CHRONOS dataset
 #file1=string(root,"s",snap_name,"/full_dtb_R_",snap_name)
 #file2=string(root,"s",snap_name,"/full_v_R_",snap_name)
 if fil == "hdf5"
 file1=string(root,"",mod,"_dtb_",snap_name,".dms")
 file2=string(root,"",mod,"_v_",snap_name,".dms")
 file3=string(root,"",mod,"_dtb_",snap_name,".dms")
 end

 if fil == "fits"
 file1=string(root,"dens",mod,".fits")
 file2=string(root,"temp",mod,".fits")
 file3=string(root,"vx",mod,".fits")
 file4=string(root,"vy",mod,".fits")
 file5=string(root,"vz",mod,".fits")
 file6=string(root,"bx",mod,".fits")
 file7=string(root,"by",mod,".fits")
 file8=string(root,"bz",mod,".fits")
 end

# file_conv=string(root,"conv/"DD0A0",snap_name,".conv2")  #file with conversion factors from enzo internal\ to cgs                                                                                            
#  file_conv=string(root,"M_512B",snap_name,".conv2")  #file with conversion factors from enzo internal\
println("assigning tracers")

#p=assign_tracers()
#p=assign_tracers_loc()
p=assign_tracers_shock()
#f=1e-2*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
#p=assign_tracers_dens(f)
#f=nothing
@everywhere tini=110
@everywhere tend=196
@everywhere tt=0
pe=Array{Float64}(ngamma,np)
pe=fill(0.0,ngamma,np)
@everywhere using Plots
pyplot()
@everywhere using LaTeXStrings


map=Array{Float64}(i2-i1+1,j2-j1+1)
map[:]=0.

for tt in tini:tend

tic()
println(tt)
step=string(tt)

if tt <10
step=string("00",tt)
end
if tt >=10 && tt <100
step=string("0",tt)
end


snap_name=step
mod=step
#file1a=string(file1,snap_name)
#file2a=string(file2,snap_name)
#file_conva=string(file_conv,snap_name,".conv2")
#a=readdlm(file_conv)
z=0. #;a[3]
cd=1.6699997e-28
cv=2e+08

 if fil == "fits"
 file1=string(root,"dens",mod,"B.fits")
 file2=string(root,"temp",mod,"B.fits")
 file3=string(root,"vx",mod,"B.fits")
 file4=string(root,"vy",mod,"B.fits")
 file5=string(root,"vz",mod,"B.fits")
 file6=string(root,"bx",mod,"B.fits")
 file7=string(root,"by",mod,"B.fits")
 file8=string(root,"bz",mod,"B.fits")
 end


#t1=calc_time()

#dt=

cb=1e6*sqrt(cd/(1+z)^3.*4*π)*cv #to have microGauss   

#println("processor number ", myid(), " is reading the following dataset:")
#println(i1," ",i2," ",j1," ",j2," ",l1," ",l2)
#if t == 1 
println(file1)

 f1=FITS(file1,"r")
 f2=FITS(file2,"r")
 f3=FITS(file3,"r")
 f4=FITS(file4,"r")
 f5=FITS(file5,"r")
 f6=FITS(file6,"r")
 f7=FITS(file7,"r")
 f8=FITS(file8,"r")

  d=read(f1[1])*cd
  t=read(f2[1])
  vx=read(f3[1])*cv
  vy=read(f4[1])*cv
  vz=read(f5[1])*cv
  bx=read(f6[1])*cb
  by=read(f7[1])*cb
  bz=read(f8[1])*cb
 
#d=cd*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
#vz=cv*h5read(file2,"z-velocity",(i1:i2,j1:j2,l1:l2))
#vy=cv*h5read(file2,"y-velocity",(i1:i2,j1:j2,l1:l2))
#vx=cv*h5read(file2,"x-velocity",(i1:i2,j1:j2,l1:l2))
#t=h5read(file1,"Temperature",(i1:i2,j1:j2,l1:l2))
#bx=cb*h5read(file3,"Bx",(i1:i2,j1:j2,l1:l2))
#by=cb*h5read(file3,"By",(i1:i2,j1:j2,l1:l2))
#bz=cb*h5read(file3,"Bz",(i1:i2,j1:j2,l1:l2))

t=convert(Array{Float64},t)
#if tt==tini 
#map=make_map(d,t,vx,vy,vz,bx,by,bz)
#imagesc((map[:,:,1]))
#file_map=string(file1,"_densX.png")
#savefig(file_map)
#end

println("assign fields")
#tic()

p=assign_tracers_fields(p,d,t,vx,vy,vz,bx,by,bz)
#p[5,:]+=(-2e8)

  map=make_mapz(convert(Array{Float32},d),convert(Array{Float32},t))

d=nothing
t=nothing

#toc()

println("move")


p=move_tracers(p,dt,dx,i1,i2,j1,j2,l1,l2)
#toc()


  p0=contour((map[:,:,2]),aspect_ratio=1,fill=true,nlevels=256,size=(1200,1200))#,c=:hsv)            
#  title!(string("z=",z))

#  scatter!(p[1,:]-i1,p[2,:]-j1,line="",label="all tracers",size=0.2)                                       
#  @inbounds for i in 1:10                                                                                  
 # scatter!(p[1,1:20]-i1,p[2,1:20]-j1,label="tracers subsample")                                            

  scatter!(p[2,:]-j1,p[1,:]-i1,label="",color="red",size=(1200,1200))
   plot(p0)
   filep1=string(root,"Bspectra_yes_reacc",tt,"_testA.png")
   savefig(filep1)
# 
 lsize=dx #kpc linear size of tracer
 
file_trac=string(file1,"_tracers_",step,"_",model,"B.dat")
 so=open(file_trac,"w")
    for pp in eachindex(p)
        write(so,convert(Float32,p[pp]))
    end
    close(so)

#[p[1,:]  p[1,:] p[1,:] p[1,:] p[1,:] p[1,:] p[1,:] p[1,:] p[1,:]
#p[1,][dm_los[:,1] rm_los[:,1] bdev[:,1] bmean[:,1] bmeanabs[:,1] bmed\
#[:,1]])

toc()


    end

