
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
#MODULES                                                                                                                           
#@everywhere using PyPlot                                                                                                          
@everywhere using FITSIO
@everywhere using Optim
@everywhere using Base
#@everywhere using Devectorize                                                                                                     

#...DOMAIN EDGES and TRACERS NUMBER                                                                                                
@everywhere i1=1
@everywhere j1=1
@everywhere l1=1
@everywhere i2=120
@everywhere j2=120
@everywhere l2=120
@everywhere nptot=40000
@everywhere np=40000
@everywhere model="test40000"
@everywhere ng1=i2-i1+1
@everywhere ng2=j2-j1+1
@everywhere ng3=l2-l1+1
#...SPECTRAL PARAMETERS                                                                                                            
#@everywhere const  g_max=50000.0
#@everywhere const  g_min=100.0
#@everywhere const  dg=100.0
@everywhere const  part=2  #...1=proton  2=electron                                     
#..constants useful to compute cooling and acceleration terms
@everywhere const norm=30 #to have all in 1e30 erg or 1e30 erg/s                                        
@everywhere const erest=0.510988946 #electron rest mass energy                                           
@everywhere const evtoerg=1.60218e-12 #                                                                  
@everywhere const  b1=1.37e-20 #...1/s                                                                   
@everywhere const  b2=1.9e-9  #...1/s                                                                    
@everywhere const  b3=7.2e-12 #...1/s                                                                                                          
@everywhere const xi=0.5 #                                                                                                                     
@everywhere const lcorr=20.0 #...in[100]kpc                                                                                                   
@everywhere const mthr=3.0
@everywhere const gyrtosec=1.0e9*3.0e7
@everywhere const fi=1e-3
#@everywhere const ngamma=(floor(Int64,(g_max-g_min)/dg))

#@everywhere const gend=ngamma
#@everywhere const gammaval=collect(g_min:dg:g_max)
@everywhere const ptot=Array{Float64}(10,nptot)
@everywhere const strategy="no"

@everywhere include("/Users/francovazza/Desktop/Julia_prog/FOKKER/bo/spectra_lab/loss_gain_speed2.jl")
@everywhere include("/Users/francovazza/Desktop/Julia_prog/FOKKER/bo/spectra_lab/advect_assign.jl")
@everywhere include("/Users/francovazza/Desktop/Julia_prog/FOKKER/bo/spectra_lab/map_making.jl")

tic()

# semilogx(gammaval,gammaval*0.0)
# axis([g_min,g_max,1,1e20])


dt=0.0027374075593099984*1.543e+15/3e7/10  #yr
sna=7
snap_name="0224"
mod="214A"
fil="fits"
#...INFO dataset
#file1=string(root,"DD0A0_dt_",snap_name)    #input files    *dt* contains gas density, dark matter \density and gas temperature + 3 mag.field components                                                
#file2=string(root,"DD0A0_v_",snap_name)   # *v* contains the 3-d components of the v-field          
 #...CHRONOS dataset
 #file1=string(root,"s",snap_name,"/full_dtb_R_",snap_name)
 #file2=string(root,"s",snap_name,"/full_v_R_",snap_name)
# file_conv=string(root,"conv/"DD0A0",snap_name,".conv2")  #file with conversion factors from enzo internal\ to cgs                                                                                            
#  file_conv=string(root,mod,"D0",snap_name,".conv2")  #file with conversion factors from enzo internal\
#println("assigning tracers")

#p=assign_tracers()
#p=assign_tracers_loc()
#f=1e-2*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
#p=assign_tracers_dens(f)
#f=nothing
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
@everywhere tini=110
@everywhere tend=196
@everywhere tt=0
#p=Array{Float64}(10,np)
#pe=Array{Float64}(ngamma,np)
#pe=fill(0.0000,ngamma,np)
@everywhere using Plots
pyplot()
@everywhere using LaTeXStrings


map=Array{Float64}(i2-i1+1,j2-j1+1)
map[:]=0.

read_trac="y"

#sor1=Array{Int64}(np)

#......MAIN LOOP WHERE THE EVOLUTION IS COMPUTED

#sor1=Array{Int64,np}
isel=Array{Int64,1}
time0=0.
for tt in tini:tend
time0+=dt
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

z=0.1 #;a[3]                                                                                                                                      
cd=1.6699997e-28
cv=2e+08

cb=1e6*sqrt(cd/(1+z)^3.*4*π)*cv #to have microGauss                                                                                              

#println("processor number ", myid(), " is reading the following dataset:")
#println(i1," ",i2," ",j1," ",j2," ",l1," ",l2)
#if t == 1 
if read_trac == "n"
d=cd*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
vz=cv*h5read(file2,"z-velocity",(i1:i2,j1:j2,l1:l2))
vy=cv*h5read(file2,"y-velocity",(i1:i2,j1:j2,l1:l2))
vx=cv*h5read(file2,"x-velocity",(i1:i2,j1:j2,l1:l2))
t=h5read(file1,"Temperature",(i1:i2,j1:j2,l1:l2))
bx=cb*h5read(file3,"Bx",(i1:i2,j1:j2,l1:l2))
by=cb*h5read(file3,"By",(i1:i2,j1:j2,l1:l2))
bz=cb*h5read(file3,"Bz",(i1:i2,j1:j2,l1:l2))

if tt==tini 
println("assign fields")

p=assign_tracers_fields(p,d,t,vx,vy,vz,bx,by,bz)
d=nothing
t=nothing
pe=Array{Float64}(ngamma,np)
 pe[:,:]=1e-20

println("move")
p=move_tracers(p,dt,dx,i1,i2,j1,j2,l1,l2)
end
end

if read_trac=="y"
if fil == "fits"
 file1=string(root,"temp",mod,"B.fits")
 filei=string(root,"dens",mod,"B.fits")
   file2=string(root,"bz",mod,"B.fits")
   file3=string(root,"bx",mod,"B.fits")
   file4=string(root,"by",mod,"B.fits")

 end
  f1=FITS(file1,"r")
  f2=FITS(file2,"r")
  f3=FITS(file3,"r")
  f4=FITS(file4,"r")
   bz=read(f2[1])*cb
   bx=read(f3[1])*cb
   by=read(f4[1])*cb

  @simd  for j in eachindex(bz)
@inbounds @fastmath  bz[j]=sqrt(bx[j]^2.+by[j]^2.+bz[j]^2.)
  end
  bx=nothing
  by=nothing

#  d=read(f1[1])*cd

 file_trac=string(filei,"_tracers_",step,"_",model,"B.dat")
  so=open(file_trac,"r")
  p=Array{Float64}(10,np)
  println(np)

 @inbounds @simd  for pp in eachindex(p)
 pa=read(so,Float32)
 p[pp]=convert(Float64,pa)
  end
@views @fastmath  p[10,:]=p[10,:]*(1+z)^2.  #...physical B-field

 map=make_mapz1(bz) #convert(Array{Float32},d),convert(Array{Float32},bz))
 
if tt==tini
 pe=Array{Float64}(ngamma,np)
 pe[:,:]=1e-30
end
end

#   tic()

 zz=convert(Float64,z)
 lsize=dx #kpc linear size of tracer                                                                                    
println(size(p))
println(size(pe))
    println("spectra")

if tt> tini
 tic()  
 end
    pe=spectra(p,pe,dt*1e-9,zz,lsize,tt)
if tt > tini
    toc()

end
    tic()

   pe_tot=Array{Float64}(ngamma)
   pe_tot=fill(1e-20,ngamma)
   gal=gammaval
 

    if (tt >= tini)    


    dmin=0.0
    dmax=0.5
    id=find(x-> (x <=dmin ),map)
    map[id]=dmin
    id=find(x-> (x >=dmax ),map)
     map[id]=dmax
#    ma=log10(map)


  p0=contour((map),aspect_ratio=1,fill=true,nlevels=256,size=(1000,1000))#,nlevels=10)   #c=:heat
  title!(string(L"|B|[$\mu G$] , t=",trunc(time0*1e-3)," kyr"))
 
  radio_trac=Array{Float64}(2,np)
  radio_trac[:]=0.
  freqf=[3e8,5.5e8,6.1e8,8.5e8,1.4e9,3e9]

   map_radio=Array{Float64}(n,n,3)
   map_radio[:]=1e-30
   map_radio[:,:,3]=3.
 

   @inbounds @simd for pp in 1:np
  

@inbounds  @simd for g in 1:ngamma
@.  pe_tot[g]=pe_tot[g]+pe[g,pp]
  end

#...synchrotron loop
   if pe[2,pp] >=1e30
    #println(p[:,pp])                                                                                                     
    #println(pe[:,pp])                                                                                                    
   spec=Vector{Float64}(6)
   spec[:]=0.

   freq=freqf[3]
   sync=synchrotron(gammaval,p[10,pp],pe[:,pp],freq)
   spec[3]=sync
   freq=freqf[5]
   sync=synchrotron(gammaval,p[10,pp],pe[:,pp],freq)
   spec[5]=sync

   y1=convert(Int32,trunc(p[2,pp]-j1))
   x1=convert(Int32,trunc(p[1,pp]-i1))
   ds=3
   if x1 <=ds+1
   x1=ds+1
   end

   if y1 <=ds+1
   y1=ds+1
   end

    if x1 >=(n-ds-1)
   x1=n-ds
   end

   if y1 >= (n-ds-1)
   y1=n-ds
   end

   for j1 in -ds:ds
@simd   for i1 in -ds:ds
@inbounds   map_radio[x1+i1,y1+j1,1]+=(spec[3]/(1.+abs(ds)))
@inbounds   map_radio[x1+i1,y1+j1,2]+=(spec[5]/(1.+abs(ds)))
   end
   end
   end  #...synchrotron loop

   end
   end

   for j in 1:n
@simd   for i in 1:n
@inbounds @views  map_radio[i,j,3]=(log10(map_radio[i,j,1])-log10(map_radio[i,j,2])) /(log10(freqf[5])-log10(freqf[3]))
   end
   end

    dmax=6e23
    dmin=1e18
    ma=map_radio[:,:,1]
    id=find(x-> (x <=dmin ),ma)
    ma[id]=dmin
    id=find(x-> (x >=dmax ),ma)
    ma[id]=dmax

  p1=contour(ma,aspect_ratio=1,fill=true,nlevels=256,size=(1000,1000))#,nlevels=10)   #c=:heat    
  title!(string("Radio Power, t=",trunc(time0*1e-3)," kyr"))

    dmin=0.5
    dmax=1.1
    ma=map_radio[:,:,3]
    id=find(x-> (x <=dmin ),ma)
    ma[id]=dmin
    id=find(x-> (x >=dmax ),ma)
    ma[id]=dmax



  p2=contour(ma,aspect_ratio=1,fill=true,nlevels=256,size=(1000,1000))#,nlevels=10)\
   #c=:heat                                                                                         
  title!(string("Radio Spectral Index, t=",trunc(time0*1e-3)," kyr"))

    p4=plot(gal,pe_tot,line=:solid,color="black",label="energy spectrum",size=(1000,1000))
    yaxis!(L"N(\gamma)",:log10,(1e38,1e51))
    xaxis!(L"\gamma",:log10,(g_min+dg,g_max))
  
  plot(p0,p1,p2,p4)
 filep1=string(root,"Bspectra_yes_reacc",tt,"_testBB_new.png")
  savefig(filep1)
#  error()
#file_trac1=string(file1,"_tracers_",step,"_",model,"_specB.dat")
# so=open(file_trac1,"w")
#    for pp in 1:#np
#    for g in 1:#ngamma
#        write(so,convert(Float32,pe[g,pp]))
#    #end
    #end
   # close(so)
  
    end
end
println("total computing time")
toc()