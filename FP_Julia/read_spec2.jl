
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
@everywhere using FITSIO
@everywhere using Optim
@everywhere using Base
@everywhere using Plots
#@everywhere using Devectorize                                                                                                     

#...DOMAIN EDGES and TRACERS NUMBER                                                                                                
@everywhere i1=1
@everywhere j1=1
@everywhere l1=1
@everywhere i2=120
@everywhere j2=120
@everywhere l2=120
@everywhere nptot=8000
@everywhere np=8000
@everywhere model="test40000"
@everywhere ng1=i2-i1+1
@everywhere ng2=j2-j1+1
@everywhere ng3=l2-l1+1
#...SPECTRAL PARAMETERS                                                                                                            
@everywhere const  g_max=50000.0
@everywhere const  g_min=100.0
@everywhere const  dg=100.0

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
@everywhere const mthr=2.0
@everywhere const gyrtosec=1.0e9*3.0e7
@everywhere const fi=1e-3
@everywhere const ngamma=(floor(Int64,(g_max-g_min)/dg))

@everywhere const gend=ngamma
const gammaval=collect(g_min:dg:g_max)
const ptot=Array{Float64}(10,nptot)
@everywhere const strategy="no"

@everywhere include("/Users/francovazza/Desktop/Julia_prog/FOKKER/bo/spectra_lab/loss_gain.jl")
@everywhere include("/Users/francovazza/Desktop/Julia_prog/FOKKER/bo/spectra_lab/advect_assign.jl")
@everywhere include("/Users/francovazza/Desktop/Julia_prog/FOKKER/bo/spectra_lab/map_making.jl")
@everywhere include("/Users/francovazza/Desktop/Julia_prog/FOKKER/bo/spectra_lab/synchrotron.jl")
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
 file1=string(root,"dens",mod,"B.fits")
 file2=string(root,"temp",mod,"B.fits")
 file3=string(root,"vx",mod,"B.fits")
 file4=string(root,"vy",mod,"B.fits")
 file5=string(root,"vz",mod,"B.fits")
 file6=string(root,"bx",mod,"B.fits")
 file7=string(root,"by",mod,"B.fits")
 file8=string(root,"bz",mod,"B.fits")
 end
@everywhere tini=190
@everywhere tend=190
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

z=0. #;a[3]                                                                                                                                      
cd=1.6699997e-28
cv=2e+08

cb=1e6*sqrt(cd/(1+z)^3.*4*π)*cv #to have microGauss                                                                                              
   file1=string(root,"dens",mod,"B.fits")
   file2=string(root,"bz",mod,"B.fits")
   file3=string(root,"bx",mod,"B.fits")
   file4=string(root,"by",mod,"B.fits")
   filei=string(root,"temp",mod,"B.fits")

  f1=FITS(file1,"r")
  f2=FITS(file2,"r")
  f3=FITS(file3,"r")
  f4=FITS(file4,"r")
  f5=FITS(filei,"r")
  d=read(f1[1])#*cd                                                                                                      
   bz=read(f2[1])*cb
   bx=read(f3[1])*cb
   by=read(f4[1])*cb
   t=read(f5[1])

@simd  for j in eachindex(bz)
 
@inbounds @fastmath  bz[j]=sqrt(bx[j]^2.+by[j]^2.+bz[j]^2.)
  end
  bx=nothing
  by=nothing
   map=make_mapz(t,convert(Array{Float32},bz))

   file1=string(root,"dens",mod,"B.fits")

 file_pos=string(file1,"_tracers_",step,"_",model,"B.dat")
 file_spec=string(filei,"_tracers_",step,"_",model,"_specB.dat")

  so1=open(file_spec,"r")

  so2=open(file_pos,"r")

  pe=Array{Float64}(ngamma,np) 
  p=Array{Float64}(10,np)

 @inbounds @simd  for pp in eachindex(p)
 pa=read(so2,Float32)
 p[pp]=convert(Float64,pa)
  end

 @inbounds @simd  for pp in eachindex(pe)
 pa=read(so1,Float32)
 pe[pp]=convert(Float64,pa)
  end
 radio_trac=Array{Float64}(2,np) 
 radio_trac[:]=0.
#  ma=map[:,:,2]
#   dmin=0.01
#   dmax=1e2
#   id=find(x-> (x <=dmin ),ma)
#    ma[id]=dmin
#    id=find(x-> (x >=dmax ),ma)
#    ma[id]=dmax
#    ma=log10(ma)
##p0=contour((ma),aspect_ratio=1,fill=true,nlevels=256,size=(1200,1200))#,nlevels=10)   #c=:heat                          
# title!(string("time=",trunc(time0*1e-3)," kyr"))

#   id=find(x-> (x >=1e5 ),pe[100,:])
#  scatter!(p[2,id]-j1,p[1,id]-i1,label="radio emitting",color="red",size=(1200,1200))

 
    pe_tot=Array{Float64}(ngamma)
    pe_tot=fill(1e-20,ngamma)
     gal=gammaval
      pol=(pe_tot+1.)
 freqf=Vector{Float64}(6)

 freqf=[3e8,5.5e8,6.1e8,8.5e8,1.4e9,3e9] 
#    p2=plot(gal,pe_tot,line=:dash,color="black",label="whole population")
#    plot!(gammaval,pe[:,1],label="",color="blue")
   map_radio=Array{Float64}(n,n,3)
   map_radio[:]=1e-30
   map_radio[:,:,3]=3.
 
   @inbounds @simd for pp in 1:np
   if pe[2,pp] >=1e5 
    #println(p[:,pp])
    #println(pe[:,pp]) 
   spec=Vector{Float64}(6)
   spec[:]=0.
#@inbounds @simd  for ff in 1:6 
   freq=freqf[3]
   sync=synchrotron(gammaval,p[10,pp],pe[:,pp],freq)
   spec[3]=sync  
   freq=freqf[5]
   sync=synchrotron(gammaval,p[10,pp],pe[:,pp],freq)
   spec[5]=sync

#   dspec=log10(spec[3]-spec[5])/log10(freqf[5]-freqf[3])
 #  radio_trac[1,pp]=spec[5]
 #  radio_trac[2,pp]=dspec  
 #  println(spec[5], " ",dspec) 
   x1=convert(Int32,trunc(p[2,pp]-j1))
   y1=convert(Int32,trunc(p[1,pp]-i1))
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
@inbounds   map_radio[x1+i1,y1+j1,1]+=(spec[3]/(1+abs(ds)))
@inbounds   map_radio[x1+i1,y1+j1,2]+=(spec[5]/(1+abs(ds)))
   end
   end
    
#   plot!(freq,sync)#,label="",color="blue",line=:dot)
#   end
#   plot!(freqf,1e30*spec,label="")

   end
   end

        for j in 1:n
@simd   for i in 1:n
@inbounds @views  map_radio[i,j,3]=(log10(map_radio[i,j,1])-log10(map_radio[i,j,2])) /(log10(freqf[5])-log10(freqf[3]))
   end
   end
 #  plot!(gammaval,pe[:,id],label="",color="red")

#  yaxis!(L"N(\gamma)",:log10,(1e2,1e28))
#  xaxis!(L"\gamma",:log10,(g_min,g_max))
#  yaxis!(L"N(\gamma)",:log10,(1e15,1e25))
#  xaxis!(L"\gamma",:log10,(3e8,3e9))
#  plot(p0,p2)

file_out=string(root,"_radio_mapB.fits")
 f = FITS(file_out, "w");
  write(f,map_radio)
  close(f)

file_out=string(root,"_radio_BfieldB.fits")
 f = FITS(file_out, "w");
  write(f,map)
  close(f)
  error()

 end
println("total computing time")
toc()