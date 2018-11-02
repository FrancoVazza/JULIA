@everywhere const h0=0.7
@everywhere const Lbox=40/0.7 #Mpc/h                                            
@everywhere const n=400 #grid size                                              
@everywhere const np=1000 #number of tracers                                    
@everywhere const dx=1000.*Lbox/n #kpc
@everywhere const pi=3.14159
@everywhere const kpc=3.086e21 #cm in one kpc
@everywhere const yr=3.154e7 #s in one yr
@everywhere using HDF5
@everywhere using PyPlot
@everywhere using Optim
const root="/Users/francovazza/Desktop/data/DATA/INFO/new/0A0/"     #main folder


#first we define all functions that are called by the main algorithm (down below). Making them @everywhere defines them on all processors.

#...DOMAIN BORDERS
@everywhere i1=200
@everywhere j1=200
@everywhere l1=200
@everywhere i2=357
@everywhere j2=357
@everywhere l2=300
@everywhere ng1=i2-i1+1
@everywhere ng2=j2-j1+1
@everywhere ng3=l2-l1+1
#...SPECTRAL PARAMETERS
@everywhere const  g_max=100000.0
@everywhere const  g_min=10.0
@everywhere const  dg=10.0
@everywhere const  part=2  #...1=proton  2=electron                                                                                             
@everywhere const  mp=1.67e-24
@everywhere const  me=0.1093897e-28   #g                                                                                                        
@everywhere const  vc=2.99792458e10   #cm s-                                                                                                    
@everywhere const  kpctocm=3.086e21
@everywhere const  cmtoMpc=3.086e24
@everywhere const  kb=1.380658e-16    #erg k-1                                                                                                  
#@everywhere const  volume=3*log10(vol_trac*kpctocm) #vol im cm^3case of                                                                        
@everywhere const norm=30 #to have all in 1e30 erg or 1e30 erg/s                                                                               
@everywhere const erest=0.510988946 #electron rest mass energy                                                                                 
@everywhere const evtoerg=1.60218e-12 #                                                                                                        
@everywhere const  b1=1.37e-20 #...1/s                                                                                                         
@everywhere const  b2=1.9e-9  #...1/s                                                                                                          
@everywhere const  b3=7.2e-12 #...1/s                                                                                                          
@everywhere const xi=0.5 #                                                                                                                     
@everywhere const lcorr=20.0 #...in[100]kpc                                                                                                   
@everywhere const mthr=1.3
@everywhere const gyrtosec=1.0e9*3.0e7
@everywhere const fi=1e-5
@everywhere const ngamma=(floor(Int64,(g_max-g_min)/dg))
println(ngamma)
@everywhere const gend=ngamma
@everywhere const gammaval=collect(g_min:dg:g_max)


@everywhere function age(idt::Float64,idg::Float64,aa1::Float64,bb::Float64,g1::Float64,q1::Float64,nn1::Float64,nn2::Float64,aa2::Float64,g2::Float64)
#@everywhere function age(idt,idg,aa1,bb,g1,q1,nn1,nn2,aa2,g2)
return @fastmath ((inv(idt+idg*(aa1+bb*g1^2.)))*(idt*q1+nn1*idt+nn2*idg*(aa2+bb*g2^2.)))
end

@everywhere function coul(cou::Float64,g1::Float64,inth::Float64)
#@everywhere function coul(cou,g1,inth)

return @fastmath cou*(1.+0.0133333*(log(g1*inth)))
end

@everywhere function bloss(p1::Int64,b2::Float64,b0::Float64,b3::Float64,zz4::Float64)
#@everywhere function bloss(p1,b2,b0,b3,zz4)
return @fastmath p1*b2*(0.666*b0^2.*1.0e-12+b3*zz4)
end

@everywhere function reacc(delta::Float64,nn::Array{Float64,1},gend::Int64,gmin_re::Int64,gam::Array{Float64,1},n_re::Array{Float64,1})
 ntemp=similar(nn)
 gtemp=similar(nn)
  @simd for j in 1:ngamma #gmin_re:ngamma                                                                   
@inbounds ntemp=nn[gmin_re:j]
@inbounds gtemp=gam[gmin_re:j]

@fastmath  @inbounds   n_re[j]=(delta+2.)*(gam[j])^(-1.*delta)*(dg*sum(dot(ntemp,gtemp)^(delta-1.)))
#nn[gmin_re:j],(gam[gmin_re:j])^(delta-1.))))
  end
  return n_re
  end


#@everywhere function evolve_spectrum(zz,v,t2,t1,nth,m,b0,ecr,nn,shock,delta_t,volume,gam,g_max,g_min,dg,ngamma,part,vshock,tacc,n_inj,q_inj,n_re)
@everywhere function evolve_spectrum(zz::Float64,v::Float64,t2::Float64,t1::Float64,nth::Float64,m::Float64,b0::Float64,ecr::Float64,nn::Array{Float64,1},shock::Int64,delta_t::Float64,volume::Float64,gam::Array{Float64,1},g_max::Float64,g_min::Float64,dg::Float64,ngamma::Int64,part::Int64,vshock::Float64,tacc::Float64,n_inj::Array{Float64,1},q_inj::Array{Float64,1},n_re::Array{Float64,1})

gend=ngamma
@fastmath m2=m^2.
@fastmath vpre=1.0e5*v
@fastmath f=(4.*m2)/(m2+3.)
@fastmath delta=2.*(m2+1.)/(m2-1.)
@fastmath zz4=(1.+zz)^4.

  ic_lose=0.0
  diff=0.0
  sh_gain=0.0

@fastmath const dt=1.0*(t2-t1)
  #losses                                                                                                          

const  b3=7.2e-12 #...1/s                                                                                          
@fastmath const   cou=nth*1.2e-12 #..1/s                                                                           
@fastmath const  idt=inv(dt)
@fastmath const  idg=inv(dg)
@fastmath const  icou=inv(cou)
@fastmath const  inth=inv(nth)
@fastmath const  cons1=2.3e29*(erest*1e-3)^(0.3333)*b0^(-0.3333)*(lcorr*0.05)^(0.66666)
@fastmath const  cons2=(f-1.)*inv(f+1.)
@fastmath const  dg2=dg*0.5

#nno=SharedArray{Float64}(ngamma)                                                                                  
 #.....SHOCK ACCELERATION OF FRESH PARTICLES                                                                       

#                                                                                                                  

  if shock >= 1 && m >= mthr   #injection of power-law spectrum of particles                                       

@fastmath    beff=sqrt(b0+3.2*(1+zz)^2.)
@fastmath    tacc=inv(vshock*inv(3000.0))*3.0e7*2.4e4*sqrt(delta*inv(beff)) #..acceleration time from Kang 2012 Eq\

@inbounds    gam_c=gam[ngamma]

@fastmath @simd   for i in 1:ngamma   #finds maximum gamma where tacc < tlosses                     
@fastmath   ic_lose=b1*((gam[1]^2.))*zz4
@inbounds   diff=(gam[i]-1.0)*cons1 # in cm^2/s                                                                             
@inbounds @fastmath   sh_gain=(inv(f)*gam[i]*(vpre)^2.*cons2*inv(3.*diff))

          if ic_lose > sh_gain
          gam_c=gam[i]
           end
               end

@fastmath   inte=(-delta+2.0)*inv(((gam_c)^(-delta+2.0)-(g_min)^(-delta+2.0)))
@fastmath   Ke=((10.^ecr)*inv(1e6*erest*evtoerg))*inv(inte)

@inbounds @fastmath @simd    for j in 1:ngamma#eachindex(n_inj)
     n_inj[j]=log10(Ke)-delta*log10(gam[j])+log10((1.-gam[j]*inv(gam_c))^(delta-2.))

        if isnan(n_inj[j]) || isinf(n_inj[j])
                n_inj[j]=-60
         end
      q_inj[j]=tacc*10.^n_inj[j]
               end

              end

#...SHOCK REACCELERATION                                                                                           
  if shock == 2
   n_re[1]=nn[1]
   gmin_re=2

 # @fastmath @simd  for j in 1:ngamma #gmin_re:ngamma                                                 
 n_re=reacc(delta,nn,gend,gmin_re,gam,n_re)

# function reacc(delta::Float64,nn::Array(Float64,2))
# @fastmath @simd  for j in 1:ngamma #gmin_re:ngamma  
# @inbounds   n_re[j]=(delta+2.)*(gam[j])^(-1.*delta)*(dg*sum(dot(nn[gmin_re:j],(gam[gmin_re:j])^(delta-1.))))
#  end
#  return n_re
   nn=n_re
 end


#setting to zero losses that are not valid either for protons (p=1) or electrons (p=2)                             
    p1=1
    p2=1
    if part==1
    p1=0
    end
    if part==2
    p2=0
    end
 #   toc()                                                                                                         

#   tic()                                                                                                          
#....INTEGRATION OF LOSSES AND INCLUSION OF FRESHLY INJECTED PARTICLES                                             
@inbounds @fastmath @simd  for j in 1:ngamma-1
      gg=gend-j
      ga=gam[gg]
       g1=ga-dg2
       g2=g1+dg

       aa1=coul(cou,g1,inth)
       aa2=coul(cou,g2,inth)
       bb=bloss(p1,b2,b0,b3,zz4)

     nn1=nn[gg]
     nn2=nn[gg+1]
     q1=q_inj[gg]
     nn[gg]=age(idt,idg,aa1,bb,g1,q1,nn1,nn2,aa2,g2)
    end

   id=find(x-> (x <= 0),nn)
   nn[id]=1e-30

 return nn
end


function eth(n::Float64,volume::Float64,t::Float64)
@fastmath ethermal=log10(n/mp)+volume+log10(1.5*kb*t)
return ethermal
end

function ecri(n::Float64,vshock::Float64,volume::Float64)
 @fastmath Ecr_inj=15.0+log10(fi)+log10(n/mp*vshock^3.0)+log10(mp)+0.6666*volume
return Ecr_inj
end

function spectra(p::Array{Float64,2},pe::Array{Float64,2},dt::Float64,zz::Float64,lsize::Float64,tt::Int64)
#function spectra(p,pe,dt,zz,lsize)

    shock=0
    n_inj=fill(0.0,ngamma)
    q_inj=fill(-60.0,ngamma)
    n_re=fill(0.0,ngamma)
    #Ecr_inj=-60. log10.(fi1)+15.0+log10(n0*vshock1^3.0)+log10(mp)+0.6666*volume
volume=3*log10(lsize*kpctocm)
@fastmath @inbounds @simd for i in 1:np
 
  ethermal=eth(p[i,7],volume,p[i,8])

#@fastmath ethermal=log10(p[i,7]/mp)+volume+log10(1.5*kb*p[i,8])  #...observed                                                                                  
 Ecr_inj=-20.0
 m=p[i,9]

@fastmath tpre=p[i,8]*(16.*m^2.)/((5.*m^2-1)*(m^2.+3))                                                                                                                  
 @fastmath cspre=1e-5*sqrt(1.666*tpre*kb*inv(mp*1.1))    #...preshock                                                                             
  vshock=0.
  ecr=0.0
  m=0.0
  v=0.
  shock=0
    
  if p[i,9] >= 2.0 
  
   m=p[i,9]
   vshock=m*cspre
    shock=1
    if pe[2,i] > 1e12
    shock=2
    end
     Ecr_inj=ecri(p[i,7],vshock,volume)
         
    ecr=Ecr_inj-norm
    v=vshock
  end

  nth=p[i,7]/mp
  t2=dt*gyrtosec
  t1=0.0
  
  nn=pe[:,i]

  delta_t=t2-t1
#  epoch+=delta_t
  tacc=0.0

  n_inj[:]=0.
  q_inj[:]=-60.
  n_re[:]=0.0
  b0=p[i,10]

  aa=evolve_spectrum(zz,v,t2,t1,nth,m,b0,ecr,nn,shock,delta_t,volume,gammaval,g_max,g_min,dg,ngamma,part,vshock,tacc,n_inj,q_inj,n_re)
  pe[:,i]=nn
  end


return pe
end



function make_map(d,t,vx,vy,vz,bx,by,bz)

map=Array{Float64}(i2-i1+1,j2-j1+1,3)

println("the 3D size of each file is", size(d))


#...main loop that produces the map 
@inbounds for i in  1:i2-i1
          for j in  1:j2-j1

@fastmath ds=sum(d[i,j,:])
@fastmath vs=sum(abs.(bz[i,j,:]))
@fastmath ts=sum(t[i,j,:])+1e4

  map[i,j,1]=ds    # map of projected gas density
  map[i,j,2]=ts
  map[i,j,3]=vs
end
end
 return map
end


function assign_tracers()

p=Array{Float64}(np,10)

p[:,1]=i1+rand(np)*(i2-i1)
p[:,2]=j1+rand(np)*(j2-j1)
p[:,3]=l1+rand(np)*(l2-l1)

return p
end


function assign_tracers_fields(p::Array{Float64,2},d::Array{Float64,3},t::Array{Float32,3},vx::Array{Float64,3},vy::Array{Float64,3},vz::Array{Float64,3},bx::Array{Float64,3},by::Array{Float64,3},bz::Array{Float64,3})
#function assign_tracers_fields(p,d,t,vx,vy,vz,bx,by,bz)

deltaT=0.0

@inbounds for i in 1:np
 x=trunc(Int,p[i,1])-i1+1
 y=trunc(Int,p[i,2])-j1+1
 z=trunc(Int,p[i,3])-l1+1
 p[i,7]=d[x,y,z]
 deltaT=p[i,8]/t[z,y,z]
 p[i,9]=0.
 if deltaT > 1   
 p[i,9]=sqrt((8.*deltaT-7.0+sqrt( (7.0-8.0*deltaT)^2.+15.))/5.0)
 end

p[i,8]=t[x,y,z]
p[i,10]=sqrt(bx[x,y,z]^2.+by[x,y,z]^2.+bz[x,y,z]^2.)
p[i,4]=vx[x,y,z]
p[i,5]=vy[x,y,z]
p[i,6]=vz[z,y,z]
end

return p
end

function move_tracers(p::Array{Float64,2},dt::Float64,dx::Float64,i1::Int64,i2::Int64,j1::Int64,j2::Int64,l1::Int64,l2::Int64)
#function move_tracers(p,dt,dx,i1,i2,j1,j2,l1,l2)
 v1=Vector{Float64}(2)
 cc=dt*yr/(dx*kpc)
@inbounds for i in 1:np
p[i,1]+=cc*p[i,4]#*yr/(dx*kpc)
p[i,2]+=cc*p[i,5]#*yr/(dx*kpc)
 p[i,3]+=cc*p[i,6]#*yr/(dx*kpc)
   v1[1]=i1
   v1[2]=p[i,1]
   p[i,1]=maximum(v1)
   v1[1]=i2
   v1[2]=p[i,1]
   p[i,1]=minimum(v1)

   v1[1]=j1
   v1[2]=p[i,2]
   p[i,2]=maximum(v1)
   v1[1]=j2
   v1[2]=p[i,2]
   p[i,2]=minimum(v1)

   v1[1]=l1
   v1[2]=p[i,3]
   p[i,3]=maximum(v1)
   v1[1]=l2
   v1[2]=p[i,3]
   p[i,3]=minimum(v1)  
end

#  println(maximum(p[:,1]),minimum(p[:,1]))
#  println(maximum(p[:,2]),minimum(p[:,2]))
#  println(maximum(p[:,3]),minimum(p[:,3])) 

return p
end




##############
#MAIN PROGRAM#
##############
tic()

using PyPlot
 clf()
# semilogx(gammaval,gammaval*0.0)
# axis([g_min,g_max,1,1e20])


dt=1e8  #yr .....fixed timestep
sna=349
snap_name=string(sna)

file1=string(root,"DD0A0_dt_",snap_name)    #input files    *dt* contains gas density, dark matter \density and gas temperature + 3 mag.field components                                                
file2=string(root,"DD0A0_v_",snap_name)   # *v* contains the 3-d components of the v-field          
file_conv=string(root,"DD0A0",snap_name,".conv2")  #file with conversion factors from enzo internal\ to cgs                                                                                            
p=assign_tracers()
@everywhere tini=1
@everywhere tend=20
pe=Array{Float64}(ngamma,np)
#pe=fill(0.0,ngamma,np)

for tt in tini:tend
 clf()
 semilogx(gammaval,gammaval*0.0)
 axis([g_min,g_max,1e1,1e30])


println(tt)

#file1a=string(file1,snap_name)
#file2a=string(file2,snap_name)
#file_conva=string(file_conv,snap_name,".conv2")
a=readdlm(file_conv)
z=a[3]
cd=a[4]
cv=a[5]
cb=1e8*1e6*sqrt(cd/(1+z)^3.*4*pi)*cv #to have microGauss   

#println("processor number ", myid(), " is reading the following dataset:")
#println(i1," ",i2," ",j1," ",j2," ",l1," ",l2)
#if t == 1 
d=cd*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
vz=cv*h5read(file2,"z-velocity",(i1:i2,j1:j2,l1:l2))
vy=cv*h5read(file2,"y-velocity",(i1:i2,j1:j2,l1:l2))
vx=cv*h5read(file2,"x-velocity",(i1:i2,j1:j2,l1:l2))
t=h5read(file1,"Temperature",(i1:i2,j1:j2,l1:l2))
bx=cb*h5read(file1,"Bx",(i1:i2,j1:j2,l1:l2))
by=cb*h5read(file1,"By",(i1:i2,j1:j2,l1:l2))
bz=cb*h5read(file1,"Bz",(i1:i2,j1:j2,l1:l2))


#if tt==tend 
#map=make_map(d,t,vx,vy,vz,bx,by,bz)
#end

#p=assign_tracers()
println("assign fields")
tic()

p=assign_tracers_fields(p,d,t,vx,vy,vz,bx,by,bz)
toc()

println("move")
tic()

p=move_tracers(p,dt,dx,i1,i2,j1,j2,l1,l2)
toc()

 zz=0.1
 lsize=dx #kpc linear size of tracer

println("spectra")
tic()
   pe=spectra(p,pe,dt*1e-9,zz,lsize,tt)
toc()
#scatter(p[:,1],p[:,2],color="blue",s=2,alpha=0.1)

println("plot")
tic()

for i in 1:100
#println(pe[1:10,i])
loglog(gammaval[1:ngamma],pe[1:ngamma,i])#,xrange=[g_min,g_max],yrange=[20,70])                                                                             
end
toc()
# axis([g_min,g_max,20,70])                                          
# println(pe[1:10,i])
# plot(gammaval,10.^pe[:,1],xrange=[g_min,g_max],yrange=[20,70])       
# axis([g_min,g_max,20,70]) 

filep1=string("/Users/francovazza/Desktop/data/DATA/INFO/new/0A0/spectra_yes_reacc",tt,".png")

savefig(filep1)

end
error()



toc()
clf()
error()
ds=map[:,:,1]
pcolor(ds, norm=matplotlib[:colors][:LogNorm](vmin=minimum(ds), vmax=maximum(ds)), cmap="PuBu_r")
xticks([])
yticks([])
colorbar()
filep1="/Users/francovazza/Desktop/data/DATA/INFO/new/0A0/map_d_slice.png"
savefig(filep1)

clf()
ds=(map[:,:,2]+1e2)
pcolor(ds, norm=matplotlib[:colors][:LogNorm](vmin=minimum(ds), vmax=maximum(ds)), cmap="autumn")
xticks([])
yticks([])
colorbar()
filep1="/Users/francovazza/Desktop/data/DATA/INFO/new/0A0/map_t_slice.png"
savefig(filep1)

clf()
ds=(1e9*map[:,:,3])
pcolor(ds, norm=matplotlib[:colors][:LogNorm](vmin=1e-5, vmax=1), cmap="jet")
xticks([])
yticks([])
colorbar()
filep1="/Users/francovazza/Desktop/data/DATA/INFO/new/0A0/map_vz_slice.png"
savefig(filep1)

toc()


