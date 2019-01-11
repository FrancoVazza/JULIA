# SIMULATION PARAMETERS                                                                                                            
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
@everywhere using HDF5
@everywhere using FITSIO
@everywhere using Optim
@everywhere using Base
#@everywhere using Devectorize                                                                                                     



function assign_tracers()

p=Array{Float64}(10,np)

p[1,:]=i1+rand(np)*(i2-i1)
p[2,:]=j1+rand(np)*(j2-j1)
p[3,:]=l1+rand(np)*(l2-l1)

return p
end

function assign_tracers_loc()

p=Array{Float64}(10,np)
p[:]=0.
xc=i1
yc=j1
zc=l1
del=10
p[1,:]=(xc+rand(np)*4.+100.) #*(i2-i1)*0.1#+150
p[2,:]=(yc+rand(np)*4.+10.) #*(j2-j1)*0.1#+100
p[3,:]=(zc+rand(np)*4.+320.) #*(l2-l1)*0.1#+50

return p
end


function assign_tracers_shock()

p=Array{Float64}(10,np)
p[:]=0.
xc=28
yc=j1 #0.15*(j2-j1)
zc=l1

p[2,:]=(yc+rand(np)*n) #*(i2-i1)*0.1#+150                                                             
p[1,:]=(xc+rand(np)*20.) #*(j2-j1)*0.1#+100                                                              
p[3,:]=(zc+rand(np)*n) #*(l2-l1)*0.1#+50                                                              

return p
end


function assign_tracers_dens(f::Array{Float64,3})
#...locate gas peak in the volume
    mad=maximum(f)
    n3=size(f) 
    ida=find(x-> (x == mad),f)
  
  ijk=ind2sub((n3[1],n3[2],n3[3]),ida)
 
  ha=(ijk[1])
  hb=(ijk[2])
  hc=(ijk[3])

  p=Array{Float64}(10,np)

  p[1,:]=rand(np)*(40.)+Float64(ha[1])+i1
  p[2,:]=rand(np)*(40.)+Float64(hb[1])+j1
  p[3,:]=rand(np)*(40.)+Float64(hc[1])+l1

return p
end
~

function assign_tracers_fields(p::Array{Float64,2},d::Array{Float64,3},t::Array{Float64,3},vx::Array{Float64,3},vy::Array{Float64,3},vz::Array{Float64,3},bx::Array{Float64,3},by::Array{Float64,3},bz::Array{Float64,3})
#function assign_tracers_fields(p,d,t,vx,vy,vz,bx,by,bz)

deltaT=0.0

@inbounds for i in 1:np

 x=trunc(Int,p[1,i])-i1+1
 y=trunc(Int,p[2,i])-j1+1
 z=trunc(Int,p[3,i])-l1+1
  p[7,i]=d[x,y,z]
 deltaT=p[8,i]/t[x,y,z]
 p[9,i]=0.
 if deltaT > 1   
 p[9,i]=sqrt((8.*deltaT-7.0+sqrt( (7.0-8.0*deltaT)^2.+15.))/5.0)
 end

 v2=[1e9,t[x,y,z]]
 p[8,i]=minimum(v2)

p[10,i]=sqrt(bx[x,y,z]^2.+by[x,y,z]^2.+bz[x,y,z]^2.)
p[4,i]=vx[x,y,z]
p[5,i]=vy[x,y,z]
p[6,i]=vz[z,y,z]

end
 
return p
end

function move_tracers(p::Array{Float64,2},dt::Float64,dx::Float64,i1::Int64,i2::Int64,j1::Int64,j2::Int64,l1::Int64,l2::Int64)
#function move_tracers(p,dt,dx,i1,i2,j1,j2,l1,l2)
 v1=Vector{Float64}(2)
 cc=dt*yr/(dx*kpc)
@inbounds @simd for i in 1:np
@fastmath  p[1,i]=p[1,i]+cc*p[4,i]#*yr/(dx*kpc)
@fastmath  p[2,i]=p[2,i]+cc*p[5,i]#*yr/(dx*kpc)
@fastmath  p[3,i]=p[3,i]+cc*p[6,i]#*yr/(dx*kpc)

   v1[1]=i1+1
   v1[2]=p[1,i]
   p[1,i]=maximum(v1)
   v1[1]=i2-1
   v1[2]=p[1,i]
   p[1,i]=minimum(v1)

   v1[1]=j1+1
   v1[2]=p[2,i]
   p[2,i]=maximum(v1)
   v1[1]=j2-1
   v1[2]=p[2,i]
   p[2,i]=minimum(v1)

   v1[1]=l1+1
   v1[2]=p[3,i]
   p[3,i]=maximum(v1)
   v1[1]=l2-1
   v1[2]=p[3,i]
   p[3,i]=minimum(v1)  
end

return p
end



