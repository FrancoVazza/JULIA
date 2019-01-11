
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
#@everywhere using PyPlot                                                                                                          
@everywhere using FITSIO
@everywhere using Optim
@everywhere using Base
#@everywhere using Devectorize                                                                                                     


function make_map(d::Array{Float64},t::Array{Float64},vx::Array{Float64},vy::Array{Float64},vz::Array{Float64},bx::Array{Float64},by::Array{Float64},bz::Array{Float64})

map=Array{Float64}(i2-i1+1,j2-j1+1,3)
map[:]=0.

ni=i2-i1+1
#...main loop that produces the map 
invn=1./ni
       @inbounds for i in  1:i2-i1
       @simd    for j in  1:j2-j1

 ds=sum(d[:,i,j])
 vs=sum(abs.(bz[:,i,j]))
 ts=sum(t[:,i,j])+1e4

  map[i,j,1]=ds*invn    # map of projected gas density
  map[i,j,2]=ts*invn
  map[i,j,3]=vs*invn
end
end
 return map
end




function make_mapx(d::Array{Float32},t::Array{Float32})#,bx::Array{Float64})

mapx=Array{Float32}(j2-j1+1,l2-l1+1,2)
mapx[:]=0.

invn=1./(i2-i1+1)
       @inbounds @simd for i in  1:i2-i1 
 @views     mapx[:,:,1]+=d[i,:,:]*invn    # map of projected gas density                                 
 @views     mapx[:,:,2]+=t[i,:,:]*invn+1e2 
#   @views mapx[:,:,3]+=abs.(bx[i,:,:])*invn

end
return mapx
end


function make_mapy(d::Array{Float32},t::Array{Float32})#,by::Array{Float64})
mapy=Array{Float32}(i2-i1+1,l2-l1+1,2)
mapy[:]=0.

  invn=1./(j2-j1+1)
       @inbounds @simd for i in  1:j2-j1
   @views  mapy[:,:,1]+=d[:,i,:]*invn    # map of projected gas density     
    @views mapy[:,:,2]+=t[:,i,:]*invn+1e2
#   @views mapy[:,:,3]+=abs.(by[:,i,:])*invn

end

 return mapy
end

function make_mapz(d::Array{Float32},t::Array{Float32})#,bz::Array{Float64})

mapz=Array{Float32}(i2-i1+1,j2-j1+1,2)
mapz[:]=0.
  
  invn=1./(l2-l1+1)
       @inbounds @simd for i in  1:l2-l1
   @views  mapz[:,:,1]+=d[:,:,i]*invn    # map of projected gas density                                   
    @views mapz[:,:,2]+=t[:,:,i]*invn#+1e2
#   @views mapz[:,:,3]+=abs.(bz[:,:,i])*invn

end

 return mapz
end



function make_mapz1(v::Array{Float64})#,bz::Array{Float64})                                              

mapz=Array{Float64}(n,n)#i2-i1+1,j2-j1+1)
mapz[:]=0.

   invn=1./(n)
       @inbounds @simd for i in  1:n
       @views  mapz[:,:]+=v[:,:,i]*invn    # map of projected gas density                                                    
end

 return mapz
end



