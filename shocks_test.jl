#in order to use any of the using <> packages one should first do Pkg.add("<package_name>") from the Julia shell

#using Devectorize    #these are possible packages for speedup
#using Optim          #their helpfulness must be tested from code to code!
                      

                      #the @everywhere below ensure that the code can run in parallel, the info is passed to all procs

 @everywhere using DistributedArrays    #necessary for parallel map making
 @everywhere using HDF5    #necessary to open hdf5 files	
 @everywhere using LaTeXStrings
# @everywhere using Optim
 @everywhere using PyPlot
#....files to read in
     
@everywhere root="/Users/francovazza/Desktop/data/DATA/INFO/new/0A0/"    

@everywhere mod=""
@everywhere snap="349"
@everywhere file1=string(root,"DD0A0_dt_",mod,snap)    #input files    *dt* contains gas density, dark matter density and gas temperature
@everywhere file2=string(root,"DD0A0_v_",mod,snap)   # *b* contains the 3-d components of the v-field,
@everywhere file3=string(root,"DD0A0_dt_",mod,snap)   # *b* contains the 3-d components of the B-field,                                                          


@everywhere using HDF5    #necessary to open hdf5 files 

@everywhere const pi=3.141592654
   #cosmological parameters and miscellanea
 @everywhere const omegab=0.0468
 @everywhere const omegam=0.308
 @everywhere const hh=0.678 
 @everywhere const z=0.0
 @everywhere const mp=1.67e-24
 @everywhere const kb=1.38e-16
 
 @everywhere const cd= 2.76e-30  #physical density 
 @everywhere const cv= 1.494e+09    #physical velocity
 @everywhere const cb=sqrt(cd*4*pi)*cv   #physical B-fild Gauss

    #enzo conversion factors

 @everywhere const lbox=40/0.7  #Mpc
 @everywhere const n=400 #ncells
 @everywhere const res=lbox/n
 @everywhere const gamma=1.66667 
 @everywhere const fac=gamma*kb/mp
 @everywhere const mthr=3.

#first we define all functions that are called by the main algorithm (down below). Making them @everywhere defines them on all processors.

@everywhere function histog(mm::Vector{Float64},hm::Vector{Float64},macx::Array{Float64})

@fastmath for i in eachindex(mm)
   tag=find( x->(x >= mm[i]),macx)
   nshock=size(tag)
   hm[i]=nshock[1]
  end
return hm
end


@everywhere  function(psi_m_all,mpsi,f)
# ............reads in the Psi(M) function of Hoeft & Bruggen
fhb="mach_psi_table.txt"
nt=13
f1=string(fhb)
ff=readdlm(f1)

n=size(v)
mpsi=v[:,1]
ff=v[:,2:n[2]]

 end




@everywhere function psi_sel(mach,t,mpsi,f)

const nt=13
const n33=244

th=[1.00e-04,3.162e-04,1.000e-03,3.162e-03,1.000e-02,3.162e-02,1.000e-01,3.162e-01,1.000e+00,3.162e+00,1.000e+01,3.162e+01,1.000e+02]
#...................selection from the HB efficiency function for Psi(M)
to=t

@inbounds @fastmath for ii in 1:nt-1 
if to > th[ii] && to < th[ii+1]
tsel=ii
end
end


@inbounds @fastmath for ii in 1:n33-1
if m >= mpsi[ii] && m < mpsi[ii+1]
msel=ii
end
end

if to < th[1] 
tsel=1
end
if to > th[nt] 
tsel=nt
end
if m < mpsi[1] 
msel=1
end

if m > mpsi[n33] 
msel=n33
end

psi=f(tsel,msel)

 return psi

end

@everywhere function shocks(d,t,vx,vy,vz,bx,by,bz)
#::Array{Float64},t::Array{Float64},vx::Array{Float64},vy::Array{Float64},vz::Array{Float64},bx::Array{Float64},by::Array{Float64},bz::Array{Float64})

#::Array{Float64},t::Array{Float64},vx::Array{Float64},vy::Array{Float64},vz::Array{Float64},bx::Array{Float64},by::Array{Float64},bz::Array{Float64},mach::Array{Float64})

n3=size(d)

 macx=Array{Float64}(n3[1],n3[2],n3[3])
 macy=Array{Float64}(n3[1],n3[2],n3[3])
 macz=Array{Float64}(n3[1],n3[2],n3[3])
 mach=Array{Float64}(n3[1],n3[2],n3[3])

 macx[:,:,:]=0.
 macy[:,:,:]=0.
 macz[:,:,:]=0.
 mach[:,:,:]=0.

  div=Array{Float64}(n3[1],n3[2],n3[3])
  div[:,:,:]=0.

 @inbounds for i in 2:n3[1]-1, j in 2:n3[2]-1, l in 2:n3[3]-1
 @fastmath div[i,j,l]=0.5*(vx[i+1,j,l]-vx[i-1,j,l]+vy[i,j+1,l]-vy[i,j-1,l]+vz[i,j,l+1]-vz[i,j,l-1])
  end
 
  tag=find( x->(x <= -1.e6), div)
  nshock=size(tag)
  ijk=ind2sub((n3[1],n3[2],n3[3]),tag)
   println(nshock, "candidate shocks") 
   v1=Vector{Float64}(2)

@simd for i in 1:nshock[1]
 # in eachindex(candidate)
    ijk=ind2sub((n3[1],n3[2],n3[3]),tag[i])
    
  a=ijk[1]
  b=ijk[2]
  c=ijk[3]
 
@inbounds  dvx=-1.*(vx[a-1,b,c]-vx[a+1,b,c])  #...we need velocity in cm/s here --> k\m/s
@inbounds  dvy=-1.*(vy[a,b-1,c]-vy[a,b+1,c])
@inbounds  dvz=-1.*(vz[a,b,c-1]-vz[a,b,c+1])

@inbounds  if dvx < 0. && t[a+1,b,c] > t[a-1,b,c]
   vvx=dvx/(sqrt(fac*t[a-1,b,c]))
   mx=(4.*abs(vvx)+sqrt(16.*vvx*vvx+36.))*0.166666    
   v1[1]=mx
   v1[2]=macx[a+1,b,c]
   macx[a+1,b,c]=maximum(v1)
   elseif dvx <0. && t[a+1,b,c] <  t[a-1,b,c] 
   mx=abs(dvx)/(sqrt(fac*t[a+1,b,c]))
   mx=(4.*mx+sqrt(16.*mx*mx+36.))*0.166666
   v1[1]=mx
   v1[2]=macx[a-1,b,c]
   macx[a-1,b,c]=maximum(v1)
  end

@inbounds if dvy < 0. && t[a,b+1,c] > t[a,b-1,c]
@fastmath   vvy=dvy/(sqrt(fac*t[a,b-1,c]))
@fastmath   my=(4.*abs(vvy)+sqrt(16.*vvy*vvy+36.))*0.166666
   v1[1]=my
   v1[2]=macy[a,b+1,c]
   macy[a,b+1,c]=maximum(v1)
   elseif dvy <0 && t[a,b+1,c] <  t[a,b-1,c]
@fastmath   my=abs(dvy)/(sqrt(fac*t[a,b+1,c]))
@fastmath   my=(4.*my+sqrt(16.*my*my+36.))*0.166666
   v1[1]=my
   v1[2]=macy[a,b-1,c]
   macy[a,b-1,c]=maximum(v1)
  end

@inbounds if dvz < 0. && t[a,b,c+1] > t[a,b,c-1]
@fastmath   vvz=dvz/(sqrt(fac*t[a,b,c-1]))
@fastmath   mz=(4.*abs(vvz)+sqrt(16.*vvz*vvz+36.))*0.166666
   v1[1]=mz
   v1[2]=macz[a,b,c+1]
   macz[a,b,c+1]=maximum(v1)
   elseif dvz <0. && t[a,b,c+1] <  t[a,b,c-1]
 @fastmath  mz=abs(dvz)/(sqrt(fac*t[a,b,c+1]))
 @fastmath  mz=(4.*mz+sqrt(16.*mz*mz+36.))*0.166666
   v1[1]=mz
   v1[2]=macz[a,b,c-1]
   macz[a,b,c-1]=maximum(v1)
  end
@inbounds @fastmath   mach[i]=sqrt(macx[i]^2. +macy[i]^2.+macz[i]^2.) 
  
 end

#@fastmath   for i in eachindex(macx)
#@inbounds   mach[i]=sqrt(macx[i]^2. +macy[i]^2.+macz[i]^2.)
#    end
    println(minimum(mach)," ",maximum(mach))

#...cleaning for shock thickness

  tag0=find(x->(x <= mthr), mach)
  mach[tag0]=0.  

  tag=find( x->(x >= mthr), mach)
  
  nshock=size(tag)
  ijk=ind2sub((n3[1],n3[2],n3[3]),tag)

 v1=Vector{Float64}(3) 
 
@simd for i in eachindex(mach[tag])
@inbounds a=ijk[1][1]
@inbounds b=ijk[2][1]
@inbounds c=ijk[3][1]

  x1=a-1
  x2=a
  x3=a+1
  y1=b-1
  y2=b
  y3=b+1
  z1=c-1
  z2=c
  z3=c+1
 
  if x1 < 1 
  x1+=1
  end  
 
  if x3 > n3[1]+1
  x3-=1
  end
  if y1 < 1
  y1+=1
  end
  if y3 > n3[2]+1  
  y3-=1
  end 
  if z1 < 1
  z1+=1
  end
  if z3 > n3[3]+1
  z3-=1
  end  

@inbounds      v1[1]=macx[x1,y2,z2]
@inbounds      v1[2]=macx[x2,y2,z2]
@inbounds      v1[3]=macx[x3,y2,z2]
      max3=maximum(v1)
    
@inbounds     if macx[x1,y2,z2] != max3 
@inbounds     macx[x1,y2,z2]=0.
    end
@inbounds     if macx[x2,y2,z2] != max3 
@inbounds     macx[x2,y2,z2]=0.
    end
@inbounds     if macx[x3,y2,z2] != max3 
@inbounds     macx[x3,y2,z2]=0.
    end

@inbounds    v1[1]=macy[x2,y1,z2]
@inbounds    v1[2]=macy[x2,y2,z2]
@inbounds    v1[3]=macy[x2,y3,z2]

     max3=maximum(v1)
@inbounds     if macy[x2,y1,z2] != max3 
@inbounds     macy[x2,y1,z2]=0.
    end   
@inbounds     if macy[x2,y2,z2] != max3
@inbounds     macy[x2,y2,z2]=0.
    end
@inbounds     if macy[x2,y3,z2] != max3
@inbounds     macy[x2,y3,z2]=0.
    end

@inbounds     v1[1]=macz[x2,y2,z1]
@inbounds     v1[1]=macz[x2,y2,z2]
@inbounds     v1[1]=macz[x2,y2,z3]

    max3=maximum(v1)

@inbounds    if macz[x2,y2,z1] != max3 
@inbounds    macz[x2,y2,z1]=0.
   end   

@inbounds     if macz[x2,y2,z2] != max3
@inbounds     macz[x2,y2,z2]=0.
    end

@inbounds     if macz[x2,y2,z3] != max3 
@inbounds     macz[x2,y2,z3]=0.
    end


   end


@simd   for i in eachindex(macx[tag])
@fastmath   macx[i]=sqrt(macx[i]^2. +macy[i]^2.+macz[i]^2.)

   end
   println(minimum(macx)," ",maximum(macx))
   clear!(:macy)
   clear!(:macz)
   clear!(:div)
   clear!(:mach)   
#  mach=abs.(div)
  return macx
 end

@everywhere function make_map(dat1)#::Array{Float64})   

n1=(size(dat1[1], 1), size(dat1[2], 1)) #size(dat1[3],1))   #this builds up the 3D structure from the dat1 array which is passed to each processor

#...here below we compute the min and max in each direction of the array.
#...notice in this implementation the domain is decomposed in rectangles which are long as the entire domain in the z direction (best for map-making)

i1=dat1[1][1]
j1=dat1[2][1]
i2=dat1[1][n1[1]]
j2=dat1[2][n1[2]]
l1=1
l2=n

#if i1 > 1
#  i1=i1-1
#end
#if j1 >	1
#  j1=j1-1
#end
#if l1 >	1
#  l1=l1-1
#end
#if i2 < n
# i2=i2+1
#end
#if j2 <	n
# j2=j2+1
#end
#if l2 <	n
# l2=l2+1
#end
println("processor number ", myid(), " is reading the following dataset:")
println(i1," ",i2," ",j1," ",j2," ",l1," ",l2)
#...each processor reads the input files file1,file2 in parallel, i.e. each processor only reads a part of the entire domain, depending on the bounds defined above.

  d=cd*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
  dm=cd*h5read(file1,"Dark_Matter_Density",(i1:i2,j1:j2,l1:l2))
  t=h5read(file1,"Temperature",(i1:i2,j1:j2,l1:l2))
  vx=cv*h5read(file2,"x-velocity",(i1:i2,j1:j2,l1:l2))
  vy=cv*h5read(file2,"y-velocity",(i1:i2,j1:j2,l1:l2))
  vz=cv*h5read(file2,"z-velocity",(i1:i2,j1:j2,l1:l2))
  bz=cb*h5read(file3,"Bz",(i1:i2,j1:j2,l1:l2))
  by=cb*h5read(file3,"By",(i1:i2,j1:j2,l1:l2))
  bx=cb*h5read(file3,"Bx",(i1:i2,j1:j2,l1:l2))

  macx=shocks(d,t,vx,vy,vz,bx,by,bz)

 m1=log10(1.5)
 m2=log10(1e4)
  nm=100
  dem=(m2-m1)/(1.0*nm)
  mm=Vector{Float64}(nm)
  hm=Vector{Float64}(nm)

  @fastmath for i in 1:nm-1
  @inbounds  mm[i]=10.^(m1+dem*i)
  end
  @fastmath for i in 1:nm-1
  iw=find((macx.>mm[i]) & (macx.<=mm[i+1]))
  nw=size(iw)
  hm[i]=nw[1]
#  iw2=find((mult_rand.>minr+bir*i) & (mult_rand.<=minr+bir*(i+1)))
#  nw=size(iw2)
#  h2[i]=nw[1]/ntot
#  xh[i]=minr+bir*i
  end

  clf()
  loglog(mm,hm,linestyle="-",color="blue")
  xlabel(L"Mach")
   ylabel(L"dN/dlog(Mach)")
  annotate(xy=[60,1.5],"double detections",color="red",size=14)
  annotate(xy=[60,1.2],"random pixels",color="blue",size=14)
  annotate(xy=[2,2.5],L"M200>10^{12} M_{\odot} halos")
  axis([1,1e4,1e3,1e6])

  
error()


# p4=histogram(macx[:,10,10],mm)
#  display(plot(p4))
#  display(

# using PyPlot
# clf()
#if semilogx(mm,hm,color="green")
 

#...this is going to be the final map
 map=Array{Float64}(n1[1],n1[2],4)

println("the 3D size of each file is", size(d))


#...main loop that produces the map 
@fastmath for i in 1:i2-i1+1, j in 1:j2-j1+1#, l in 1:n

@inbounds     d1=d[i,j,:]
@inbounds      ds=sum(d[i,j,:])
@inbounds        ts=sum(t[i,j,:])
@inbounds          b1=bz[i,j,:]
@inbounds           m1=sum(macx[i,j,:])/n
#          ma=sum(dot(d1,m1))
@inbounds           rm1=sum(dot(d1,b1))
         
@inbounds         map[i,j,1]=ds/n  #density map
        map[i,j,2]=ts  #temperature map
       map[i,j,3]=rm1 #RM map
      map[i,j,4]=m1
 end
 
 return map
end



#...MAIN PROGRAM
tic()   #...this starts a timer, closed by toc() below, useful to check performance
#@everywhere using Plots
 map=DArray(make_map,(n,n,4))  #...this builds up a Distributed Array (see Julia's documentation), i.e. an array made of chunks, each passed to a different processor. The distributed array is directly passed to the make_map function defined above

 map_total=convert(Array,map)  #....this step assembles a unique array by combining the different pieces of the distributed array 

println("The minimum and maximum values in the map are ",minimum(map_total[:,:,1])," ",maximum(map_total[:,:,1]))


  ma1=(map_total[:,:,1])
  ma2=(map_total[:,:,2])
  ma3=abs.(map_total[:,:,3])  
  ma4=(map_total[:,:,4])
#  for i in eachindex(ma4)
#  ma4[i]=(ma4[i])/ma1[i]
#  end
  
 toc()

@everywhere using FITSIO
 filep1=string(root,"map_dens4",mod,snap,".fits")
 f = FITS(filep1, "w");
 write(f,ma1)
 close(f)

 filep2=string(root,"map_temp4",mod,snap,".fits")
 f = FITS(filep2, "w");
 write(f,ma2)
 close(f)

 filep3=string(root,"map_RM4",mod,snap,".fits")
 f = FITS(filep3, "w");
 write(f,ma3)
 close(f)

filep4=string(root,"map_div4",mod,snap,".fits")
 f = FITS(filep4, "w");
 write(f,ma4)
 close(f)







