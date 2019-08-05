#......basic routine to compute complexity from 2 snapshots of an ENZO simulation
#......the data are suppposed to be regular 3D datasets, with conversion factors externally provided to convert internal units into cgs
#......@Franco Vazza, august 2019

using StatsBase
using HDF
using FITSIO

sna1=100
sna2=sna1+ 5

root=""  #root folder for input filts

rooto=root  #...folder for output files 

#.....1st dataset (t)
file1=string(root,)   #....file with gas density, temperature and dark matter density
file2=string(root,)   #....file with 3d velocities
file3=string(root,)   #....file with 3d magnetic field

#.....2nd dataset (t+dt)
file1b=string(root,)
file2b=string(root,)
file3b=string(root,)

#.....definition of constants
 const tlag=10   #...time step separation between snapshots to be used 
 const n =400    #...1D grid size
 const dx=100.   #...spatial resolution in kpc
 const vol=3.*log10(3.085e21*dx)  #....cell volume
 const mp=1.67e-24
 const kb=1.38e-16

#....ENZO conversion file (necessary to transform from code to cgs units)
#....1st dataset (t)
filec=string(root,) 
    a=readdlm(filec)
    z=a[3]  #redshift                                                        
    cd=a[4] #conversion factor for density                                      
    cv=a[5] #conversion factor for velocity                                     

#....2nd dataset (t+dt)
 filec=string(root,)
    a=readdlm(filec)
    zb=a[3]     #redshift                                                      
    cdb=a[4] #conversion factor for density                                    
    cvb=a[5] #c                                     

 zed=z
 cb=sqrt(cd*4*pi)*cv   #..in Gauss                                              
 cd=cd/float(1+zed)^3.
 cv=cv/float(1+zed)
 cb=cb/float(1+zed)^2.
#...conversion second snap                                                      
 zed2=zb
 cbb=sqrt(cdb*4*pi)*cvb   #..in Gauss                                          
 cdb=cdb/float(1+zed2)^3.
 cvb=cvb/float(1+zed2)
 cbb=cbb/float(1+zed2)^2.


# FUNCTION TO COMPUTE TRANSITION PROBABILITY MATRIX


function prob(ni::Int64,v1::Array{Int64},probi::Array{Float64},idi::Array{Int64})
 to=0.
 @inbounds @simd for j in 1:ni
  jj=v1[idi[j]]
 @fastmath  probi[jj]+=1.0
 @fastmath  to+=1.0
  end
@fastmath  probi/=to
return probi
end




#....FUNCTION TO COMPUTE STATISTICAL COMPLEXITY

 function compute_complex(n)

     nbin=200    #....number of logaritmic energy bins
     #...define adjustable edges of subset to read from the HDF5 file
     i1=1
     j1=1
     i2=i1+n-1
     j2=j1+n-1
     l1=1
     l2=l1+n-1
   #....read 1st dataset (t)     
     d=h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
     vz=h5read(file2,"z-velocity",(i1:i2,j1:j2,l1:l2))
     vy=h5read(file2,"y-velocity",(i1:i2,j1:j2,l1:l2))
     vx=h5read(file2,"x-velocity",(i1:i2,j1:j2,l1:l2))
     t=h5read(file1,"Temperature",(i1:i2,j1:j2,l1:l2))
     bz=h5read(file3,"Bz",(i1:i2,j1:j2,l1:l2))
     by=h5read(file3,"By",(i1:i2,j1:j2,l1:l2))
     bx=h5read(file3,"Bx",(i1:i2,j1:j2,l1:l2))

     #.....arrays of Energy fields for the two datasets 
     eb0=Array{Float32}(n,n,n)
     ed0=Array{Float32}(n,n,n)
     ev0=Array{Float32}(n,n,n)
     eb1=Array{Float32}(n,n,n)
     ed1=Array{Float32}(n,n,n)
     ev1=Array{Float32}(n,n,n)
     logcd=log10(cd)
     logcv=log10(cv)
     logmp=log10(kb/mp)
     

 #....conversion factors from ENZO code units into cgs
     c1=1e31*0.5/0.6*cd*cv^2.
     c2=1e31*0.6*cd*kb/mp
     c3=1e41*cb^2./(8.*pi)
#.....the energy fields are computed 
     @inbounds @simd    for i in eachindex(d)
         di=d[i]
         ti=t[i]
         @fastmath   vi=vx[i]^2.+vy[i]^2.+vz[i]^2.
         @fastmath   bi=bx[i]^2.+by[i]^2.+bz[i]^2.
         @fastmath  ev0[i]=c1*di*vi
         @fastmath  ed0[i]=c2*di*ti
         @fastmath  eb0[i]=c3*(bi)
     end

    @views ev0[:]=log10.(ev0[:])
    @views ed0[:]=log10.(ed0[:])
    @views eb0[:]=log10.(eb0[:])
    @views d[:]=log10.(d[:])
    @views t[:]=log10.(t[:])


   #...read 2nd dataset (t+dt)
     d=h5read(file1b,"Density",(i1:i2,j1:j2,l1:l2))
     vz=h5read(file2b,"z-velocity",(i1:i2,j1:j2,l1:l2))
     vy=h5read(file2b,"y-velocity",(i1:i2,j1:j2,l1:l2))
     vx=h5read(file2b,"x-velocity",(i1:i2,j1:j2,l1:l2))
     t=h5read(file1b,"Temperature",(i1:i2,j1:j2,l1:l2))
     bz=h5read(file3b,"Bz",(i1:i2,j1:j2,l1:l2))
     by=h5read(file3b,"By",(i1:i2,j1:j2,l1:l2))
     bx=h5read(file3b,"Bx",(i1:i2,j1:j2,l1:l2))

     #....conversion factors from ENZO code units into cgs
     c1=1e31*0.5/0.6*cdb*cvb^2.
     c2=1e31/0.6*cdb*kb/mp
     c3=1e41*cbb^2./(8.*pi)
     @inbounds @simd    for i in eachindex(d)
         di=d[i]
         ti=t[i]
         @fastmath vi=vx[i]^2.+vy[i]^2.+vz[i]^2.
         @fastmath bi=bx[i]^2.+by[i]^2.+bz[i]^2.
         @fastmath  ev1[i]=c1*di*vi
         @fastmath  ed1[i]=c2*di*ti
         @fastmath  eb1[i]=c3*(bi)
     end
     
    @views ev1[:]=log10.(ev1[:])
    @views ed1[:]=log10.(ed1[:])
    @views eb1[:]=log10.(eb1[:])
    @views d[:]=log10.(d[:])
    @views t[:]=log10.(t[:])


     minE1=minimum(ed1)
     minE2=minimum(ev1)
     minE3=minimum(eb1)
     mino1=minimum(ed0)
     if minE1 >= mino1
         minE1=mino1
     end
     mino1=minimum(ev0)
     if minE2 >= mino1
         minE2=mino1
     end
     mino1=minimum(eb0)
     if minE3 >= mino1
         minE3=mino1
     end
   
      maxE1=maximum(ed1)
      maxE2=maximum(ev1)
      maxE3=maximum(eb1)
      bin1=minE1
      bin2=minE2
      bin3=minE3
    
     nbin=200
     bin1=(maxE1-minE1)/nbin
     bin2=(maxE2-minE2)/nbin
     bin3=(maxE3-minE3)/nbin

     clear!(:bx)
     clear!(:by)
     clear!(:bz)
     clear!(:vx)
     clear!(:vy)
     clear!(:vz)

     #.....arrays of statistical complexity
     info1=Array{Float64}(n,n,n)
     info2=Array{Float64}(n,n,n)
     info3=Array{Float64}(n,n,n)
     info1[:]=1e-30
     info2[:]=1e-30
     info3[:]=1e-30

     #.....auxiliary files to map each energy field into energy bins
    a0=Array{Int64}(n,n,n)
    a1=Array{Int64}(n,n,n)
    b0=Array{Int64}(n,n,n)
    b1=Array{Int64}(n,n,n)
    c0=Array{Int64}(n,n,n)
    c1=Array{Int64}(n,n,n)
 
     @fastmath @inbounds @simd    for i in eachindex(ed0)  
     a0[i]=Int64(round((ed0[i]-minE1)/bin1))
     a1[i]=Int64(round((ed1[i]-minE1)/bin1))
     b0[i]=Int64(round((ev0[i]-minE2)/bin2))
     b1[i]=Int64(round((ev1[i]-minE2)/bin2))
     c0[i]=Int64(round((eb0[i]-minE3)/bin3))
     c1[i]=Int64(round((eb1[i]-minE3)/bin3))
     end


     xbt=Vector{Float64}(nbin)
     xbb=Vector{Float64}(nbin)
     xbv=Vector{Float64}(nbin)

     #....matrixes of transition probabilities 
     probt=Array{Float64}(nbin,nbin)
     probv=Array{Float64}(nbin,nbin)
     probb=Array{Float64}(nbin,nbin)

     probt[:]=1e-30
     probb[:]=1e-30
     probv[:]=1e-30
     

     @inbounds @simd for i in 1:nbin-1

     id1=find(x-> (x <= i+1 && x > i),a0)
     nn=size(a0[id1])
     probt[i,:]=0. 
     probt[i,:]=prob(nn[1],a1,probt[i,:],id1)

     id1=find(x-> (x <= i+1 && x > i),b0)
     nn=size(b0[id1])
     probv[i,:]=0.
     probv[i,:]=prob(nn[1],b1,probv[i,:],id1)

     id1=find(x-> (x <= i+1 && x > i),c0)
     nn=size(c0[id1])
     probb[i,:]=0.
     probb[i,:]=prob(nn[1],c1,probb[i,:],id1)
     end  


     #....output transition probability matrix are written as FITS files
     filep=string(rooto,)
     f = FITS(filep, "w")
     write(f,probt[:,:])
     close(f)

     
     filep=string(rooto,)
     f = FITS(filep, "w")
     write(f,probv[:,:])
     close(f)


     filep=string(rooto,)
     f = FITS(filep, "w")
     write(f,probb[:,:])
     close(f)

     
    dlog2=1./(log2(2.))

     @inbounds @simd for l in eachindex(b0)
     @fastmath     info1[l]=-1.*log2(probt[(a0[l]),(a1[l])])*dlog2
     @fastmath     info2[l]=-1.*log2(probv[(b0[l]),(b1[l])])*dlog2
     @fastmath     info3[l]=-1.*log2(probb[(c0[l]),(c1[l])])*dlog2
     end


     #....the statistical complexity in 3D is written as FITS files
     filep2=string(rooto,)
     f = FITS(filep2, "w")
     write(f,info1[:,:,:])
     close(f)

     filep2=string(rooto,)
     f = FITS(filep2, "w")
     write(f,info2[:,:,:])
     close(f)


     filep2=string(rooto,)
     f = FITS(filep2, "w")
     write(f,info3[:,:,:])
     close(f)

info1=nothing
info2=nothing
info3=nothing

     return 
end



#...MAIN PROGRAM

     compute_complex(n)


