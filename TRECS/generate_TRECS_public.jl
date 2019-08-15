#.......test code to produce lightcones for TRECS simulations

using PyPlot
using PyCall
using LaTeXStrings
using Distributions

  lbox=200.
  n0=2400
  root="/Users/francovazza/Desktop/data/DATA/CHRONOS++/200Mpc/2400/halos/"
  snaps=["7","6","5"]
  #...the 200^2 Mpc^2 FOV is divided into 4x4 tiles
  #...at z=0.025, 200 Mpc -> 130 degrees

  dl=10.68     #..Mpc
  fov=5 #..degrees
  adist0=50. #...Mpc distance at z=0.01

#...sequence of redshif snapshots as in TRECS catalog
  zz=[0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00,1.20,1.40,1.60,1.80,2.00,2.20,2.40,2.60,2.80,3.00,3.20,3.40,3.60,3.80,4.00,4.20,4.40,4.60,4.80,5.00,5.20,5.40,5.60,5.80,6.00,6.20,6.40,6.60,6.80,7.00,7.20,7.40,7.60,7.80,8.00,8.20,8.40,8.60,8.80,9.00,9.20,9.40,9.60,9.80,10.00]

#.....file with pre-computed luminosity distance
filed=string(root,"distl_z.dat")
a=readdlm(filed)
distl=a[:,2]
a=nothing

#.....file with halo mass function from TRECS as a sanity check
filed=string(root,"TRECS_halos.dat")
trecs=readdlm(filed)
ntrecs=22

  nt=15
file_dat1=string(root,"snapshot_slices_new2.txt")     
sol=open(file_dat1,"w")
writedlm(sol,["      snap        red           xmin           xmax            dy  "])

file_sim=string(root,"snapshot_slices_sim.txt")
so3=open(file_sim,"w")
writedlm(so3,["      snap        red      xmin       xmax      ymin     ymax     zmin     zmax    "])


  da_dist=1757 #Mpc
  dl_dist=15818 #Mpc
  dl_comov=5400 #...Mpc

  initial_dist=44 #...Mpc initial distance (consistent with zmin)
  zmax=2.16     #...16x 200 Mpc = 3200 Mpc (comoving) ->  zmax~1.0
  zmin=0.01
  dz=(zmax-zmin)/(nt)  #...a bit simplistic, but fair enough for halos...
  fi=1
  nima=512
  ima_gal=Array{Float64,3}(nima+1,nima+1,nt)
  ima_gal[:,:,:]=0.
  lbox2=lbox

#....arrays for halos mass function
 nh=22
 xh=Array{Float64}(nh)
 massh=Array{Float64}(nh,ntrecs)
 massh[:]=0.
 for i in 1:nh
 xh[i]=9.5+0.25*i
 end


clf()
#.....loop over snapshots
ib=1
 x0=Array{Float64}(1)
 y0=Array{Float64}(1)
 z0=Array{Float64}(1)
 logm0=Array{Float64}(1)
 nbig0=0
#.....MAIN LOOP OVER REDSHIFT BINS
for ii in 1:nt-1 

file_gal=string(root,"cone_5x5_z",zz[ii],"txt.dat")
so2=open(file_gal,"w")


  ngal_trecs=sum(trecs[ib,3:24])  #...number of galaxies in the corresponding TRECS redshif bin

if zz[ii] >=0 && zz[ii] <=0.5
    file=string(root,"gal_cat_radio7.dat")
end

if zz[ii] >0.5 && zz[ii] <=1.0
    file=string(root,"gal_cat_radio6.dat")
end

if zz[ii] > 1.0 
    file=string(root,"gal_cat_radio5.dat")
end

    a=readdlm(file)
    close(file)

    logm=a[:,1]
    x=a[:,2] 
    y=a[:,3]
    z=a[:,4]
    lp1=a[:,5]
    lp2=a[:,6]
 
    a=nothing
    x0=x
    y0=y
    z0=z
    logm0=logm
           
        zz_pre=zz[ii]
        if ii==1
        zz_pre=0.
        end
        zz_post=zz[ii+1]
        zz_mean=(zz_post+zz_pre)*0.5
        dzz=zz_post-zz_pre
        zeds=zz_mean
        println("doing tile ",fi,"with mean redshift z=",zeds)
        ldist=distl[ii]
        adist=ldist/(1.+zz[ii]^2.)

     writedlm(sol,[ii zz[ii] distl[ii] distl[ii+1] fov ])
        
        dla=dl
        dla=dl*(adist/(adist0))

        x1=0.
        x2=dla
        y1=0.
        y2=dla


     writedlm(so3,[ii zz[ii] x1 x2 y1 y2 distl[ii] distl[ii+1] fov ])

        println(x2," ",y2, "Mpc")

    #....computing whether a replica of the galaxy distribution (assuming periodicity) is needed
    nbig=convert(Int32,trunc(x2/(lbox2+0.1)))
    if nbig>=1 
     println("replicating volume ",nbig, "times" )
  @inbounds   for l in 1:nbig       

  @views        append!(x,x0+l*lbox2) 
  @views        append!(x,x0+l*lbox2)
          append!(x,x0)

          append!(y,y0+l*lbox2)
          append!(y,y0)
          append!(y,y0+l*lbox2)

          append!(z,z0)
          append!(z,z0)
          append!(z,z0)
          append!(logm,logm0)
          append!(logm,logm0)
          append!(logm,logm0)
        end
       nbig0=nbig
        end

        ngg=size(x)

        ngal=0
    leap=1
    if zz[ii]==0.3 || zz[ii]==0.25
    leap=3
    end
    if zz[ii]>0.3 
    leap=2
    end
    if zz[ii]>=0.6
    leap=4
    end

  
     @inbounds  for i in 1:leap:ngg[1]   #....loop over the (augmeted by replication, or not) galaxy distribution 
     
        if x[i]>=x1 && x[i]<=x2 && y[i]>=y1 && y[i]<=y2 
        ngal+=1
      
      #....the following is the attempt of fixing the overaboundance of >1e14 halos in my input catalog, compared to TRECS, by introducing an ad-hoc exponential cut.
         if logm[i] > 14.
          bor=rand(1)
         logm[i]= logm[i]*1./(1.+exp(-5*bor[1]))
         end
      @fastmath    ix=convert(Int32,trunc((logm[i]-9.5)/0.25))

        if ix <=1
        ix=1
        end
        if ix>nh
        ix=nh
        end

   #....this is to produce mass functions within the same redshift bins of the ones produced by Torrence

       ib=convert(Int32,trunc(zz[ii]/0.1))+1
       massh[ix,ib]+=1.
       
      ntotal_bin=sum(massh[:,ib])

            xo=(x[i]-x1)/(x2-x1)
            yo=(y[i]-y1)/(y2-y1)
            xu=convert(Int64,trunc((nima*xo)))
            yu=convert(Int64,trunc((nima*yo)))

       xs=convert(Float64,fov*(xo-0.5))
       ys=convert(Float64,fov*(yo-0.5))
       zs=convert(Float64,zz_pre+dzz*z[i]/lbox)

    #....WRITING OF THE redshift_slice log
      writedlm(so2,[logm[i] zs xs ys])
    

    #......MAP MAKING    
        if xu>nima
         xu=nima
        end
       if yu>nima
        yu=nima
       end      
       if xu <=1 
       xu=1
       end
       if yu <=1
       yu=1
       end

         ima_gal[xu,yu,ii]+=logm[i]

       end
       end    
      println("ngalaxies=",ngal)
       

 #.....compare halo mass function with TRECS

         semilogy(xh,massh[:,ib],alpha=0.7,marker="o")   #....VAZZA
          plot(xh,trecs[ib,3:24],alpha=0.5)              #....TRECS
        axis("tight")
   ylabel(L"N_{\rm gal}",size=15)
   xlabel(L"log10(M_{\odot})",size=15)
   axis([9.,15., 1, 1e6])
   annotate(xy=[13,1e4],string("z=",zz[ii]),color="black",size=14,alpha=1)
   savefig(string(root,"mass_function_z",zz[ii],".png"))
   clf()
 
   close(so2)
        end        
  

#....write galaxy map on a fits file 
     using FITSIO
  
 filep4=string(root,"map_gal.fits")
 f = FITS(filep4, "w");
 write(f,ima_gal)
 close(f)

 close(sol)
 close(so3)
