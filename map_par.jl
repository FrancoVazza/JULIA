#......SIMPLE PARALLEL (WITHIN NODE) MAP MAKING TOOL IN Julia (WORKING ON HDF5 FILES)
#......tested up to 64 cores with a 1024 grid on Piz Daint
#......on a Macbook 2,6 GHz Intel Core i7 it takes ~1s for a 320^3 grid using 8 cores
#......time approximately scales as 1/nrprocs


@everywhere using DistributedArrays      
@everywhere root="/Users/francovazza/Desktop/data/DATA/CHRONOS/new_clusters/E5A/athena/"     #main folder
@everywhere file1=string(root,"Density.h5")    #input files in HDF5 format -gas density
@everywhere file2=string(root,"Temperature.h5")    #input files in HDF5 format  - temperature
@everywhere file3=string(root,"Bz.h5")   # z-component of the B-field,

@everywhere using HDF5    #necessary to open hdf5 files 


#first we define all functions that are called by the main algorithm (down below). Making them @everywhere defines them on all processors.


@everywhere function make_map(dat1)   

n1=(size(dat1[1], 1), size(dat1[2], 1),n) #size(dat1[3],1))   #this builds up the 3D structure from the dat1 array which is passed to each processor

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
 t=h5read(file2,"Temperature",(i1:i2,j1:j2,l1:l2))   
 d=h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))
 bz=h5read(file3,"Bz",(i1:i2,j1:j2,l1:l2))
#...this is going to be the final map

 map=Array{Float64}(n1[1],n1[2],3)
 ni=inv(n)

#...main loop that produces the map 
 @inbounds for i in 1:i2-i1+1, j in 1:j2-j1+1#, l in 1:n

 @fastmath ds=sum(d[i,j,:])
 @fastmath ts=sum(t[i,j,:])
 @fastmath rms=sum(dot(d[i,j,:],bz[i,j,:]))
 map[i,j,1]=ds*ni  # map of projected gas density
 map[i,j,2]=ts*ni  #map of temperature
 map[i,j,3]=rms #map of RM
 end
 return map
 end



#...MAIN PROGRAM
tic()   #...this starts a timer, closed by toc() below, useful to check performance

@everywhere n=320   #...the 1D size of input grid in the HDF5 file. 

map=DArray(make_map,(n,n,3))  #...this builds up a Distributed Array (see Julia's documentation), i.e. an array made of chunks, each passed to a different processor. The distributed array is directly passed to the make_map function defined above

map_total=convert(Array,map)  #....this step assembles a unique array by combining the different pieces of the distributed array 

println("The minimum and maximum values in the map are ",minimum(map_total[:,:,1])," ",maximum(map_total[:,:,1]))

toc()
#...very basic plotting routine using Winston
#...better alternative are available, i.e. using the Pyplot package

using Winston
filep=string(root,"map1_dens.png")
imagesc(map_total[:,:,1])  
savefig(filep) #...the image is saved into this file

filep=string(root,"map1_temp.png")
imagesc(map_total[:,:,2])
savefig(filep) #...the image is saved into this file      
filep=string(root,"map1_rm.png")
imagesc(map_total[:,:,3])
savefig(filep) #...the image is saved into this file      

using FITSIO
filep=string(root,"map_total.fits")
 f = FITS(filep, "w");
 write(f,map_total)
 close(f)



