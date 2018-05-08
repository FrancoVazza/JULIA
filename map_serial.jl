#SIMPLE SERIAL MAP MAKING FOR HDF5 FILES

#in order to use any of the using <> packages one should first do Pkg.add("<package_name>") from the Julia shell                                                     

#using Devectorize    #these are possible packages for speedup                                                                                                       
#using Optim          #their helpfulness must be tested from code to code!                                                                                           


root="/homea/hhh42/hhh420/DATA/KP/E18B/256/new/nested/"  #main folder                                                                                                
file2=string(root,"5mdmd001R_dt_015")    #input files    *dt* contains gas density, dark matter density and gas temperature                                          
file1=string(root,"5mdmd001R_b_015")   # *b* contains the 3-d components of the B-field,                                                                             
using HDF5    #necessary to open hdf5 files                                                                                                                          

const convd=1 #conversion factors from ENZO code units to cgs                                                                                                        
const convv=1
const convb=1


#first we define all functions that are called by the main algorithm (down below).                                                                                   

function make_map(n)

#...here below we compute the min and max in each direction of the array.                                                                                            
#...this way of deriving the array size may seem obscure, but it will make sense in the parallel version (due to domain decomposition).                              
i1=1
j1=1
i2=n
j2=n
l1=1
l2=n
println("processor number ", myid(), " is reading the following dataset:")
println(i1," ",i2," ",j1," ",j2," ",l1," ",l2)

bz=h5read(file1,"Bz",(i1:i2,j1:j2,l1:l2))
d=h5read(file2,"Density",(i1:i2,j1:j2,l1:l2))

#...this is going to be the final map                                                                                                                                
map=Array{Float64}(n,n,2)

println("the 3D size of each file is", size(d))

#...main loop that produces the map                                                                                                                                  
for i in 1:i2-i1+1, j in 1:j2-j1+1#, l in 1:n                                                                                                                        

ds=sum(d[i,j,:])
bs=sum(bz[i,j,:])

map[i,j,1]=ds    # map of projected gas density                                                                                                                      
map[i,j,2]=bs*ds  #map of RM (gas density * Bz)                                                                                                                      
end
return map
end



#...MAIN PROGRAM                                                                                                                                                     
tic()   #...this starts a timer, closed by toc() below, useful to check performance                                                                                  

@everywhere const n=160   #...the 1D size of input grid in the HDF5 file.                                                                                            

map_total=make_map(n)  #...main function that produces the map                                                                                                       


println("The minimum and maximum values in the map are ",minimum(map_total[:,:,1])," ",maximum(map_total[:,:,1]))

ma1=map_total[:,:,1]
ma2=map_total[:,:,2]
toc()
#...very basic plotting routine using Winston                                                                                                                        
#...better alternative are available, i.e. using the Pyplot package                                                                                                  

using Winston
filep1=string(root,"map_Dens_serial.png")
imagesc(ma1)
savefig(filep1) #...the image is saved into this file                                                                                                                

filep2=string(root,"map_RM_serial.png")
imagesc(ma2)
savefig(filep2)

