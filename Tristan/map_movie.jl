
#...MAIN PROGRAM                                                                                                                                             
#@everywhere using DistributedArrays                                                                                                                        
@everywhereusing HDF5
@everywhereusing FITSIO
 
sna="0"
s=1
ns=200   #...number of snapshots
n2=64     #2nd dimension size
n1=1024       #1st dimension size                                                                                                        
movie=Array{Float64}(n1,n2,ns)                                                                                                     
movie[:]=0.
 for s in 1:ns
tic()
@everywhere    sna=string(s)
    if s < 100 && s > 9
@everywhere    sna=string("0",s)
    end
    if s < 10
@everywhere    sna=string("00",s)
    end
 
println("reducing snapshot ",sna)
 
@everywhere root="/marconi_scratch/userexternal/fvazza00/Tristan/run1/output/"
 
@everywhere file1=string(root,"flds.tot.",sna)    #input files    *dt* conta\ins gas density, dark matter density and gas temperature                      \
  
f=h5read(file1,"bx",(1:n1,1:n2,1))
 
println("The minimum and maximum values in the map are ",minimum(f[:,:,1])," ",maximum(f[:,:,1]))
movie[:,:,s]=f[:,:,1]
  
toc()
end
 
filep1=string(file1,"_bx_movie.fits")
f=FITS(filep1,"w")
write(f,movie)
close(f)
