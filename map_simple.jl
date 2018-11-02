#this program reads in an hdf5 file containing 3D physical fields (written in code units)
#converts all fields into cgs units
#produce a map of quantities in various formats (png,fits)

#...list of Julia packages (if not present, install them via Pkg.add("name") in the command line


@everywhere using HDF5      #...support to read HDF5 files
@everywhere using FITSIO    #...support to read fits files

# optional packages to optimize performance (they may speed up the code or not, depending on whethe the code is already well written or not!)
#using Devectorize  #...
#using Optim        #...

root="/home/PERSONALE/franco.vazza2/Desktop/DATA/INFO/new/20/"  #folder  with data

nam="DD20"  #...name of model
snap=372    #...number of snapshot
file1=string(root,nam,"_dt_",snap)    #input files - (Density,Dark_Matter_Density,Temperature,Bx,By,Bz)
file2=string(root,nam,"_v_",snap)    #input files - (x-velocity,y-velocity,z-velocity)


file_conv=string(root,nam,snap,".conv2")  #file with conversion factors from Enzo to cgs

#useful constants (can also go in a "module")
#@everywhere is useful when switching to parallel computing
@everywhere const pi=3.141592654
#cosmological parameters and miscellanea
@everywhere const omegab=0.0468
@everywhere const omegam=0.308
@everywhere const hh=0.678 
@everywhere const mp=1.67e-24
@everywhere const kb=1.38e-16
#.....functions



#...MAP MAKING OF A SINGLE FIELD ALONG Z AXIS
function mapsz(field::Array{Float64},lim1::Int64,lim2::Int64)   #function to make the maps of project fields 
    
    n=size(field)
    proj=Array{Float64}(n[2],n[2])
        
    @inbounds for i=1:n[2], j=1:n[2] 
        @fastmath proj[i,j]=sum(field[i,j,lim1:lim2])   #projection along z-axis
    end
    @fastmath proj=proj/(n[2])
    return proj
end


#...MAP MAKING OF A SINGLE FIELD ALONG Z AXIS
function slicez(field::Array{Float64},lim1::Int64)   #function to take slices 
    n=size(field)
    slice=Array{Float64}(n[2],n[2])
        @fastmath slice=field[:,:,lim1]   #slice
    return slice
    #..PS yes, this routine is very stupid!
end

#..MODULE OF 3D VECTOR
function mod_calc(vx::Array{Float64},vy::Array{Float64},vz::Array{Float64})
    mod=similar(vx)


    @inbounds for i in eachindex(mod)
        @fastmath mod[i]=sqrt(vx[i]^2.+vy[i]^2.+vz[i]^2.)
    end
    return mod
end



#.................
#MAIN PROGRAM
const ng=400   #1D size of 3D grid

const i1=1     #edges of subbox, in case we do not want to read the entire grid
const i2=128
const j1=128
const j2=255
const l1=1
const l2=128

#....read conversion factor file
conv=readdlm(file_conv,'\t',Float64,'\n')
@everywhere const cd=conv[4]  #...from code density to g/cm^3
@everywhere const cv=conv[5]  #...from code velocity to cm/s
@everywhere const cb=sqrt(cd*4*pi)*cv  #...from code b-field to Gauss
@everywhere const z=conv[3]  #...redshift
 
#...read Enzo datasets
d=cd*h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))   #read various data
dm=cd*h5read(file1,"Dark_Matter_Density",(i1:i2,j1:j2,l1:l2))
t=h5read(file1,"Temperature",(i1:i2,j1:j2,l1:l2))*1.0e0

bx=cb*h5read(file1,"Bx",(i1:i2,j1:j2,l1:l2))   #read magnetic field data
by=cb*h5read(file1,"Bz",(i1:i2,j1:j2,l1:l2))
bz=cb*h5read(file1,"Bz",(i1:i2,j1:j2,l1:l2)) 

println("the size of the last dataset is", size(bz))

vx=cv*h5read(file2,"x-velocity")   #...notice here we read the entire dataset


println("the size of the last dataset is", size(vx))

tic()   #...routines that gives timing informations time=toc-tic

bb=mod_calc(bx,bx,bz)

toc()

lim1=l1
lim2=l2
tic()
map_bfield=mapsz(bb,lim1,lim2)
toc()

p=d.*t   #...notice Julia notation for the element-wise multiplication of matrices!

map_pressure=mapsz(p ,lim1,lim2)

slice_vx=slicez(vx,288)

filep1=string(file2,"MAP_Bfield.fits")
f=FITS(filep1,"w")
write(f,map_bfield)
close(f)

filep1=string(file2,"MAP_Pressure.fits")
f=FITS(filep1,"w")
write(f,map_pressure)
close(f)


filep1=string(file2,"SLICE_velocity.fits")
f=FITS(filep1,"w")
write(f,slice_vx)
close(f)

using Winston
Winston.colormap("blues")  #....other color maps I've found for Winston "greens" "jet" ... else?

Winston.imagesc(slice_vx)

 filep1=string(file2,"SLICE_velocity.png")

savefig(filep1)

error()


#...all is below should work, but gives problem on some machines
#...see here for installation https://github.com/JuliaPy/PyPlot.jl
#using PyPlot
#figure()
#subplot(211)
#pcolor(slice_vx, norm=matplotlib[:colors][:LogNorm](vmin=minimum(slice_vx), vmax=maximum(slice_vx)), cmap="RdBu")
#xticks([])
#yticks([])
#colorbar()

#subplot(212)
#pcolor(map_pressure, norm=matplotlib[:colors][:LogNorm](vmin=minimum(map_pressure), vmax=maximum(map_pressure)), cmap="RdBu")
#xticks([])
#yticks([])
#colorbar()
#filep1=string(file2,"SLICE_velocity.png")
#savefig(filep1)

#show()
