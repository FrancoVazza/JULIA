
fold="/Users/francovazza/Desktop/data/DATA/CHRONOS/new_clusters/E18B/snap/"

 file1="Bx.fits"  #input files
 file2="By.fits"
 file3="Bz.fits"

 
using FITSIO


 fbx=FITS(string(fold,file1),"r")
 fby=FITS(string(fold,file2),"r")
 fbz=FITS(string(fold,file3),"r")
  bx=read(fbx[1],100:299,100:299,100:299)    #focus on some innermost region
  by=read(fby[1],100:299,100:299,100:299)
  bz=read(fbz[1],100:299,100:299,100:299)

  bxr=similar(bx) 
  byr=similar(by)
  bzr=similar(bz)
using CoordinateTransformations
using StaticArrays
using Rotations
#for aa in 1:50    #...to produce movies

ang1=0.0
ang2=0.0
ang3=-π/2.*70./90.
 rot=LinearMap(RotX(ang3))
 x = SVector(100, 100, 100)
 rot_around_x = recenter(rot, x)
  n=200
  for l in 1:n
  for j in 1:n         
  @simd for i in 1:n
  @inbounds       p3=[bx[i,j,l],by[i,j,l],bz[i,j,l]]
                  x3=[i,j,l]
  
        pnew=rot(p3)
        xnew=rot_around_x(x3)   
    x1=convert(Int32,trunc(xnew[1]))
      x2=convert(Int32,trunc(xnew[2]))
      x3=convert(Int32,trunc(xnew[3]))
      if x1 < 1
      x1=1
      end
       if x2 < 1
      x2=1
      end
 if x3 < 1
      x3=1
      end
 if x1 > n
      x1=n
      end
if x2 > n
      x2=n
      end
if x3 > n
      x3=n
      end
  @inbounds     bxr[x1,x2,x3]=pnew[1]
  @inbounds     byr[x1,x2,x3]=pnew[2]
  @inbounds     bzr[x1,x2,x3]=pnew[3]
       end
       end
       end
  bx=nothing
  by=nothing
  bz=nothing

  #...maps of field along the LOS 
 mapx=Array{Float64}(n,n)
 mapy=Array{Float64}(n,n)
 mapz=Array{Float64}(n,n)
mapz[:]=0.
mapy[:]=0.
mapx[:]=0.
@inbounds @simd for i in 1:n
  @views mapx[:,:]+=bxr[:,:,i]
  @views mapz[:,:]+=bzr[:,:,i]
  @views mapy[:,:]+=byr[:,:,i]
  end

  for i in eachindex(mapx)
  mapx[i]=sqrt(mapx[i]^2.+mapy[i]^2.+mapz[i]^2.)
  end


file_out=string(fold,"RMx_test2.fits")

f = FITS(file_out, "w");
  write(f,mapx)
close(f)

    error()
file_out=string(fold,"RMz_test2.fits")
f = FITS(file_out, "w");
  write(f,mapz)
close(f)

file_out=string(fold,"RMy_test2.fits")
f = FITS(file_out, "w");
  write(f,mapy)
close(f)

error()
#...this is to write out the rotated 3D fileds as fits files.
file_out1=string(fold,"Bx_rot_test.fits")
file_out2=string(fold,"By_rot_test.fits")
file_out3=string(fold,"Bz_rot_test.fits")

 f = FITS(file_out1, "w");
  write(f,bxr)
  close(f)

  
f = FITS(file_out2, "w");
  write(f,byr)
  close(f)

f = FITS(file_out3, "w");
  write(f,bzr)
  close(f)
