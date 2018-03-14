#...THIS IS THE SIMPLIFIED JULIA VERSIONE OF THE "MIRO" CODE TO PRODUCE MOCK RM MAPS FROM 3D MAGNETIC FIELDSGENERATED IN FOURIER SPACE. SEE https://github.com/FrancoVazza/MIRO for the IDL version and more comments.


#using Optim
#using ParallelAccelerator
#using Devectorize
using FFTViews
using DistributedArrays

@everywhere using DistributedArrays
@everywhere const n = 512 #grid size
@everywhere const index = 1.666 #1D power spectrum (negative index)
@everywhere const pig = 3.141592653
@everywhere const beta  = 0.75
@everywhere const rc =290.0 #kpc
@everywhere const ncentre =3.44e-3  #central core density
@everywhere const centre=[n*0.5,n*0.5,n*0.5] #center of cluster in box units
@everywhere const res=10.0 #cell resolution in kpc
@everywhere const bnorm=4.7 #central B-field in muG
@everywhere const kmin=1  #lowest k for power spectrum
@everywhere const kmax=Int(floor(n*0.5)) #largest k for power spectrum

@everywhere function pk(index,kk)
return @fastmath kk.^(-index-2)
end


@everywhere function normk(l,zc,j,yc,i,xc)
return @fastmath sqrt((l-zc)*(l-zc)+(j-yc)*(j-yc)+(i-xc)*(i-xc))
end

@everywhere function rayleigh(sigma,x)
return @fastmath sigma*sqrt(-2.*log(x))
end

@everywhere function cross(i1,i2,v1,v2)
 return @fastmath i1*v1-i2*v2
end

#@everywhere function cross!(c,a,b)
#    c[1] = a[2]*b[3]-a[3]*b[2]
#    c[2] = a[3]*b[1]-a[1]*b[3]
#    c[3] = a[1]*b[2]-a[2]*b[1]
#end

@everywhere function miro(I)

  d=(size(I[1], 1), size(I[2], 1),size(I[3],1))
  bkz=Array(Complex64,d)


  imin=I[1][1]
  jmin=I[2][1]
  lmin=I[3][1]

 @inbounds @fastmath for i=I[1], j=I[2],l=I[3]

  kk=normk(l,centre[1],j,centre[2],i,centre[3])   
 
  am=0
  if kk >= kmin && kk<kmax


  am=rayleigh(sqrt(pk(index,kk)),rand())


  ph1=rand()*2.*pi
  ph2=rand()
  ph3=rand()*pi

  c2=cos(ph2)
  s2=sin(ph2)
  c3=cos(ph3)
  s3=sin(ph3)
  c1=cos(ph1)
  s1=sin(ph1)

   ampx=abs(am*s2*c3)
   ampy=abs(am*s2*s3)
   ampz=abs(am*c2)
cc=Complex64(c1,s1)
   akx=ampx.*cc
   aky=ampy.*cc
   akz=ampz.*cc

c11=real(cross(j,l,akz,aky))
c21=imag(cross(l,i,akx,akz))

   in1=i-imin+1
   in2=j-jmin+1
    in3=l-lmin+1
 bkz[in1,in2,in3]=Complex64(c11,c21)
 end


end

return bkz

end

tic()
Bdistr=DArray(miro,(n,n,n))

bkz=convert(Array,Bdistr)

FFTW.set_num_threads(nprocs())

bz=real(ifft(bkz))


toc()
using Winston
#using Plots
#using Images
filep6="/Users/francovazza/Desktop/Julia_prog/JUliAN_MIRO/Bfield_cut_p4.png"
zc=Int(floor(centre[3]))
ima=Array{Float64}(n,n)
ima=log10.(abs.(bz[:,:,zc]))

mi=minimum(ima)
ma=maximum(ima)
ima=256.0*(ima-mi)/(ma-mi)
println(sizeof(ima))
println(typeof(ima))
imagesc(ima[:,:,1])

savefig(filep6)
