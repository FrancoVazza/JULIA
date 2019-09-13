###########################################################
##### main caller of X-ray emissivity functions
##....

@everywhere function Xray_em(d::Array{Float64,3},t::Array{Float64,3},lx1::Array{Float64,3},lx2::Array{Float64,3}) 


 n3=size(d)
 filex="temp_tab_0.3_2.0kev_Julia_03.dat"
 file_spec=string(main,"/Xem/",filex)
# using do means the file is closed automatically
# in the same way "with" does in python
aa=readdlm(file_spec) 
te=aa[:,1]  #...temperature
emiss1=aa[:,2] #...free free
emiss2=aa[:,3]  #...line emission
nt2=size(te)
nt=nt2[1]

 lx1[:]=0.
 lx2[:]=0.
  de=0.88*d/(mp) #....gas density transformed into 1/cm^3                                                   
 # volx=10.^(3*log10(res*3.08e21)-40) #.....we must multiply everything for 10^40 at the end                                          
 volx=(res*3.08e24)^3.                                                             
                               
  dt=0.050000714
  Tmax=1e9
  Tmin=105925

@inbounds @simd  for i in eachindex(d)
   if t[i] > Tmin && t[i] < Tmax
   tem=convert(Int64,trunc(log10(t[i]/Tmin)/dt))
   if tem > nt 
   tem=nt
   end
   if tem < 1
   tem=1
   end
#  it=find(x->(x >= 1e5),t)
   @fastmath ded=(d[i]*0.88/mp )^2.
@fastmath   lx1[i]=(ded*emiss1[tem])*volx
@fastmath   lx2[i]=(ded*emiss2[tem])*volx
   end
   end

  return lx1,lx2                                                                                                       
                            

end
 
