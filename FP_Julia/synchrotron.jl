#CONSTANTS                                                                                                                         
@everywhere const kpc=3.086e21 #cm in one kpc                                                            
@everywhere const yr=3.154e7 #s in one yr                                                                
@everywhere const  mp=1.67e-24
@everywhere const  me=0.1093897e-28   #g                                                                
@everywhere const  vc=2.99792458e10   #cm s-                                                             
@everywhere const  kpctocm=3.086e21
@everywhere const  cmtoMpc=3.086e24
@everywhere const  kb=1.380658e-16    #erg k-1                                                                                     
#MODULES                                                                                                                           
@everywhere using Optim
@everywhere using Base
#@everywhere using Devectorize                                                                                                     
@everywhere const  part=2  #...1=proton  2=electron                                     
#..constants useful to compute cooling and acceleration terms
@everywhere const norm=30 #to have all in 1e30 erg or 1e30 erg/s                                        
@everywhere const erest=0.510988946 #electron rest mass energy                                           
@everywhere const evtoerg=1.60218e-12 #                                                                  
@everywhere const  b1=1.37e-20 #...1/s                                                                   
@everywhere const  b2=1.9e-9  #...1/s                                                                    
@everywhere const  b3=7.2e-12 #...1/s 
@everywhere const re=2.81e-13  #...electron radius

@everywhere const  keverg=1.602e-9 #
@everywhere const  e=4.8e-10 #....esu
@everywhere const gyrtosec=1.0e9*3.0e7

@everywhere function bessel_integ1(vv::Float64)
    totf=0.0
    dx=0.2e-3
    ymin=convert(Int64,trunc(vv/dx))
#    totf=Array{Float64}(1)
#    totf=Array{Float64}(1)
  #  totf[1]=0. 

@simd     for y in ymin:2e5  #....K53 modified bessel function                                
            @fastmath  argb=y*dx#(1.+y*dx
    #if argb >= 30       
     @fastmath  bbes=besselk(1.666,argb)
      
             totf+=(vv*bbes*dx)
    # end 
            end

   return totf
   end

@everywhere function bessel_integ2(vv::Float64)
#   ak1=[-0.97947838884478688,-0.83333239129525072,0.15541796026816246]
#   ak2=[-4.69247165562628882e-2,-0.70055018056462881,1.03876297841949544e-2]
   ak1=[-1.0194198041210243,+0.28011396300530672,-7.71058491739234908]
   ak2=[-15.761577796582387]

    h1=0.
    h2=0.
    vi=0.
   for k in 1:3
   vi=vv^(1./k)
   h1+=(ak1[k]*vi)
   end
   h2=ak2[1]*vv
#   end
   ffn=exp(h1)+(1.-1.*exp(h2))
   println(ffn)
   return ffn

end


@everywhere function bessel_integ3(vv::Float64)
#   ak1=[-0.97947838884478688,-0.83333239129525072,0.15541796026816246]         
#   ak2=[-4.69247165562628882e-2,-0.70055018056462881,1.03876297841949544e-2]   
   a1=-0.14602
   a2=-0.36648
   a3=9.69376e-2
   b1=-0.20250
   p=3.
   F1=π*2.^1.666/(sqrt(3.)*gamma(0.333))
   κp=2.*F1/(p-0.333)
   Cp=2.*((p+1.)*0.5)/(p+1.)*gamma(p*0.25+19./12.)*gamma(p*0.25-1./12.)
   ffn=κp*vv^0.333/(1.+(κp/Cp)*vv^((p-0.333)*0.5))

   return ffn

end

@everywhere function synchrotron(gammaval::Vector{Float64},b0::Float64,pe::Vector{Float64},freq::Float64)

psync=0.0
#dp=Array{Float64,ngamma}
g2=ngamma-1
g1=1
  
@fastmath  c1=3.*e*0.5*2.^0.5/(4*π*me*vc)
@fastmath  c2=1.732*e^3./(me*vc^2.)
  b0*=1e-6
 #.,,,,compute syncrothron emissivity
     nang=10
     dth=1./nang

    to=0.5*π*dth*(collect(1:1:nang))
    integ=0.0
    freqc=0.0
  @inbounds  @simd for gg in g1:g2 #...gamma factors
               gu=gammaval[gg]
  @fastmath    freqc=c1*b0*(gu*gu)
  @fastmath    vv=freq/(freqc)
  @inbounds   @simd  for ang in 1:nang #....angles
        thi=to[ang]
        sinthi=sin(thi)

#        ffn1=bessel_integ1(vv)     
        ffn=bessel_integ3(vv)
       

#        totf=0.0
        
#  @inbounds    @simd    for y in 1:100  #....K53 modified bessel function
#            @fastmath  argb=(1.+y*0.005)*vv
#             bbes=besselk(1.666,argb)
#             totf+=bbes*0.25
#             end
       
#          ffn=vv*totf


         @fastmath   integ+=pe[gg]*ffn*sinthi*sinthi*dth
         end
       @fastmath    psync+=integ*(b0*c2)
    end


 return psync
end






