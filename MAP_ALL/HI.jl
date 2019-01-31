###########################################################
##### main caller of HI emissivity functions
##....based on Horii+2017 MNRAS

@everywhere function HI_compute(d::Array{Float64,3},t::Array{Float64,3},dHI::Array{Float64,3}) 

 tic()
 n3=size(d)

    local  Ts=fill(0.0,n3[1],n3[2],n3[3])
@everywhere Tcmb=3.27*(1+z)#->OK
    local ν10=1.420405751e9 #Hz->OK 
    local να=2.47e15 #Hz natural frequency of HI (include redshift here?)
    local Tstar=hp*ν10/kb #->OK 

#    local TA=t #6000.0 #K
    local fα=0.4162 #5.75e-12 #0.4162#   5.75e-12 #oscillator strength !!double check
    local mu=0.59
    local mue=1.18
    local Jα=1e-22 #Lyman- background mean intensity Jα,#....double check Haardt & Madau 2012  (L
    local A10=2.85e-15   #1/s Einstein coefficient for transition -> OK

   

    local xHI=1.0
if z <= 0.5 
Jα=9.59e-21
end
if z > 0.5 && z <= 1.0 
Jα=2.89e-20
end
if z > 1.0 && z <= 2.0 
Jα=7.25e-20
end
if z > 2.0  
Jα=6.31e-20
end


@inbounds @simd for i in eachindex(t)
if t[i]>=1e3 && t[i]<=1e9
  @fastmath  np=d[i]/(mu*mp)
  @fastmath  xHI=dHI[i]/d[i]
  @fastmath  ne=d[i]/(mue*mp)
  @fastmath  nHI=dHI[i]/(mu*mp)   #xHI*np

 # if t[i] <= 1e3
#@fastmath   γHI=3.1e-11*t[i]^0.357*exp(-32./t[i])  #....cm^3/s
#  end
 # if t[i] >1e3
@fastmath  γHI=4.*3.1e-11*(t[i]/3.)^0.357*exp(-32./(t[i]/3.))
 # end 


#   if t[i]<=1e4 
#@fastmath    γe=exp(-9.607+0.5*log(t[i])*exp((-(log(t[i]))^4.5)/1800))        #...double check if log or log10!!
#! (it seems log in the article, eq.9)

#    end
#    if t[i] > 1e4
@fastmath     γe=exp(-9.607+0.5*log(1e4)*exp((-(log(1e4))^4.5)/1800))       
#@fastmath γe=10.^(-9.607+0.5*log10(1e4)*exp((-(log10(1e4))^4.5)/1800.))
#    end

    γp=3.2*γHI  
  
    Cp=np*γp     #....collisional de-excitation rate by protons 
    Ce=ne*γe     #....collisiona de-excitation by electrons    
    CHI=nHI*γHI   #....collisional de-excitation rate by HI 

    
@fastmath    chiα=fα*(pi*qe^2./(me*vc)) #...
@fastmath    γ=(Hz*να)/(nHI*vc*chiα)   #...Sobolev parameter=inverse of Gunn-Peterson optical depth
@fastmath    Sα=exp(-0.803*t[i]^(-0.666)*(1e-6/γ)^0.3333)    #...scattering amplitude factor eq.12
@fastmath    xa=1.81e11/(1+z)*Sα*Jα/(hp*να)  #coupling coefficient for Wouthuysen-Field (eq.11)

@fastmath  xc=(Tstar/(A10*Tcmb))*(CHI+Cp+Ce)  #...collisional coupling coeff. Eq.6

@fastmath @.   Ts[i]=(1.+xc+xa)/(1./Tcmb+xc/t[i]+xa/t[i])    #spin Temp. eq.5
 #         println(1./Tcmb," ",xc/t[i]," ",xa/t[i])
#@fastmath @. Ts[i]= (3./(32.*pi))*A10*hp*vc^3./(kb*ν10^2.)*(1./(1.+z))*(nHI/Hz)*(1.-Tcmb/Ts[i])
#           println("1=",Ts[i])
    @fastmath @. Ts[i]=9*nHI*(d[i]/(rhocr*0.133))*(1.+z)^0.5*(1.-Tcmb/Ts[i])
#           println("2=",Ts[i])
end
end

    return Ts
end


