@everywhere function radio(d::Float64,t::Float64,m::Float64,bx::Float64,by::Float64,bz::Float64,nu::Float64,pradio1::Float64,psi::Float64)

    xsi0=1e-4
    pradio1=0.
 const   logs=(res/(1.+z))^2.  #...pysical surface of shocked cell 
   
     m2=m*m
    @fastmath fm1=(m2+3.)/(4.*m2)
    @fastmath fm2= (16.*m2)/((5.*m2-1.)*(m2+3.))
        d1=d   #..post-shock density
        t1=t  #..post-shock temperature
        bf=sqrt(bx^2.+by^2.+bz^2.)*1e6       
        s=2.*(m2+1.)/(m2-1.)  #...injection spectral index                                         
    xsi=1e-2*psi#_sel(m,t1/keV,ff,mpsi)                         
    if xsi <=2e-5
    xsi=2e-5
    end
#     xsi=1e-4  #...electron acceleration efficiency (here fixed)
    sradio=(s-1.)*0.5 +0.5 
@fastmath    pradio=6.4e34*logs*(d1/(mp*1e-4))*xsi/(5e-2)*(t1/(7.*1.16e7))^1.5*(nu/(1.4))^(-sradio)
@fastmath    pradio1=pradio*bf^(1.+sradio)/(3.2^2.*(1.+z)^4.+bf^2.)
    return pradio1
    end
