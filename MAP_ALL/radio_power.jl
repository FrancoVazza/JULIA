@everywhere function radio(d::Float64,t::Float64,m::Float64,bx::Float64,by::Float64,bz::Float64,nu::Float64,pradio1::Float64,pradio2::Float64,psi::Float64)

    xsi0=1e-4
    pradio1=0.
    pradio2=0.
 const   logs=(res/(1.+z))^2.  #...pysical surface of shocked cell 
   
     m2=m*m
    @fastmath fm1=(m2+3.)/(4.*m2)
    @fastmath fm2= (16.*m2)/((5.*m2-1.)*(m2+3.))
        d1=d*fm1   #..pre-shock density
        t1=t*fm2   #..pre-shock temperature
        bf=sqrt(bx^2.+by^2.+bz^2.)*1e6       
        s=2.*(m2+1.)/(m2-1.)  #...injection spectral index                                         
    xsi=1e-2*psi#_sel(m,t1/keV,ff,mpsi)                         
    if xsi <=2e-5
    xsi=2e-5
    end
#     xsi=1e-4  #...electron acceleration efficiency (here fixed)
    sradio=(s-1.)*0.5 +0.5 
    pradio=6.4e34*logs*(d1/(mp*1e-4))*xsi/(5e-2)*(t1/(7.*1.16e7))^1.5*(nu/(1.4))^(-sradio)
    pradio1=pradio*bf^(1.+sradio)/(3.2^2.*(1.+z)^4.+bf^2.)
    pradio2=pradio1*xsi0/xsi
    println(pradio1," ",pradio2)
    return pradio1,pradio2
    end

@everywhere function fpradio(d::Float64,t::Float64,mach::Float64,bx::Float64,by::Float64,bz::Float64,nu::Float64,pradio::Float64) 
    xsi=1e-4
    
    pradio=0.
    logs=(res/(1.+z))^2.  #...pysical surface of shocked cell 
    println("here")
    println(size(d))
    
    #@inbounds  @simd  for i in eachindex(d)
        
        m2=mach*mach
    fm1=(m2+3.)/(4.*m2)
    fm2= (16.*m2)/((5.*m2-1.)*(m2+3.))
        d=d*fm1   #..pre-shock density
        t=t*fm2   #..pre-shock temperature
   bf=sqrt(bx^2.+by^2.+bz^2.)*1e6       
   s=2.*(m2+1.)/(m2-1.)  #...injection spectral index                                         
                                                                                                                                        xsi=psi_sel(m,t[i]/keV)   
   sradio=(s-1.)*0.5 +0.5 
   pradio=6.4e34*logs*(d*mw_pm/(1e-4))*xsi/(5e-2)*(t/(7.*keV))^1.5*(nu/(1.4))^(-sradio)
   pradio=pradio*bf^(1.+sradio)/(3.2^2.+bf^2.)
   println(pradio)

    return pradio

end

@everywhere function fpradio2(d::Vector{Float64},t::Vector{Float64},mach::Vector{Float64},bx::Vector{Float64},by::Vector{Float64},bz::Vector{Float64},nu::Float64,pradio::Vector{Float64}) 
    xsi=1e-4
    pradio=similar(d)
    pradio[:]=0.
    logs=(res/(1.+z))^2.  #...pysical surface of shocked cell 
    println("here")
    println(size(d))
    error()
    @inbounds  @simd  for i in eachindex(d)
        m=mach[i]
        m2=m*m
    fm1=(m2+3.)/(4.*m2)
    fm2= (16.*m2)/((5.*m2-1.)*(m2+3.))
        d[i]=d[i]*fm1   #..pre-shock density
        t[i]=t[i]*fm2   #..pre-shock temperature
   bf=sqrt(bx[i]^2.+by[i]^2.+bz[i]^2.)*1e6       
   s=2.*(m2+1.)/(m2-1.)  #...injection spectral index                                         
                                                                                                                                        xsi=psi_sel(m,t[i]/keV)    
   sradio=(s-1.)*0.5 +0.5 
   pradio[i]=6.4e14*logs*(d[i]*mw_pm/(1e-4))*xsi/(5e-2)*(t[i]/(7.*1.16e7))^1.5*(nu/(1.4))^(-sradio)
   pradio[i]=pradio[i]*bf^(1.+sradio)/(3.2^2.+bf^2.)
   println(pradio[i])
    end 
   stop


    return pradio

end



@everywhere  function fsel(psi_m_all,mpsi)
# ............reads in the Psi(M) function of Hoeft & Bruggen
fhb="mach_psi_table.txt"
nt=13
f1=string(fhb)
ff=readdlm(f1)

n=size(v)
mpsi=v[:,1]
ff=v[:,2:n[2]]

    return ff
 end




@everywhere function psi_sel(m,t,ff,mpsi)

const nt=13
const n33=244

th=[1.00e-04,3.162e-04,1.000e-03,3.162e-03,1.000e-02,3.162e-02,1.000e-01,3.162e-01,1.000e+00,3.162e+00,1.000e+01,3.162e+01,1.000e+02]
#...................selection from the HB efficiency function for Psi(M)
to=t

 for ii in 1:nt-1 
if to > th[ii] && to < th[ii+1]
tsel=ii
end
end


for ii in 1:n33-1
if m >= mpsi[ii] && m < mpsi[ii+1]
msel=ii
end
end

if to < th[1] 
tsel=1
end
if to > th[nt] 
tsel=nt
end
if m < mpsi[1] 
msel=1
end

if m > mpsi[n33] 
msel=n33
end
#(psi_m_all,mpsi)
#    psi=
psi=ff[tsel,msel]# fsel(tsel,msel)
println(tsel," ",msel," ",psi)
 return psi

end
