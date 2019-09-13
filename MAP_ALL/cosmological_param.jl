 #cosmological parameters and miscellanea
@everywhere const Ωb=0.0468
@everywhere const Ωm=0.308
@everywhere const ΩΛ=1.-Ωm
@everywhere const hh=0.678 
@everywhere const  Hz=hh*100.*sqrt(Ωm*(1+z)^3.*ΩΛ)
@everywhere const lbox=100#200 #187/hh #40/hh #200
@everywhere const n=2400
@everywhere const res=lbox/n
@everywhere const G=6.67385e-8
@everywhere const rhocr=10.^(2.*log10(hh*1e5)-2.*log10(Mpctocm)+log10(3./(8.*pi*G)))
