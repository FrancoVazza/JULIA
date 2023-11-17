#...SPECTRAL PARAMETERS
@everywhere using SpecialFunctions
 @everywhere const  p_max=log10(1e6)
 @everywhere const  p_min=log10(2)
 @everywhere const  dp=0.1 #...was 20 for best runs
 @everywhere const  part=2  #...1=proton  2=electron
#..@everywhere constants useful to compute cooling and acceleration terms
 @everywhere const norm=30 #to have all in 1e30 erg or 1e30 erg/s
 @everywhere const evtoerg=1.60218e-12 #
 @everywhere const c1=1.3e-12 #...1/s
 @everywhere const  b1=1.37e-20 #...1/s
 @everywhere const  b2=1.9e-9  #...1/s
 @everywhere const  b3=7.2e-12 #...1/s
 @everywhere const xi=0.5 #
 @everywhere const lcorr=20.0 #...inkpc
 @everywhere const mthr=2.0 #...lower limit for shock injection of CRe - efficiency very uncertain here! 
 @everywhere const gyrtosec=1.0e9*3.0e7
 @everywhere const fi=1e-30
 @everywhere const np=1+(floor(Int64,(p_max-p_min)/dp))
 @everywhere const pend=np
 @everywhere const pval=collect(p_min:dp:p_max)

# @everywhere idp=similar(pval)
#@everywhere for i in 1:pval-1
# idp[i]=10^(pval[i+1])-10^(pval[i])
# end
 #@everywhere idp=similar(pval)
#other costants
 @everywhere const kpc=3.086e21 #cm in one kpc
 @everywhere const yr=3.154e7 #s in one yr
 @everywhere const  mp=1.67e-24
 @everywhere const  me=9.193897e-28   #g
 @everywhere const qe=4.8032e-10 #...stattcoulomb (needed in CGS! )
 @everywhere const  msol=1.988e33 #g
 @everywhere const  vc=2.99792458e10   #cm s-
 @everywhere const  erest=me*vc^2.
 @everywhere const  kpctocm=3.086e21
 @everywhere const  cmtoMpc=3.086e24
 @everywhere const  kb=1.380658e-16    #erg k-1
  @everywhere const  cj1=3*qe*0.5*2^0.5/(4*π*me*vc)
  @everywhere const  cj2=1.732*qe^3/(me*vc^2.)
  @everywhere const nang=12
  @everywhere const  dth=1/nang
  @everywhere const to=(0.5*π*dth*(collect(1:1:nang)))
  @everywhere const sinthi=sin.(to)

    @everywhere const  a1=-0.14602
    @everywhere const  a2=-0.36648
    @everywhere const  a3=9.69376e-2
    @everywhere const  bs1=-0.20250
    @everywhere const  F1=π*2^1.666/(sqrt(3.)*gamma(0.333))

    @everywhere const vc2me=vc^2*me
