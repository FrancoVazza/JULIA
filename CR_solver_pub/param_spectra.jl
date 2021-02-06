#...SPECTRAL PARAMETERS

     const  p_max=11000.
     const  p_min=1.0
     const  dp=20.0
     const  part=2  #...1=proton  2=electron
#..    constants useful to compute cooling and acceleration terms
     const evtoerg=1.60218e-12 #
     const c1=1.3e-12 #...1/s
     const  b1=1.37e-20 #...1/s
     const  b2=1.9e-9  #...1/s
     const  b3=7.2e-12 #...1/s
     const xi=0.5 #
     const lcorr=20.0 #...inkpc
     const mthr=2.25
     const gyrtosec=1.0e9*3.0e7
     const fi=1e-3

     const np=(floor(Int64,(p_max-p_min)/dp))
     const pend=np
     const pval=collect(p_min:dp:p_max)
#other costants
     const kpc=3.086e21 #cm in one kpc
     const yr=3.154e7 #s in one yr
     const  mp=1.67e-24
     const  me=0.1093897e-28   #g
     const  msol=1.988e33 #g
     const  vc=2.99792458e10   #cm s-
     const  erest=me*vc^2.
     const  kpctocm=3.086e21
     const  cmtoMpc=3.086e24
     const  kb=1.380658e-16    #erg k-1
