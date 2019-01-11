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
@everywhere using HDF5
@everywhere using FITSIO
#@everywhere using Optim
@everywhere using Base
@everywhere using Devectorize                                                                                                     
@everywhere const  g_max=50000.0
@everywhere const  g_min=50.0
@everywhere const  dg=100.0
@everywhere const  part=2  #...1=proton  2=electron                                     
#..constants useful to compute cooling and acceleration terms
@everywhere const norm=30 #to have all in 1e30 erg or 1e30 erg/s                                        
@everywhere const erest=0.510988946 #electron rest mass energy                                           
@everywhere const evtoerg=1.60218e-12 #                                                                  
@everywhere const  b1=1.37e-20 #...1/s                                                                   
@everywhere const  b2=1.9e-9  #...1/s                                                                    
@everywhere const  b3=7.2e-12 #...1/s                                                                                                          
@everywhere const xi=0.5 #                                                                                                                     
@everywhere const lcorr=20.0 #...in[100]kpc                                                                                                   
@everywhere const mthr=2.0
@everywhere const gyrtosec=1.0e9*3.0e7
@everywhere const fi=1e-3
@everywhere const ngamma=(floor(Int64,(g_max-g_min)/dg))

@everywhere const gend=ngamma
@everywhere const gammaval=collect(g_min:dg:g_max)
pe=Array{Float64}

@everywhere function age(idt::Float64,idg::Float64,aa1::Float64,bb::Float64,g1::Float64,q1::Float64,nn1::Float64,nn2::Float64,aa2::Float64,g2::Float64)
#@everywhere function age(idt,idg,aa1,bb,g1,q1,nn1,nn2,aa2,g2)
return  ((inv(idt+idg*(aa1+bb*g1^2.)))*(idt*q1+nn1*idt+nn2*idg*(aa2+bb*g2^2.)))
end

@everywhere function coul(cou::Float64,g1::Float64,inth::Float64)
return cou*(1.+0.0133333*(log(g1*inth)))
end

@everywhere function bloss(p1::Int64,b2::Float64,b0::Float64,b3::Float64,zz4::Float64)

return  p1*b2*(0.666*b0^2.*1.0e-12+b3*zz4)
end

@everywhere function reacc(delta::Float64,nn::Array{Float64,1},gend::Int64,gmin_re::Int64,gam::Array{Float64,1},n_re::Array{Float64,1})
 ntemp=similar(nn)
 gtemp=similar(nn)
  @simd for j in 1:ngamma #gmin_re:ngamma                                                                   
@inbounds @views ntemp=nn[gmin_re:j]
@inbounds @views gtemp=gam[gmin_re:j]

@fastmath  @inbounds   n_re[j]=(delta+2.)*(gam[j])^(-1.*delta)*(dg*sum(dot(ntemp,gtemp)^(delta-1.)))
#nn[gmin_re:j],(gam[gmin_re:j])^(delta-1.))))
  end
  return n_re
  end


@everywhere function reacc2(delta::Float64,nn::Array{Float64,1},gend::Int64,gmin_re::Int64,gam::Array{Float64,1},n_re::Array{Float64,1})
 ntemp=similar(nn)
 gtemp=similar(nn)
 nre=0.0
  @simd for j in 1:ngamma #gmin_re:ngamma                                                                                                                                                                  
@inbounds @views ntemp=nn[gmin_re:j]
@inbounds @views gtemp=gam[gmin_re:j]
     nre=0.0
    @simd     for i in gmin_re:j
@. @fastmath  @inbounds   nre+=(delta+2.)*(gam[j])^(-1.*delta)*(dg*ntemp[i]*gtemp[i])^(delta-1.)

   end
   n_re[j]=nre  
end
  
  return n_re
  end


@everywhere function intej(delta::Float64,gam_c::Float64)
return  (-delta+2.0)*inv(((gam_c)^(-delta+2.0)-(g_min)^(-delta+2.0)))
end


@everywhere function kej(ecr::Float64,inte::Float64)
 return  (ecr)*inv(1e6*erest*evtoerg)*inv(inte)
end


@everywhere function ninji(Ke::Float64,gamj::Float64,delta::Float64,gam_c::Float64)
return  Ke*(gamj^(-delta))*(1.-gamj/gam_c)^(delta-2.)
end

@everywhere function tacci(vshock::Float64,delta::Float64,beff::Float64)
return  inv(vshock*inv(3000.0))*3.0e7*2.4e4*sqrt(delta*inv(beff)) #..acceleration time from Kang 2012 Eq\                                                                                        
end



@everywhere function sh_gaini(f::Float64,gami::Float64,vpre::Float64,cons2::Float64,diff::Float64)
return (inv(f)*gami*(vpre)^2.*cons2*inv(3.*diff))
end

@everywhere function evolve_spectrum(zz::Float64,v::Float64,t2::Float64,t1::Float64,nth::Float64,m::Float64,b0::Float64,ecr::Float64,nn::Array{Float64,1},shock::Int64,delta_t::Float64,volume::Float64,gam::Array{Float64,1},g_max::Float64,g_min::Float64,dg::Float64,ngamma::Int64,part::Int64,vshock::Float64,tacc::Float64,n_inj::Array{Float64,1},q_inj::Array{Float64,1},n_re::Array{Float64,1})

gend=ngamma
@fastmath m2=m^2.
@fastmath vpre=1.0e5*v
@fastmath f=(4.*m2)/(m2+3.)
@fastmath delta=2.*(m2+1.)/(m2-1.)
@fastmath zz4=(1.+zz)^4.

  ic_lose=0.0
  diff=0.0
  sh_gain=0.0

@fastmath const dt=1.0*(t2-t1)
  #losses                                                                                                          

const  b3=7.2e-12 #...1/s                                                                                          
@fastmath const   cou=nth*1.2e-12 #..1/s                                                                           
@fastmath const  idt=inv(dt)
@fastmath const  idg=inv(dg)
@fastmath const  icou=inv(cou)
@fastmath const  inth=inv(nth)
@fastmath const  cons1=2.3e29*(erest*1e-3)^(0.3333)*b0^(-0.3333)*(lcorr*0.05)^(0.66666)
@fastmath const  cons2=(f-1.)*inv(f+1.)
@fastmath const  dg2=dg*0.5

#nno=SharedArray{Float64}(ngamma)                                                                                  
 #.....SHOCK ACCELERATION OF FRESH PARTICLES                                                                       

 #                                                                                                                  

  if shock >= 1 && m >= mthr   #injection of power-law spectrum of particles                                       

   beff=sqrt(b0+3.2*(1+zz)^2.)
   tacc=tacci(vshock,delta,beff)

@inbounds    gam_c=gam[ngamma]

@inbounds @simd   for i in 1:ngamma   #finds maximum gamma where tacc < tlosses                     
            gami=gam[i]
            ic_lose=b1*((gami^2.))*zz4
            diff=(gami-1.0)*cons1 # in cm^2/s                                                                             
            sh_gain=sh_gaini(f,gami,vpre,cons2,diff)

          if ic_lose > sh_gain
          gam_c=gami
           end
               end

 inte=intej(delta,gam_c) 
  Ke=kej(ecr,inte)


@inbounds @simd    for j in eachindex(n_inj)
        n_inj[j]=ninji(Ke,gam[j],delta,gam_c)

        if isnan(n_inj[j]) || isinf(n_inj[j])
                n_inj[j]=1e-60
         end
       q_inj[j]=tacc*n_inj[j]
               end

#   id=find(x-> (x <= 1e-60),n_inj)
#   n_inj[id]=1e-60
#   q_inj[id]=0.

           end



#...SHOCK REACCELERATION                                                                                           
  if shock == 2
   n_re[1]=nn[1]
   gmin_re=2
   n_re=reacc(delta,nn,gend,gmin_re,gam,n_re)
   nn=n_re
 end


#setting to zero losses that are not valid either for protons (p=1) or electrons (p=2)                             
    p1=1
    p2=1
    if part==1
    p1=0
    end
    if part==2
    p2=0
    end
 #   toc()                                                                                                         

#   tic()                                                                                                          
#....INTEGRATION OF LOSSES AND INCLUSION OF FRESHLY INJECTED PARTICLES                                             
@inbounds @simd  for j in 1:ngamma-1
      gg=gend-j
      ga=gam[gg]
      g1=ga-dg2
      g2=g1+dg

       aa1=coul(cou,g1,inth)
       aa2=coul(cou,g2,inth)
       bb=bloss(p1,b2,b0,b3,zz4)

     nn1=nn[gg]
     nn2=nn[gg+1]
     q1=q_inj[gg]
     nn[gg]=age(idt,idg,aa1,bb,g1,q1,nn1,nn2,aa2,g2)
 end

   id=find(x-> (x <= 1e-10),nn)
   nn[id]=1e-10

 return nn
end


function eth(n::Float64,volume::Float64,t::Float64)
@fastmath ethermal=(n/mp)*volume*(1.5*kb*t)
return ethermal
end

function ecri(n::Float64,vshock::Float64,volume::Float64)
@fastmath  Ecr_inj=1e15*(fi*(n/mp)*vshock^3.0)*mp*volume^0.6666
return Ecr_inj
end

function tprei(ti::Float64,m2::Float64)
 @fastmath  tpre=ti*(16.*m2)/((5.*m2-1.)*(m2+3.))
  return tpre
end

function csprei(tpre::Float64)
@fastmath cspre=1e-5*sqrt(1.666*tpre*kb/(mp*1.1))
 return cspre
end

function spectra(p::Array{Float64,2},pe::Array{Float64,2},dt::Float64,zz::Float64,lsize::Float64,tt::Int64)
#function spectra(p,pe,dt,zz,lsize)
#print("1")
#tic()
    shock=0
    n_inj=fill(0.0,ngamma)
    q_inj=fill(1e-60,ngamma)
    n_re=fill(0.0,ngamma)

    #Ecr_inj=-60. log10.(fi1)+15.0+log10(n0*vshock1^3.0)+log10(mp)+0.6666*volume
   volume=(lsize*kpctocm)^3.
#toc()
 @inbounds @simd for i in 1:np
#print("2")
#tic() 
  ethermal=eth(p[7,i],volume,p[8,i])
#toc()
#@fastmath ethermal=log10(p[i,7]/mp)+volume+log10(1.5*kb*p[i,8])  #...observed                                                                                  
#print("3")
#tic()
 Ecr_inj=1e-20
 m=p[9,i]
 m2=m*m
  tpre=tprei(p[8,i],m2)
#p[8,i]*(16.*m2)/((5.*m2-1.)*(m2+3.))                                                                                                                  
  cspre=csprei(tpre)

#1e-5*sqrt(1.666*tpre*kb*inv(mp*1.1))    #...preshock                                                                              

  vshock=0.
  ecr=0.0
  m=0.0
  v=0.
  shock=0
#toc()
#   if tt==tini
#   p[9,i]=0.
#   m=p[9,i]
#   ecr=0.001
#   tpre=1.e7
#   cspre=1e-5*sqrt(1.666*tpre*kb*inv(mp*1.1))    
#   v=m*cspre
#   vshock=v
#   shock=1
    
#   end
#print("4")
#tic()  
   if p[9,i] >= 1.5 && p[8,i] >=1e6 && tt>= tini+1 

       m=p[9,i]
       vshock=m*cspre
       shock=1
#     if m<=3.0
#     shock=2
#     end
     Ecr_inj=ecri(p[7,i],vshock,volume)
         
    ecr=Ecr_inj/norm
    v=vshock
  end

#toc()

  nth=p[7,i]/mp
  t2=dt*gyrtosec
  t1=0.0
  
  nn=pe[:,i]

  delta_t=t2-t1

  tacc=0.0
#print("5")
#tic()
  n_inj[:]=0.
  q_inj[:]=1e-60
  n_re[:]=0.0
  b0=p[10,i]
  aa=evolve_spectrum(zz,v,t2,t1,nth,m,b0,ecr,nn,shock,delta_t,volume,gammaval,g_max,g_min,dg,ngamma,part,vshock,tacc,n_inj,q_inj,n_re)
   pe[:,i]=nn
# toc()
e
  end


return pe
end

