@everywhere using LinearAlgebra


#.....combination of loss and gains
@everywhere  function age(idt::Float64,idg::Float64,aa1::Float64,bb::Float64,g1::Float64,q1::Float64,nn1::Float64,nn2::Float64,aa2::Float64,g2::Float64,cc::Float64,dd::Float64)

return @fastmath ((inv(idt+idg*(+cc+aa1 + bb*g1^2)))*(q1+ nn1*(dd)+nn1*idt+nn2*idg*(+cc+aa2+bb*g2^2)))

end

#....advection terms compression/rarefaction terms via Div(v)
@everywhere  function adv(div::Float64)
 return @fastmath 8.11e16*div/(7720.0*3e13)
end

#....turbulent reacceleration term (~Fermi II) - only systematic acceleration
@everywhere  function asa(turb::Float64,n1::Float64,scale::Float64,bf::Float64)
 return @fastmath sqrt(n1*1e3)*(turb*1e-7)^3/(3e13*1.25e5*scale*bf/0.5)
end

#...Coulomb losses
  @everywhere  function coul(cou::Float64,g1::Float64,inth::Float64)

return @fastmath cou*(1+0.0133333*(log(g1*inth)))
end

#....JP 1973 model
@everywhere  function bloss(part1::Int64,b0::Float64,zz4::Float64)
return  @fastmath part1*b2*(b0^2. +zz4*b3)
end


#.....shock re-acceleration
@everywhere  function reacc3(delta::Float64,nn::Array{Float64,1},pend::Int64,pmin_re::Int64,pval::Array{Float64,1},n_re::Array{Float64,1})
pval_re=similar(pval)
for i in eachindex(pval)
pval_re[i]=10^pval[i]
end
n_re[1]=nn[1]
pmin_re=1

@inbounds    for p in pmin_re:np-2
nnp=0.
@inbounds @simd    for p2 in pmin_re:p
  dp2=(pval_re[p2+1]-pval_re[p2])

@fastmath    nnp+=(delta+2.)*(pval_re[p])^(-delta)*((nn[p2]*dp2*(pval_re[p2])^(delta-1.)))#; as in Eq.6 of Kang & Ryu 2011
end
n_re[p]=nnp
end
 return n_re
end


@everywhere  function evolve_spectrum(zz::Float64,v::Float64,t2::Float64,t1::Float64,nth::Float64,m::Float64,b0::Float64,tpost::Float64,ecr::Float64,nn::Array{Float64,1},shock::Int64,delta_t::Float64,volume::Float64,pval::Array{Float64,1},p_max::Float64,p_min::Float64,dp::Float64,np::Int64,part::Int64,vshock::Float64,tacc::Float64,n_inj::Array{Float64,1},q_inj::Array{Float64,1},n_re::Array{Float64,1},div::Float64,curl::Float64)

local pend=np
local @fastmath m2=m^2.
local @fastmath m4=inv(m^4.)
local @fastmath vpre=1.0e5*v
local @fastmath f=(4. *m2)/(m2+ 3.)
local @fastmath  vpost=vpre/(f)
local @fastmath delta=2. *(m2+ 1.)/(m2-1.)
local @fastmath zz4=(1. +zz)^4.

local   ic_lose=0.0
local   diff=0.0
local   sh_gain=0.0

local @fastmath dt=1.0*(t2-t1)
  #losses
@fastmath  local   cou=nth*c1 #..1/s   from Sarazin 1999
@fastmath  local  idt=inv(dt)

@fastmath  local  icou=inv(cou)
@fastmath  local  inth=inv(nth)
@fastmath  local  cons1=2.3e29*(erest*1e-3)^(0.3333)*b0^(-0.3333)*(lcorr*0.05)^(0.66666)
@fastmath  local  cons2=(f-1.)*inv(f+1.)

 #
 n_inj[:].=0.
 q_inj[:].=0.

  if shock >= 1 && m >= mthr   #injection of power-law spectrum of particles
  local @fastmath   beff=sqrt(1e6*b0+3.2*(1+zz)^2.)

      local  p_cut=10^pval[np]

      local        @fastmath        pth=sqrt(2*tpost*kb*mp)   #...thermal momentum of protons
      local        @fastmath        pthe=sqrt(2*tpost*kb*me)  #...thermal momentum of electrons
      local                         eB=0.23                   #...Kang et al. 2011
      local        @fastmath        xinj=(1.17*mp*vpost/(pth))*(1+1.07/(eB))*(m/(3))^(0.1)  #...Eq.7 in Pizke+13 -

      local                p_inj=pthe*xinj              #....injection momentum for electrons
      local                p_injp=pth*xinj              #....injection momentum for electrons

  pmin_inj=pval[1]


 #....computing the normalisation of injected Cosmic Ray electrons by shocks (thermal leakage model from DSA)
 if m >=2
 local @fastmath   eta=5.46*m4-9.78*(m-1)*m4+4.17*(m-1)^2*m4-0.334*(m-1)^3*m4+0.57*(m-1)^4*m4 #Kang & Jones 2007 efficiency
end
 if m <2
 local @fastmath  eta=1.96e-3*(m2-1)
 end
# eta=1e-3
 local @fastmath   q=4*m2/(m2-1)
 local @fastmath   Kpe=(mp/me)^(0.5*(q-3))   #...proton to electron ratio From Kang+11, eq.10
 local @fastmath   Kep=inv(Kpe)              #...electron to proton ratio to be used below

 @inbounds @fastmath   for i in 1:np   #finds maximum gamma where tacc < tlosses
 local @fastmath    gam1=sqrt(1+(10. ^pval[i])^2)
 local @fastmath    ic_lose=b1*((gam1^2))*zz4
 local @fastmath             diff=(gam1-1.0)*cons1 # in cm^2/s
 local @fastmath    sh_gain=(inv(f)*gam1*(vpre)^2. *cons2*inv(3. *diff))

          if ic_lose > sh_gain
          p_cut=10. ^pval[i]
          end
          end

 deltaB=delta
 if deltaB >=2.99
 deltaB=2.99
 end
 if deltaB<=2.1
 deltaB=2.1
 end
 local   @fastmath  ibeta=beta_inc(0.5*(deltaB-2),0.5*(3-deltaB),1/(1+p_cut^2))
 local   @fastmath  inteB=(erest/(delta-1))*(0.5*ibeta[1]+(p_cut^(1-delta))*(sqrt(1+p_cut^2)-1))
 local   @fastmath  Ke= (Kep*eta*ecr)/(vpost*inteB)
 ipmin_inj=1

  @inbounds @simd    for j in ipmin_inj:np-1
    @fastmath    n_inj[j]=Ke*(10. ^pval[j])^(-delta)*(1-(10. ^pval[j])/p_cut)^(delta-2)
    if isnan(n_inj[j]) || isinf(n_inj[j])
       n_inj[j]=0.
         end
        q_inj[j]=n_inj[j]
          end
           end


#...SHOCK REACCELERATION
  if shock == 2
   local   pmin_re=1
 n_re=reacc3(delta,nn,pend,pmin_re,pval,n_re)
@views  nn[:]=n_re[:]
 end


#setting to zero losses that are not valid either for protons (p=1) or electrons (p=2)
    part1=1
    part2=1
    if part==1
    part1=0
    end
    if part==2
    part2=0
    end

#....INTEGRATION OF LOSSES AND INCLUSION OF FRESHLY INJECTED PARTICLES
local      np2=np-1
@inbounds @simd  for j in 1:np2
      gg=pend-j

local       @fastmath g1=sqrt(1+(10. ^pval[gg])^2)
local       @fastmath g2=sqrt(1+(10. ^pval[gg+1])^2)

local        a1=coul(cou,g1,inth)
local        a2=coul(cou,g2,inth)

local        bb=bloss(part1,b0,zz4)
local        cc=adv(div)

local        turb=scale*1.08e21*curl   #...in cm/s

      dd=asa(turb,nth,1e-3*scale,b0*1e6)
       if turb<=1e1 || turb >=1e8  || shock >=1   #velocity and shock limiters to prevent spurious acceleration from post-shock velocity
       dd=0.
       end

        nn1=nn[gg]
        nn2=nn[gg+1]
        q1=q_inj[gg]
        idg=1/((g2-g1))
          nn0=nn[gg]
          nn[gg]=age(idt,idg,a1,bb,g1,q1,nn1,nn2,a2,g2,cc,dd)

end


 return nn
end

@everywhere  function eth(n::Float64,volume::Float64,t::Float64)   #...thermal energy of tracers
local  ethermal=n/mp*10^volume*1.5*kb*t
return ethermal
end

@everywhere  function ecri(n::Float64,vshock::Float64,volume::Float64)  #...injection of energy from shocks
local   Ecr_inj=0.5*n*vshock^3.0*10^volume
  return Ecr_inj
end

@everywhere  function spectra(tra::Array{Float64,2},pev::Array{Float64,2},dt0::Float64,zz::Float64,lsize::Float64,tt::Int64)
  last_tsub=4 #...subcycles
  for tsub in 1:last_tsub
  dt=dt0*0.25

  local   ntras=size(tra)
  local   ntra=ntras[2]

    shock=0
local     n_inj=fill(0.0,np)
local     q_inj=fill(0.0,np)
local     n_re=fill(0.0,np)
local     n_trac_v = 100
local    volume=3*log10(lsize*kpctocm)   #..users must customise their choice here: at the moment the volume associated with each tracer, relevant to compute the shock injection of CRe, is a generic cell size (side lsize)


  time_total=0.
 @inbounds @simd for i in 1:ntra

local   ethermal=eth(tra[7,i],volume,tra[8,i])

 Ecr_inj=0.
local  m=tra[9,i]

local            tmin=1e4
@fastmath  tpre=tra[8,i]*(16. *m^2)/((5. *m^2-1)*(m^2+3.))
@fastmath  tpost=tra[8,i]
local      v2=[tpre,tmin]
local      tpre=maximum(v2)
local @fastmath   cspre=1e-5*sqrt(1.666*tpre*kb*inv(mp*1.1))    #...preshock

local    vshock=0.
local    ecr=0.0
local    v=0.
local    shock=0

   if m >= mthr && tsub == last_tsub   #...shock injection only done at the last subcycling step

       vshock=m*cspre
       shock=1

       if m<=5.0
       shock=2
       end


     Ecr_inj=ecri(tra[7,i],vshock*1e5,volume)   #...electrons injected by new shocks

     ecr=Ecr_inj
      v=vshock

   end

local    nth=tra[7,i]/mp
local    t2=dt*gyrtosec
local    t1=0.0

  nn=pev[:,i]

local    delta_t=t2-t1

  tacc=0.0

  n_inj[:].=0.
  q_inj[:].=0.
  n_re[:].=0.
  b0=tra[10,i]
  div=tra[11,i]
  curl=tra[12,i]

  tnow=time()
  aa=evolve_spectrum(zz,v,t2,t1,nth,m,b0,tpost,ecr,nn,shock,delta_t,volume,pval,p_max,p_min,dp,np,part,vshock,tacc,n_inj,q_inj,n_re,div,curl)
  tnow2=time()
  time_total+=(tnow2-tnow)
  pev[:,i].=nn

end
end

return pev[:,:]
end
