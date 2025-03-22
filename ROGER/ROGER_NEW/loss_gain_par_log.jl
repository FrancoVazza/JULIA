@everywhere using LinearAlgebra
#.....combination of loss and gains
@everywhere  function age(idt::Float64,idg::Float64,aa1::Float64,bb::Float64,g1::Float64,q1::Float64,nn1::Float64,nn2::Float64,aa2::Float64,g2::Float64,cc::Float64,dd::Float64)
  return @fastmath q1+((inv(idt+idg*(+cc+aa1 + dd+bb*g1^2)))*(nn1*(dd)+nn1*idt+nn2*idg*(+cc+aa2+bb*g2^2)))
end

#....advection terms compression/rarefaction terms via Div(v)
@everywhere  function adv(div::Float64)
  return @fastmath 8.11e16*div/(7720.0*3e13)
end

#....turbulent reacceleration term (~Fermi II) - only systematic acceleration
@everywhere  function asa(turb::Float64,n1::Float64,scale_turbo::Float64,bf::Float64)
  return @fastmath sqrt(n1/1e-3)*(turb*1e-7)^3/(1.25e5*1e6*yr*scale_turbo*bf/0.5)
end


@everywhere  function coul_nr(cou::Float64,g1::Float64,inth::Float64)
  betav=(1.0-1.0/g1^2)   #...this is sqrt(beta) but beta^2 in the next step 
  return @fastmath   (1/betav)*cou*(1+(log(g1*inth))/75.0)

end

#...Bremsstrahlung losses
@everywhere  function bremss(bs::Float64,g1::Float64)
  return @fastmath bs*g1*(log(2*g1)-0.33333)
end
#....JP 1973 model
@everywhere  function bloss(part1::Int64,b0::Float64,zz4::Float64)
  return  @fastmath part1*(b_syn*b0^2. +zz4*b_IC)
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


#.....FP solver based on Donnert, J. & Brunetti, G. 2014, MNRAS, 443, 3564
@everywhere  function evolve_spectrum_tridiag(zz::Float64,v::Float64,t2::Float64,t1::Float64,nth::Float64,m::Float64,b0::Float64,tpost::Float64,ecr::Float64,nn::Array{Float64,1},shock::Int64,delta_t::Float64,volume::Float64,pval::Array{Float64,1},p_max::Float64,p_min::Float64,dp::Float64,np::Int64,part::Int64,tacc::Float64,n_inj::Array{Float64,1},Qinj::Array{Float64,1},n_re::Array{Float64,1},div::Float64,curl::Float64,posts::Int64,volume_tr::Float64) #,model)

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

  local @fastmath dt=(t2-t1)

  #losses
  @fastmath  local  cou=nth*c1
  @fastmath  local  idt=inv(dt)
  @fastmath  local  bs=nth*cbs*(Hfrac+3*Hefrac)  #...constant for Bremsstrahlung losses, from Donnert+14
  @fastmath  local  icou=inv(cou)
  @fastmath  local  inth=inv(nth)
  @fastmath  local  cons1=2.3e29*(erest*1e-3)^(0.3333)*b0^(-0.3333)*(lcorr*0.05)^(0.66666)
  @fastmath  local  cons2=(f-1.)*inv(f+1.)

  #
  n_inj[:].=0.
  Qinj[:].=0.

  if shock == 1 && m >= mthr   #injection of power-law spectrum of particles

    local @fastmath   beff=sqrt(1e6*b0+3.2*(1+zz)^2.)

    local  p_cut=10^pval[np]

    local        @fastmath        pth=sqrt(2*tpost*kb*mp)   #...thermal momentum of protons
    local        @fastmath        pthe=sqrt(2*tpost*kb*me)  #...thermal momentum of electrons
    local                         eB=0.23                   #...Kang et al. 2011
    local        @fastmath        xinj=(1.17*mp*vpost/(pth))*(1+1.07/(eB))*(m/(3))^(0.1)  #...Eq.7 in Pizke+13 - xinj ~2.5.3-5
    local                p_inje=pthe*xinj*sqrt(mp/me) #....injection momentum for electrons #...Kang 2021, eq.4
    local                p_injp=pth*xinj              #....injection momentum for protons

    local                pmin_inje=p_inje/(me*vc)   #..injection momentum of electrons normalised to me*vc to be inserted in our spectrum
    local                ipmin_inje=1
    if pmin_inje > 10^p_min
    local                ipmin_inje=convert(Int64,trunc((log10(pmin_inje)-p_min)/dp)) #...number of array element of momenta to start injection
    end 


      local @fastmath   q=4*m2/(m2-1)
      local @fastmath   Kpe=(mp/me)^(0.5*(q-3))   #...proton to electron ratio From Kang+11, eq.10
      local @fastmath   Kep=inv(Kpe)              #...electron to proton ratio to be used below

      @inbounds @fastmath   for i in 1:np   #finds maximum gamma where tacc < tlosses based on DSA model
      local @fastmath    gam1=sqrt(1+(10. ^pval[i])^2)
      local @fastmath    ic_lose=b_syn*((gam1^2))*zz4
      local @fastmath    diff=(gam1-1.0)*cons1 # in cm^2/s
      local @fastmath    sh_gain=(inv(f)*gam1*(vpre)^2. *cons2*inv(3. *diff))
      if ic_lose > sh_gain
        p_cut=10. ^pval[i]
      end
      end

    local    ξinj=3.5   #...injection momentum 
    local    ξ=(4.0/(1.77))*((ξinj^3.0)/(q-3))*exp(-ξinj^2)
    local   @fastmath  Ke=Kep*ξ*nth*volume_tr

      @inbounds @simd    for j in ipmin_inje:np-2
      n_inj[j]=Ke*(10.0^pval[j])^(-delta)*(1-(10.0^pval[j])/p_cut)^(delta-2)
      if isnan(n_inj[j]) || isinf(n_inj[j]) || Ke==0
        n_inj[j]=0.
      end
      Qinj[j]=n_inj[j]    
      end
     end

  #...SHOCK REACCELERATION
  if shock == 2
    local   pmin_re=ipmin_inje    #1
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
  local      np2=np
  H=similar(nn)
  Escape=similar(nn)
  Wp=similar(nn)
  Wm=similar(nn)
  Dpp=similar(nn)
  E=similar(nn)
  F=similar(nn)

  H.=0.0
  Wp.=0.0
  Wm.=0.0
  Dpp.=0.0
  Escape.=1e30
  F.=0.0

  @inbounds @simd  for j in 2:np
    gg=j

    local       g1=sqrt(1+(10. ^pval[gg])^2)
    local       pv=10^pval[gg]
    idg=1/(g1*dp)
    local        coulomb_loss=coul(cou,g1,inth)     #...already has gamma dependence inside
    local        bs_loss=bremss(bs,g1)              #...already has gamma dependence inside
    local        syn_IC_loss=bloss(part1,b0,zz4)    #...must be multiplied for gamma^2
    local        adv_term=adv(div)                  #... must be multiplied for gamma
    local        turb=scale_turbo*kpc*curl   #...in cm/s
    if posts ==1   #velocity and shock limiters to prevent spurious acceleration from post-shock velocity
        turb=0.0
      end
      
      local        dd=asa(turb,nth,1e-3*scale_turbo,b0*1e6) #...scale in Mpc, B in muG

    if 1/dd <=dt
      dd=1/dt
    end

    local        Dpp[j]=pv^2*0.25*dd
 

    H[j]=syn_IC_loss*pv^2+pv*adv_term-4.0*Dpp[j]/pv+coulomb_loss+bs_loss
    ww=H[j]/Dpp[j]/idg
    if ww >=700
      ww=700
    end
    if  abs(ww)<=1e-8
      ww=sign(ww)*1e-8
    end
    expw=exp(ww)
    Wp[j]=1.0/ (1.0 - 1.0/expw)
    Wm[j]= 1.0/(expw-1.0)
  end


  Wp[np]=0.0
  Wm[np]=0.0
  Wp[1]=0.0
  Wm[1]=0.0
  H[np]=0.0
  Dpp[np]=0.0


  @inbounds @simd for j in 2:np

    local       @fastmath g1=10. ^pval[j]
    idg=1/(g1*dp)
    A=  -H[j]*Wp[j]*idg/idt
    C = -H[j-1]*Wm[j-1]*idg/idt
    B = 1.0+(H[j-1]*Wp[j-1]+ H[j]*Wm[j])*idg/(idt)#+1.0/(idt*Escape[j])   #...escape term in momentum space is commented out 
    O = nn[j] + Qinj[j]
        E[j] = (A)/(B - C*E[j-1])
    if abs(E[j]) <=1e-10#
      E[j]=sign(E[j])*1e-10
    end
    F[j]=(O - C*F[j-1])/(B - C*E[j-1])

  end

  E[np-1:np].=0.0
  F[np-1:np].=0.0

  @inbounds  @simd for j in 1:np-2
    i=pend-j
    nn0=nn[i]
    nn[i]=F[i]-E[i]*nn[i+1]
    if nn[i]<=1e-20 || nn[i]>=1e70
      nn[i]=1e-20
    end

  end

  return nn
end


@everywhere  function eth(n::Float64,volume::Float64,t::Float64)   #...thermal energy of tracers
  local  ethermal=n/mp*volume*1.5*kb*t
  return ethermal
end

@everywhere  function ecri(n::Float64,vshock::Float64,volume::Float64,eth::Float64)  #...injection of energy from shocks
  local   Ecr_inj=0.5*n*vshock^3.0*volume
  if      Ecr_inj>=eth
    Ecr_inj=eth
  end
  return Ecr_inj
end

@everywhere  function spectra(tra::Array{Float64,2},pev::Array{Float64,2},dt0::Float64,zz::Float64,lsize::Float64,tt::Int64,psync::Array{Float64,2},model::String,maxp::Int64)

 ntras=size(tra)
 ntra=ntras[2]

#....dynamic setting of timestep: the timestep is the minimum between tcool at the maximum momentum maxp, and the turbulent re-acceleration timescale
#....cooling time
 γref=10^pval[maxp]
 maxb=maximum(tra[10,:])
 btime=1/(γref*bloss(2,maxb,(1+zz)^4))/(3e7*1e9)  #....cooling time in Gyr

 #....minimum reacceleratioe time - we take the maximum turbulence, density and minimum B-field of these tracers
 maxt=maximum(tra[12,:])
 minb=minimum(tra[10,:])
 maxn=maximum(tra[7,:]/mp)
 maxt*=scale_turbo*dx
 dd=asa(maxt,maxn,1e-3*scale_turbo,minb*1e6)

 turb_time=1/dd/(3e7*1e9)  #....reacceleration time in Gyr
 time2=[btime,turb_time]
 dt_min=minimum(time2)

 if dt_min<=dt0/30   #...100 best
  dt_min=dt0/30
 end

 dt=dt_min        #...the subcylcing time is set
 last_tsub=convert(Int64,trunc(dt0/dt)) #...number of time subcyles
 dtlast=dt0-last_tsub*dt #...residual timestep to match dt0

 #..if m>mthr, a tracer must be shocked withih the timestep
 #...since we have last_tsub subcyling timesteps, we must select a random one to shock the tracer
 #...turbulent re-acceleration is de-activated for a time corresponding to the shock crossing time for the dx scale
 tsub_shock=last_tsub*rand(ntra)  #.....array of randomly drawn times (within dt) for each tracer to be shocked

 tt=0

@inbounds  for tsub in 1:last_tsub+1
tt+=1
  if tt==last_tsub+1
  dt=dtlast
  end

    shock=0
    local     n_inj=fill(0.0,np)
    local     q_inj=fill(0.0,np)
    local     n_re=fill(0.0,np)
    local    volume=3*log10(lsize*kpctocm/(1+z))   #..users must customise their choice here: at the moment the volume associated with each tracer, relevant to compute the shock injection of CRe, is a generic cell size (side lsize)
    

    time_total=0.

    @inbounds @simd for i in 1:ntra
      itsub_shock=convert(Int64,trunc(tsub_shock[i]))  #...integer definiting the subcycle of shock time for this tracer
      posts=0
      volume_tr=tra[15,i]*(dx*kpctocm/(1+z))^3 #volume#*(Dmin/tra[7,i])*((1+z)/(1+z_in))^3  #...the volume of each tracer is corrected by the gas compression/expansion relative to the initial density at injection
      local   ethermal=eth(tra[7,i],volume_tr,tra[8,i])
      Ecr_inj=0.
      local  m=tra[9,i]
      local            tmin=1e4
      @fastmath  tpre=tra[8,i]*(16. *m^2)/((5. *m^2-1)*(m^2+3.))
      @fastmath  tpost=tra[8,i]
      local      v2=[tpre,tmin]
      local      tpre=maximum(v2)
      local @fastmath   cspre=1e-5*sqrt(1.666*tpre*kb*inv(mp*1.1))    #...preshock
      local @fastmath   cspost=1e-5*sqrt(1.666*tpost*kb*inv(mp*1.1))    #...postshock

      local    vshock=0.
      local    ecr=0.0
      local    v=0.
      local    shock=0
      local    tcross=(dx*kpc)/(1e5*vshock)/(yr*1e9)/(1+zz) #...crossing time of element dx/(1+zz) in Gyr. if (time-tsub_shock)<tcross the particle is being shocked and turbulence is deactivated


      if m >= mthr && tsub == itsub_shock  #...shock injection only done at the last subcycling step
        posts=1
        if m<minj #...also shock re-acceleration is activated
          shock=2
        end
        if m>minj #...also shock re-acceleration is activated
          shock=1
        end
        if m>10.0
          m=10.0
        end
        vshock=m*cspre
        Ecr_inj=ecri(tra[7,i],vshock*1e5,volume_tr,ethermal)   #...electrons injected by new shocks
        ecr=Ecr_inj
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
      turb=scale_turbo*kpc*curl   #


 if model=="model1"   #...shocks injected new CRe 
 shock=1
   if turb >= cspost*1e5*0.5
  ratiot=cspost*1e5*0.5/turb
  curl*=ratiot
   end
  end


 if model=="model2"  #...the effect of shock injection is neglected 
  shock=0
  if turb >= cspost*1e5*0.5
  ratiot=cspost*1e5*0.5/turb
  curl*=ratiot
   end
 end

  if model=="model3"   #...as model 1  (but the input population was different )
  shock=1
  if turb >= cspost*1e5*0.5
    ratiot=cspost*1e5*0.5/turb
    curl*=ratiot
   end
end
 
        aa=evolve_spectrum_tridiag(zz,vshock,t2,t1,nth,m,b0,tpost,ecr,nn,shock,delta_t,volume_tr,pval,p_max,p_min,dp,np,part,tacc,n_inj,q_inj,n_re,div,curl,posts,volume_tr)
      tnow2=time()
      time_total+=(tnow2-tnow)
      pev[:,i].=nn

    end
  end

  return pev[:,:]
end

@everywhere function sync_wrapper(pval::Array{Float64},b0::Array{Float64},nn::Array{Float64,2},freqf::Array{Float64},psync::Array{Float64,2})
  local   ntras=size(b0)
  @inbounds for i in 1:ntras[1]
    @inbounds @simd    for f in 1:nfreq[1]
      sync=synchrotronK(pval,b0[i],nn[:,i],freqf[f])
      psync[i,f]=sync
    end
  end
  return psync
end
