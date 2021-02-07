  using LinearAlgebra
#.....combination of loss and gain terms

      function age(idt::Float64,idg::Float64,aa1::Float64,bb::Float64,g1::Float64,q1::Float64,nn1::Float64,nn2::Float64,aa2::Float64,g2::Float64,cc::Float64)

      return @fastmath ((inv(idt+idg*(+cc+aa1 + bb*g1^2)))*(q1+nn1*idt+nn2*idg*(+cc+aa2+bb*g2^2)))

      end

#....advection terms compression/rarefaction terms via Div(v)
      function adv(div::Float64)

      return @fastmath 8.11e16*div/(7720*3e13)

    end

#...Coulomb losses
      function coul(cou::Float64,g1::Float64,inth::Float64)

        return @fastmath cou*(1+0.0133333*(log(g1*inth)))#,@fastmath cou*(1+0.0133333*(log(g2*inth)))

      end

#....Synchtrotron + Inverse Compton Losses (Jaffe-Perola 1973 model)
      function bloss(part1::Int64,b0::Float64,zz4::Float64)

        return  @fastmath part1*b2*(0.666*b0^2. +zz4*b3)

      end

#....shock Re-acceleration
      function reacc3(delta::Float64,nn::Array{Float64,1},pend::Int64,pmin_re::Int64,pval::Array{Float64,1},n_re::Array{Float64,1})

        n_re[1]=nn[1]
        pmin_re=1

        @inbounds    for p in pmin_re:np-2
          nnp=0.
          @inbounds @simd    for p2 in pmin_re:p
            @fastmath    nnp+=(delta+2.)*pval[p]^(-delta)*((nn[p2]*dp*(pval[p2])^(delta-1.)))#; as in Eq.6 of Kang & Ryu 2011

                          end
                          n_re[p]=nnp
                        end
                        return n_re
                      end


   function evolve_spectrum(zz::Float64,v::Float64,t2::Float64,t1::Float64,nth::Float64,m::Float64,b0::Float64,tpost::Float64,ecr::Float64,nn::Array{Float64,1},shock::Int64,delta_t::Float64,volume::Float64,pval::Array{Float64,1},p_max::Float64,p_min::Float64,dp::Float64,np::Int64,part::Int64,vshock::Float64,tacc::Float64,n_inj::Array{Float64,1},q_inj::Array{Float64,1},n_re::Array{Float64,1},div::Float64,curl::Float64)

        pend=np
        #..useful quantitites
        @fastmath m2=m^2.
        @fastmath m4=inv(m^4.)
        @fastmath vpre=1.0e5*v
        @fastmath f=(4. *m2)/(m2+ 3.)
        @fastmath  vpost=vpre/(f)
        @fastmath delta=2. *(m2+ 1.)/(m2-1.)
        @fastmath zz4=(1. +zz)^4.

          ic_lose=0.0
          diff=0.0
          sh_gain=0.0

        @fastmath dt=1.0*(t2-t1)
  #losses
    @fastmath       cou=nth*c1 #..1/s   from Sarazin 1999
    @fastmath      idt=inv(dt)
    @fastmath      idp=inv(dp)
    @fastmath      icou=inv(cou)
    @fastmath      inth=inv(nth)
    @fastmath      cons1=2.3e29*(erest*1e-3)^(0.3333)*b0^(-0.3333)*(lcorr*0.05)^(0.66666)
    @fastmath      cons2=(f-1.)*inv(f+1.)
    @fastmath      dp2=dp*0.5

 #.....SHOCK ACCELERATION

 n_inj[:].=0.
 q_inj[:].=0.

  if shock >= 1 && m >= mthr   #injection of power-law spectrum of particles
     @fastmath   beff=sqrt(1e6*b0+3.2*(1+zz)^2.)

          p_cut=pval[np]

                @fastmath        pth=sqrt(2*tpost*kb*mp)   #...thermal momentum of protons
                @fastmath        pthe=sqrt(2*tpost*kb*me)  #...thermal momentum of electrons
                                 eB=0.23                   #...Kang et al. 2011
                @fastmath        xinj=(1.17*mp*vpost/(pth))*(1+1.07/(eB))*(m/(3))^(0.1)  #...Eq.7 in Pizke+13 -

                        p_inj=pthe*xinj              #....injection momentum for electrons


 #....computing the normalisation of injected Cosmic Ray electrons by shocks (thermal leakage model from DSA)
    if m >=2
       @fastmath   eta=5.46*m4-9.78*(m-1)*m4+4.17*(m-1)^2*m4-0.334*(m-1)^3*m4+0.57*(m-1)^4*m4 #Kang & Jones 2007 efficiency
     end
     if m <2
       @fastmath  eta=1.96e-3*(m2-1)
     end
      @fastmath   q=4*m2/(m2-1)
      @fastmath   Kpe=(mp/me)^(0.5*(q-3))   #...proton to electron ratio From Kang+11, eq.10
      @fastmath   Kep=inv(Kpe)              #...electron to proton ratio to be used below

   @inbounds @fastmath @simd   for i in 1:np   #finds maximum gamma where tacc < tlosses
        @fastmath    gam1=sqrt(1+pval[i]^2)
        @fastmath    ic_lose=b1*((gam1^2.))*zz4
        @fastmath             diff=(gam1-1.0)*cons1 # in cm^2/s
        @fastmath    sh_gain=(inv(f)*gam1*(vpre)^2. *cons2*inv(3. *diff))

          if ic_lose > sh_gain
          p_cut=pval[i]
           end
               end

               deltaB=delta
     @fastmath  ibeta=beta_inc(0.5*(deltaB-2),0.5*(3-deltaB),1/(1+p_cut^2))
     @fastmath  inteB=(erest/(delta-1))*0.5*ibeta[1]+(p_cut^(1-delta))*(sqrt(1+p_cut^2)-1)
     @fastmath  Ke= (Kep*eta*ecr)/(vpost*inteB)
     @inbounds @simd    for j in 1:np

       @fastmath    n_inj[j]=Ke*pval[j]^(-delta)*(1-pval[j]/p_cut)^(delta-2) 

         q_inj[j]=n_inj[j]

               end

           end



#...SHOCK REACCELERATION
    if shock == 2

        pmin_re=1

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
        np2=np-1
        @inbounds @simd for j in 1:np-1

      end
      @inbounds @simd  for j in 1:np2
      gg=pend-j

         @fastmath g1=sqrt(1+(pval[gg])^2)
         @fastmath g2=sqrt(1+(pval[gg]+dp2)^2)

          a1=coul(cou,g1,inth)
          a2=coul(cou,g2,inth)

          bb=bloss(part1,b0,zz4)
          cc=adv(div)

        nn1=nn[gg]
        nn2=nn[gg+1]
        q1=q_inj[gg]
        idg=1/(2*(g2-g1))

        nn[gg]=age(idt,idg,a1,bb,g1,q1,nn1,nn2,a2,g2,cc)

end


 return nn
end

#...thermal energy
   function eth(n::Float64,volume::Float64,t::Float64)
    ethermal=n/mp*10^volume*1.5*kb*t
    return ethermal
  end

#...cosmic ray energy flux
    function ecri(n::Float64,vshock::Float64,volume::Float64)
     Ecr_inj=0.5*n*vshock^3*10^volume
     return Ecr_inj
   end

   function spectra(tra::Array{Float64,2},pe::Array{Float64,2},dt::Float64,zz::Float64,lsize::Float64,tt::Int64,psync::Array{Float64,2})

         shock=0
         n_inj=fill(0.0,np)
         q_inj=fill(0.0,np)
         n_re=fill(0.0,np)
         volume=3*log10(lsize*kpctocm)   #volume associated to each tracers

  @inbounds @simd for i in 1:ntra


    Ecr_inj=0.
    m=tra[9,i]

    @fastmath  tpre=tra[8,i]*(16. *m^2)/((5. *m^2-1)*(m^2+3.))
    @fastmath  tpost=tra[8,i]
        v2=[tpre,tmin]
        tpre=minimum(v2)
        @fastmath   cspre=1e-5*sqrt(1.666*tpre*kb*inv(mp*1.1))    #...preshock

      vshock=0.
      ecr=0.0
      v=0.
      shock=0

      if m >= 1.6

       vshock=m*cspre
       shock=1

       if m<=3.5
       shock=2
       end

     Ecr_inj=ecri(tra[7,i],vshock*1e5,volume)
     ecr=Ecr_inj
      v=vshock

   end

      nth=tra[7,i]/mp
      t2=dt*gyrtosec
      t1=0.0

  nn=pe[:,i]

      delta_t=t2-t1

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
  pe[:,i].=nn


 @inbounds @simd    for f in 1:nfreq[1]
   sync=synchrotron(pval,1e6*b0,pe[:,i],freqf[f])
   psync[:,f].=sync
 end
end
#  end


return pe,psync
end
