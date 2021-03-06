using Statistics
using Plots
using Random
using LaTeXStrings

gr()
Plots.GRBackend()

  root="/Users/francovazza/Desktop/"

  #...PARAMETERS
       pdead=0.01  #....death probability for an infected
       pinf0=0.40    #....probability to infect another one (in a day)
       day_max=38
       t_incub=5.5
       t_death=20
       R0=pinf0*t_incub   #...never explicitly used be, it is estimated to be 2.5 for Covid19.
       nm=4      #...number of random realisation
       day_intervention=17  #....day at which the government takes action and the infection rate is reduced (used only in scenario ss=1 )


  #...daily official data from Italian Gov.
  infreal=[97689,92472,86498,80539,74386,69176,63927,59138,53578,47021,41035,35713,31506,27980,24747,21157,17660,15113,12462,10149,9172,7375,5883,4636,3858,3089,2502,2036,1694,1128,888,650,400,322,229,157,79,16]
  deadreal=[10779,10023,9134,8165,7503,6820,6077,5476,4825,4032,3405,2978,2503,2158,1809,1441,1266,1016,827,631,463,366,233,197,148,107,79,52,34,29,21,17,12,10,7,3,2,1]
  dreal=   [38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14, 13, 12,11,10,9,8,7,6,5,4,3,2,1 ]  #21 febbraio = day1

 println(size(infreal))
 println(size(dreal))
  println(size(deadreal))

@inbounds     for ss in 0:1  #...loop over 2 possible scenarios 0=no intervention, 1=intervention which reduces pinf0 starting from a given day

  model = Array{Float64}(undef, nm,day_max,3)   #....array with global statistics for each model

   @inbounds   for mm in 1:nm #....loop over random realisation

     println("doing trial number=",mm)
     pinf=pinf0
     n_dead=0    #...initial counter of victims
     n_inf0=400  #...initial counter of infected (day=1). This is uknown - but it sets the overall normalisation and can be calibrated after the first few days
     n_diag=0    #...number counter of  diagnosed

    age=Array{Int64}(undef,n_inf0)    #...age (starting from t=1) of infection for each infected
    alive=Array{Int64}(undef,n_inf0)  #...1=alive, 0=dead
    status=Array{Int64}(undef,n_inf0) #...1=diagnosed (hence isolated/treated/dead), 0=undiagnosed

    alive[:].=1
    age[:].=0
    status[:].=0

    @inbounds for dd in 1:day_max    #...loop over days
    pdead=0.05
    if dd >=12
    pdead=0.06   #....ad-hoc increased lethality, as seems to be required by data since
    end

    if dd >=17#
    pdead=0.08   #....ad-hoc increased lethality, as seems to be required by data since
    end



    if ss==1 &&  dd >= day_intervention   #....in scenario ss=1, we can model a reduced contagion rate here
    pinf=pinf0*0.44        #...test change of infectivity after day - just a wild guess
    end
    nninfo=0



    rng = MersenneTwister(mm)     #...we first generate random set of numbers for the Monte Carlo
    a=rand!(rng,zeros(n_inf0))
    tt=randn!(rng,zeros(n_inf0))


    n_new=0
    t1=time()
    @inbounds  @simd  for i in 1:n_inf0  #loop over people already infected on this day
    if alive[i]==1 && age[i]>=0         #...how many people are infected and alive

    @fastmath  dt=t_incub+0.5*t_incub*tt[i]   #...individual incubation time: it's  a normal distribution peaked at t_incub, with 1sigma=t_incub/2
    @fastmath  dt2=t_death+0.5*t_death*tt[i]   #...individual incubation time: it's  a normal distribution peaked at t_incub, with 1sigma=t_incub/2

    if a[i] <=pdead && age[i]>=dt2#...patient is removed: the time is negative, the status is dead (0) but there likely is a diagnosis (1)
    age[i]=-1   #
    alive[i]=0  #dead
    status[i]=1 #diagnosed
    else


     age[i]+=1   #...patient is infected for one day
     alive[i]=1  #alive
     status[i]=0 #undiagnosed
     time=dt
     if age[i]>=time
     status[i]=1  #after an incubation time, the patient becomes diagnosed/isolated and cannot infect anymore
     end
     end
     end

    #...now we compute the number of new infections
     if alive[i]==1 && status[i]==0       #...identifies alive & undiagnosed patients
     if a[i] <=pinf                       #a new infection
     n_new+=1
     append!(alive,[1])                   #...all arrays are increased by the latest infection
     append!(age,[0])
     append!(status,[0])


     end
     end
 end #end of loop over infected
     t2=time()
     n_inf0+=n_new  #..updated counter of infected

   #....update daily statistics

   id=findall(x-> x==0,alive)
   n=size(id)
   ndead=n[1]  #....total number of casualties

    id=findall(x-> x==1,alive)
    n=size(id)
    ninfected=n[1]  #....total number of infected

    id=findall(x-> x==1,status)
    n=size(id)
    ndiagnosed=n[1]  #....total number of diagnosed

    println("day=",dd," ",ndead," ",n_inf0," ",ndiagnosed," CPU t=",t2-t1)
    model[mm,dd,1]=ninfected
    model[mm,dd,2]=ndead
    model[mm,dd,3]=ndiagnosed

    end #....days

  end  #...random realizations
#

#.....PLOTTING STUFF
      plo=Array{Float64}(undef,day_max,3)      #....mean
      splo1=Array{Float64}(undef,day_max,3)    #...mean-3sigma
      splo2=Array{Float64}(undef,day_max,3)    #...mean+3sigma
      day=Array{Float64}(undef,day_max)

@inbounds   @simd for dd=1:day_max
    plo[dd,1]=mean(model[:,dd,1])
    plo[dd,2]=mean(model[:,dd,2])
    plo[dd,3]=mean(model[:,dd,3])
    splo1[dd,1]=plo[dd,1]-3*std(model[:,dd,1])
    splo1[dd,2]=plo[dd,2]-3*std(model[:,dd,2])
    splo1[dd,3]=plo[dd,3]-3*std(model[:,dd,3])

    splo2[dd,1]=plo[dd,1]+3*std(model[:,dd,1])
    splo2[dd,2]=plo[dd,2]+3*std(model[:,dd,2])
    splo2[dd,3]=plo[dd,3]+3*std(model[:,dd,3])
    day[dd]=1.0*dd
    end


if ss==0
        plot(day,plo,label = ["Infected (model)" "Dead(model)" "Diagnosed(model)"],color=["blue" "red" "green"],lw=1,legend=:bottomright)
end
if ss==1
       plot!(day,plo,label = ["Infected (model+lockdown)"  "Dead(model+lockdown)" "Diagnosed(model+lockdown)"],color=["blue" "red" "green"],lw=3,ls=[:dash :dash :dash])
end


    plot!(day,splo1[:,1],color="blue",label=nothing)#"model-3stdev")
    plot!(day,splo2[:,1],color="blue",label=nothing)#"model+3stdev")
    plot!(day,splo1[:,2],color="red",label=nothing)#
    plot!(day,splo2[:,2],color="red",label=nothing)#
    plot!(day,splo1[:,3],color="green",label=nothing)#
    plot!(day,splo2[:,3],color="green",label=nothing)#
  end

#    plot(day,)

#    axis([1,day_max,1,1e5])
    plot!(xlabel="days",ylabel="number of people")
    plot!(xlims=(1,day_max*1.4),ylims=(1,2e5))
    yaxis!(:log10)

    plot!(dreal,infreal,seriestype = :scatter,color="green",label="Infected (REAL)")
    plot!(dreal,deadreal,seriestype = :scatter,color="red",label="Dead(REAL)")
  #  println(plo[day_max,2]," ",splo1[day_max,2]," ",splo2[day_max,2])

 #....scenario
     savefig("/Users/francovazza/Desktop/Julia_prog/fig_23marchB.png")
