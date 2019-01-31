
##############################################################################################
#...HISTOGRAM CALCULATION, IF NEEDED

@everywhere function histog(mm::Vector{Float64},hm::Vector{Float64},macx::Array{Float64}) 

@inbounds @simd for i in eachindex(mm)
   tag=find( x->(x >= mm[i]),macx)
    nshock=size(tag)
   hm[i]=nshock[1]
  end
return hm
end




#######################################################################################
 #....3D DIVERGENCE OF VELOCITY FIELD
 @everywhere function div_compute(vx::Array{Float64,3},vy::Array{Float64,3},vz::Array{Float64,3},n3::Tuple{Int64,Int64,Int64}) 

    div=fill(0.0,n3[1],n3[2],n3[3])
    for l in 2:n3[3]-1
    for j in 2:n3[2]-1 
 @simd    for i in 2:n3[1]-1
 @inbounds  div[i,j,l]=0.5*(vx[i+1,j,l]-vx[i-1,j,l]+vy[i,j+1,l]-vy[i,j-1,l]+vz[i,j,l+1]-vz[i,j,l-1])
    end
    end
    end

     #...takes care of boundaries (also relevant for multi-core runs)
  
    l=1
    for j in 2:n3[2]-1 
 @simd    for i in 2:n3[1]-1
 @inbounds  div[i,j,l]=0.5*(vx[i+1,j,l]-vx[i-1,j,l]+vy[i,j+1,l]-vy[i,j-1,l])+(vz[i,j,l+1]-vz[i,j,l])
    end
    end
    
 l=n3[3]
    for j in 2:n3[2]-1 
 @simd    for i in 2:n3[1]-1
 @inbounds  div[i,j,l]=0.5*(vx[i+1,j,l]-vx[i-1,j,l]+vy[i,j+1,l]-vy[i,j-1,l])+(vz[i,j,l]-vz[i,j,l-1])
    end
    end

     
    j=1
    for l in 2:n3[3]-1 
 @simd    for i in 2:n3[1]-1
 @inbounds  div[i,j,l]=0.5*(vx[i+1,j,l]-vx[i-1,j,l]+vz[i,j,l+1]-vz[i,j,l-1])+(vy[i,j+1,l]-vy[i,j,l])
    end
    end
    
   j=n3[2]
    for l in 2:n3[3]-1 
 @simd    for i in 2:n3[1]-1
 @inbounds  div[i,j,l]=0.5*(vx[i+1,j,l]-vx[i-1,j,l]+vz[i,j,l+1]-vz[i,j,l-1])+(vy[i,j,l]-vy[i,j-1,l])
    end
    end

     
      
    i=1
      for l in 2:n3[3]-1 
 @simd    for j in 2:n3[2]-1
 @inbounds  div[i,j,l]=0.5*(vy[i,j+1,l]-vy[i,j-1,l])+(vz[i,j,l+1]-vz[i,j,l-1])+(vx[i+1,j,l]-vx[i,j,l])
    end
    end

     i=n3[1]
      for l in 2:n3[3]-1 
 @simd    for j in 2:n3[2]-1
 @inbounds  div[i,j,l]=0.5*(vy[i,j+1,l]-vy[i,j-1,l])+(vz[i,j,l+1]-vz[i,j,l-1])+(vx[i,j,l]-vx[i-1,j,l])
    end
    end



     
     
  return div        

 end

########################################################################
#.....COMPUTES NORM OF MACH VECTOR
@everywhere function machm(vvx::Float64)

 @fastmath    mach=(4.*abs(vvx)+sqrt(16.*vvx*vvx+36.))*0.166666  
 return mach
    end

#######################################################################
  #...VELOCITY JUMP MACH NUMBER
  #....see Vazza+2009 MNRAS for details

@everywhere    function mach_compute(ijk::Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}},tag::Array{Int64,1},nshock::Tuple{Int64},n3::Tuple{Int64,Int64,Int64},vx::Array{Float64,3},vy::Array{Float64,3},vz::Array{Float64,3},macx::Array{Float64,3},macy::Array{Float64,3},macz::Array{Float64,3},t::Array{Float64,3},fac::Float64)   #...VELOCITY JUMP MACH NUMBER
    #....see Vazza+2009 MNRAS for details
    
    mach=similar(macx)

    #----first, Mach numbers are computed for all tagged cells (based on 3D divergence)
    @inbounds @simd   for i in 1:nshock[1]
        ijk=ind2sub((n3[1],n3[2],n3[3]),tag[i])

        local     a=ijk[1]
        local     b=ijk[2]
        local     c=ijk[3]

        #...the following steps to take care of boundaries (i.e. multi-core runs)    
 #...Option A: nice and compact and without IFs...however it runs SLOWER than uglier option B. 
#    cc2=[c+1,n3[3]]
#    bb2=[b+1,n3[2]]
#    aa2=[a+1,n3[1]]
#    cc1=[1,c-1]
#    bb1=[1,b-1]
#    aa1=[1,a-1]
#    c2=minimum(cc2)
#    b2=minimum(bb2)
#    a2=minimum(aa2)
#    c1=maximum(cc1)
#    b1=maximum(bb1)    
#    a1=maximum(aa1)

 #..Option B: ugly and with IFs, but FASTER 
    
        local  a1=a-1
        local  a2=a+1
        local  b1=b-1
        local  b2=b+1
        local  c1=c-1
        local  c2=c+1
        
        if a1 < 1 
            a1=1
        end  
 
        if a2 > n3[1]
            a2=n3[1]
        end
        if b1 < 1
            b1=1
        end
        if b2 > n3[2] 
            b2=n3[2]
        end 
        if c1 < 1
            c1=1
        end
        if c2 > n3[3]
            c2=n3[3]
        end  
        #......end of Option B

        
        @fastmath dvx=-1.*(vx[a1,b,c]-vx[a2,b,c]) 
        @fastmath dvy=-1.*(vy[a,b1,c]-vy[a,b2,c])
        @fastmath dvz=-1.*(vz[a,b,c1]-vz[a,b,c2])
        
        local t2=t[a2,b,c]
        local t1=t[a1,b,c]
        if dvx < 0. && t2 > t1 
            local  vvx=dvx/(sqrt(fac*t1))
            local  mx=machm(vvx)
            
            local v1=[mx,macx[a2,b,c]]
            macx[a2,b,c]=maximum(v1)
        elseif dvx <0. && t2 < t1 
               vvx=abs(dvx)/(sqrt(fac*t2))
             mx=machm(vvx)
            
            v1=[mx,macx[a1,b,c]]
            macx[a1,b,c]=maximum(v1)
        end
        t2=t[a,b2,c]
        t1=t[a,b1,c]
        if dvy < 0. && t2 > t1 
            local  vvy=dvy/(sqrt(fac*t1))
            local my=machm(vvy)

            v1=[my,macy[a,b2,c]]
            
            macy[a,b2,c]=maximum(v1)
        elseif dvy <0 && t2 < t1 
                 vvy=abs(dvy)/(sqrt(fac*t2))
               my=machm(vvy) 
            
            v1=[my,macy[a,b1,c]]
            macy[a,b1,c]=maximum(v1)
        end
        t2=t[a,b,c2]
        t1=t[a,b,c1]
        if dvz < 0. && t2 > t1 
            local   vvz=dvz/(sqrt(fac*t1))
            local mz=machm(vvz)
            
            v1=[mz,macz[a,b,c2]]
            
            macz[a,b,c2]=maximum(v1)
            
        elseif dvz <0. && t2 < t1 
              vvz=abs(dvz)/(sqrt(fac*t2))
              mz=machm(vvz)

            v1=[mz,macz[a,b,c1]]
     macz[a,b,c1]=maximum(v1)
        end
        
    end


    @inbounds  @simd  for i in eachindex(macx)
        @fastmath  mach[i]=sqrt(macx[i]^2. +macy[i]^2.+macz[i]^2.) 
    end
    
    mach= mach_clean(mach,macx,macy,macz,n3)  #...cleaning of shocks in case of multiple shocks in the same patch

    return mach
end





##################################################################################################
#### #...cleaning of shocks in case of multiple shocks in the same patch

@everywhere function mach_clean(mach::Array{Float64,3},macx::Array{Float64,3},macy::Array{Float64,3},macz::Array{Float64,3},n3::Tuple{Int64,Int64,Int64})   #SHOCK CLEANER 
    #.........shocks can be identifed across multiple cells because of numerical smearing
    #.........in this routine we retain only the maximum along the shock normal
    println(mthr)
    tag0=find( x->(x >= mthr), mach)

    nshock=size(tag0)
    println("number of cells with M>Mthr=", nshock)
    ijk0=ind2sub((n3[1],n3[2],n3[3]),tag0)


    local v1=Vector{Float64}(3) 

      
@inbounds    @simd for i in eachindex(tag0)

        local a=ijk0[1][i]
        local b=ijk0[2][i]
        local c=ijk0[3][i]

        x2=a
        y2=b
        z2=c
    
#...the following steps to take care of boundaries (i.e. multi-core runs)    

    #...Option A: nice and compact and without IFs...however it runs SLOWER than uglier option B. 
    #local    cc2=[c+1,n3[3]]
#local    bb2=[b+1,n3[2]]
#local    aa2=[a+1,n3[1]]
#local    cc1=[1,c-1]
#local    bb1=[1,b-1]
#local    aa1=[1,a-1]
#local    z3=minimum(cc2)
#local    y3=minimum(bb2)
#local    x3=minimum(aa2)
#local    z1=maximum(cc1)
#local    y1=maximum(bb1)    
#local    x1=maximum(aa1)

     #..Option B: ugly and with IFs, but FASTER 
    
        x1=a-1
        x3=a+1
        y1=b-1
        y3=b+1
        z1=c-1
        z3=c+1
 
        if x1 < 1 
            x1+=1
        end  
 
        if x3 >= n3[1]+1
            x3-=1
        end
        if y1 < 1
            y1+=1
        end
        if y3 >= n3[2]+1  
            y3-=1
        end 
        if z1 < 1
            z1+=1
        end
        if z3 >= n3[3]+1
            z3-=1
        end  
    #......end of Option B
    
        v1[1]=macx[x1,y2,z2]
        v1[2]=macx[x2,y2,z2]
        v1[3]=macx[x3,y2,z2]
        max3=maximum(v1)

        if macx[x1,y2,z2] != max3
            macx[x1,y2,z2]=0.
        end
        if macx[x2,y2,z2] != max3 
            macx[x2,y2,z2]=0.
        end
        if macx[x3,y2,z2] != max3 
            macx[x3,y2,z2]=0.
        end
        
        v1[1]=macy[x2,y1,z2]
        v1[2]=macy[x2,y2,z2]
        v1[3]=macy[x2,y3,z2]
        
        max3=maximum(v1)
        if macy[x2,y1,z2] != max3 
            macy[x2,y1,z2]=0.
        end   
        if macy[x2,y2,z2] != max3
            macy[x2,y2,z2]=0.
        end
        if macy[x2,y3,z2] != max3
            macy[x2,y3,z2]=0.
        end
        
        v1[1]=macz[x2,y2,z1]
        v1[2]=macz[x2,y2,z2]
        v1[3]=macz[x2,y2,z3]
        
        max3=maximum(v1)

        if macz[x2,y2,z1] != max3 
            macz[x2,y2,z1]=0.
        end   

        if macz[x2,y2,z2] != max3
            macz[x2,y2,z2]=0.
        end

        if macz[x2,y2,z3] != max3 
            macz[x2,y2,z3]=0.
        end
 @fastmath  mach[a,b,c]=sqrt(macx[a,b,c]^2. +macy[a,b,c]^2.+macz[a,b,c]^2.)
      
   end
   ijk0=nothing   
   tag0=nothing
  return mach

                    end
                    
###########################################################
##### main caller of shock functions
@everywhere function shocks(d::Array{Float64,3},t::Array{Float64,3},vx::Array{Float64,3},vy::Array{Float64,3},vz::Array{Float64,3},macx::Array{Float64,3},macy::Array{Float64,3},macz::Array{Float64,3},mach::Array{Float64,3}) 

 tic()
 n3=size(d)

#    local  macx=fill(0.0,n3[1],n3[2],n3[3])

#    local  macy=similar(macx)

#    local  macz=similar(macx)

#   local mach=similar(macx) #   mach=fill(0.0,n3[1],n3[2],n3[3])   
     toc()
    println("allocated")
    tic()
    div=div_compute(vx,vy,vz,n3)   #computes 3D divergence 
    toc()
    println("div")
    tic()
    tag=find( x->(x <= -1.e6), div)  #selection of negative divergence cells
    nshock=size(tag)
    toc()
    println("tag")
    tic()
    ijk=ind2sub((n3[1],n3[2],n3[3]),tag)
    println(nshock, "candidate shocks") 
    toc()
    println("index")
    tic()
    mach=mach_compute(ijk,tag,nshock,n3,vx,vy,vz,macx,macy,macz,t,fac)   #...computes Mach based on velocity jumps

    toc()
   println("mach computed")
    println(minimum(mach)," ",maximum(mach))

    clear!(:macy)
    clear!(:macz)
    clear!(:div)
    clear!(:macx)
    clear!(:ijk)
    clear!(:tag)   
    return mach
end

