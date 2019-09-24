############################################
## functions to compute turbulent velocity fields
###############################################


#.....single scale approach using spline
@everywhere function turb1(vx::Array{Float64,3},vy::Array{Float64,3},vz::Array{Float64,3},sigma::Array{Float64,3})


#for i in 1:i2-i1,j in 1:j2-j1, l in 1:l2-l1 

#@fastmath vsmx= std(vx[a1:a2,b1:b2,c1:c2])
#@fastmath vsmy= std(vy[a1:a2,b1:b2,c1:c2])
#@fastmath vsmz= std(vz[a1:a2,b1:b2,c1:c2])

#sigma=minvalue[vsmx,vsmy,vsmz]
#end

#@fastmath vtx=vx-interpolate(vx, BSpline(Cubic(Line())), OnCell())
#@fastmath vty=vy-interpolate(vy, BSpline(Cubic(Line())), OnCell())
#@fastmath vtz=vz-interpolate(vz, BSpline(Cubic(Line())), OnCell())

@fastmath vtx=vx-interpolate(vx, BSpline(Quadratic(Reflect())), OnCell())
@fastmath vty=vy-interpolate(vy, BSpline(Quadratic(Reflect())), OnCell())
@fastmath vtz=vz-interpolate(vz, BSpline(Quadratic(Reflect())), OnCell())


@fastmath sigma=(vtx*vtx+vty*vty+vyz*vtz)^0.5
vtx=nothing
vty=nothing
vtz=nothing

return sigma
end






#.....single scale approach using local velocity dispersion                                       
@everywhere function turb2(vx::Array{Float64,3},vy::Array{Float64,3},vz::Array{Float64,3},sigma::Array{Float64,3})

n3=size(vx)
scale=4 #..number of cells in 1D for the kernel

#I know it's written in an awful way, thank u
@inbounds for l in 1:n3[3]
for j in 1:n3[2]
@simd for i in 1:n3[1]
a1=i-scale
a2=i+scale
b1=j-scale
b2=j+scale
c1=l-scale
c2=l+scale
if a1 < 1
a1=1
end
if a2 > n3[1]
a2=n3[1]
end

if b1 <1
b1=1
end
if b2 >n3[2]
b2=n3[2]
end

    if c1 <1
c1=1
end
    if c2 >n3[3]
c2=n3[3]
end

@fastmath vtx= std(vx[a1:a2,b1:b2,c1:c2])                                     
@fastmath vty= std(vy[a1:a2,b1:b2,c1:c2])                                     
@fastmath vtz= std(vz[a1:a2,b1:b2,c1:c2])                                     

@fastmath sigma[i,j,l]=(vtx*vtx+vty*vty+vtz*vtz)^0.5

end
end
end

vtx=nothing
vty=nothing
vtz=nothing

return sigma
end
