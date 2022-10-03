
@everywhere  function total_spectra(pval::Array{Float64},pe::Array{Float64},pe_tot::Array{Float64},nn1::Int64,nn2::Int64,x::Array{Float64},y::Array{Float64},z::Array{Float64})

   @inbounds @simd   for i in nn1:nn2
      @inbounds @simd for g in 1:np
    pe_tot[g]+=(pe[g,i]*10^pval[g])
    end
    end

return pe_tot
end

@everywhere function map2d(pradio::Array{Float64,2},x::Array{Float64},y::Array{Float64},z::Array{Float64},nn1::Int64,nn2::Int64,freqf::Array{Float64},map::Array{Float64,3},xc::Float64,pee::Array{Float64,2})

   @inbounds @simd   for i in nn1:nn2

      @fastmath    x1=convert(Int32,trunc(xc*x[i]))
      @fastmath    y1=convert(Int32,trunc(xc*y[i]))
      @fastmath    z1=convert(Int32,trunc(xc*z[i]))
   if x1 <=3
   x1=3
   end
   if y1 <=3
   y1=3
   end
   if z1 <=3
   z1=3
   end
   if x1 >=n
   x1=n
   end
   if y1 >= n
   y1=n
   end
   if z1 >= n-3
   z1=n-3
   end

   @inbounds   map[x1,y1,1]+=pradio[i,1]
   @inbounds   map[x1,y1,2]+=pradio[i,2]
   @inbounds   map[x1,y1,3]+=pradio[i,3]
   @inbounds   map[x1,y1,4]+=pradio[i,4]
   @inbounds   map[x1,y1,5]+=pradio[i,5]
   @inbounds   map[x1,y1,6]+=pradio[i,6]
   @inbounds   map[x1,y1,7]+=pradio[i,7]
   @inbounds   map[x1,y1,8]+=pee[2,i]
   @inbounds   map[x1,y1,9]+=pee[10,i]
   @inbounds   map[x1,y1,10]+=pee[100,i]
end
return map
end
