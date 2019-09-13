###########################################################
##### SZ decrement (ref: Brunetti+04, JKAS)
##....

@everywhere function SZ(d::Array{Float64,3},t::Array{Float64,3},freqSZ::Float64,sz_dec::Array{Float64,3}) 

#freqSZ in Hz
#cell in Mpc
#d in part/cm^3
#t in K

const term=(kb^2.*(tcmb*(1.+z))*sigmat/(me*vc^2.))*(freqSZ/vc)^2.*res^3.*cmtoMpc*1e26/(mp*mu)
@inbounds @simd for i in eachindex(d)
sz_dec[i]= d[i]*t[i]*term 
end

return sz_dec   #....Jy/std
#Y_dec*4.25e10/float(3600) ;..gives decrement in Jy to arcmin^2
end
 
