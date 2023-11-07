using FITSIO
using Statistics
using SpecialFunctions

const freqf = ["0.1", "0.35", "0.7", "1.4", "1.76"]
const freqv = [0.1, 0.35, 0.7, 1.4, 1.76]  #...GHz
const nf=size(freqf)
const df=0.05  #GHz
const freq_min=0.1 #GHz
const freq_max=2.0 #GHz
const ng=2400 #...pixels in each image




function interpIQU_λ(Vin::Array{Float64,2},Vin2::Array{Float64,2},λ2::Float64,λ1::Float64,λ::Float64)
 return @fastmath @views Vin+(λ1-λ)*(Vin2-Vin)/(λ1-λ2)
end

function interpIQU_single(Vin::Float64,Vin2::Float64,f2::Float64,f1::Float64,freq::Float64)
 return @fastmath @views 10^(Vin+(freq-f1)*(Vin2-Vin)/(f2-f1))
end

function interpIQU2(Vin::Array{Float64,2},wei1::Float64,Vin2::Array{Float64,2},wei2::Float64,f2::Float64,f1::Float64)

 return @fastmath @views Vin2+ (Vin2-Vin)*(wei1)/(f2-f1)
end

nfreq=convert(Int64,trunc((freq_max-freq_min)/df+1))

I_in=Array{Float64}(undef,ng,ng,nf[1])
Q_in=Array{Float64}(undef,ng,ng,nf[1])
U_in=Array{Float64}(undef,ng,ng,nf[1])


I_out=Array{Float64}(undef,ng,ng,nfreq)
Q_out=Array{Float64}(undef,ng,ng,nfreq)
U_out=Array{Float64}(undef,ng,ng,nfreq)

println(nfreq)
    fold = "/Users/francovazza/Desktop/data/DATA/CHRONOS++/100Mpc/2400/snap/SKA/highz/"
    fold_out="/Users/francovazza/Dropbox/DATA/SKA_DC/"
    snap = ["0.025", "0.046", "0.069", "0.0925", "0.116", "0.141", "0.165", "0.189", "0.215", "0.244"]
    z = [0.025, 0.046, 0.069, 0.0925, 0.116, 0.141, 0.165, 0.189, 0.215, 0.244]


    nz=size(z)

for s in 1:nz[1] #...loop over snapshots

for ff in 1:nf[1] #...loop over available frequencies
    file_stokes=string(fold,"IQUf_y-los_nu_",freqf[ff],"_M2.5_z",snap[s],"_zinter.fits")
    f1=FITS(string(file_stokes),"r")
    map=read(f1[1])
    close(f1)


        I_in[:,:,ff]=map[:,:,1]
        Q_in[:,:,ff]=map[:,:,2]/6.5
        U_in[:,:,ff]=map[:,:,3]/6.5

    end

@inbounds for fi in 1:nfreq
    fis=string(fi)
    if fi <10
    fis=string("00",fi)
    end
    if fi >=10 && fi <100
    fis=string("0",fi)
    end

    println(fi)
 @fastmath freq=freq_min+df*(fi-1)

@inbounds @simd for i in 1:nf[1]-1
if freqv[i] <= freq && freqv[i+1]>=freq
f1=freqv[i]
f2=freqv[i+1]

    λ1=1/f1
    λ2=1/f2
    λ=1/freq
 I_out[:,:,fi]= interpIQU_λ(I_in[:,:,i],I_in[:,:,i+1],λ2,λ1,λ)
 Q_out[:,:,fi]= interpIQU_λ(Q_in[:,:,i],Q_in[:,:,i+1],λ2,λ1,λ)
 U_out[:,:,fi]= interpIQU_λ(U_in[:,:,i],U_in[:,:,i+1],λ2,λ1,λ)
end


if freq>freqv[5]
f1=freqv[4]
f2=freqv[5]
λ1=1/f1
λ2=1/f2
λ=1/freq
I_out[:,:,fi]= interpIQU_λ(I_in[:,:,4],I_in[:,:,5],λ2,λ1,λ)
Q_out[:,:,fi]= interpIQU_λ(Q_in[:,:,4],Q_in[:,:,5],λ2,λ1,λ)
U_out[:,:,fi]= interpIQU_λ(U_in[:,:,4],U_in[:,:,5],λ2,λ1,λ)
end

end
end
    filep1=string(fold_out,"I_y-los_nu_all_freq_M2.5_z",snap[s],"_zinter.fits")
    f = FITS(filep1, "w")
    write(f,I_out)
    close(f)
    filep1=string(fold_out,"Q_y-los_nu_all_freq_M2.5_z",snap[s],"_zinter.fits")
    f = FITS(filep1, "w")
    write(f,Q_out)
    close(f)
    filep1=string(fold_out,"U_y-los_nu_all_freq_M2.5_z",snap[s],"_zinter.fits")
    f = FITS(filep1, "w")
    write(f,U_out)
    close(f)


    end
