using FITSIO
using Statistics
using SpecialFunctions

# Constants for frequency and related values
const freqf = ["0.1", "0.35", "0.7", "1.4", "1.76"]
const freqv = [0.1, 0.35, 0.7, 1.4, 1.76]  # ...GHz
const nf = size(freqf)
const df = 0.025  # GHz
const freq_min = 0.1  # GHz
const freq_max = 2.0  # GHz
const ng = 2400  # ...pixels in each image

# Function for interpolation
function interpIQU(Vin::Array{Float64, 2}, wei1::Float64, Vin2::Array{Float64, 2}, wei2::Float64, f2::Float64, f1::Float64)
    return @fastmath @views (Vin * wei1 + Vin2 * wei2) / (f2 + f1)
end

# Function for another type of interpolation
function interpIQU2(Vin::Array{Float64, 2}, wei1::Float64, Vin2::Array{Float64, 2}, wei2::Float64, f2::Float64, f1::Float64)
    return @fastmath @views Vin2 + (Vin2 - Vin) * (wei1) / (f2 - f1)
end

nfreq = convert(Int64, trunc((freq_max - freq_min) / df + 1))

I_in = Array{Float64}(undef, ng, ng, nf[1])
Q_in = Array{Float64}(undef, ng, ng, nf[1])
U_in = Array{Float64}(undef, ng, ng, nf[1])

I_out = Array{Float64}(undef, ng, ng, nfreq)
Q_out = Array{Float64}(undef, ng, ng, nfreq)
U_out = Array{Float64}(undef, ng, ng, nfreq)

println(nfreq)

fold = "/Users/francovazza/Desktop/data/DATA/CHRONOS++/100Mpc/2400/snap/SKA/highz/"
snap = ["0.025", "0.046", "0.069", "0.0925", "0.116", "0.141", "0.165", "0.189", "0.215", "0.244"]
z = [0.025, 0.046, 0.069, 0.0925, 0.116, 0.141, 0.165, 0.189, 0.215, 0.244]
nz = size(z)

for s in 1:nz[1] # ...loop over snapshots

    for ff in 1:nf[1] # ...loop over available frequencies
        file_stokes = string(fold, "IQUf_y-los_nu_", freqf[ff], "_M2.5_z", snap[s], "_zinter.fits")
        f1 = FITS(string(file_stokes), "r")
        map = read(f1[1])
        close(f1)

        I_in[:, :, ff] = map[:, :, 1]
        Q_in[:, :, ff] = map[:, :, 2]
        U_in[:, :, ff] = map[:, :, 3]
    end

    @inbounds for fi in 1:nfreq
        fis = string(fi)
        if fi < 10
            fis = string("00", fi)
        end
        if fi >= 10 && fi < 100
            fis = string("0", fi)
        end

        println(fi)
        @fastmath freq = freq_min + df * (fi - 1)

        @inbounds @simd for i in 1:nf[1] - 1
            if freqv[i] <= freq && freqv[i + 1] >= freq
                f1 = freqv[i]
                f2 = freqv[i + 1]
                wei1 = (f2 - freq)
                wei2 = (freq - f1)
                I_out[:, :, fi] = interpIQU(I_in[:, :, i], wei1, I_in[:, :, i + 1], wei2, f2, f1)
                Q_out[:, :, fi] = interpIQU(Q_in[:, :, i], wei1, Q_in[:, :, i + 1], wei2, f2, f1)
                U_out[:, :, fi] = interpIQU(U_in[:, :, i], wei1, U_in[:, :, i + 1], wei2, f2, f1)
            end
            if freq > freqv[5]
                f1 = freqv[4]
                f2 = freqv[5]
                wei1 = (freq - f2)
                wei2 = (freq - f1)
                I_out[:, :, fi] = interpIQU2(I_in[:, :, 4], wei1, I_in[:, :, 5], wei2, f2, f1)
                Q_out[:, :, fi] = interpIQU2(Q_in[:, :, 4], wei1, Q_in[:, :, 5], wei2, f2, f1)
                U_out[:, :, fi] = interpIQU2(U_in[:, :, 4], wei1, U_in[:, :, 5], wei2, f2, f1)
            end
        end
    end

    filep1 = string(fold, "IQUf_y-los_nu_all_freq_M2.5_z", snap[s], "_zinter.fits")
    f = FITS(filep1, "w")
    write(f, I_out)
    write(f, Q_out)
    write(f, U_out)
    close(f)
end
