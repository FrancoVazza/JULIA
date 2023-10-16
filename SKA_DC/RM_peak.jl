using FITSIO
using HDF5

# Constants for frequency and related values
const freqf = ["0.1", "0.35", "0.7", "1.4", "1.76"]
const freqv = [0.1, 0.35, 0.7, 1.4, 1.76]  # ...GHz
const nf = size(freqf)
const df = 0.025  # GHz
const freq_min = 0.1  # GHz
const freq_max = 2.0  # GHz
const ng = 2400  # ...pixels in each image
const res = 0.0416  # ...kpc
const mu = 0.59
const mp = 1.67e-24  # proton mass

nfreq = convert(Int64, trunc((freq_max - freq_min) / df + 1))

println(nfreq)

# Define file paths
fold_sim = "/home/PERSONALE/franco.vazza2/Desktop/DATA/CHRONOS++/100Mpc/2400/"
fold_ska = "/snap/SKA_DEEP/"
fold_snap = "/conv/"
mod = "D"

# Redshifts and corresponding snapshots
zs = ["0.025", "0.046", "0.069", "0.0925", "0.116", "0.141", "0.165", "0.189", "0.215", "0.244"]
snap = [188, 188, 188, 188, 188, 166, 166, 166, 166, 166]
z = [0.025, 0.046, 0.069, 0.0925, 0.116, 0.141, 0.165, 0.189, 0.215, 0.244]
nz = size(z)

# Loop over snapshots
for s in 1:nz[1]
    # Construct the filename for density and B-field
    filec = string(fold_sim, "/conv/DD0", snap[s], ".conv2")
    println(filec)
    a = readdlm(filec)
    zed = a[3]  # redshift
    cd = a[4]  # conversion factor for density
    cv = a[5]  # conversion factor for velocity
    cb = sqrt(cd * 4.0 * pi) * cv * (1 + zed)^2.0  # b-field in Gauss

    # Read RM data from FITS file
    file_rm = string(fold_sim, fold_ska, "rmy_simp_cool_full_z", zs[s], "_zinter.fits")
    f1 = FITS(string(file_rm), "r")
    rm = read(f1[1])
    close(f1)
    rm .= abs.(rm)

    rm_max = maximum(rm)

    # Find peak location in RM data
    iw = find(x -> (x == rm_max), rm)
    ii = ind2sub((ng, ng), iw)

    x = ii[1]
    y = ii[2]
    println(x, " ", y)

    # Read HDF5 files for density and Bz
    file1 = string(fold_sim, "/snap/full_dens_", mod, "_s", snap[s])
    file3 = string(fold_sim, "snap/full_bz_", mod, "_s", snap[s])

    d = cd * h5read(file1, "Density", (x[1], y[1], :))
    bz = cb * h5read(file3, "Bz", (x[1], y[1], :))
    d = convert(Array{Float64, 3}, d)
    bz = convert(Array{Float64, 3}, bz)

    # Compute and process RM data
    rmd = abs.(d .* bz) * (812. * res / (1e-6 * 1e-3 * mu * mp))
    rmd_maximum = maximum(rmd)
    println(size(rmd))

    iw = find(x -> (x == rmd_maximum), rmd)
    ipeak = ind2sub((ng, ng), iw)
    println(ipeak[1])
    println(d[ipeak[1]])
    println(bz[ipeak[1]])
    println(rmd[ipeak[1]])

    # Write results to file
    file_out = string(fold_sim, "RM_peak_", zs[s], ".dat")
    writedlm(file_out, [x y ipeak[1] d[ipeak[1]] bz[ipeak[1]] rmd[ipeak[1]]])
end
