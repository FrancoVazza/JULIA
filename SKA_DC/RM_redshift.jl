# Import necessary packages
using FITSIO
using Statistics
using SpecialFunctions
using Random
using DelimitedFiles

# Constants and parameters
Lbox = 85.0
h = 0.678
nboxb = [10, 8, 10, 7, 10, 10, 3]  # Number of boxes to generate for each redshift interval

zin = [2.441, 1.84, 1.455, 1.065, 0.835, 0.555, 0.312]
zfin = [1.84, 1.455, 1.065, 0.835, 0.555, 0.312, 0.244]

zed = ["004", "004", "004", "005", "005", "006", "006"]
los = ["X", "Y", "Z"]

nz = size(zfin)

# Path for file storage
fold_ska = "/Users/francovazza/Dropbox/DATA/SKA_DC/z0.244-2.5_85Mpc/"

const ng = 1024  # Pixels in each image

# Loop over redshift intervals
for s in 1:nz[1]
    dz = (zfin[s] - zin[s]) / nboxb[s]

    # Initialize variables
    bi = 0
    z0 = zfin[s]
    
    for b in 1:nboxb[s]
        z0 -= dz
        zs = string(convert(Int64, trunc(z0 * 10000)) / 10000)
        println(zs)

        # Random extraction of LOS
        seedp = s * b
        rng = MersenneTwister(seedp)
        aa = rand!(rng, zeros(1))
        whi = convert(Int64, 1 + trunc(aa[1] * 3))
        println(whi)
        println(los[whi])

        # Construct file path for RM map
        file_rm = string(fold_ska, "RM_map_allR_", zed[s], "_1024_", los[whi], "_SKA_FULL.fits")
        f1 = FITS(string(file_rm), "r")
        rm = read(f1[1])
        close(f1)

        # Cosmological 1/(1+z)^2 correction to RM
        rm ./= (1 + z0)^2.0

        # Read data from files
        filec = string(fold_ska, "RM_peak_", zed[s], "_LOS", los[whi], ".dat")
        rmp = readdlm(filec)  # x, y, z

        # Write data to FITS file
        filep1 = string(fold_ska, "rmy_simp_cool_full_IQUf_y-los_nu_all_freq_M2.5_z", zs, "_zinter.fits")
        f = FITS(filep1, "w")
        write(f, rm)
        close(f)

        # Write data to output file
        file_out = string(fold_ska, "RM_peak_z", zs, ".dat")
        writedlm(file_out, [rmp[1] rmp[2] rmp[3] rmp[4] rmp[5] rmp[6]])
    end
end
