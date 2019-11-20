using FITSIO
using DelimitedFiles: readdlm, writedlm

# Simulation parameters
lbox = 100.   # 1D box size in Mpc
n0 = 2400     # 1D box size in cells
res = lbox / n0 # spatial resolution

snaps = ["188","166","156"]   #...set of available snapshots

#....observable-like parameters (all assumes concordance LCDM from Cosmology Calculator
fov = 1.      #..field of view [degrees]
dl = 14.      #..1D size of the fov at z=0.01 [Mpc]
adist0 = 50.  #..angular distance at z=0.01 [Mpc] distance at z=0.01
nt = 8 #number of snapshot slices
nima = 2400
ima_gal = zeros(nima + 1, nima + 1, nt)
lbox2 = lbox
#...sequence of redshift snapshots as in TRECS catalog
zz = vcat([0.01, 0.02], range(0.05, stop=1, step=0.05), range(1.2, stop=10, step=0.2))

#.....file with pre-computed luminosity distances for TRECS redshfits
filed = string(pwd(), "/distl_z.dat")
a = readdlm(filed)
distl = a[:,2]
a = nothing

#.....loop over snapshots
box_number = 1   #....this parameter varies the centering of the FOV in the sky model (see below)

#...there arrays will contain the 2D positions and power of the final sky model
x0 = Array{Float64}(undef, 1)  #...x,y positions in degrees
y0 = Array{Float64}(undef, 1)
logm0 = Array{Float64}(undef, 1)   #.....power
nbig0 = 0


#...loop over nt snapshots until redshift zz[nt]
for ii in 1:nt - 1

    ncloud = 1  #...ncloud sets whether we need to distribute the info of one pixel in the input sky model into ncloud^2 smaller pixels, as a function of the final resolution
    #....the following choices for ncloud are referred to a 15" resolution and were just hardcoded
    if zz[ii] == 0.01   #....referred to a 15" resolution observation
        ncloud = 18
    end
    if zz[ii] == 0.02
        ncloud = 12
    end
    if zz[ii] == 0.05
        ncloud = 6
    end
    if zz[ii] == 0.1
        ncloud = 2
    end
    #...for larger redshift the input pixel size is < than the beam size

    #...output file with information on pixels for the MWA map making
    file_gal = string(pwd(), "/cone_1x1_z", zz[ii], "txt_web", box_number, ".dat")
    so2 = open(file_gal, "w")


#...in principle, here we choose which snapshots slice and projection along the LOS to use in order to minimize repeating structures. In this case however just one LOS and redshift, as example

    if zz[ii] >= 0 && zz[ii] <= 0.1
        los = "Y"
        filefits = string(pwd(), "/map_allD_", snaps[1], "_1024_", los, "newHB2.fits")
        pthr = 1e23   #lower limit on pixels to use (i.e. pixels with Flux<pthr are not used) in [erg/s/Hz]
    end

    if zz[ii] > 0.1 && zz[ii] <= 0.2
        los = "Y"
        filefits = string(pwd(), "/map_allD_", snaps[1], "_1024_", los, "newHB2.fits")
        pthr = 1e23
    end

    if zz[ii] > 0.2
        los = "Y"
        filefits = string(pwd(), "/map_allD_", snaps[1], "_1024_", los, "newHB2.fits")
        pthr = 1e24
    end


#......here we read from the appropriate sky model with the cosmic web emission

    imaf = FITS(filefits, "r")
    ima = read(imaf[1])
    close(imaf)

    iradio = 9  #...plane in the FITS file which contains the radio emission
    ispec = 5   #...weighted distribution of Mach number -> spectral indices
    radio = ima[:,:,iradio]
    mach = ima[:,:,ispec]
    ima = nothing
    ν_in = 1400.    #input frequency in the sky model [MHz]
    ν_out = 155.     #desired output frequency  [MHz]

    for i in eachindex(mach)  #...loops over cells with a potential signal
    #..cures for too small or too large Mach numbers
        if mach[i] <= 2.0
            mach[i] = 2.0
        end
        if mach[i] >= 10.
            mach[i] = 10.
        end

        @fastmath  α = 0.5 * (mach[i]^2 + 1) / (mach[i]^2 - 1) + 0.5   #...radio spectral index
        @fastmath  radio[i] = (radio[i] * ((ν_in / ν_out)^(α))) #....correcting the sky model for

    end

    id = findall(radio .> pthr)    #find index of all pixels brighter than p_thr [erg/s/Hz]

    ng = size(id)
    x = Array{Float64}(undef, ng[1])
    y = Array{Float64}(undef, ng[1])
    logm = Array{Float64}(undef, ng[1])

    for jj in eachindex(id)
        ij = CartesianIndices((n0, n0))[id[jj]]
        # ij = ind2sub((n0, n0), id[jj])   #...this converts id index into x,y information
        x[jj] = ij[1] * res
        y[jj] = ij[2] * res
        logm[jj] = radio[id[jj]]
    end
    #....these temp.arrays are needed if we need to produce replicas of the structures to fill high-z parts of the lightcone
    x0 = x
    y0 = y
    logm0 = logm

    #...sets redshift edges of each slice
    zz_pre = zz[ii]
    if ii == 1
        zz_pre = 0.
    end
    zz_post = zz[ii + 1]
    zz_mean = (zz_post + zz_pre) * 0.5
    dzz = zz_post - zz_pre
    zeds = zz_mean

    println("doing tile with mean redshift z=", zeds)
    ldist = distl[ii]      #...luminosity distance
    adist = ldist / (1 + zz[ii]^2)   #..angular size distance

    dla = dl * (adist / (adist0))      #...proportion between the side in Mpc in the FOV at z=0.01 and in all subsequent redshifts
    shift = 20                     #...arbitrary displacement to point our FOV in different positions as a function of box_number

    x1 = 1 + shift * box_number        #...some fuzzy way of randomly changing the sky model edges along the sequence (all meant to avoid the overlap of the same structures along the LOS)
    x2 = x1 + dla
    y1 = 90 - shift * box_number
    y2 = y1 + dla

#....computing whether a replica of the galaxy distribution (assuming periodicity) is needed
    nbig = convert(Int32, trunc(x2 / (lbox2 + 0.1)))   #...if x2>lbox, our input sky model is not large enough, and we need to replicate it nbig times
    if nbig >= 1
        println("replicating volume ", nbig, "times")
        for l in 1:nbig
            x00 = y0
            y00 = x0
            @views        append!(x, y00 .+ (l * lbox2))
            @views        append!(x, y00 .+ (l * lbox2))

        #....swapping of x and y wrt to the previous (ii-1) step, to minimize repeating patterns
        #....the following sequence is tasselate the entire 2D extension of the FOV
            append!(x, y00)
            append!(y, x00 .+ (l * lbox2))
            append!(y, x00)
            append!(y, x00 .+ (l * lbox2))

            append!(logm, logm0)
            append!(logm, logm0)
            append!(logm, logm0)
        end
        nbig0 = nbig
    end

    ngg = size(x)

    for i in 1:ngg[1]   #....loop over the (augmeted by replication, or not) distribution of pixels

        if x[i] >= x1 && x[i] <= x2 && y[i] >= y1 && y[i] <= y2
            pflux = logm[i] / (ncloud^2.)  #...weights the pixel contribution for the number of ncloud^2 resampling that were ncessary
            for ll in 1:ncloud
                for jj in 1:ncloud
        #...we fill the area in a square-like tasselation + random displacement
                    @fastmath          xi = x[i] - res * 0.5 + ll * res / ncloud + 0.75 * res * rand()
                    @fastmath          yi = y[i] - res * 0.5 + jj * res / ncloud + 0.75 * res * rand()

                    xo = (xi - x1) / (x2 - x1)
                    yo = (yi - y1) / (y2 - y1)
                    xu = convert(Int64, trunc((nima * xo)))
                    yu = convert(Int64, trunc((nima * yo)))

                    xs = convert(Float32, fov * (xo - 0.5))
                    ys = convert(Float32, fov * (yo - 0.5))
                    logmo = convert(Float32, pflux)
                    writedlm(so2, [logmo zz[ii] xs ys])  #...write on disk all (galaxy-like) data for this redshift slice


        #....additional map making as sanity check
                    if xu > nima
                        xu = nima
                    end
                    if yu > nima
                        yu = nima
                    end
                    if xu <= 1
                        xu = 1
                    end
                    if yu <= 1
                        yu = 1
                    end

                    ima_gal[xu,yu,ii] += ((pflux) / (4 * pi * ldist^2 * 1e23 * 35.1^2)) #generate a map in ...Jy/arcsec^2)

                end
            end
        end
    end


    close(so2)
end

filep4 = string(pwd(), "/map_web", box_number, ".fits")
f = FITS(filep4, "w");
write(f, ima_gal)
close(f)
