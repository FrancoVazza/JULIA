#....simple map making routine for colors and contours out of a sequence of
#...  snapshots from an ENZO simulations, originally saved in FITS format

# SOME SIMULATION PARAMETERS AND CONSTANTS
 const n=640 #grid size
    min_pradio=1e29 #.....min and max values for final map 2
    max_pradio=1e33
    min_x=1e41     #.....min and max values for final map 1
    max_x=1e44
    res=0.008    #...resolution in Megaparsec
    lbox=res*n   #...side of the box in physical units
    #...DOMAIN EDGES FOR MAP MAKING
     j1=150
     i1=150
     j2=520
     i2=520
     sin=105  #...initial snapshot
     sfin=193 #...final snapshot

 const root="/DATA/ESCAPE/0/0_1/tracers_all_radio/" #..folder for input data 1
 const rootx="/DATA/ESCAPE/0/0_1/D/"   #...folder for input data 2


#MODULES
   using PyPlot
   using PyCall
   using FITSIO


##############
#MAIN PROGRAM#
##############

for mo in sin:sfin

  println("doing snapshot=",mo)
  #....two sets of FITS images to extract maps for the final plot
  file1=string("map_radio_tracers0_1_",mo,"CSTY_momenta.fits")
  file2=string("Dbh56HDD_dt_",mo,"_dtb_proj_proj5_sliY_lx.fits")
  #....reading FITS files
  file_fits=string(root,file1)
  f1=FITS(file_fits,"r")

  file_fits=string(rootx,file2)
  f2=FITS(file_fits,"r")
  map1=read(f1[1])
  map2=read(f2[1])
  maps_radio=(transpose(map1[:,:,2]))
  maps_x=(transpose(map2[:,:,7]))

 maps_radio=convert(Array{Float64,2},maps_radio)
 maps_x=convert(Array{Float64,2},maps_x)


#...first map : X-ray in colors with radio contours in overlay
        pcolormesh(1e-3*(j1*dx:dx:j2*dx),1e-3*(i1*dx:dx:i2*dx),mapd2[i1:i2,j1:j2],norm=matplotlib[:colors][:LogNorm](vmin=min_x,vmax=max_x),cmap="Oranges")
            axis("scaled")
            title(string(tag[model1]))
            xticks(fontsize=9)
            yticks(fontsize=9);
            cbar=colorbar()
            cbar[:set_label](L"F_X[erg/s/pixel]",fontsize=13)
        contour(1e-3*(j1*dx:dx:j2*dx),1e-3*(i1*dx:dx:i2*dx),sqrt.(mapd1[i1:i2,j1:j2]),colors="black",linewidths=0.4)
             xlabel(string("x ",L"[\mathrm{Mpc}]"),fontsize=12)
             ylabel(string("y ",L"[\mathrm{Mpc}]"),fontsize=12)
             filep1=string(root,file,"_Xrayradio_over",mo,".png")
             savefig(filep1,dpi=400)
#...second map : radio in colors with X-ray  contours in overlay
             clf()

        pcolormesh(1e-3*(j1*dx:dx:j2*dx),1e-3*(i1*dx:dx:i2*dx),mapd1[j1:j2,i1:i2],norm=matplotlib[:colors][:LogNorm](vmin=min_pradio),cmap="Blues")
              axis("scaled")
              title(string(tag[model2]))
              xticks(fontsize=9)
              yticks(fontsize=9);
              cbar=colorbar()
              cbar[:set_label](L"F_R(140MHz)[erg/s/Hz/pixel]",fontsize=13)
        contour(1e-3*(j1*dx:dx:j2*dx),1e-3*(i1*dx:dx:i2*dx),sqrt.(mapd2[i1:i2,j1:j2]),colors="red",linewidths=0.4,levels=20)
               xlabel(string("x ",L"[\mathrm{Mpc}]"),fontsize=12)
               ylabel(string("y ",L"[\mathrm{Mpc}]"),fontsize=12)
              filep1=string(root,file,"_radiooverX_",mo,"_",model_list[model2],".png")
              savefig(filep1,dpi=400)
             clf()

end
