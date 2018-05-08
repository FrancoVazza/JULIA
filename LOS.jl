root="/Users/francovazza/Dropbox/DATA/FRB/"

file1=string(root,"200Mpc_R_1Dbeam_z1_00",i,"_new.dat")

file3=string(root,"200Mpc_cqhhR_1Dbeam_z1_00",i,"_new.dat")
a1=readdlm(file1)   #each file is a N x M tabulated file, a1[:,1] is redshift, a1[:,2] density etc...                                                                
a3=readdlm(file3)


#bz=a1[:,6]  #bfield                                                                                                                                                 
#x=a1[:,1]  #x                                                                                                                                                       
#d=a1[:,2] #density                                                                                                                                                  

using PyPlot
figure(width=600,height=400)
p=FramedPlot(ylable="B[\\mu G]",xlabel="\\rho /<\\rho>",title="FRB")
add(p,Curve(a1[1:1000,1],a1[1:1000,6],color="red"))
file_out=string(root,"line",i,"_new.png")

savefig(file_out)


