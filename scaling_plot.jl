fileo="/Users/francovazza/Desktop/articles/proposals/weak.png"
fileo2="/Users/francovazza/Desktop/articles/proposals/strong.png"

np=[32.0,64.0,128.0,256.0,512.0,1024.0]
grid1=[322,406,512,644,812,1024]
time=[99.29,103.,107.15,134.88,318.,1846.]
tperf=[99.29,99.29,99.29,99.29,99.29]

np2=[5.0,40.,320.]
time2=[50.,38.,67.]
grid2=[864.*0.5,864.,864.*2.]
using PyPlot,Colors

 clf()
rc("font", family="arial")
 loglog(np[1:5],time[1:5],color="#0f87bf",linestyle="-",marker="o")
 plot(np[1:5],tperf[1:5],color="red",linestyle=":")

axis("tight")
ylabel("Wall Clock Time [s]",size=15)
xlabel("# cores",size=15)
title("weak scaling on Jureca",size=18)
 axis([20,600,50,1000])

  savefig(fileo)

clf()

np=[7.0,14.0,28.0,56.0,112.0]
time=5.*[2.38e1,1.411546e+01,7.623709e+00,4.500006e+00,2.90]/(1.75)

tperfect=[time[1],time[1],time[1],time[1],time[1]]

for i in 1:5
tperfect[i]=tperfect[i]/(np[i]/7)
end
#tperfect=[789.000,      394.500,      197.250,      98.6250,      49.3125]                                                                                          

rc("font", family="arial")
 loglog(np,time,color="#0f87bf",linestyle="-",marker="o")
 plot(np,tperfect,color="red",linestyle=":")
ylabel("Wall Clock Time [s]",size=15)
xlabel("# cores",size=15)
title("strong scaling on Jureca",size=18)
 axis([2,200,4,1e2])

  savefig(fileo2)


