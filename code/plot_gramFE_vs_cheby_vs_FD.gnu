set logscale y
plot 'gamer_output/Gaussian_Dens_000001_N1024_FCm12' u 1:(abs($2-$3)) w l,\
     'gamer_output/Gaussian_Dens_000001_N4096' u 1:(abs($2-$3)),\
     'gamer_output/Gaussian_Dens_000001_N1024_FCm10' u 1:(abs($2-$3)) w l,\
     'output/40.dat' u 1:(abs($2-$3)) w lp title 'Cheby N=144'    