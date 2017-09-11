## DWD evolution analysis

This FORTRAN code finds the ejecta/circumstellar/binary masses plus their energies in 10<sup>46</sup> ergs. The code was developed for double white dwarf formation, and it assumes

* Similar masses in their code, i.e. the RG + Companion has similar He cores
* It also assumes that the He core are modelled with a point mass
* Also assumes that the separation between the stars are given by the distance between the 2 point masses

This codes will generate few DAT files, which contains the evolution for each file, the files are as follows,

* conservation.dat contains conservation quantity checks, cols: time [day], E\_tot [10^46 erg], J [code units]
* ejecta.dat contains the ejecta quantities, cols: time [day], Munb [Msun], E\_kin^unb [10^46 erg], E\_int^unb [10^46 erg], E\_pot^unb [10^46 erg], # of particles 
* binary.dat contains binary quantities, cols: time [day], Mbin [Msun], E\_kin^bin [10^46 erg], E\_int^bin [10^46 erg], E\_pot^bin [10^46 erg], # of particles
* ang\_mom.dat contains angular momentum quantities, cols: t [day], lz, lz\_bin, lz\_cir, lz\_unb in units of g cm^2 s^-1 10^51  
* circumbinary.dat contains the circumbinary matter, cols: time [day], Mcir [Msun], E\_kin^cir [10^46 erg], E\_int^cir [10^46 erg], E\_pot^cir [10^46 erg], # of particles 
* orbit.dat cotains orbital values, col: Nout, time [day], a [Rsun], E\_tot^bound, E\_tot^unb, x [Rsun of RG],y [Rsun of RG],z [Rsun of RG], vx [code units of RG], vy [code units of RG], vz [code units of RG], x [Rsun of WD comp], y [Rsun of WD comp], z [Rsun of WD comp],vx [code units of WD comp], vy [code units of WD comp], vz [code units of WD comp], h [smoothing length of RG in Rsun], h [smoothing length of WD comp in Rsun], m [mass of RG He core in Msun], m [mass of WD comp in Msun].

