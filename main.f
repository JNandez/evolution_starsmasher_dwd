      program para_of_WD
      include 'common_sph_var.h'
      integer i,ni,nf,nd
      real*8 s1,s2,s1a,s2a
      write(*,*)'Give the initial value of the dump-file,final,delta'
      read(*,*)ni,nf,nd
      write(*,100)'Units(t [day],E [erg],v [km/s], lunit):',tunit,eunit,vunit,lunit
 100  format(A,ES15.6,ES15.6,ES15.6,ES15.6)
      open(50,file='conservation.dat')
      open(51,file='ejecta.dat')
      open(52,file='binary.dat')
      open(53,file='ang_mom.dat')
      open(55,file='circumbinary.dat')
      open(54,file='orbit.dat')
      do i=ni,nf,nd
        call time_elapsed(s1a)
        call read_file(i)
        call time_elapsed(s2a)
        !write(*,*)'It takes ',s2a-s1a, ' sec to read out*.dat'
        call classification
        call ejecta
        !call orbital_parameter
        write(*,*)'It takes ',s2a-s1a, ' sec to read out*.dat'
      enddo
      close(50)
      close(51)
      close(52)
      close(53)
      close(54)
      end
