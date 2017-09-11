      subroutine read_file(ni)
      include 'common_sph_var.h'
      integer ni,i,ncheck
      real*8 dummy
      logical it_exists,it_existsb
      character*13 filename
      write(filename,100)ni
 100  format('out',I4.4,'.sph')
 101  format('out',I5.5,'.sph')
      inquire(file=filename,exist=it_exists)
      if(it_exists) then 
          open(10,file=filename,form='unformatted')
          print*,'Reading:',filename
      else
          write(filename,101)ni
          inquire(file=filename,exist=it_existsb)
          if (it_existsb) then 
               open(10,file=filename,form='unformatted')
               print*,'Reading:',filename
          else
               write(*,*)'Check if the file exists'
               stop
          endif
      endif
      read(10)ntot,nnopt,hmin,hmax,sep0,tf,dtout,nout,nit,t,
     &          nav,alpha,beta,tskipahead,ngr,nrelax,trelax,dt,omega2
     
      if (nout.eq.0) then 
          write(*,102)ntot,nnopt,hmin,hmax,sep0,tf,dtout,nout,nit,t,
     &          nav,alpha,beta,tskipahead,ngr,nrelax,trelax,dt,omega2
 102  format(I6,I6,ES15.6,ES15.6,ES15.6,ES15.6,ES15.6,I6,I6,ES15.6,
     &      I6,ES15.6,ES15.6,ES15.6,I6,I6,ES15.6,ES15.6,ES15.6)
      endif
 
      do i=1,ntot
         read(10) x(i),y(i),z(i),m(i),h(i),rho(i),vx(i),vy(i),vz(i),
     &            vxdot(i),vydot(i),vzdot(i),u(i),udot(i),grpot(i),
     &            dummy,cc(i),divv(i)
         mu(i)=dummy/mp
      enddo
      read(10) ncheck
      if (ncheck.ne.ntot) then 
          write(*,*)'The file is corrupt'
          stop
      endif
      write(*,*)filename,' read!'
      close(10)
      end
