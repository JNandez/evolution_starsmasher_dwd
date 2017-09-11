      subroutine ejecta
      include 'common_sph_var.h'
      integer i
C Bound material
      real*8 lxb,lyb,lzb,mb
      real*8 ekb,epb,eib
C Ejecta
      real*8 ep,munb,ekunb,epunb,eiunb
      real*8 lxunb,lyunb,lzunb 
C Total angular momentum and energy
      real*8 lx,ly,lz,etot
C Velocities
      real*8 vtan,vrad,r
C Circumbinary disk
      real*8 mc,lxc,lyc,lzc,ekc,epc,eic,rmin,ravg
C Number of particles
      integer cmb,cmunb,cmc
 
      real*8 a,m1,m2,etotbound,etotunb
      common/sep_max/a,m1,m2   

      cmb=0
      cmunb=0
      cmc=0
      munb=0.d0
      ekunb=0.d0
      epunb=0.d0
      eiunb=0.d0
      ekb=0.d0
      epb=0.d0
      eib=0.d0
      etot = 0.d0
      lx=0.d0
      ly=0.d0
      lz=0.d0
      lxunb=0.d0
      lyunb=0.d0
      lzunb=0.d0
      lxb=0.d0
      lyb=0.d0
      lzb=0.d0
      mb=0.d0
      lxc=0.d0
      lyc=0.d0
      lzc=0.d0
      mc=0.d0
      ekc=0.d0
      epc=0.d0
      eic=0.d0
      rmin=1d30
      ravg=0.d0
      do i=1,ntot
C Get the total energy of a particle
c         ep=0.5d0*(vx(i)**2+vy(i)**2+vz(i)**2)+grpot(i)+u(i)
C Get the tangential and radial velocity
         r = sqrt(x(i)**2+y(i)**2+z(i)**2)
         vrad = x(i)*vx(i)+y(i)*vy(i)+z(i)*vz(i)
         vrad = vrad/r
         vtan = sqrt((y(i)*vz(i)-z(i)*vy(i))**2
     &              +(z(i)*vx(i)-x(i)*vz(i))**2
     &              +(x(i)*vy(i)-y(i)*vx(i))**2)
         vtan = vtan/r
C Check if the total energy of a particle is positive, if so, get ejecta
         if (id(i).eq.4) then
            munb=munb+m(i)
            ekunb=ekunb+0.5d0*m(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
            eiunb=eiunb+m(i)*u(i)
            epunb=epunb+m(i)*grpot(i)
            lxunb = lxunb + m(i)*(y(i)*vz(i)-z(i)*vy(i))
            lyunb = lyunb + m(i)*(z(i)*vx(i)-x(i)*vz(i))
            lzunb = lzunb + m(i)*(x(i)*vy(i)-y(i)*vx(i))
            cmunb = cmunb + 1
         elseif(id(i).le.2) then 
            !write(70,*)i,x(i),y(i),z(i),vrad,vtan,vx(i),vy(i),vz(i)
            mb = mb + m(i)
            ekb=ekb+0.5d0*m(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
            eib=eib+m(i)*u(i)
            epb=epb+m(i)*grpot(i)
            lxb = lxb + m(i)*(y(i)*vz(i)-z(i)*vy(i))
            lyb = lyb + m(i)*(z(i)*vx(i)-x(i)*vz(i))
            lzb = lzb + m(i)*(x(i)*vy(i)-y(i)*vx(i))
            cmb = cmb + 1
         else
            mc = mc + m(i)
            ekc=ekc+0.5d0*m(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
            eic=eic+m(i)*u(i)
            epc=epc+m(i)*grpot(i)
            lxc = lxc + m(i)*(y(i)*vz(i)-z(i)*vy(i))
            lyc = lyc + m(i)*(z(i)*vx(i)-x(i)*vz(i))
            lzc = lzc + m(i)*(x(i)*vy(i)-y(i)*vx(i))
            cmc = cmc + 1
            rmin = min(rmin,sqrt((x(1)-x(i))**2+(y(1)-y(i))**2+(z(1)-z(i))**2 ))
            ravg = ravg +sqrt((x(1)-x(i))**2+(y(1)-y(i))**2+(z(1)-z(i))**2 )!*m(i)
         endif
C Get the total angular momentum and energy
         etot = etot + 0.5d0*m(i)*(vx(i)**2+vy(i)**2+vz(i)**2)+0.5d0*m(i)*grpot(i)+m(i)*u(i)
         lx = lx + m(i)*(y(i)*vz(i)-z(i)*vy(i))
         ly = ly + m(i)*(z(i)*vx(i)-x(i)*vz(i))
         lz = lz + m(i)*(x(i)*vy(i)-y(i)*vx(i))
      enddo
      write(50,*)t*tunit,etot*eunit,sqrt(lx**2+ly**2+lz**2)
      write(51,*)t*tunit,munb,ekunb*eunit,eiunb*eunit,epunb*eunit,cmunb
      write(52,*)t*tunit,mb,ekb*eunit,eib*eunit,epb*eunit,cmb
      write(55,*)t*tunit,mc,ekc*eunit,eic*eunit,epc*eunit,cmc
      write(53,*)t*tunit,lz*lunit,lzb*lunit,lzc*lunit,lzunb*lunit
cyc,
c     &      lzc,sqrt(lx**2+ly**2+lz**2)*lunit,
c     &      sqrt(lxunb**2+lyunb**2+lzunb**2)*lunit,
c     &      sqrt(lxb**2+lyb**2+lzb**2)*lunit,
c     &      sqrt(lxc**2+lyc**2+lzc**2)*lunit
      etotbound=(ekb+eib+0.5*epb+ekc+eic+0.5*epc)*eunit
      etotunb=(ekunb+eiunb+0.5*epunb)*eunit
      write(54,*)nout,t*tunit,a,etotbound,etotunb,x(1),y(1),z(1),vx(1),vy(1),vz(1),
     &    x(ntot),y(ntot),z(ntot),vx(ntot),vy(ntot),vz(ntot),h(1),h(ntot),m(1),m(ntot)
      write(*,*)'Distance circumbinary min = ',rmin,ravg/cmc

      end
