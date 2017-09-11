      subroutine classification
      include 'common_sph_var.h'
      integer i
      real*8 etest,sep,phi1,xL1
      real*8 xL1p,yL1p,theta,zL1p
      real*8 rcmi,r1,r2
      integer idmax1,idmax2 
      real*8 rhomax1,rhomax2
      real*8 r1i,r2i,m1,m2,c1,c2
      real*8 epot1i,epot2i,cent 
      real*8 cmmbin,cmxbin,cmybin,cmzbin
      real*8 cmvxbin,cmvybin,cmvzbin
      real*8 x1,x2,y1,y2,z1,z2,vx1,vx2,vy1,vy2,vz1,vz2
      real*8 Ixx,Iyy,Izz,Ixy,Iyz,Ixz
      real*8 lx12,ly12,lz12,omegax,omegay,omegaz
      real*8 bin,binp
      integer count1,count2
      common/sep_max/sep,m1,m2
      common/rho_max/idmax1,idmax2,rhomax1,rhomax2
      common/l_mom/lx12,ly12,lz12
      common/cofm_bin/cmmbin,cmxbin,cmybin,cmzbin,cmvxbin,cmvybin,cmvzbin 
      common/bin_par/x1,x2,y1,y2,z1,z2,vx1,vx2,vy1,vy2,vz1,vz2

      rhomax1=-1d30
      rhomax2=-1d30 
      cmmunb = 0.d0
      cmxunb = 0.d0
      cmyunb = 0.d0
      cmzunb = 0.d0 
      cmvxunb = 0.d0
      cmvyunb = 0.d0
      cmvzunb = 0.d0
      cmmb = 0.d0
      cmxb = 0.d0
      cmyb = 0.d0
      cmzb = 0.d0
      cmvxb = 0.d0
      cmvyb = 0.d0
      cmvzb = 0.d0
      do i=1,ntot
         etest = 0.5d0*m(i)*(vx(i)**2+vy(i)**2+vz(i)**2)+m(i)*grpot(i)+m(i)*u(i)
         if (etest.gt.0.d0) then 
             id(i) = 4
             cmmunb = cmmunb + m(i)
             cmxunb = cmxunb + x(i)*m(i)
             cmyunb = cmyunb + y(i)*m(i)
             cmzunb = cmzunb + z(i)*m(i)
             cmvxunb = cmvxunb + vx(i)*m(i)
             cmvyunb = cmvyunb + vy(i)*m(i)
             cmvzunb = cmvzunb + vz(i)*m(i)
c             write(74,*)x(i),y(i),z(i),rho(i)*rhounit
         else
             if (cc(i).eq.cc(1)) then 
                id(i) = 1
                if (u(i).eq.0.d0) then
                   write(*,*)'Core point in ',1
                   idmax1 = i
                else
                   if(rhomax1.gt.rho(i)) then
                     idmax1 = i
                     rhomax1 = rho(i)
                   endif
                endif
             else
                id(i) = 2
                if (u(i).eq.0.d0) then
                   write(*,*)'Core point in ',2
                   idmax2 = i
                else
                   if(rhomax2.gt.rho(i)) then
                      idmax2 = i
                      rhomax2 = rho(i)
                   endif
                endif
             endif
             cmmb = cmmb + m(i)
             cmxb = cmxb + x(i)*m(i)
             cmyb = cmyb + y(i)*m(i)
             cmzb = cmzb + z(i)*m(i)
             cmvxb = cmvxb + vx(i)*m(i)
             cmvyb = cmvyb + vy(i)*m(i)
             cmvzb = cmvzb + vz(i)*m(i)        
         endif
      enddo
      if (cmmunb.gt.0.d0) then 
         cmxunb = cmxunb / cmmunb
         cmyunb = cmyunb / cmmunb
         cmzunb = cmzunb / cmmunb
         cmvxunb = cmvxunb / cmmunb
         cmvyunb = cmvyunb / cmmunb
         cmvzunb = cmvzunb / cmmunb
      endif
      if (cmmb.gt.0.d0) then
         cmxb = cmxb / cmmb
         cmyb = cmyb / cmmb
         cmzb = cmzb / cmmb
         cmvxb = cmvxb / cmmb
         cmvyb = cmvyb / cmmb
         cmvzb = cmvzb / cmmb
      endif

      write(*,*) t*tunit,cmmunb,cmxunb,cmyunb,cmzunb,cmvxunb,cmvyunb,cmvzunb, 'unbmat'
      write(*,*) t*tunit,cmmb,cmxb,cmyb,cmzb,cmvxb,cmvyb,cmvzb,
     &           "boundmat"
      sep = sqrt((x(idmax1)-x(idmax2))**2+(y(idmax1)-y(idmax2))**2
     &          +(z(idmax1)-z(idmax2))**2)
      write(70,*)nout,t*tunit,sep
      c1 = m(idmax1)
      c2 = m(idmax2)
      do i=1,ntot
         rcmi = sqrt((x(i)-cmxb)**2+(y(i)-cmyb)**2+(z(i)-cmzb)**2)
c         write(72,*)sqrt(x(i)**2+y(i)**2+z(i)**2),rcmi
         if (rcmi.lt.sep .and. id(i).le.2) then
            if((i.ne.idmax1).and.(i.ne.idmax2)) then 
              r1i = sqrt((x(i)-x(idmax1))**2+(y(i)-y(idmax1))**2+(z(i)-z(idmax1))**2)
              r2i = sqrt((x(i)-x(idmax2))**2+(y(i)-y(idmax2))**2+(z(i)-z(idmax2))**2)
              !write(72,*)x(i),y(i),z(i),rho(i)*rhounit,-m(idmax1)*m(i)/r1i,-m(idmax2)*m(i)/r2i
              epot1i = -m(idmax1)*m(i)/r1i
              epot2i = -m(idmax2)*m(i)/r2i
              if(epot1i.lt.epot2i) then 
                 id(i) = 1
              else
                 id(i) = 2
c                 write(72,*)x(i),y(i),z(i),rho(i)*rhounit,-m(idmax1)*m(i)/r1i,-m(idmax2)*m(i)/r2i
              endif
            endif
         else
            if (id(i).ne.4 .and. (i.ne.idmax1 .and. i.ne.idmax2)) then  
               id(i) = 3
c               write(73,*)x(i),y(i),z(i),rho(i)*rhounit
            endif
         endif
      enddo
      
      m1 = 0.d0
      m2 = 0.d0
      cmmbin = 0.d0
      cmxbin = 0.d0
      cmybin = 0.d0
      cmzbin = 0.d0
      cmvxbin = 0.d0
      cmvybin = 0.d0
      cmvzbin = 0.d0
      Ixx = 0.d0
      Iyy = 0.d0
      Izz = 0.d0
      Ixy = 0.d0
      Iyz = 0.d0
      Ixz = 0.d0
      lx12 = 0.d0
      ly12 = 0.d0
      lz12 = 0.d0
      do i=1,ntot
         if (id(i).eq.1) then 
            m1 = m1 + m(i)
            cmmbin = cmmbin + m(i)
            cmxbin = cmxbin + x(i)*m(i)
            cmybin = cmybin + y(i)*m(i)
            cmzbin = cmzbin + z(i)*m(i)
            cmvxbin = cmvxbin + vx(i)*m(i)
            cmvybin = cmvybin + vy(i)*m(i)
            cmvzbin = cmvzbin + vz(i)*m(i)
c            write(75,*)x(i),y(i),z(i),rho(i)*rhounit
         elseif (id(i).eq.2) then  
            m2 = m2 + m(i)
            cmmbin = cmmbin + m(i)
            cmxbin = cmxbin + x(i)*m(i)
            cmybin = cmybin + y(i)*m(i)
            cmzbin = cmzbin + z(i)*m(i)
            cmvxbin = cmvxbin + vx(i)*m(i)
            cmvybin = cmvybin + vy(i)*m(i)
            cmvzbin = cmvzbin + vz(i)*m(i)
c            write(76,*)x(i),y(i),z(i),rho(i)*rhounit
         endif
         if (id(i).le.2) then 
            Ixx = Ixx + m(i)*(y(i)**2+z(i)**2)
            Iyy = Iyy + m(i)*(x(i)**2+z(i)**2)
            Izz = Izz + m(i)*(x(i)**2+y(i)**2)
            Ixy = Ixy - m(i)*x(i)*y(i)
            Iyz = Iyz - m(i)*z(i)*y(i)
            Ixz = Ixz - m(i)*x(i)*z(i)
            lx12 = lx12 + m(i)*(y(i)*vz(i)-z(i)*vy(i))
            ly12 = ly12 + m(i)*(z(i)*vx(i)-x(i)*vz(i))
            lz12 = lz12 + m(i)*(x(i)*vy(i)-y(i)*vx(i))
         endif
      enddo
      if (cmmbin.ne.0.d0) then 
          cmxbin = cmxbin / cmmbin
          cmybin = cmybin / cmmbin
          cmzbin = cmzbin / cmmbin
          cmvxbin = cmvxbin / cmmbin
          cmvybin = cmvybin / cmmbin
          cmvzbin = cmvzbin / cmmbin
      endif
      x1 = x(idmax1)
      x2 = x(idmax2)
      y1 = y(idmax1)
      y2 = y(idmax2)
      z1 = z(idmax1)
      z2 = z(idmax2)  
      vx1 = vx(idmax1)
      vx2 = vx(idmax2)    
      vy1 = vy(idmax1)
      vy2 = vy(idmax2)
      vz1 = vz(idmax1)
      vz2 = vz(idmax2)

      write(*,*)'I tensor',Ixx,Iyy,Izz,Ixy,Ixz,Iyz
      write(*,*)'L vec',lx12,ly12,lz12
      call omega_vec(Ixx,Iyy,Izz,Ixy,Ixz,Iyz,lx12,ly12,lz12,omegax,omegay,omegaz)
      write(*,*)'Omega (x,y,z)',omegax/tunit,omegay/tunit,omegaz/tunit
      write(*,*)t*tunit,m1,m2,cmmbin,cmxbin,cmybin,cmzbin,
     &          cmvxbin,cmvybin,cmvzbin

      call phi(phi1,xL1,m1,m2,sep)
      write(*,*)'m1,m2,sep',m1,m2,sep 
      write(*,*)'Potential,xL1: ',phi1,xL1
      if ((y2.gt.0.d0).and.(x2.gt.0.d0)) then 
          theta = atan(y2/x2)
      elseif ((y2.gt.0.d0).and.(x2.lt.0.d0)) then
          theta = pi + atan(y2/x2)
      elseif ((y2.lt.0.d0).and.(x2.lt.0.d0)) then
          theta = pi + atan(y2/x2)
      elseif ((y2.lt.0.d0).and.(x2.gt.0.d0)) then
          theta = 2.d0*pi + atan(y2/x2)
      endif
      xL1p = xL1*cos(theta) + cmxbin
      yL1p = xL1*sin(theta) + cmybin
      zL1p = cmzbin
      write(*,*)'xL1 prime (x,y)',xL1p,yL1p,theta
      write(*,*)'x1 (x,y,z)(vx,vy,vz)',x1,y1,z1,vx1,vy1,vz1
      write(*,*)'x2 (x,y,z)(vx,vy,vz)',x2,y2,z2,vx2,vy2,vz2
      r1 = sqrt((x1-xL1p)**2+(y1-yL1p)**2+(z1-zL1p)**2)
      r2 = sqrt((x2-xL1p)**2+(y2-yL1p)**2+(z2-zL1p)**2)
      write(*,*)'Spheres 1,2',r1,r2

CC Recalculate I, L, omega, and centre of mass for the particles in 
C the sphere r1 and r2 which gives a better approximation of the
C particles around m1 and m2
      m1 = 0.d0
      m2 = 0.d0
      cmmbin = 0.d0
      cmxbin = 0.d0
      cmybin = 0.d0
      cmzbin = 0.d0
      cmvxbin = 0.d0
      cmvybin = 0.d0
      cmvzbin = 0.d0
      Ixx = 0.d0
      Iyy = 0.d0
      Izz = 0.d0
      Ixy = 0.d0
      Iyz = 0.d0
      Ixz = 0.d0
      lx12 = 0.d0
      ly12 = 0.d0
      lz12 = 0.d0
      do i=1,ntot
       if(i.ne.idmax1 .and. i.ne.idmax2) then 
         r1i = sqrt((x(i)-x1)**2+(y(i)-y1)**2+(z(i)-z1)**2)
         r2i = sqrt((x(i)-x2)**2+(y(i)-y2)**2+(z(i)-z2)**2)     
         etest = 0.5d0*m(i)*(vx(i)**2+vy(i)**2+vz(i)**2)+m(i)*grpot(i)+m(i)*u(i)
         if (r1i.lt.r1 .and. etest.lt.0.d0) then 
             id(i)=1
         elseif(r2i.lt.r2 .and. etest.lt.0.d0) then 
             id(i)=2
         elseif(etest .lt. 0.d0) then 
             id(i)=3
         else
             id(i)=4
         endif
       endif
       if (id(i).eq.1) then
            m1 = m1 + m(i)
            cmmbin = cmmbin + m(i)
            cmxbin = cmxbin + x(i)*m(i)
            cmybin = cmybin + y(i)*m(i)
            cmzbin = cmzbin + z(i)*m(i)
            cmvxbin = cmvxbin + vx(i)*m(i)
            cmvybin = cmvybin + vy(i)*m(i)
            cmvzbin = cmvzbin + vz(i)*m(i)
         elseif (id(i).eq.2) then
            m2 = m2 + m(i)
            cmmbin = cmmbin + m(i)
            cmxbin = cmxbin + x(i)*m(i)
            cmybin = cmybin + y(i)*m(i)
            cmzbin = cmzbin + z(i)*m(i)
            cmvxbin = cmvxbin + vx(i)*m(i)
            cmvybin = cmvybin + vy(i)*m(i)
            cmvzbin = cmvzbin + vz(i)*m(i)
       endif
      enddo
      
      if (cmmbin.ne.0.d0) then
          cmxbin = cmxbin / cmmbin
          cmybin = cmybin / cmmbin
          cmzbin = cmzbin / cmmbin
          cmvxbin = cmvxbin / cmmbin
          cmvybin = cmvybin / cmmbin
          cmvzbin = cmvzbin / cmmbin
      endif
      count1=0
      count2=0  
      do i=1,ntot
         r1i = sqrt( (x(i)-cmxbin)**2+(y(i)-cmybin)**2+(z(i)-cmzbin)**2)
         if(r1i .le. 10.0d0) then
             count1 = count1 + 1
         endif
         if (id(i).le.2) then
            Ixx = Ixx + m(i)*((y(i)-cmybin)**2+(z(i)-cmzbin)**2)
            Iyy = Iyy + m(i)*((x(i)-cmxbin)**2+(z(i)-cmzbin)**2)
            Izz = Izz + m(i)*((x(i)-cmxbin)**2+(y(i)-cmybin)**2)
            Ixy = Ixy - m(i)*(x(i)-cmzbin)*(y(i)-cmybin)
            Iyz = Iyz - m(i)*(z(i)-cmzbin)*(y(i)-cmybin)
            Ixz = Ixz - m(i)*(x(i)-cmxbin)*(z(i)-cmzbin)
            lx12 = lx12 + m(i)*((y(i)-cmybin)*(vz(i)-cmvzbin)-(z(i)-cmzbin)*(vy(i)-cmvybin))
            ly12 = ly12 + m(i)*((z(i)-cmzbin)*(vx(i)-cmvxbin)-(x(i)-cmxbin)*(vz(i)-cmvzbin))
            lz12 = lz12 + m(i)*((x(i)-cmxbin)*(vy(i)-cmvybin)-(y(i)-cmybin)*(vx(i)-cmvxbin))
            !write(75,*)x(i)-cmxbin,y(i)-cmybin,z(i)-cmzbin,rho(i)*rhounit
         endif
         if(id(i).eq.3) then
           !write(76,*)x(i)-cmxbin,y(i)-cmybin,z(i)-cmzbin,rho(i)*rhounit
           if(r1i .gt. 10.0d0) then
             count2 = count2 + 1
           endif
         endif
      enddo
      write(*,*) 'Particles within 10Rsun ',count1
      write(*,*) 'Particles more than 10Rsun ',count2
      write(*,*) 'Ekin (2points) ', (0.5*m(1)*(vx(1)**2+vy(1)**2+vz(1)**2)+
     &                             0.5*m(ntot)*(vx(ntot)**2+vy(ntot)**2+vz(ntot)**2))*eunit
      write(*,*) 'Epot (2points) ', (0.5*m(1)*grpot(1)+0.5*m(ntot)*grpot(ntot))*eunit
      write(*,*) 'Eint (2points) ', (m(1)*u(1)+m(ntot)*u(ntot))**eunit


      write(*,*)'I tensor',Ixx,Iyy,Izz,Ixy,Ixz,Iyz
      write(*,*)'L vec',lx12,ly12,lz12
      call omega_vec(Ixx,Iyy,Izz,Ixy,Ixz,Iyz,lx12,ly12,lz12,omegax,omegay,omegaz)
      write(*,*)'Omega (x,y,z)',omegax/tunit,omegay/tunit,omegaz/tunit
      write(*,*)t*tunit,m1,m2,cmmbin,cmxbin,cmybin,cmzbin,
     &          cmvxbin*vunit,cmvybin*vunit,cmvzbin*vunit,' bin_evol'
      write(*,*)t*tunit,2.d0*pi/(sqrt(omegax**2+omegay**2+omegaz**2))*tunit,
     &       'Period'
      write(*,*)t*tunit,m1,x1-cmxbin,y1-cmybin,z1-cmzbin,
     &     (vx1-cmvxbin)*vunit,(vy1-cmvybin)*vunit,(vz1-cmvzbin)*vunit,
     &       ' star_1'
      write(*,*)t*tunit,m2,x2-cmxbin,y2-cmybin,z2-cmzbin,
     &     (vx2-cmvxbin)*vunit,(vy2-cmvybin)*vunit,(vz2-cmvzbin)*vunit,
     &       ' star_2'
      write(*,*)t*tunit,sep,m1-m(1),m2-m(ntot),cmmbin-m(1)-m(ntot),' N'

      end
