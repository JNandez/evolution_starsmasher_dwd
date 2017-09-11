      implicit none
c Define the maximum number of particles
      integer nmax
      parameter(nmax=300000)
c Definition of variables
      integer ntot,nnopt,nout,nit,nav,ngr,nrelax,cc(nmax)
      real*8 hmin,hmax,sep0,tf,dtout,t,alpha,beta,
     &       tskipahead,trelax,dt,omega2
      real*8 x(nmax),y(nmax),z(nmax),m(nmax),h(nmax),rho(nmax),vx(nmax),vy(nmax),vz(nmax),
     &       vxdot(nmax),vydot(nmax),vzdot(nmax),u(nmax),udot(nmax),grpot(nmax),
     &       mu(nmax),divv(nmax) 
      integer id(nmax)
      real*8 cmxb,cmyb,cmzb,cmvxb,cmvyb,cmvzb,cmmb
      real*8 cmxunb,cmyunb,cmzunb,cmvxunb,cmvyunb,cmvzunb,cmmunb

c Common variables
      common/nconst/ntot,nnopt,nout,nit,nav,ngr,nrelax
      common/init_real/hmin,hmax,sep0,tf,dtout,t,alpha,beta,
     &       tskipahead,trelax,dt,omega2      
      common/dynamics/m,x,y,z,vx,vy,vz,vxdot,vydot,vzdot,h,grpot,divv
      common/thermodynamics/rho,u,udot,mu,cc
      common/id_tracks/id
      common/bound_ma/cmxb,cmyb,cmzb,cmvxb,cmvyb,cmvzb,cmmb
      common/unb_ma/cmxunb,cmyunb,cmzunb,cmvxunb,cmvyunb,cmvzunb,cmmunb

c Constants (note that all constant are in cgs units)
      real*8 pi,mp,G,Msun,Rsun,eunit,tunit
      real*8 tunits,lunit,rhounit,vunit
      parameter(pi=3.1415926535897932384626d0)
      parameter(mp=1.67262158d-24,G=6.67390d-8)
      parameter(Msun=1.9891d33,Rsun=6.9599d10)
      parameter(eunit=G*Msun**2/Rsun*1.0d-46) 
      parameter(tunits=sqrt(Rsun**3/G/Msun))
      parameter(tunit=tunits/60.d0/60.d0/24.d0) !in days
      parameter(lunit=Msun*Rsun**2/tunits*1.0d-51)
      parameter(rhounit=Msun/Rsun**3)
      parameter(vunit=sqrt(G*Msun/Rsun)*1e-5)
