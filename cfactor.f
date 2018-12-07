
!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * cfactor.f *                                   galprop package * 4/14/2000 
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|


      REAL function EMISS(r)
c***********************************************************************
c calculation of emissivity (EMISS) of the Galaxy as a function of the 
c distance from its center (r, kpc)
c***********************************************************************
       !parameter(z=0)
       REAL r,z
c       REAL isrf_energy_density ! added 10 Jan 2000
   
c      EMISS = exp(-r/3.0)
c      EMISS = 1.d0           ! for test1
c      EMISS = exp(-r/0.2)    ! for test2
    
      z=0.0
      call isrf_energy_density(r,z,emiss)
c      print*," EMISS #1 ->", r,emiss

      return
      end

      real*8 function CFactor(kbg,E0,E,PLindex,RG,rho,xi,z,RS)
c***********************************************************************
c                            *** I.Moskalenko, version of 31.03.1998 ***
c calculation of the correction factor for
c INVERSE COMPTON scattering in an ANISOTROPIC photon field
c    INPUT:
c kbg =-1 for isotropic background photon field (constant); 
c     = 0 for opaque disk photon field (~cos(theta));
c     = 1 for transparent disk photon field (~1/cos(theta)); 
c E0 is the energy of a background photon (in electron rest mass units mc^2);
c E  is the energy of scattered photon (in electron rest mass units mc^2);
c PLindex is the power-law index of the electron spectrum (>0);
c RG is the Galactic radius (kpc);
c (rho,xi,z) are the cylindrical coordinates of the electron position
c   in respect to the Galactic center (rho,z in kpc, 0< xi (radians) <2Pi 
c   is polar angle between the electron position and observer direction);
c RS is the Solar distance from the Galactic center (kpc).
c    internal OUTPUT:
c SPEC - the differential energy spectrum of gamma-rays (phot/ccm/sec/energy)
c    per one electron in ccm divided by Pi*r0^2*c;
c DENS - number density of background photons at the electron position
c    (the emissivity of unit area of the Galactic plane is equal to 1).
c    REFERENCES:
c I.V.Moskalenko & A.W.Strong 2000, ApJ 528, 357
c***********************************************************************
      implicit real*8 (a-h,o-z)
      real EPS1

      CFactor = 0.
      EPS1 = 5.e-3        ! relative accuracy in integration over gamma
      ig = 8
      dg = 1.2**(1./ig)
      Gmax = 1.d13
      gam = 1.001 *(E+dsqrt(E/E0+E**2))/2.d0
      NG = dlog(Gmax/gam)/dlog(dg) +1
      S = 0.d0
      SFJ = 0.d0
      if (5.*gam .gt. Gmax) return
c## photon number density calculation ##
      call AIC(0,kbg,E0,E,gam,RG,rho,0.d0,z,RS, SPEC,DENS)
c## integration over the Lorentz factor ##
      do n=1,NG
         ddg = dlog(dg)
         call AIC(1,kbg,E0,E,gam,RG,rho,xi,z,RS, SPEC,DENS)
         S1   = SPEC            /gam**PLindex *gam*ddg  !!!!
         S = S +S1
         SFJ1 = FJONES(gam,E0,E)/gam**PLindex *gam*ddg  !!!!
         SFJ = SFJ + SFJ1
         if(dabs(S1/S) .lt. EPS1 .and. dabs(SFJ1/SFJ) .lt. EPS1) goto 10
         if(n/ig*ig .eq. n .and. dg. lt. 1.15) dg = dg**2
         gam = gam*dg
      enddo
   10 CFactor = S/SFJ
      return
      end



      subroutine AIC(key1,kbg1,E0,E,gamma,RG,rho1,xi,z1,RS, SPEC,DENS1)
c***********************************************************************
c                            *** I.Moskalenko, version of 31.03.1998 ***
c INVERSE COMPTON scattering in an ANISOTROPIC photon field
c    INPUT:
c key1 = 0 for call first (photon number density calculation); 
c      = 1 for subsequent calls (if rho1,z1,RS,kbg are not changed);
c kbg1 =-1 for isotropic background photon field (constant); 
c      = 0 for opaque disk photon field (~cos(theta));
c      = 1 for transparent disk photon field (~1/cos(theta)); 
c E0 is the energy of a background photon (in electron rest mass units mc^2);
c E  is the energy of scattered photon (in electron rest mass units mc^2);
c gamma is the electron Lorentz factor;
c RG is the Galactic radius (kpc);
c (rho1,xi,z1) are the cylindrical coordinates of the electron position
c   in respect to the Galactic center (rho1,z1 in kpc, 0< xi <2Pi - 
c   polar angle between the electron position and observer direction);
c RS is the Solar distance from the Galactic center (kpc).
c    OUTPUT:
c SPEC - the differential energy spectrum of gamma-rays (phot/ccm/sec/energy)
c    per one electron in ccm;
c DENS - number density of background photons in the electron position
c    (the emissivity of unit area of the Galactic plane is equal to 1).
c    REFERENCES:
c I.V.Moskalenko & A.W.Strong 2000, ApJ 528, 357
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common/xyz/R,rho,z,akappa,psi,costheta /gamma/E1,E4,gam,A2,B2
     #      /ZETA/DENS,key,kbg
      external FTHETA
!Gulli20070821
!$omp threadprivate(/xyz/,/gamma/,/ZETA/)

      EPS = 1.d-2
      E1 = E0
      E4 = E
      gam = gamma
      R = RG
      rho = rho1
      z = z1
      key = key1
      kbg = kbg1

      SPEC = 0.d0
      DENS1=DENS
      Pi=3.141592653589793d0

c## limits on cos(theta) from the geometry ## eqs.(19),(20)
      c1 = z/dsqrt(z*z+(rho+R)**2)
      c2 = z/dsqrt(z*z+(rho-R)**2)
      if(rho .lt. R)    c2 = 1.d0
      if(c2  .gt. 1.d0) c2 = 1.d0 !IMOS20060420

c## calculation of the photon density (essentially normalization) ##
      if(key .eq. 0) then
c## absolute accuracy estimation ##
         AEPS = 0.1*EPS*FTHETA((c2+c1)/2.d0)*(c2-c1)
         call SIM1(c1,c2,(c2-c1)/30.,EPS,AEPS,FTHETA,DENS) !## eq.(8)
         DENS1=DENS
c#         beta = 1.d0-1.d0/(2.d0*gam**2)
c#         if(gam .lt. 1.d5) beta = dsqrt(1.d0-1.d0/gam**2)
c#         A1 = E1*gam*(1.d0+beta*CZmin)
c#         B1 = E1*2.d0*gam   !commented IMOS20060405
         return
      endif
c## check if E4 lower then maximal energy ##
      if(E4/gam .ge. 4.d0*E1*gam/(1.d0+4.d0*E1*gam)) return

c#      B2 = B1*E4/(E1**2*(gam-E4))-1.d0  !commented IMOS20060405
c#      A2 = A1*E4/(E1**2*(gam-E4))-1.d0
      A2 = E4/E1/gam/(gam-E4)/2.d0-1.d0   ! zeta0 calculation eq.(13)

      a = dsqrt((RS-rho*dcos(xi))**2 +(rho*dsin(xi))**2) ! observer distance
      akappa = Pi-datan(a/z)                             ! obs.'s polar angle 
      b = 1.d0
      if(a .ne. 0.d0) b = (rho-RS*dcos(xi))/a   ! cosine angle between a and rho
      if(dabs(b) .gt. 1.d0) b = dsign(1.d0,b)
      psi = Pi-dacos(b)                         ! angle between a and rho
      if(dsin(xi) .gt. 0.d0) psi = Pi+dacos(b)

c## limits on theta from kinematics ## eq.(15)
      theta1 = Pi-akappa+dacos(A2)
      theta2 = Pi-akappa-dacos(A2)
      if(theta1 .gt. Pi)   theta1 = Pi
      if(theta2 .lt. 0.d0) theta2 = 0.d0

c## choose the intersection ##
      c1 = dmax1(c1,dcos(theta1))
      c2 = dmin1(c2,dcos(theta2))
      if(c1 .ge. c2) return

c## absolute accuracy estimation ##
      AEPS = 0.1*EPS*FTHETA((c2+c1)/2.d0)*(c2-c1)
c## integration over cos(theta) ##
      call SIM1((c2+c1)/2.d0,c1,(c2-c1)/50.,EPS,AEPS,FTHETA,SPEC1)
      call SIM1((c2+c1)/2.d0,c2,(c2-c1)/50.,EPS,AEPS,FTHETA,SPEC2)
      SPEC = (SPEC2-SPEC1)/E1/(gam-E4)**2  /DENS !## eq.(8)

      return
      end

      real*8 function FTHETA(costheta1)
c***********************************************************************
c calculation of the integral over fi
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common/xyz/R,rho,z,akappa,psi,costheta /gamma/E1,E4,gam,A2,B2
     #      /ZETA/DENS,key,kbg
      external FFI
!Gulli20070821
!$omp threadprivate(/xyz/,/gamma/,/ZETA/)

      EPS = 5.d-3
      FTHETA = 0.d0
      Pi=3.141592653589793d0
      costheta = costheta1
      sintheta = dsqrt(1.d0-costheta**2)
      tantheta=0.
      if(costheta .ne. 0.d0) tantheta = sintheta/costheta !IMOS20060420
      RRzz = R*R-rho**2-z*z*tantheta**2
      cosfi = 1.d0
c## limits on cos(fi) from the geometry ## eqs.(20)-(22)
      if(rho*z*tantheta .ne. 0.d0) then
         cosfi = RRzz/(2.d0*rho*z*tantheta)
         if(dabs(cosfi) .gt. 1.d0) cosfi = dsign(1.d0,cosfi)
      else
         cosfi = dsign(1.d0,RRzz)
      endif
      if(cosfi .eq. -1.d0) return
      fi1 = dacos(cosfi)
      fi2 = 2.d0*Pi-dacos(cosfi)

c## weight function ##
      WEIGHT = 1.d0                            ! "isotropic" disk
      if(kbg. eq. 0) WEIGHT = costheta         ! opaque disk
      if(kbg. eq. 1) then                      ! transparent disk
         WEIGHT = 1.d8                                 !IMOS20060420
         if(costheta .ne. 0.d0) WEIGHT = 1.d0/costheta !IMOS20060420
      endif
c## calculation of the photon density ##
      if(key .eq. 0) then
c## absolute accuracy estimation ##
         AEPS = 0.1*EPS*(fi2-fi1)*FFI((fi2+fi1)/2.d0)
         call SIM2(fi2,fi1,(fi2-fi1)/30.,EPS,AEPS,FFI,FTHETA)
         FTHETA = -FTHETA *WEIGHT
         return
      endif

c## limits on cos(fi) from kinematics ##
      sinkappa = dsin(akappa)
      if(sinkappa .lt. 0.d0) sinkappa = 0.d0
      if(sinkappa*sintheta .ne. 0.d0) then
         cosfipsi=(A2+dcos(akappa)*costheta)/(sinkappa*sintheta)
         if(dabs(cosfipsi) .gt. 1.d0) cosfipsi = dsign(1.d0,cosfipsi)
      else
         cosfipsi = dsign(1.d0,A2+dcos(akappa)*costheta)
      endif
      if(cosfipsi .eq. 1.d0) return
      fi11 = psi-dacos(cosfipsi)
      fi12 = psi+dacos(cosfipsi)
c      write(*,*) "I",fi1,fi2,fi11,fi12

      cf1  = dcos(fi1)
      cf11 = dcos(fi11)
      cf12 = dcos(fi12)
      sf11 = dsin(fi11)
      sf12 = dsin(fi12)
      a = 0.d0
      b = 0.d0
      c = 0.d0
      d = 0.d0
c## check if there is the intersection ##
      if((cf1-cf11)*(cf1-cf12) .gt. 0.d0) then
         if(cf11 .lt. cf1) then
            a = fi11
            b = fi12
            if( (sf11 .lt. 0.d0 .and. sf12 .gt. 0.d0)
     #         .or. (fi12-fi11 .gt. Pi .and. sf11*sf12 .gt. 0.d0)
     #        ) then
               a = fi1
               b = fi12
               c = fi11
               d = fi2
            endif
         else
            a = fi1
            b = fi2
            if( (fi12-fi11 .lt. Pi .and. sf11*sf12 .gt. 0.d0)
     #         .or. (sf11*sf12 .lt. 0.d0 .and. sf11 .lt. 0.d0)
     #        ) return
         endif
      else
         if(cf11 .gt. cf1) then
            a = fi1
            b = fi12
         else
            a = fi11
            b = fi2
         endif
      endif
c      write(*,*) "II",fi1,fi2,fi11,fi12,a,b,c,d

c## 1st intersection; reduce angles to a standard range ##
      if(dsin(a) .ge. 0.d0) then 
         a = dacos(dcos(a))
      else
         a = 2.d0*Pi-dacos(dcos(a))
      endif
      if(dsin(b) .ge. 0.d0) then 
         b = dacos(dcos(b))
      else
         b = 2.d0*Pi-dacos(dcos(b))
      endif
      if(b .lt. a) b = b+2.d0*Pi
c## absolute accuracy estimation ##
      AEPS = 0.1*EPS*(b-a)*FFI((b+a)/2.d0)
c## integration over fi over the 1st intersection ##
      call SIM2(b,a,(b-a)/30.,EPS,AEPS,FFI,FTHETA)

c## 2st intersection; reduce angles to a standard range ##
      if(dsin(c) .ge. 0.d0) then 
         c = dacos(dcos(c))
      else
         c = 2.d0*Pi-dacos(dcos(c))
      endif
      if(dsin(d) .ge. 0.d0) then 
         d = dacos(dcos(d))
      else
         d = 2.d0*Pi-dacos(dcos(d))
      endif
      if(d .lt. c) d = d+2.d0*Pi
c      write(*,*) "III",a,b,c,d
      FTHETA1 = 0.d0
c## integration over fi over the 2nd intersection ##
      if(c .ne. d) then
         AEPS = 0.1*EPS*(d-c)*FFI((c+d)/2.d0)
         call SIM2(d,c,(d-c)/30.,EPS,AEPS,FFI,FTHETA1)
      endif
      FTHETA = -(FTHETA+FTHETA1) *WEIGHT
      return
      end

      real*8 function FFI(fi)
c***********************************************************************
c calculation of the integrand 
c***********************************************************************
      implicit real*8 (a-h,o-z)
      real RRR,EMISS
      common/xyz/R,rho,z,akappa,psi,costheta /gamma/E1,E4,gam,A2,B2
     #      /ZETA/DENS,key,kbg
!Gulli20070821
!$omp threadprivate(/xyz/,/gamma/,/ZETA/)

      sintheta=dsqrt(1.d0-costheta**2)
      coszeta =-dcos(akappa)*costheta+dsin(akappa)*sintheta*dcos(fi-psi)
c#      if(coszeta .lt. A2-1.d-7) return

      ztantheta = 0.
      if(costheta .ne. 0.) ztantheta = z*sintheta/costheta !IMOS20060420
c## radius via angles  ## eq.(24)
      RRR = dsqrt(rho**2+ztantheta**2+2.d0*rho*ztantheta*dcos(fi))

c      print*," FFI #1 ->", RRR
c## should be the integration over all -Z_h < z1 < z ##
      FFI = EMISS(RRR)                            ! emissivity calc.
c      print*," FFI #2 ->", FFI

c## calculation of the photon density ##
      if(key .eq. 0) return

c## calculation of the integrand ## eq.(8)
      E2 = E1*gam*(1.d0+coszeta)
      FFI = FFI *( 2.d0-2.d0*E4/gam*(1.d0/E2+2.d0)
     #   +(E4/gam)**2*(1.d0/E2**2+2.d0/E2+3.d0) -(E4/gam)**3 )
c      print*," FFI #3 ->", FFI
      return
      end

      real*8 function FJONES(gam,E1,E4)
c***********************************************************************
c                            *** I.Moskalenko, version of 26.03.1998 ***
c calculation of the spectrum from INVERSE COMPTON scattering
c of ISOTROPIC background photons off isotropic electrons
c    INPUT:
c gam is the electron Lorentz factor;
c E1 is the energy of a background photon (in electron rest mass units mc^2);
c E4 is the energy of scattered photon (in electron rest mass units mc^2);
c    OUTPUT:
c FJONES - the differential energy spectrum of gamma-rays
c    (phot ccm/sec/energy) per one electron in ccm divided by Pi*r0^2*c;
c    REFERENCES:
c F.C.Jones 1968, Phys.Rev. 167, 1159
c***********************************************************************
      implicit real*8 (a-h,o-z)

      FJONES = 0.d0
      if(E4/gam .gt. 4.d0*E1*gam/(1.d0+4.d0*E1*gam)) return
      q = E4/(4.d0*E1*gam**2)*(1.d0+4.d0*E1*gam)
      if(4.d0*E1*gam .lt. 1.d10) q = E4/(4.d0*E1*gam*(gam-E4))
      FJONES = 2.d0*q*dlog(q)+(1.d0+2.d0*q)*(1.d0-q)
     #   +(1.d0-q)/2.d0*(4.d0*E1*gam*q)**2/(1.d0+4.d0*E1*gam*q)
      FJONES = FJONES*2.d0/E1/gam**2  ! in units Pi*r0^2*c  eq.(12)

      return
      end

 
