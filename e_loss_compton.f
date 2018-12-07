
!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * e_loss_compton.f *                            galprop package * 4/14/2000 
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

      function e_loss_compton(w,gam)
c***********************************************************************
c                            *** I.Moskalenko, version of 25.02.1998 ***
c ELECTRON/POSITRON energy losses due to the inverse Compton, dE/dt. 
c    INPUT: 
c w - photon energy (in mc^2 units);
c gam - electron/positron Lorentz factor;
c The normalization of photon distrib. should be Int{dw f(w)}=number density.
c    OUTPUT: 
c FCOMPT the energy losses, normalization [dE/dt]/Ne [MeV cm^3/c] 
c    REFERENCES:
c Jones, F.C. 1965, Phys. Rev. 137, B1306
c Moskalenko, I.V., Jourdain E. 1997, A&A 325, 401
c Strong A.W., Moskalenko I.V. 1998, ApJ 509, 212
c***********************************************************************
      implicit real*8 (A-H,M,O-Z), integer (i-l,n)
      FC(A)=gam*((A+6.d0+3.d0/A)*dlog(2.d0*A+1.d0)-(22.d0/3.d0*A**3
     +      +24.d0*A**2+18.d0*A+4.d0)/(2.d0*A+1.d0)**2)
     -      -w*((A+31.d0/6.d0+5.d0/A+3.d0/2.d0/A**2)*dlog(2.d0*A+1.d0)
     -      -(22.d0/3.d0*A**3+28.d0*A**2+103.d0/3.d0*A+17.d0+3.d0/A)
     /      /(2.d0*A+1.d0)**2)
      POW2(w)=                        !! expansion series !!
     +  +56.d0/25.d0*w**5*gambet*gam*(7.d0*gambet**4                    ! w5
     +     +5.d0*gam**2*(-5.d0+7.d0*gam**2)
     +     +5.d0*gambet**2*(-5.d0+14.d0*gam**2))
     -  -16.d0/45.d0*w**4*gambet*gam**2*(-48.d0 +63.d0*gam**2           ! w4
     +     +bet2*(-16.d0+63.d0*gam**2))
     +  +16.d0/9.d0*w**3*gambet*gam*(-3.d0+(3.d0+bet2)*gam**2)          ! w3

      Pi = 3.14159265359d0
      R02 = 2.81794d-13**2 ! classical electron radius square (cm^2)
      mc2 = 0.511   ! MeV, electron rest mass
      cc = 2.99792d10   ! cm/s, light speed

      gambet = dsqrt(gam**2-1.d0)
      bet2 = 1.d0-1.d0/gam**2
c normalization of the gamma-ray spectra is int{dw f(w)*w**2} = Ngam,
c number density.
      FCOMPT = 1.d0/w/w
   
c $v$ calcs of Int_A^B {dx ln(2x+1)/x} = Li2(-2A)-Li2(-2B)
      A = w*(gam-gambet)
      if(gam .gt. 1.d5) A = w/2.d0/gam
      B = w*(gam+gambet)
      C = 0.2d0
c	     A = 0.d0
c	     B = 0.3d0
      if((A-C)*(B-C).gt.0.d0) then
         if(B.lt.C) C = B
         if(A.gt.C) C = A
      endif
      sum = 0.d0
      if(A.lt.C) then
         j = -13.*alog(10.)/dlog(2.*C)+3
         up   = C*2.d0/(j+1)**2
         down = A*2.d0/(j+1)**2
         do k = j,1,-1
            up   = 2.d0*C*(1.d0/k**2-up)
            down = 2.d0*A*(1.d0/k**2-down)
         enddo
         sum = up-down
      endif
      if(C.lt.B) then
         up   = 2.d0*B+1.d0
         down = 2.d0*C+1.d0
         j = 15.*alog(10.)/dlog(down)+3
         sum1 = 1.d0/(j+1)**2/up
         z    = 1.d0/(j+1)**2/down
         do k = j,1,-1
            sum1 = (1.d0/k**2+sum1)/up
	    z	 = (1.d0/k**2+z)/down
         enddo
         up   = dlog(2.d0*B+1.d0)
         down = dlog(2.d0*C+1.d0)
	 sum=sum+up*(dlog(2.d0*B)-up/2.d0)-sum1
     -	     -(down*(dlog(2.d0*C)-down/2.d0)-z)
      endif
c #!attention! u(b)=1.6449340668482+y*(dlog(2.d0*b)-y/2.d0)-v; a=0, b>0.2#
c $^$
c $$      call SIM1(A,B,(B-A)/1.d3,1.d-5,1.d-18,FL,sum)
      FBA = FC(B)-FC(A)-(2.d0*gam-w)*sum
      if(gam*(w*gam)**3 .le. 4.d-5) FBA = POW2(w)
      e_loss_compton = FCOMPT*FBA  *Pi*R02*mc2*cc    /2.d0/gam/gambet
      return
      end

