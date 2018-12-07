
!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * synchrotron.f *                               galprop package * 4/14/2000 
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

          real*8 function SYNCHROTRON(gamma,Esynch,B)
c**********************************************************************
c         *** I.Moskalenko (MPE,Garching) *** version of 12.11.1997 ***
c PURPOSE:
c Calculates the synchrotron emissivity of an electron integrated
c over the emission angle for an isotropic pitch angle distribution.
c UNITS: erg/s/Hz;
c        for units MeV/s/Hz remove a comment "c[0]";
c        for transformation to units 1/MeV/s remove comments "c[0]","c[2]".
c REFERENCES:
c [1] Ginzburg 1979 ,"Theor. Phys. & Astrophys."
c [2] Ghisellini, Guilbert, Svensson 1988, ApJ 334, L5 (eq.[9])
c [3] Strong A.W., Moskalenko I.V., Reimer O. 2000, ApJ 537, 763
c INPUT parameters:
c    gamma   is the electron Lorentz factor;
c    Esynch (Hz)   is the frequency of synchrotron photons;
c       if you prefer that Esynch to be the energy (MeV) of synchro-photons,
c       then remove a comment "c[1]";
c    B (gauss=g^1/2 cm^-1/2 s^-1)   is the (total) magnetic field strength.
c**********************************************************************
      implicit real*8 (a-h,m-z)
      Pi = 3.14159265d0
      MeV_to_erg = 1.60217733d-6
      h  = 6.6260755d-27 /MeV_to_erg                ! MeV*s; Planck constant
      emc = 4.80321d-10/9.10939d-28/2.99792d10      ! (cm/g)^1/2; emc =e/mc
      e2c = 2.81794d-13*0.511*MeV_to_erg/2.99792d10 ! erg*s; =e^2/c =r0*mc^2/c
c[0]      e2c = 2.81794d-13*0.511/2.99792d10   ! MeV*s; =e^2/c= r0*mc^2/c
      nuB = emc*B/2.d0/Pi                           ! Hz; frequency

      x=Esynch/nuB/gamma**2/3.d0           ! if [Esynch] = Hz
c[1]      x=(Esynch/h)/nuB/gamma**2/3.d0       ! if [Esynch] = MeV
      BSK1 = BSKR3(x,1)                    ! Bessel K{1/3}
      BSK4 = BSKR3(x,-2)+2.d0/3.d0/x*BSK1  ! Bessel K{4/3}
      SYNCHROTRON=x*x*(BSK4*BSK1-3.d0/5.d0*x*(BSK4+BSK1)*(BSK4-BSK1))
     2   *4.d0*Pi*dsqrt(3.d0)*e2c*nuB  !<- units erg/s/Hz (or MeV/s/Hz if [0])
c[2]     3   /Esynch/h                     !<- transformation to units 1/MeV/s
      return
      end
                                       
      
      real*8 FUNCTION BSIKR3(X,NU)
c**********************************************************************
c Entries BSIR3(X,NU) and BSKR3 calculate the modified Bessel functions
c   I{nu/3}(X) and K{nu/3}(X) for real arguments X>0 and NU=-2,-1,1,2.
c Entries EBSIR3(X,NU) and EBSKR3 calculate exp(-X)*I{nu/3}(X) and
c   exp(X)*K{nu/3}(X), correspondingly.
c The value x=0 is permitted for the functions I if NU>0.
c The functions K are even with respect to NU.       ** CERNLIB C340 **
c**********************************************************************
      implicit real*8 (a-h,o-z)
      LOGICAL L,E,NUL,NOE

      ENTRY EBSIR3(X,NU)

      E=.TRUE.
      GO TO 7

      ENTRY BSIR3(X,NU)

      E=.FALSE.
    7 L=.TRUE.
      N=IABS(NU)
      IF(N .NE. 1 .AND. N .NE. 2) GO TO 10
      IF(X .LE. 0.0d0) GO TO 11
      NUL=NU .LT. 0
      IF(X .GT. 8.0d0) GO TO 1
      IF(NUL) N=N+2
    6 Y=0.0625d0*X**2-2.0d0

      GO TO (21,22,23,24), N

   21 A=        0.00000 00000 00001d0
      B=Y*A  +  0.00000 00000 00038d0
      A=Y*B-A+  0.00000 00000 01998d0
      B=Y*A-B+  0.00000 00000 90184d0
      A=Y*B-A+  0.00000 00034 97936d0
      B=Y*A-B+  0.00000 01152 07185d0
      A=Y*B-A+  0.00000 31776 13643d0
      B=Y*A-B+  0.00007 22075 54088d0
      A=Y*B-A+  0.00132 57975 71063d0
      B=Y*A-B+  0.01921 26868 62831d0
      A=Y*B-A+  0.21351 04236 00374d0
      B=Y*A-B+  1.75552 58978 34804d0
      A=Y*B-A+ 10.20622 18597 97179d0
      B=Y*A-B+ 39.59331 03435 21063d0
      A=Y*B-A+ 95.23190 27040 92164d0
      A=Y*A-B+130.34733 09768 33054d0
      BSIKR3=0.5d0*(A-B)*X**0.333333333333333d0
      GO TO 25

   22 B=       0.00000 00000 00012d0
      A=Y*B  + 0.00000 00000 00650d0
      B=Y*A-B+ 0.00000 00000 30014d0
      A=Y*B-A+ 0.00000 00011 92919d0
      B=Y*A-B+ 0.00000 00403 31192d0
      A=Y*B-A+ 0.00000 11441 45898d0
      B=Y*A-B+ 0.00002 68011 17197d0
      A=Y*B-A+ 0.00050 85533 07618d0
      B=Y*A-B+ 0.00763 78220 14940d0
      A=Y*B-A+ 0.08824 17620 77464d0
      B=Y*A-B+ 0.75671 34574 71448d0
      A=Y*B-A+ 4.60121 27076 45103d0
      B=Y*A-B+18.68656 02563 49831d0
      A=Y*B-A+46.83793 90742 83087d0
      A=Y*A-B+65.54958 48060 70052d0
      BSIKR3=0.5d0*(A-B)*X**0.666666666666666d0
      GO TO 25

   23 A=        0.00000 00000 00007d0
      B=Y*A  +  0.00000 00000 00371d0
      A=Y*B-A+  0.00000 00000 18463d0
      B=Y*A-B+  0.00000 00007 94846d0
      A=Y*B-A+  0.00000 00293 09948d0
      B=Y*A-B+  0.00000 09143 62383d0
      A=Y*B-A+  0.00002 37873 65337d0
      B=Y*A-B+  0.00050 74091 61720d0
      A=Y*B-A+  0.00869 83095 97421d0
      B=Y*A-B+  0.11697 76116 93127d0
      A=Y*B-A+  1.19852 10167 30412d0
      B=Y*A-B+  9.02594 18911 35279d0
      A=Y*B-A+ 47.81292 47970 09021d0
      B=Y*A-B+168.94690 54665 51934d0
      A=Y*B-A+374.90419 63977 61772d0
      A=Y*A-B+493.52893 55220 53900d0
      BSIKR3=0.5d0*(A-B)/X**0.333333333333333d0
      GO TO 25

   24 A=        0.00000 00000 00020d0
      B=Y*A  +  0.00000 00000 01141d0
      A=Y*B-A+  0.00000 00000 55474d0
      B=Y*A-B+  0.00000 00023 30625d0
      A=Y*B-A+  0.00000 00837 21455d0
      B=Y*A-B+  0.00000 25392 92566d0
      A=Y*B-A+  0.00006 40820 86137d0
      B=Y*A-B+  0.00132 26199 57421d0
      A=Y*B-A+  0.02187 48334 95706d0
      B=Y*A-B+  0.28291 28347 05747d0
      A=Y*B-A+  2.77810 11375 55156d0
      B=Y*A-B+ 19.98542 28723 93609d0
      A=Y*B-A+100.90075 85218 30026d0
      B=Y*A-B+340.10272 33771 79942d0
      A=Y*B-A+726.03216 50664 95924d0
      A=Y*A-B+939.90625 59898 40844d0
      BSIKR3=0.5d0*(A-B)/X**0.666666666666666d0

   25 IF(L .AND. E) BSIKR3=dEXP(-X)*BSIKR3
      IF(L) RETURN

      IF(N .GT. 2) GO TO 2
      F=BSIKR3
      N=N+2
      GO TO (21,22,23,24), N
    2 BSIKR3=1.813799364234218d0*(BSIKR3-F)
      IF(E) BSIKR3=dEXP(X)*BSIKR3
      RETURN

    1 Y=32.0d0/X-2.0d0

      GO TO (31,32), N

   31 B=     -0.00000 00000 00002d0*Y
      A=Y*B  +0.00000 00000 00014d0
      B=Y*A-B+0.00000 00000 00016d0
      A=Y*B-A-0.00000 00000 00091d0
      B=Y*A-B-0.00000 00000 00265d0
      A=Y*B-A+0.00000 00000 00298d0
      B=Y*A-B+0.00000 00000 03074d0
      A=Y*B-A+0.00000 00000 05374d0
      B=Y*A-B-0.00000 00000 14443d0
      A=Y*B-A-0.00000 00001 16415d0
      B=Y*A-B-0.00000 00003 75617d0
      A=Y*B-A-0.00000 00006 11915d0
      B=Y*A-B+0.00000 00012 91332d0
      A=Y*B-A+0.00000 00227 82111d0
      B=Y*A-B+0.00000 02520 73238d0
      A=Y*B-A+0.00000 37262 16111d0
      B=Y*A-B+0.00009 08034 04815d0
      A=Y*B-A+0.00467 34791 99874d0
      A=Y*A-B+2.00917 23421 86414d0
      GO TO 33

   32 B=      0.00000 00000 00001d0
      A=Y*B  -0.00000 00000 00001d0
      B=Y*A-B-0.00000 00000 00005d0
      A=Y*B-A+0.00000 00000 00005d0
      B=Y*A-B+0.00000 00000 00040d0
      A=Y*B-A+0.00000 00000 00012d0
      B=Y*A-B-0.00000 00000 00312d0
      A=Y*B-A-0.00000 00000 00677d0
      B=Y*A-B+0.00000 00000 01331d0
      A=Y*B-A+0.00000 00000 10138d0
      B=Y*A-B+0.00000 00000 19007d0
      A=Y*B-A-0.00000 00000 39372d0
      B=Y*A-B-0.00000 00004 20439d0
      A=Y*B-A-0.00000 00019 23009d0
      B=Y*A-B-0.00000 00075 81590d0
      A=Y*B-A-0.00000 00365 71574d0
      B=Y*A-B-0.00000 02916 95418d0
      A=Y*B-A-0.00000 41406 57716d0
      B=Y*A-B-0.00010 60188 22352d0
      A=Y*B-A-0.00646 71526 00616d0
      A=Y*A-B+1.98726 99734 33850d0

   33 BSIKR3=0.199471140200717d0*(A-B)/dSQRT(X)
      NOE=.NOT. E
      IF(NUL .OR. NOE) V=dEXP(-X)
      IF(NUL) GO TO 5
      IF(NOE) BSIKR3=BSIKR3/V
      RETURN

    5 Y=20.0d0/X-2.0d0

      GO TO (41,42), N

   41 A=     -0.00000 00000 00001d0
      B=Y*A  +0.00000 00000 00006d0
      A=Y*B-A-0.00000 00000 00040d0
      B=Y*A-B+0.00000 00000 00295d0
      A=Y*B-A-0.00000 00000 02327d0
      B=Y*A-B+0.00000 00000 19868d0
      A=Y*B-A-0.00000 00001 85898d0
      B=Y*A-B+0.00000 00019 39555d0
      A=Y*B-A-0.00000 00231 23232d0
      B=Y*A-B+0.00000 03265 50333d0
      A=Y*B-A-0.00000 57870 60592d0
      B=Y*A-B+0.00014 30095 80961d0
      A=Y*B-A-0.00631 44392 60799d0
      A=Y*A-B+1.98707 28245 52186d0
      GO TO 43

   42 A=      0.00000 00000 00001d0
      B=Y*A  -0.00000 00000 00006d0
      A=Y*B-A+0.00000 00000 00042d0
      B=Y*A-B-0.00000 00000 00308d0
      A=Y*B-A+0.00000 00000 02441d0
      B=Y*A-B-0.00000 00000 20946d0
      A=Y*B-A+0.00000 00001 97224d0
      B=Y*A-B-0.00000 00020 74924d0
      A=Y*B-A+0.00000 00250 24412d0
      B=Y*A-B-0.00000 03595 19190d0
      A=Y*B-A+0.00000 65547 92550d0
      B=Y*A-B-0.00017 13895 98262d0
      A=Y*B-A+0.00897 12068 42484d0
      A=Y*A-B+2.01829 90761 45578d0

   43 F=(A-B)/dSQRT(X)
      IF(L) GO TO 4
      BSIKR3=0.626657068657750d0*F
      IF(E) RETURN
      BSIKR3=dEXP(-X)*BSIKR3
      RETURN

    4 BSIKR3=BSIKR3+0.345494149471336d0*F*V**2
      IF(E) RETURN
      BSIKR3=BSIKR3/V
      RETURN

      ENTRY EBSKR3(X,NU)

      E=.TRUE.
      GO TO 8

      ENTRY BSKR3(X,NU)

      E=.FALSE.
    8 L=.FALSE.
      N=IABS(NU)
      IF(N .NE. 1 .AND. N .NE. 2) GO TO 10
      IF(X .LE. 0.0d0) GO TO 11
      IF(X .GE. 5.0d0) GO TO 5
      GO TO 6

   10 BSIKR3=0.d0
      PRINT 12,NU
      RETURN
   11 BSIKR3=0.d0
      IF(.NOT. L) GO TO 9
      IF(NU .GT. 0 .AND. X .EQ. 0.0) RETURN
    9 PRINT 13,X
      RETURN
   12 FORMAT(1X,30HBSIKR3 ... ILLEGAL VALUE NU = ,I15)
   13 FORMAT(1X,24HBSIKR3 ... NEGATIVE X = ,E15.4)

      END
