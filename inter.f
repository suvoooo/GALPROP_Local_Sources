
!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * inter.f *                                     galprop package * 4/14/2000 
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

      SUBROUTINE INTER(XX,FF,NN,X0,F1)
C********************************************************************
C* INTERPOLATION BY THE CUBIC POLINOM, REAL*4 version               *
C* XX,FF  ARE THE ARGUMENT AND FUNCTION ARRAYS (j=1..NN);           *
C* X0  IS THE ARGUMENT AT WHICH THE FUNCTION IS TO BE EVALUATED;    *
C* F1  ON EXIT, IS THE EVALUATED VALUE OF THE FUNCTION.             *
C********************************************************************
      IMPLICIT real (A-H,O-Z)
      IMPLICIT integer (I-N)
      DIMENSION XX(NN),FF(NN)
      IN=0
      N1=NN-1
      DO 1 M=1,N1
         IF(X0.GE.XX(M).AND.X0.LE.XX(M+1)) GOTO 2  
1           CONTINUE
2     IN=M
      IF(IN.EQ.1.OR.X0.LT.XX(1)) IN=2
      IF(IN.EQ.N1.OR.X0.GT.XX(NN)) IN=NN-2
      X1=XX(IN-1)
      X2=XX(IN)
      X3=XX(IN+1)
      X4=XX(IN+2)
      Y1=FF(IN-1)
      Y2=FF(IN)
      Y3=FF(IN+1)
      Y4=FF(IN+2)
      X01=X0-X1             
      X02=X0-X2       
      X03=X0-X3       
      X04=X0-X4       
      X12=X1-X2       
      X13=X1-X3       
      X14=X1-X4       
      X23=X2-X3       
      X24=X2-X4       
      X34=X3-X4
      X1= X02*X03*X04/X12/X13/X14
      X2=-X01*X03*X04/X12/X23/X24
      X3= X01*X02*X04/X13/X23/X34
      X4=-X01*X02*X03/X14/X24/X34
      F1= Y1*X1+Y2*X2+Y3*X3+Y4*X4
      RETURN
      END

      SUBROUTINE DINTER(XX,FF,NN,X0,F1)
C********************************************************************
C* INTERPOLATION BY CUBIC POLINOM, REAL*8 version                   *
C* XX,FF  ARE THE ARGUMENT AND FUNCTION ARRAYS (j=1..NN);           *
C* X0  IS THE ARGUMENT AT WHICH THE FUNCTION IS TO BE EVALUATED;    *
C* F1  ON EXIT, IS THE EVALUATED VALUE OF THE FUNCTION.             *
C********************************************************************
      IMPLICIT real*8  (A-H,O-Z)
      IMPLICIT integer (I-N)
      DIMENSION XX(NN),FF(NN)
      IN=0
      N1=NN-1
      DO 1 M=1,N1
         IF(X0.GE.XX(M).AND.X0.LE.XX(M+1)) GOTO 2  
1           CONTINUE
2     IN=M
      IF(IN.EQ.1.OR.X0.LT.XX(1)) IN=2
      IF(IN.EQ.N1.OR.X0.GT.XX(NN)) IN=NN-2
      X1=XX(IN-1)
      X2=XX(IN)
      X3=XX(IN+1)
      X4=XX(IN+2)
      Y1=FF(IN-1)
      Y2=FF(IN)
      Y3=FF(IN+1)
      Y4=FF(IN+2)
      X01=X0-X1             
      X02=X0-X2       
      X03=X0-X3       
      X04=X0-X4       
      X12=X1-X2       
      X13=X1-X3       
      X14=X1-X4       
      X23=X2-X3       
      X24=X2-X4       
      X34=X3-X4
      X1= X02*X03*X04/X12/X13/X14
      X2=-X01*X03*X04/X12/X23/X24
      X3= X01*X02*X04/X13/X23/X34
      X4=-X01*X02*X03/X14/X24/X34
      F1= Y1*X1+Y2*X2+Y3*X3+Y4*X4
      RETURN
      END

