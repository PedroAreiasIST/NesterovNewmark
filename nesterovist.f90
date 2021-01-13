MODULE smsutility
END MODULE smsutility
!-----------------------------------------------
!*** Areias 2020
!*** Nesterov-based optimizer
!*** with restart and automated step calculation
!-----------------------------------------------
MODULE NVALS
  INTEGER,SAVE::NFE=0
END MODULE NVALS

!-------------------------------
!*** golden section, trump style
!-------------------------------
SUBROUTINE goldensection(FUNC,N,Q,PHI,G,DIR,STP)
  IMPLICIT REAL(8)(A-H,O-Z)
  EXTERNAL::FUNC
  REAL(8),DIMENSION(N)::G,Q,GASH,DIR
!*** MITER=20
!*** GIVES THE SMALLEST
  INTEGER::MITER=500
  STP0=50.0d00*STP
  IF(ABS(STP0).LE.1.0D-30)THEN
     STP0=1.0D00
     MITER=20
  END IF
  gr=1.61803d00!0.5d00*(SQRT(5.0d00)+1.0d00)
  a=0.0d00
  b=stp0
  phia=phi
  GASH=Q+a*DIR
  CALL func(n,gash,phia)
  GASH=Q+b*DIR
  CALL func(n,gash,phib)
  c=b-(b-a)/gr
  d=a+(b-a)/gr
  DO I=0,MITER
     GASH=Q+C*DIR
     CALL FUNC(N,GASH,PHIC)
     GASH=Q+D*DIR
     CALL FUNC(N,GASH,PHID)     
     IF(PHIc.LT.phid)THEN
        b=d
     ELSE
        a=c
     END IF
     c=b-(b-a)/gr
     d=a+(b-a)/gr
  END DO
  stp=0.5d00*(a+b)
END SUBROUTINE GOLDENSECTION

!---------
!*** alpha
!---------
REAL(8) FUNCTION ALPHADET(N,Q0,Q1,PHI0,PHI1,G0,G1)
  IMPLICIT REAL(8)(A-H,O-Z)
  REAL(8),DIMENSION(N)::Q0,Q1,G0,G1
  ALPHADET=DOTP(N,Q1-Q0,Q1-Q0)/(2.0D00*(PHI0-PHI1+DOTP(N,Q1-Q0,G1)))
END FUNCTION ALPHADET

!---------------------------
!*** FIRST STEP
!*** STEP SIZE DETERMINATION
!---------------------------
REAL(8) FUNCTION FIRSTSTEP(FUNC,PROC,N,Q)
  IMPLICIT REAL(8)(A-H,O-Z)
  EXTERNAL::PROC,FUNC
  REAL(8),DIMENSION(N)::G,Q,DIR
  CALL PROC(N,Q,PHI,G)
  DIR=-G
  STP=0.0D00
  CALL SUFFICIENTLINESEARCH(FUNC,N,Q,PHI,G,DIR,STP)
  FIRSTSTEP=STP
END FUNCTION FIRSTSTEP

!-----------------------------------
!*** SUFFICIENTLY DECREASE BACKTRACK
!*** ALGORITHM NOCEDAL AND WRIGHT
!-----------------------------------
SUBROUTINE SUFFICIENTLINESEARCH(FUNC,N,Q,PHI,G,DIR,STP)
  IMPLICIT REAL(8)(A-H,O-Z)
  EXTERNAL::FUNC
  REAL(8),DIMENSION(N)::G,Q,GASH,DIR
!*** MITER=20
!*** GIVES THE SMALLEST
!*** STEP AS 9.53674E-7
  INTEGER::MITER=10
  STP0=STP
  IF(ABS(STP0).LE.1.0D-30)THEN
     STP0=1.0D00
     MITER=20
  END IF
  C=1.0D-4
  GP=MIN(0.0d00,DOTP(N,G,DIR))
  DO I=0,MITER
     STP=4.0D00*STP0*0.5D00**I
     GASH=Q+STP*DIR
     CALL FUNC(N,GASH,PHII)
     IF(PHII.LE.PHI+C*STP*GP)EXIT
  END DO
END SUBROUTINE SUFFICIENTLINESEARCH

!---------------
!*** nonmonotone
!---------------
SUBROUTINE nonmonotone(FUNC,N,Q,PHI0,PHI1,PHI2,GPHI0,GPHI1,GPHI2,DIR,STP)
  IMPLICIT REAL(8)(A-H,O-Z)
  EXTERNAL::FUNC
  REAL(8),DIMENSION(N)::GPHI0,GPHI1,GPHI2,DIR,Q,GASH
!*** MITER=20
!*** GIVES THE SMALLEST
!*** STEP AS 9.53674E-7
  INTEGER::MITER=10
  STP0=STP
  IF(ABS(STP0).LE.1.0D-30)THEN
     STP0=1.0D00
     MITER=20
  ELSE
     STP0=4.0D00*STP0
  END IF
  C=1.0D-4
  FMAX=MAX(PHI0,PHI1,PHI2)
  GP=DOTP(N,GPHI1,DIR)
  DO I=0,MITER
     STP=STP0*0.5D00**I
     GASH=Q+STP*DIR
     CALL FUNC(N,GASH,PHII)
     IF(PHII.LE.FMAX+C*STP*GP)EXIT
  END DO
END SUBROUTINE NONMONOTONE

!-----------------------
!*** UPDATE THE SOLUTION
!-----------------------
SUBROUTINE UPDASTEP(N,IALGO,Q0,Q1,Q2,K,C,S,LAMBDA,PROC,PHI,GPHI)
  IMPLICIT REAL(8)(A-H,O-Z)
  EXTERNAL::PROC
  REAL(8)::LAMBDA,S,PHI
  REAL(8),DIMENSION(N)::Q0,Q1,Q2,GPHI,Q1TEMP
!--------------------
!*** TEMPORARY VALUES
!--------------------
  H=S*(2.0D00*LAMBDA)**(-1.0D00/2.0D00)
  IF(K.NE.0)THEN
     T=H*K
     CEQ=C/T
  END IF
!------------------------
!*** SELECT THE ALGORITHM
!------------------------
  IF(K.LE.0)THEN
     Q2=Q1-(H**2.0d00)*GPHI
  ELSE
     SELECT CASE(IALGO)
     CASE(1)
        CONST1=(2.0D00-H*CEQ)/(2.0D00+H*CEQ)
        CONST2=(2.0D00*H**2.0D00)/(2.0D00+H*CEQ)
        Q2=Q1+CONST1*(Q1-Q0)-CONST2*GPHI
     CASE(2)
        CONST1=(2.0D00-H*CEQ)/(2.0D00+H*CEQ)
        CONST2=(2.0D00*H**2.0D00)/(2.0D00+H*CEQ)        
        Q1TEMP=Q1+CONST1*(Q1-Q0)
        CALL PROC(N,Q1TEMP,PHI,GPHI)
        Q2=Q1TEMP-CONST2*GPHI
     CASE(3)
        CONST1=(1.0D00-H*CEQ)
        CONST2=H**2
        Q2=Q1+CONST1*(Q1-Q0)-CONST2*GPHI
     CASE(4)
        CONST1=(1.0D00-H*CEQ)
        CONST2=H**2
        Q1TEMP=Q1+CONST1*(Q1-Q0)
        CALL PROC(N,Q1TEMP,PHI,GPHI)
        Q2=Q1TEMP-CONST2*GPHI
     END SELECT
  END IF
END SUBROUTINE UPDASTEP

!------------------
!*** EUCLIDEAN NORM
!------------------
REAL(8) FUNCTION EUCLIDEANNORM(N,G)
  IMPLICIT REAL(8)(A-H,O-Z)
  REAL(8),DIMENSION(N)::G
  EUCLIDEANNORM=0.0D00
  RMAX=0.0D00
  DO I=1,N
     RMAX=MAX(RMAX,ABS(G(I)))
  END DO
  IF(ABS(RMAX).GT.1.0D-30)THEN
     RMAX1=1.0D00/RMAX
     DO I=1,N
        T=G(I)*RMAX1
        EUCLIDEANNORM=EUCLIDEANNORM+T*T
     END DO
     EUCLIDEANNORM=SQRT(EUCLIDEANNORM)*RMAX
  ELSE
     EUCLIDEANNORM=0.0D00
  END IF
END FUNCTION EUCLIDEANNORM

!---------------
!*** DOT PRODUCT
!---------------
REAL(8) FUNCTION DOTP(N,A,B)
  IMPLICIT REAL(8)(A-H,O-Z)
  REAL(8),DIMENSION(*)::A,B  
  DOTP=0.0D00
  DO I=1,N
     DOTP=DOTP+A(I)*B(I)
  END DO
END FUNCTION DOTP

!-------------
!*** OPTIMIZER
!------------- 
SUBROUTINE OPTIMIZER(N,NITER,IALGO,S,Q,C,LAMBDA,FUNC,PROC,HISTPHI,HISTGRAD,tol)
  IMPLICIT REAL(8)(A-H,O-Z)
  REAL(8)::LAMBDA
  LOGICAL::CONST
  REAL(8),DIMENSION(N)::Q,Q0,Q1,Q2,GPHI,GPHI0,GPHI1,GPHI2,DIR
  REAL(8),DIMENSION(NITER)::HISTPHI,HISTGRAD
  LOGICAL::CLARGE,CSMALL
  EXTERNAL::FUNC,PROC
!------------------
!*** INITIAL VALUES
!------------------
  Q0=Q
  Q1=Q
  Q2=Q
  K=0
  CONST=.FALSE.
  IF(S.GT.1.0D-20)CONST=.TRUE.
  PHI1=0.0D00
  PHI2=0.0D00
  GPHI1=0.0D00
  GPHI2=0.0D00
  DIR=0.0D00
  DO ITER=0,NITER-1
     PHI0=PHI1
     PHI1=PHI2
     GPHI0=GPHI1
     GPHI1=GPHI2
!-------------------------------------
!*** EVALUATES FUNCTION AND DERIVATIVE
!-------------------------------------
     CALL PROC(N,Q2,PHI2,GPHI2)
     IF(ITER.EQ.0)THEN
        PHI0=PHI2
        PHI1=PHI2
        GPHI0=GPHI2
        GPHI1=GPHI2
     END IF
!---------------------------------------
!*** NOW IT'S ALL SYNCRONIZED, 0 1 AND 2
!---------------------------------------     
!----------------------
!*** DETERMINES S AND K
!----------------------
     CALL SIZEDETERMINATIONANDRESTART(CONST,LAMBDA,C,FUNC,PROC,K,ITER,S,N,Q0,Q1,Q2,PHI0,PHI1,PHI2,GPHI0,GPHI1,GPHI2)
!-----------------
!*** DETERMINES Q2
!-----------------
     Q0=Q1
     Q1=Q2
     CALL UPDASTEP(N,IALGO,Q0,Q1,Q2,K,C,S,LAMBDA,PROC,PHI2,GPHI2)
!--------------------------------
!*** UPDATE STEP BETWEEN RESTARTS
!--------------------------------
     K=K+1
!------------------
!*** STORE FUNCTION
!------------------
     HISTPHI(ITER+1)=PHI2
     HISTGRAD(ITER+1)=EUCLIDEANNORM(N,GPHI2)
!-------------------------
!*** CHECK FOR CONVERGENCE
!-------------------------
     IF(EUCLIDEANNORM(N,GPHI2).LE.TOL)EXIT
  END DO
END SUBROUTINE OPTIMIZER

!----------------
!*** CG OPTIMIZER
!----------------
SUBROUTINE CG(N,NITER,IALGO,S,Q,C,LAMBDA,FUNC,PROC,HISTPHI,HISTGRAD,TOL)
  IMPLICIT REAL(8)(A-H,O-Z)
  REAL(8)::LAMBDA
  LOGICAL::CONST
  REAL(8),DIMENSION(N)::Q,Q0,Q1,Q2,GPHI,DX,dxo,dir,gphio
  REAL(8),DIMENSION(NITER)::HISTGRAD,HISTPHI
  LOGICAL::CLARGE,CSMALL
  EXTERNAL::FUNC,PROC
!------------------
!*** INITIAL VALUES
!------------------
  DIR=0.0D00
  STP=0.0D00
  dx=0.d00
  DO ITER=0,NITER-1
     GPHIO=GPHI
     CALL PROC(N,Q,PHI,GPHI)
     DXO=DX
     DX=-GPHI
     IF(ITER.EQ.0)THEN
        BETA=0.0D00
     ELSE
        BETA=MAX(0.0D00,DOTP(N,DX,DX-DXO)/DOTP(N,DXO,DXO))
     END IF
     DIR=DX+BETA*DIR
     CALL SUFFICIENTLINESEARCH(FUNC,N,Q,PHI,GPHI,DIR,STP)
     Q=Q+STP*DIR
     HISTGRAD(ITER+1)=EUCLIDEANNORM(N,GPHI)
     HISTPHI(ITER+1)=PHI
     IF(EUCLIDEANNORM(N,GPHI).LE.TOL)EXIT
  END DO
END SUBROUTINE CG

SUBROUTINE NST(N,NITER,IALGO,S,Q,C,LAMBDA,FUNC,PROC,HISTGRAD,TOL)
  IMPLICIT REAL(8)(A-H,O-Z)
  REAL(8)::LAMBDA
  LOGICAL::CONST
  REAL(8),DIMENSION(N)::diro,Q,Q0,Q1,Q2,GPHI,DX,dxo,dir,gphio
  REAL(8),DIMENSION(NITER)::HISTGRAD
  LOGICAL::CLARGE,CSMALL
  EXTERNAL::FUNC,PROC
!------------------
!*** INITIAL VALUES
!------------------
  DIR=0.0D00
  dx=0.0d00
  DO ITER=0,NITER-1
     CALL PROC(N,Q,PHI,GPHI)
     diro=dir
     dir=q-s*gphi
     q=dir+3.0d00*(dir-diro)/(iter+2)
     WRITE(*,*)"ITER,GPHI1",ITER,EUCLIDEANNORM(N,GPHI)
     HISTGRAD(ITER+1)=EUCLIDEANNORM(N,GPHI) 
     IF(EUCLIDEANNORM(N,GPHI).LE.TOL)EXIT
  END DO
END SUBROUTINE NST

!---------------------------------------------
!*** DETERMINATION OF STEP SIZES, RESTART, ETC
!---------------------------------------------
SUBROUTINE SIZEDETERMINATIONANDRESTART(CONST,LAMBDA,C,FUNC,PROC,K,ITER,S,N,Q0,Q1,Q2,PHI0,PHI1,PHI2,GPHI0,GPHI1,GPHI2)
  IMPLICIT REAL(8)(A-H,O-Z)
  LOGICAL::CONST
  REAL(8)::S,LAMBDA,PHI0,PHI1,PHI2
  REAL(8),DIMENSION(N)::Q0,Q1,Q2,GPHI0,GPHI1,GTEMP,GPHI2,DIR,qtemp
  LOGICAL::C1,C2,C3,C4,C5,C6,C7,C8
  EXTERNAL::PROC,FUNC
  IF(ITER.EQ.0.AND..NOT.CONST)THEN
     S=SQRT(2.0D00*LAMBDA)*SQRT(FIRSTSTEP(FUNC,PROC,N,Q2))
  ELSE
     IF(K.GE.2)THEN
        GTEMP=Q2-Q0
        C5=SIGN(1.0D00,DOTP(N,GPHI1,GTEMP)).GT.0.0D00        
        GTEMP=Q2-Q0
        C6=SIGN(1.0D00,DOTP(N,GPHI2,GTEMP)).GT.0.0D00
        C7=(PHI2.GT.1.1d00*PHI0).OR.(PHI2.GT.1.1D00*PHI1).OR.(PHI1.GT.1.1D00*PHI0)
        C8=EUCLIDEANNORM(N,GPHI2).GT.1.1d00*EUCLIDEANNORM(N,GPHI0)
        IF(C5.OR.C6.OR.C7.OR.C8)THEN
           IF(.NOT.CONST)THEN
              WRITE(*,*)"C5,C6,C7,C8",C5,C6,C7,C8
              H1=S**2/(2.0D00*LAMBDA)
              DIR=-GPHI1
              CALL SUFFICIENTLINESEARCH(FUNC,N,Q1,PHI1,GPHI1,DIR,H1)
              S1=SQRT(2.0D00*LAMBDA*H1)
              H2=ALPHADET(N,Q0,Q2,PHI0,PHI2,Gphi0,Gphi2)
              IF(H2.LE.0.0D00)H2=1.0D30
              H3=ALPHADET(N,Q0,Q1,PHI0,PHI1,Gphi0,Gphi1)
              IF(H3.LE.0.0D00)H3=1.0D30
              H4=ALPHADET(N,Q1,Q2,PHI1,PHI2,Gphi1,Gphi2)
              IF(H4.LE.0.0D00)H4=1.0D30
              HMED=(H2+H3+H4)/3.0D00
              SMED=SQRT(2.0D00*LAMBDA*HMED)
              IF(K.GE.4)THEN
                 S=S1 !MIN(S1,SMED)
              ELSE
                 S=S1
              END IF
           END IF
           K=1
        END IF
     END IF
  END IF
END SUBROUTINE SIZEDETERMINATIONANDRESTART

SUBROUTINE ROSENBROCK(N,X,F,DF)
  USE nvals
  IMPLICIT REAL(8)(A-H,O-Z)
  REAL(8),DIMENSION(N)::DF,X
  NFE=NFE+1
  CALL minimizationtestfunctions(x,f,df)
END SUBROUTINE ROSENBROCK

SUBROUTINE ROSENBROCKF(N,X,F)
  USE nvals
  IMPLICIT REAL(8)(A-H,O-Z)
  REAL(8),DIMENSION(N)::X,df
  NFE=NFE+1
  CALL minimizationtestfunctions(x,f,df)
END SUBROUTINE ROSENBROCKF

PROGRAM TESTE
  USE nvals
  IMPLICIT REAL(8)(A-H,O-Z)
  EXTERNAL ROSENBROCK, ROSENBROCKF
  REAL(8)::LAMBDA
  REAL(8),DIMENSION(40)::X,XO
  REAL(8),DIMENSION(1000)::HISTPHI,HISTGRAD

  CALL initialconditions(xo,n)
  NITER=500
  IALGO=1

  C=3.0D00
  LAMBDA=1.0d2
  S=0.00d00
  TOL=1.0D-3
  histgrad=0.0d00
  histphi=0.0d00
  X=XO
  
  CALL OPTIMIZER(N,NITER,IALGO,S,X,C,LAMBDA,ROSENBROCKF,ROSENBROCK,HISTPHI,HISTGRAD,TOL)  

  OPEN(12,FILE="nesterov.txt",STATUS="UNKNOWN")
  WRITE(*,*)"NFE Nesterov",nfe
  DO I=1,NITER
     TEMP=HISTGRAD(I)
     WRITE(12,"(I6,2E15.4)")I,HISTPHI(I),HISTGRAD(I)
  END DO
  CLOSE(12)

  S=0.00d00
  histgrad=0.0d00
  histphi=0.0d00
  X=XO
  NFE=0
  CALL CG(n,NITER,IALGO,S,X,C,LAMBDA,ROSENBROCKF,ROSENBROCK,HISTPHI,HISTGRAD,TOL)

  OPEN(12,FILE="cg.txt",STATUS="UNKNOWN")  
  DO I=1,NITER
     TEMP=HISTGRAD(I)
     WRITE(12,"(I6,2E15.4)")I,HISTPHI(I),HISTGRAD(I)
  END DO
  WRITE(*,*)"NFE cg",nfe
  CLOSE(12)
END PROGRAM TESTE

!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           22 Nov 20 17:17:03 *
!**************************************************************
! User     : Full professional version
! Notebook : minimizationtestfunctions
! Evaluation time                 : 1 s     Mode  : Optimal
! Number of formulae              : 10      Method: Automatic
! Subroutine                      : minimizationtestfunctions size: 188
! Total size of Mathematica  code : 188 subexpressions
! Total size of Fortran code      : 563 bytes

!******************* S U B R O U T I N E **********************
SUBROUTINE initialconditions(q,n)
  IMPLICIT REAL(8)(a-h,o-z)
  REAL(8),DIMENSION(*)::q
  ijob=1
  SELECT CASE(ijob)
!*** ROSENBROCK
  CASE(1)     
     n=2
     q(1)=-2.0d00 !-1.2d00
     q(2)=-2.0d00 !1.0d00
!*** FREUDENSTEIN AND ROTH
  CASE(2)
     n=2
     q(1)=0.5d00
     q(2)=-2.0d00
!*** BEALE
  CASE(3)
     n=2
     q(1)=1.0D00
     q(2)=1.0d00
!*** POWELL SINGULAR
  CASE(4)
     n=4
     q(1)=3.0d00
     q(2)=-1.0d00
     q(3)=0.0d00
     q(4)=1.0d00
!*** Wood     
  CASE(5)
     n=4
     q(1)=-3.0d00
     q(2)=-1.0d00
     q(3)=-3.0d00
     q(4)=-1.0d00
!*** Brown and Dennis
  CASE(6)
     n=4
     q(1)=25.0d00
     q(2)=5.0d00
     q(3)=-5.0d00
     q(4)=-1.0d00
!*** Biggs Exp 6
  CASE(7)
     n=6
     q(1)=1.0d00
     q(2)=2.0d00
     q(3)=1.0d00
     q(4)=1.0d00
     q(5)=1.0d00
     q(6)=1.0d00
!*** trigonometric function
  CASE(8)
     n=6
     q(1:n)=1.0d0/n
  END SELECT
END SUBROUTINE initialconditions


SUBROUTINE minimizationtestfunctions(q,f,df)
  USE SMSUtility
  IMPLICIT NONE
  INTEGER::ijob,n
  DOUBLE PRECISION v(50000),q(*),f,df(*)
  ijob=1
  SELECT CASE(ijob)
  CASE(1)! rosenbrock
     v(11)=-q(1)**2+q(2)
     v(10)=1d0-q(1)
     f=(v(10)*v(10))+100d0*(v(11)*v(11))
     df(1)=(-2d0)*v(10)-400d0*q(1)*v(11)
     df(2)=200d0*v(11)
  CASE(2)! Freudenstein and Roth
     v(15)=2d0*q(2)
     v(14)=-(((-5d0)+q(2))*q(2))
     v(16)=(-2d0)+v(14)
     v(13)=q(2)*(1d0+q(2))
     v(18)=(-14d0)+v(13)
     v(8)=(-29d0)+q(1)+q(2)*v(18)
     v(19)=2d0*v(8)
     v(7)=(-13d0)+q(1)+q(2)*v(16)
     v(17)=2d0*v(7)
     f=(v(7)*v(7))+(v(8)*v(8))
     df(1)=v(17)+v(19)
     df(2)=(q(2)*(5d0-v(15))+v(16))*v(17)+(q(2)*(1d0+v(15))+v(18))*v(19)
  CASE(3)! beale
     v(17)=1d0-q(2)
     v(16)=1d0-q(2)**3
     v(15)=q(2)**2
     v(12)=0.2625d1-q(1)*v(16)
     v(9)=1d0-v(15)
     v(10)=0.225d1-q(1)*v(9)
     v(7)=0.15d1-q(1)*v(17)
     v(18)=2d0*v(7)
     f=(v(10)*v(10))+(v(12)*v(12))+(v(7)*v(7))
     df(1)=(-2d0)*v(12)*v(16)-v(17)*v(18)-2d0*v(10)*v(9)
     df(2)=q(1)*(4d0*q(2)*v(10)+6d0*v(12)*v(15)+v(18))
  CASE(5)! wood
     v(25)=q(2)-q(4)
     v(24)=(-2d0)+q(2)+q(4)
     v(23)=-q(3)**2+q(4)
     v(22)=1d0-q(3)
     v(21)=-q(1)**2+q(2)
     v(20)=1d0-q(1)
     f=(v(20)*v(20))+100d0*(v(21)*v(21))+(v(22)*v(22))+90d0*(v(23)*v(23))+10d0*(v(24)*v(24))+(v(25)*v(25))/10d0
     v(14)=v(25)/5d0
     v(16)=20d0*v(24)
     df(1)=(-2d0)*v(20)-400d0*q(1)*v(21)
     df(2)=v(14)+v(16)+200d0*v(21)
     df(3)=(-2d0)*v(22)-360d0*q(3)*v(23)
     df(4)=-v(14)+v(16)+180d0*v(23)
  CASE(6)! brown and dennis
     v(35)=q(3)-dcos(0.2d0)+q(4)*dsin(0.2d0)
     v(34)=q(3)-dcos(0.4d0)+q(4)*dsin(0.4d0)
     v(33)=q(3)-dcos(0.6d0)+q(4)*dsin(0.6d0)
     v(32)=q(3)-dcos(0.8d0)+q(4)*dsin(0.8d0)
     v(31)=(-0.12214027581601698d1)+q(1)+0.2d0*q(2)
     v(30)=(-0.14918246976412702d1)+q(1)+0.4d0*q(2)
     v(29)=(-0.1822118800390509d1)+q(1)+0.6d0*q(2)
     v(28)=(-0.22255409284924677d1)+q(1)+0.8d0*q(2)
     v(18)=(v(28)*v(28))+(v(32)*v(32))
     v(36)=v(18)*v(28)
     v(16)=(v(29)*v(29))+(v(33)*v(33))
     v(37)=v(16)*v(29)
     v(14)=(v(30)*v(30))+(v(34)*v(34))
     v(38)=v(14)*v(30)
     v(12)=(v(31)*v(31))+(v(35)*v(35))
     v(39)=v(12)*v(31)
     f=(v(12)*v(12))+(v(14)*v(14))+(v(16)*v(16))+(v(18)*v(18))
     v(20)=4d0*v(12)*v(35)
     v(22)=4d0*v(14)*v(34)
     v(24)=4d0*v(16)*v(33)
     v(26)=4d0*v(18)*v(32)
     df(1)=4d0*(v(36)+v(37)+v(38)+v(39))
     df(2)=0.32d1*v(36)+0.24d1*v(37)+0.16d1*v(38)+0.8d0*v(39)
     df(3)=v(20)+v(22)+v(24)+v(26)
     df(4)=v(20)*dsin(0.2d0)+v(22)*dsin(0.4d0)+v(24)*dsin(0.6d0)+v(26)*dsin(0.8d0)
  CASE(7)!! Biggs Exp6
     v(57)=dexp((-0.6000000000000001d0)*q(5))
     v(56)=dexp((-0.5d0)*q(5))
     v(55)=dexp((-0.4d0)*q(5))
     v(54)=dexp((-0.30000000000000004d0)*q(5))
     v(53)=dexp((-0.2d0)*q(5))
     v(52)=dexp((-0.1d0)*q(5))
     v(51)=dexp((-0.6000000000000001d0)*q(2))
     v(50)=dexp((-0.5d0)*q(2))
     v(49)=dexp((-0.4d0)*q(2))
     v(48)=dexp((-0.30000000000000004d0)*q(2))
     v(47)=dexp((-0.2d0)*q(2))
     v(46)=dexp((-0.1d0)*q(2))
     v(45)=dexp((-0.6000000000000001d0)*q(1))
     v(44)=dexp((-0.5d0)*q(1))
     v(43)=dexp((-0.4d0)*q(1))
     v(42)=dexp((-0.30000000000000004d0)*q(1))
     v(41)=dexp((-0.2d0)*q(1))
     v(40)=dexp((-0.1d0)*q(1))
     v(26)=(-0.10764003502856656d1)+q(3)*v(40)-q(4)*v(46)+q(6)*v(52)
     v(70)=2d0*v(26)
     v(64)=v(26)*v(46)
     v(58)=(-0.2d0)*v(26)
     v(24)=(-0.1490041229246583d1)+q(3)*v(41)-q(4)*v(47)+q(6)*v(53)
     v(71)=2d0*v(24)
     v(65)=v(24)*v(47)
     v(59)=(-0.4d0)*v(24)
     v(22)=(-0.13954655145790046d1)+q(3)*v(42)-q(4)*v(48)+q(6)*v(54)
     v(72)=2d0*v(22)
     v(66)=v(22)*v(48)
     v(60)=(-0.6000000000000001d0)*v(22)
     v(20)=(-0.11844314055759346d1)+q(3)*v(43)-q(4)*v(49)+q(6)*v(55)
     v(73)=2d0*v(20)
     v(67)=v(20)*v(49)
     v(61)=(-0.8d0)*v(20)
     v(18)=(-0.9788467744270442d0)+q(3)*v(44)-q(4)*v(50)+q(6)*v(56)
     v(75)=v(18)*v(56)
     v(68)=v(18)*v(50)
     v(62)=v(18)*v(44)
     v(16)=(-0.8085717350789321d0)+q(3)*v(45)-q(4)*v(51)+q(6)*v(57)
     v(74)=2d0*v(16)
     v(69)=v(16)*v(51)
     v(63)=(-0.12000000000000002d1)*v(16)
     f=(v(16)*v(16))+(v(18)*v(18))+(v(20)*v(20))+(v(22)*v(22))+(v(24)*v(24))+(v(26)*v(26))
     df(1)=q(3)*(v(40)*v(58)+v(41)*v(59)+v(42)*v(60)+v(43)*v(61)-v(62)+v(45)*v(63))
     df(2)=q(4)*(0.2d0*v(64)+0.4d0*v(65)+0.6000000000000001d0*v(66)+0.8d0*v(67)+v(68)+0.12000000000000002d1*v(69))
     df(3)=2d0*v(62)+v(40)*v(70)+v(41)*v(71)+v(42)*v(72)+v(43)*v(73)+v(45)*v(74)
     df(4)=(-2d0)*(v(64)+v(65)+v(66)+v(67)+v(68)+v(69))
     df(5)=q(6)*(v(52)*v(58)+v(53)*v(59)+v(54)*v(60)+v(55)*v(61)+v(57)*v(63)-v(75))
     df(6)=v(52)*v(70)+v(53)*v(71)+v(54)*v(72)+v(55)*v(73)+v(57)*v(74)+2d0*v(75)
  CASE(8) !! Trigonometric function
     v(47)=dcos(q(1))
     v(46)=dcos(q(6))
     v(45)=dsin(q(6))
     v(44)=dcos(q(5))
     v(43)=dsin(q(5))
     v(42)=dcos(q(4))
     v(41)=dsin(q(4))
     v(40)=dcos(q(3))
     v(39)=dsin(q(3))
     v(38)=dcos(q(2))
     v(37)=dsin(q(2))
     v(36)=dsin(q(1))
     v(15)=-v(38)-v(40)-v(42)-v(44)-v(46)
     v(23)=7d0+v(15)-v(36)-2d0*v(47)
     v(50)=2d0*v(23)
     v(20)=6d0+v(15)-v(47)
     v(29)=6d0+v(20)-v(45)-6d0*v(46)
     v(49)=2d0*v(29)
     v(52)=v(49)+v(50)
     v(28)=5d0+v(20)-v(43)-5d0*v(44)
     v(54)=2d0*v(28)
     v(55)=v(52)+v(54)
     v(27)=4d0+v(20)-v(41)-4d0*v(42)
     v(57)=2d0*v(27)
     v(59)=v(55)+v(57)
     v(26)=3d0+v(20)-v(39)-3d0*v(40)
     v(56)=2d0*v(26)
     v(25)=2d0+v(20)-v(37)-2d0*v(38)
     v(58)=2d0*v(25)
     v(53)=v(56)+v(58)
     v(51)=v(53)+v(57)
     v(48)=v(51)+v(54)
     f=(v(23)*v(23))+(v(25)*v(25))+(v(26)*v(26))+(v(27)*v(27))+(v(28)*v(28))+(v(29)*v(29))
     df(1)=v(36)*(v(48)+v(49))+(2d0*v(36)-v(47))*v(50)
     df(2)=(3d0*v(37)-v(38))*v(58)+v(37)*(v(56)+v(59))
     df(3)=(4d0*v(39)-v(40))*v(56)+v(39)*(v(58)+v(59))
     df(4)=v(41)*(v(53)+v(55))+(5d0*v(41)-v(42))*v(57)
     df(5)=v(43)*(v(51)+v(52))+(6d0*v(43)-v(44))*v(54)
     df(6)=(7d0*v(45)-v(46))*v(49)+v(45)*(v(48)+v(50))
  END SELECT
END SUBROUTINE minimizationtestfunctions
