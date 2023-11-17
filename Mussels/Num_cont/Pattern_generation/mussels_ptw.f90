!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   enz :    A two-cell, one-substrate enzyme model 
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION m,mp,a,s,sp,delta,eta,alpha,beta,theta,c,d,nu

       m=U(1)
       mp=U(2)
       a=U(3) 
	   s=U(4)
	   sp=U(5)
       alpha=PAR(1)
       beta=PAR(2)
       d=PAR(3)
       nu=PAR(4)
       c=PAR(5)
	   delta=PAR(6)
	   eta=PAR(7)
	   theta=PAR(8)
    
       F(1) = mp
       F(2) = m - delta*a*m*(s+eta)/(s+1) - c*mp
       F(3) = 1/(nu+c)*(-1+alpha*a+beta*a*m*(s+eta)/(s+1))
	   F(4) = sp
	   F(5) = 1/d*(-m+theta*s-c*sp)

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION m,mp,a,s,sp,delta,eta,alpha,beta,theta,c,d,nu,meq,aeq,seq,temp

 
	   
	  OPEN(14,FILE='../parameters.dat',STATUS='OLD')
      REWIND(14)
      READ(14,*)PAR(1),PAR(2),PAR(7),PAR(8),PAR(5),PAR(3),PAR(4)
      CLOSE(14)
      
	  !parameters_start = [alpha, beta, delta, eta, theta, c, d, nu];
	   
	   PAR(6)=400.0
	   alpha=PAR(1)
       beta=PAR(2)
       d=PAR(3)
       nu=PAR(4)
       c=PAR(5)
	   delta=PAR(6)
	   eta=PAR(7)
	   theta=PAR(8)
       
	  temp = delta - alpha - beta*theta*eta
	   seq = 1/(2*beta*theta)*(temp + sqrt((temp)**2 -4*beta*theta*(alpha-delta*eta)))
       U(1)=theta*seq
       U(2)=0 
	   U(3)=(seq+1)/(delta*(seq+eta))
       U(4)=seq
	   U(5)=0

      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS(NDIM,U,PAR)
	  IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

      DOUBLE PRECISION, EXTERNAL :: GETP,GETU2
      INTEGER NDX,NCOL,NTST
	  
	  
	  PAR(15)=GETP('NRM',1,U)
	  PAR(16)=GETP('NRM',3,U)
	 ! PAR(17)=GETP('NRM',5,U)
	  PAR(18)=GETP('MIN',3,U)
	  PAR(19)=GETP('MAX',3,U)
	  PAR(20)=GETP('MIN',1,U)
	  PAR(21)=GETP('MAX',1,U)
      END SUBROUTINE PVLS
