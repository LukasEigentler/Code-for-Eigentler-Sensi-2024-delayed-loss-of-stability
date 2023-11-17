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

      DOUBLE PRECISION p1,p1p,w,wp,A,B1,c,d,nu

       
       p1=U(1)
       p1p=U(2)
	   
       w=U(3)
       wp=U(4)
       
       A=PAR(1)
       B1=PAR(2)
   
       d=PAR(3)
      
       
       nu=PAR(4)
       c=PAR(5)
    
	   
	   
       
     
		
       F(1) = p1p
       F(2) = -p1**2*w+B1*p1-c*p1p
       
     
       F(3) = wp
       F(4) = 1/d*(-A+w*p1**2+w-(nu+c)*wp)

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION A,B1,nu,c,d,u1eq,weq

 
	   
	  OPEN(14,FILE='../parameters.dat',STATUS='OLD')
      REWIND(14)
      READ(14,*)PAR(2),PAR(5),PAR(3),PAR(4)
      CLOSE(14)
      
	  !parameters_start = [ B1, c, d,  nu];
	   
	   PAR(1)=4.5
	   !PAR(2)=0.45
	   !PAR(3)=500.0
	   !PAR(4)=182.5
	   !PAR(5)=0.2
	   A=PAR(1)
       B1=PAR(2)
	   d=PAR(3)
       nu=PAR(4)
       c=PAR(5)
       
 
	   
	
	  
	   u1eq = (A+sqrt(A**2-4*B1*(B1)))/(2*(B1))
       U(1)=u1eq
       !U(1) = 0
	
      U(2)=0
	   
	   
	U(3)=A/(1+u1eq**2)
	
	
	
	U(4)=0

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
