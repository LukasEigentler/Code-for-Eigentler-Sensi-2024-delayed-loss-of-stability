      ! ------------------ Equations ----------------------------------
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION dhdw,A,B1,c,d,nu,FF,H,B2,Dp,chi, gamma, FF2, H2, BB2, D2, k1
      DOUBLE PRECISION p1,p1p,p2,p2p,w,wp,al1,al2,al3,al4,al5,al6,be1,be2,be3,be4,be5,be6
      DOUBLE PRECISION reallambda,imaglambda,L,dfdu1,dfdu2,dfdw,dgdu1,dgdu2,dgdw,dhdu1,dhdu2
       p1 = U(1)
       p1p = U(2)
       w = U(3)
       wp = U(4)
       
       al1 = U(5)
       al2 = U(6)
       al3 = U(7)
       al4 = U(8)

       
       be1 = U(9)
       be2 = U(10)
       be3 = U(11)
       be4 = U(12)

	   
       
       A = PAR(1) 
       d = PAR(3)
       nu = PAR(6)
       c = PAR(7)
       reallambda = PAR(8)
       imaglambda = PAR(9)
       gamma = PAR(5)
       L = PAR(2)
	   B1 = PAR(12)

       
       dfdu1 = 2*w*p1-B1 
   
       dfdw = p1**2
       dhdu1 = -2*w*p1

       dhdw = -1 - p1**2
     
		
       F(1) = L*p1p
       F(2) = L*(-p1**2*w+B1*p1-c*p1p) 
       F(3) = L*wp
       F(4) = (L/d)*(-A+w*p1**2+w-(nu+c)*wp)
       
       F(5) = L*al2 + gamma*be1
       F(6) = L*reallambda*al1 -L*imaglambda*be1 - L*dfdu1*al1 - L*c*al2  - L*dfdw*al3 + gamma*be2
       
       F(7) = L*al4 + gamma*be3
       F(8) = (L/d)*reallambda*al3 - (L/d)*imaglambda*be3 - (L/d)*dhdu1*al1 - (L/d)*dhdw*al3 - &
	   (L/d)*(c+nu)*al4 + gamma*be4
       
       F(9) = L*be2 - gamma*al1
       F(10) = L*reallambda*be1 + L*imaglambda*al1 - L*dfdu1*be1 - L*c*be2 - L*dfdw*be3 - gamma*al2
       
       F(11) = L*be4 - gamma*al3
       F(12) = (L/d)*reallambda*be3 + (L/d)*imaglambda*al3 - (L/d)*dhdu1*be1 - (L/d)*dhdw*be3 - &
	   (L/d)*(c+nu)*be4 - gamma*al4

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      
      DOUBLE PRECISION	A,B1,d,chi,nu,c,reallambda,imaglambda,gamma,FF,H,B2,DP, FF2, H2, BB2, D2, k1
	  
	   !     A = PAR(1) 
       !d = PAR(3)
       !chi = PAR(4)
       !nu = PAR(6)
       !c = PAR(7)
       !reallambda = PAR(8)
       !imaglambda = PAR(9)
       !gamma = PAR(5)
       !L = PAR(2)
	   !B1 = PAR(12)
	   !FF2 = PAR(13)
	   !H2 = PAR(14)
	   !BB2 = PAR(15)
	   !D2 = PAR(16)
	   
 
      PAR(5) = 0.0

 
	  
	  OPEN(17,FILE='eigenvalues.dat',STATUS='OLD')
      REWIND(17)
      READ(17,*)PAR(8),PAR(9)
      CLOSE(17)
	  
	  OPEN(18,FILE='../parameters.dat',STATUS='OLD')
      REWIND(18)
      READ(18,*)PAR(12),PAR(7),PAR(3),PAR(6)
      CLOSE(18)
      
	  !parameters_start = [B1, b2, c, chi, d, D2, F2, H2, nu];
	  
	  OPEN(19,FILE='stab_data/period_rain.dat',STATUS='OLD')
      REWIND(19)
      READ(19,*)PAR(2), PAR(1), PAR(7)
      CLOSE(19)
      
	   
 
      
      END SUBROUTINE STPNT

      ! ---------------------- Boundary conditions --------------------------------------
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM,ICP(*),NBC,IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*),U0(NDIM),U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

      INTEGER j
	   do j=1,NDIM
			FB(j)=U0(j) - U1(j)
	   end do
       

      END SUBROUTINE BCND

      ! ------------- Integral contraints -------------------------
      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
      DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)
      INTEGER j
      ! Phase condition on V (U(1) to U(6) in AUTO)
       FI(1) = 0.0
       do j = 1,NDIM/3
	    FI(1) = FI(1) + UPOLD(j) * (UOLD(j) - U(j))
       end do
       
       ! Normalisation of eigenvector w 
       FI(2) = -1.0
       do j = NDIM/3 + 1, 2*NDIM/3
	    FI(2) = FI(2) + U(j)*U(j) + U(j+NDIM/3)*U(j+NDIM/3)
       end do
       
       ! Phase shift condition on eigenvector w
       FI(3) = 0.0
       do j = NDIM/3 + 1, 2*NDIM/3
	    FI(3) = FI(3) + UOLD(j)*U(j+NDIM/3) - UOLD(j+NDIM/3)*U(j)
       end do

      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS	  
      END SUBROUTINE PVLS
