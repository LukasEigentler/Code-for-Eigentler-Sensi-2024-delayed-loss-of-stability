      ! ------------------ Equations ----------------------------------
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION dhdw,A,B1,c,d,nu,FF,H,B2,Dp,chi, gamma, FF2, H2, BB2, D2,k1
      DOUBLE PRECISION p1,p1p,p2,p2p,w,wp,al1,al2,al3,al4,al5,al6,be1,be2,be3,be4,be5,be6
      DOUBLE PRECISION reallambda,imaglambda,L,dfdu1,dfdu2,dfdw,dgdu1,dgdu2,dgdw,dhdu1,dhdu2
	  DOUBLE PRECISION be1gam, be2gam, be3gam, be4gam, be5gam, be6gam, al1gam2, al2gam2, al3gam2
	  DOUBLE PRECISION al4gam2, al5gam2, al6gam2, imaglambda_gamma, reallambda_gamma2, dummy
       p1 = U(1)
       p1p = U(2)  
       w = U(3)
       wp = U(4)
       
       al1 = U(5)
       al2 = U(6)
       al5 = U(7)
       al6 = U(8)
       
       be1 = U(9)
       be2 = U(10)
       be5 = U(11)
       be6 = U(12)
	   
	   be1gam = U(13)
	   be2gam = U(14)
	   be5gam = U(15)
	   be6gam = U(16)
	   
	   al1gam2 = U(17)
	   al2gam2 = U(18)
	   al5gam2 = U(19)
	   al6gam2 = U(20)

	   
       
       A = PAR(1) 
       d = PAR(3)
       nu = PAR(6)
       c = PAR(7)
       reallambda = PAR(8)
       imaglambda = PAR(9)
       gamma = PAR(5)
       L = PAR(2)
	   B1 = PAR(12)
	   
	   
	   imaglambda_gamma = PAR(17)
	   reallambda_gamma2 = PAR(18)
	   dummy = PAR(19)
	   
	   
	   
       
       dfdu1 = 2*w*p1-B1
       
       dfdw = p1**2
      
       dhdu1 = -2*w*p1
  
       dhdw = -1 - p1**2
     
		
       F(1) = L*p1p
       F(2) = L*(-p1**2*w+B1*p1-c*p1p ) 
       
       F(3) = L*wp
       F(4) = (L/d)*(-A+w*p1**2+w-(nu+c)*wp)
       
       F(5) = L*al2 + gamma*be1
       F(6) = L*reallambda*al1 -L*imaglambda*be1 - L*dfdu1*al1 - L*c*al2 - L*dfdw*al5 + gamma*be2
       F(7) = L*al6 + gamma*be5
       F(8) = (L/d)*reallambda*al5 - (L/d)*imaglambda*be5 - (L/d)*dhdu1*al1 - (L/d)*dhdw*al5 - &
	   (L/d)*(c+nu)*al6 + gamma*be6
       
       F(9) = L*be2 - gamma*al1
       F(10) = L*reallambda*be1 + L*imaglambda*al1 - L*dfdu1*be1 - L*c*be2 - L*dfdw*be5 - gamma*al2
       F(11) = L*be6 - gamma*al5
       F(12) = (L/d)*reallambda*be5 + (L/d)*imaglambda*al5 - (L/d)*dhdu1*be1  - (L/d)*dhdw*be5 - &
	   (L/d)*(c+nu)*be6 - gamma*al6
	   
	   F(13) = L*(be2gam - 1/L*al1)
	   F(14) = L*(-dfdu1*be1gam - c*be2gam  - dfdw*be5gam + imaglambda_gamma*al1 - 1/L*al2)
	   F(15) = L*(be6gam - 1/L*al5)
	   F(16) = L*(1/d*(-dhdu1*be1gam  - dhdw*be5gam - (c+nu)*be6gam) + 1/d*imaglambda_gamma*al5 - 1/L*al6)
	   
	   F(17) = L*(al2gam2 + 2/L*be1gam)
	   F(18) = L*(-dfdu1*al1gam2 - c*al2gam2  - dfdw*al5gam2 - 2*imaglambda_gamma*be1gam + 2/L*be2gam + & 
	   reallambda_gamma2*al1)
	   F(19) = L*(al6gam2 + 2/L*be5gam)
	   F(20) = L*(1/d*(-dhdu1*al1gam2  - dhdw*al5gam2 - (c+nu)*al6gam2) - 2/d*imaglambda_gamma*be5gam + & 
	   2/L*be6gam + reallambda_gamma2/d*al5)
	   
	   

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      
      DOUBLE PRECISION	A,B1,d,chi,nu,c,reallambda,imaglambda,gamma,FF,H,B2,DP, FF2, H2, BB2, D2,k1
	  
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
	   
 
      PAR(5) = 0.0d0
	  PAR(17) = 0.0d0
	  PAR(18) = 0.0d0
	  PAR(19) = 0.0d0
	  !PAR(8) = 0.0d0
	  !PAR(9) = 0.0d0

 
	  
	  OPEN(30,FILE='eigenvalues.dat',STATUS='OLD')
      REWIND(30)
      READ(30,*)PAR(8),PAR(9)
      CLOSE(30)
	  
	  OPEN(31,FILE='../../parameters.dat',STATUS='OLD')
      REWIND(31)
      READ(31,*)PAR(12),PAR(7),PAR(3),PAR(6)
      CLOSE(31)
      
	  !parameters_start = [B1, b2, c, chi, d, D2, F2, H2, nu];
	  
	  OPEN(32,FILE='period_rain.dat',STATUS='OLD')
      REWIND(32)
      READ(32,*)PAR(2), PAR(1)
      CLOSE(32)
      
	   
      FF = 1-chi*(1-FF2)
      H = 1-chi*(1-H2)
      B2 = B1-chi*(B1-0.004)
      Dp = 1-chi*(1-D2)
      
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
       do j = 1,NDIM/5
	    FI(1) = FI(1) + UPOLD(j) * (UOLD(j) - U(j))
       end do
       
       ! Normalisation of eigenvector w 
       FI(2) = -1.0
       do j = NDIM/5 + 1, 2*NDIM/5
	    FI(2) = FI(2) + U(j)*U(j) + U(j+NDIM/5)*U(j+NDIM/5)
       end do
       
       ! Phase shift condition on eigenvector w
       FI(3) = 0.0
       do j = NDIM/5 + 1, 2*NDIM/5
	    FI(3) = FI(3) + UOLD(j)*U(j+NDIM/5) - UOLD(j+NDIM/5)*U(j)
       end do
	   
	   ! uniqueness eigenfunction
	   FI(4) = 0.0
	   FI(5) = 0.0
	   do j = 1,NDIM/5
		FI(4) = FI(4) + U(NDIM/5+j)*U(3*NDIM/5+j)
		FI(5) = FI(5) + U(NDIM/5+j)*U(4*NDIM/5+j)
	   end do

      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS	  
      END SUBROUTINE PVLS
