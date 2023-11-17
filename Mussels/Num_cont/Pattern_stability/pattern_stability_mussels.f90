      ! ------------------ Equations ----------------------------------
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION dhds,alpha,beta,delta,eta,theta,c,d,nu,gamma,ff,g,h
      DOUBLE PRECISION m,n,a,s,v,al1,al2,al3,al4,al5,al6,be1,be2,be3,be4,be5,be6
      DOUBLE PRECISION relam,imlam,L,dfda,dfdm,dfds,dgda,dgdm,dgds,dhdm,dhda
       m=U(1)
       n=U(2)
       a=U(3) 
	   s=U(4)
	   v=U(5)
	   
       
       
       al1 = U(6)
       al1 = U(6)
       al2 = U(7)
       al3 = U(8)
       al4 = U(9)
	   al5 = U(10)

       
       be1 = U(11)
       be2 = U(12)
       be3 = U(13)
       be4 = U(14)
	   be5 = U(15)

	   
	   alpha=PAR(1)
       beta=PAR(2)
       d=PAR(3)
       nu=PAR(4)
       c=PAR(5)
	   delta=PAR(6)
	   eta=PAR(7)
	   theta=PAR(8)
       
      
       relam = PAR(9)
       imlam = PAR(10)
       gamma = PAR(12)
       L = PAR(13)
	   
 
	   
       dfdm = (a*delta*(eta+s))/(s+1.0D0)-1.0D0
	   dfda = (delta*m*(eta+s))/(s+1.0D0)
	   dfds = (a*delta*m)/(s+1.0D0)-a*delta*m*(eta+s)*1.0D0/(s+1.0D0)**2
	   dgdm = -(a*beta*(eta+s))/(s+1.0D0)
	   dgda = -alpha-(beta*m*(eta+s))/(s+1.0D0)
	   dgds = -(a*beta*m)/(s+1.0D0)+a*beta*m*(eta+s)*1.0D0/(s+1.0D0)**2
	   dhdm = 1.0D0
	   dhda = 0.0D0
	   dhds = -theta
 
     
		
       F(1) = L*n
       F(2) = L*(m-c*n-(a*delta*m*(eta+s))/(s+1.0D0))
	   F(3) = (L*(a*alpha+(a*beta*m*(eta+s))/(s+1.0D0)-1.0D0))/(c+nu)
       F(4) = L*v
       F(5) = -(L*(m+c*v-s*theta))/D
       
       F(6) = L*al2+be1*gamma
       F(7) = be2*gamma-L*al2*c-L*al3*dfda-L*al1*dfdm-L*al4*dfds &
     -L*be1*imlam+L*al1*relam
	   
	   F(8) = be3*gamma-(L*al3*dgda)/(c+nu)-(L*al1*dgdm)/(c+&
     nu)-(L*al4*dgds)/(c+nu)+(L*al3*relam)/(c+nu)-(L*be3*imlam)/(c+nu)
       
       F(9) = L*al5+be4*gamma
       F(10) = be5*gamma-(L*al5*c)/D-(L*al3*dhda)/D-(L*al1*dhdm)/D- &
	   (L*al4*dhds)/D-(L*be4*imlam)/D+(L*al4*relam)/D
       
       F(11) = L*be2-al1*gamma
       F(12) = -al2*gamma-L*be2*c-L*be3*dfda-L*be1*dfdm-L*be4*dfds+&
	   L*al1*imlam+L*be1*relam
	   
	   F(13) = -al3*gamma-(L*be3*dgda)/(c+nu)-(L*be1*dgdm)/(c+nu)-&
	   (L*be4*dgds)/(c+nu)+(L*al3*imlam)/(c+nu)+(L*be3*relam)/(c+nu)
       
       F(14) = L*be5-al4*gamma
       F(15) = -al5*gamma-(L*be5*c)/D-(L*be3*dhda)/D-(L*be1*dhdm)/D-&
	   (L*be4*dhds)/D+(L*al4*imlam)/D+(L*be4*relam)/D

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      
      
	  
	
	   
 
      PAR(12) = 0.0 !gamma

 
	  
	  OPEN(17,FILE='eigenvalues.dat',STATUS='OLD')
      REWIND(17)
      READ(17,*)PAR(9),PAR(10)
      CLOSE(17)
	  
	  OPEN(18,FILE='../parameters.dat',STATUS='OLD')
      REWIND(18)
      READ(18,*)PAR(1),PAR(2),PAR(7),PAR(8),PAR(50),PAR(3),PAR(4)
      CLOSE(18)
      
	  ![alpha, beta, eta, theta, c, d, nu];
	  
	  OPEN(19,FILE='stab_data/period_delta.dat',STATUS='OLD')
      REWIND(19)
      READ(19,*)PAR(13), PAR(6), PAR(5)
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
