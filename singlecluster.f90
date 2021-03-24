
SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

! Evaluates the algebraic equations or ODE right hand side

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      U      :   State variables
!      ICP    :   Array indicating the free parameter(s)
!      PAR    :   Equation parameters

! Values to be returned :
!      F      :   ODE right hand side values

! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      INTEGER I, J, K
      INTEGER NE, NI
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION G, Mee, Mei, Mie, Mii
      DOUBLE PRECISION S
      ! REAL, dimension(NDIM, NDIM) :: H

      G   = PAR(1)
      Mee = PAR(2)
      Mei = PAR(3)
      Mie = PAR(4)
      Mii = PAR(5)

      NE = NDIM * 0.8
      NI = NDIM - NE

      ! make matrix
      DO I = 1,NE
            DO J = 1,NE
                  H(I,J) = Mee
            END DO
            DO J = NE+1,NDIM
                  H(I,J) = Mei
            END DO
            H(I,I) = 0
      END DO
      DO I = NE+1,NDIM
            DO J = 1,NE
                  H(I,J) = Mie
            END DO
            DO J = NE+1,NDIM
                  H(I,J) = Mii
            END DO
            H(I,I) = 0
      END DO

      ! DO I=1,NE
      !       S = 0
      !       DO J = 1,NE
      !             IF (J /= I) THEN
      !                   S = S + Mee*TANH( G * U(J) )
      !                   ! PRINT *,J,I
      !             END IF
      !       END DO
      !       DO J = NE+1,NDIM
      !             S = S + Mei*TANH( G * U(J) )
      !       END DO
      !       S = S / SQRT(REAL(NDIM))
      !       F(I) = -U(I) + S 
      ! END DO

      ! DO I=NE+1,NDIM
      !       S = 0
      !       DO J = 1,NE
      !             S = S + Mie*TANH( G * U(J) )
      !       END DO
      !       DO J = NE+1,NDIM
      !             IF (J /= I) THEN
      !                   S = S + Mii*TANH( G * U(J) )
      !             END IF
      !       END DO
      !       S = S / SQRT(REAL(NDIM))
      !       F(I) = -U(I) + S 
      ! END DO

      DO I=1,NDIM
            S = 0
            DO J = 1,NDIM
                  S = S + H(I,J)*TANH( G * U(J) )
            END DO
            S = S / SQRT(REAL(NDIM))
            F(I) = -U(I) + S 
      END DO

END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

! Input arguments :
!      NDIM   :   Dimension of the ODE system 

! Values to be returned :
!      U      :   A starting solution vector
!      PAR    :   The corresponding equation-parameter values

      IMPLICIT NONE
      INTEGER NDIM
      INTEGER LOFFSET, ROFFSET
      INTEGER I
      DOUBLE PRECISION U(NDIM), PAR(*), T
      DOUBLE PRECISION D, K, W, PHI
      DOUBLE PRECISION G, Mee, Mei, Mie, Mii

! Initialize the equation parameters
      G = 5
      Mee = 0.25
      Mei = 0.25
      Mie = -1
      Mii = -1

      ! initialize to 0
      DO I = 1,NDIM
          U(I) = 0
      END DO
      U(NDIM) = 0.1204
      U(NDIM-1) = 0.1204
      U(NDIM-2) = -0.1204
      U(NDIM-3) = -0.1204
      
      PAR(1) = G
      PAR(2) = Mee
      PAR(3) = Mei
      PAR(4) = Mie
      PAR(5) = Mii

END SUBROUTINE STPNT

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! The following subroutines are not used here,
! but they must be supplied as dummy routines

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------