
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
    INTEGER NE, NI, Nc, p
    DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(NDIM,NDIM), DFDP(NDIM,*)
    DOUBLE PRECISION G, Mee, Mei, Mie, Mii
    DOUBLE PRECISION frac, S
    DOUBLE PRECISION H(NDIM,NDIM)

    G   = PAR(1)
    Mee = PAR(2)
    Mei = PAR(3)
    Mie = PAR(4)
    Mii = PAR(5)

    frac = 0.8
    Nc = 4
    p = FLOOR( NDIM * frac / Nc)
    NE = Nc * p
    NI = NDIM - NE

    ! make matrix
    DO I = 1, NDIM
        DO J = 1, NDIM
            H(I,J) = 0
        END DO
    END DO

    IF (Nc == 1) THEN
        DO I = 1,NE
            DO J = 1,NE
                H(I,J) = Mee
            END DO
            DO J = NE+1,NDIM
                H(I,J) = Mei
            END DO
            H(I,I) = 0
        END DO
    ELSE
        DO I = 1,Nc
            DO J = 1,p
                DO K = 1,p
                    H( (I-1)*p + J, (I-1)*p + K) = Mee
                END DO
                H( (I-1)*p + J, (I-1)*p + J) = 0
            END DO
        END DO
        DO I = 1,NE
            DO J = NE+1,NDIM
                H(I,J) = Mei
            END DO
        END DO
    END IF

    DO I = NE+1,NDIM
        DO J = 1,NE
                H(I,J) = Mie
        END DO
        DO J = NE+1,NDIM
                H(I,J) = Mii
        END DO
        H(I,I) = 0
    END DO

    DO I=1,NDIM
        S = 0
        DO J = 1,NDIM
            S = S + H(I,J)*TANH( G * U(J) ) / SQRT(REAL(NDIM))
        END DO
        F(I) = -U(I) + S 
    END DO

    ! ! Jacobian
    IF(IJAC.EQ.0) RETURN 
    DO I = 1,NDIM
        DO J = 1, NDIM
            IF (I == J) THEN
                DFDU(I,J) = -1
            ELSE
                DFDU(I, J) = (G / SQRT(REAL(NDIM))) * H(I,J) / ( COSH(G*U(J))**2 ) 
            END IF
        END DO
    END DO

    IF(IJAC.EQ.1) RETURN 
    DO I=1,NDIM
        S = 0
        DO J = 1,NDIM
            S = S + H(I,J) * U(J) / ( COSH(G*U(J))**2 )
        END DO
        DFDP(I,1) = S / SQRT(REAL(NDIM))
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
    DOUBLE PRECISION G, Mee, Mei, Mie, Mii
    DOUBLE PRECISION frac, a, Gstar, Xstart
    INTEGER NE, NI, Nc, p

! Initialize the equation parameters
    frac = 0.8
    a = frac/(1-frac)

    Nc = 4
    p = FLOOR( NDIM * frac / Nc)
    NE = Nc * p
    NI = NDIM - NE

    Mei = -1
    Mii = -1
    Mie = -((1-frac)/frac)*Mii
    Mee = -((1-frac)/frac)*Mei
    ! Mee = -Nc*((1-frac)/frac)*Mei

    ! initialize to 0
    DO I = 1,NDIM
        U(I) = 0
    END DO

    ! IF (Nc == 1) THEN
        Gstar = SQRT( REAL(NDIM) ) / (a*Mee)
        G = Gstar + 0.1;
        Xstart = SQRT( 3*(G - Gstar) / (Gstar**3) )
        ! G = 4;
        DO I = 1, NI/2
            U(NE + I) = Xstart
            U(NE + NI/2 + I) = -Xstart
        END DO
    ! ELSE
    !     ! Gstar = SQRT( REAL(NDIM) / ( (p-1)*Mee) )
    !     G = 4
    !     ! Xstart = SQRT( 3*(G - Gstar) / (Gstar**3) ) 
    !     ! Xstart = 0.6642
    !     ! DO I = 1, NE/2
    !     !     U(I) = Xstart
    !     !     U(NE/2 + I) = -Xstart
    !     ! END DO
    !     ! Xstart = 0.1204
    !     ! DO I = 1, NI/2
    !     !     U(NE + I) = Xstart
    !     !     U(NE + NI/2 + I) = -Xstart
    !     ! END DO
    ! END IF

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