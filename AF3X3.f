      SUBROUTINE AF3X3(RDIST,CmVAL,RE,RE3,RE6,C6adj,C9adj,ULR,
     1 DEIGM1,DEIGR,MXMLR)

      REAL*8           H(3,3),DM1(3,3),DR(3,3),Q(3,3),CmVAL(MXMLR)
      REAL*8           DEIGM1(1,1),DEIGR(1,1)		
      DOUBLE PRECISION W(3)
      DOUBLE PRECISION EIGVEC(3,1)	
      REAL*8           RDIST,RE,RE3,RE6,C6adj,C9adj,DELTAE,ULR 
      REAL*8           RESID(3,1)

      REAL*8           M1,M2

      INTEGER          I,J,L,K
      INTEGER          ISTATE,MXMLR
      
      DOUBLE PRECISION MODULUS
      REAL*8           Z

      M1 = CmVAL(1)
      M2 = CmVAL(3)
      DELTAE = CmVAL(2)

c WRITE(25,*) 'Variables = "r", "U(r)","U(r)-U(r)^2/(4De)" ' 
c WRITE(25,*) 'zone T = "U(r)"'
 
	    DM1(1,1)=-1.D0/(3.D0*RDIST**3)
	    DM1(1,2)=SQRT(2.D0)/(3.D0*RDIST**3)
	    DM1(1,3)=SQRT(6.D0)/(6.D0*RDIST**3)
	    DM1(2,1)=DM1(1,2)
	    DM1(2,2)=-2.D0/(3.D0*RDIST**3)
	    DM1(2,3)=SQRT(3.D0)/(6.D0*RDIST**3)
	    DM1(3,1)=DM1(1,3)
	    DM1(3,2)=DM1(2,3)
	    DM1(3,3)=0.D0

          DR(1,1)=M1/RDIST**4
          DR(1,2)=-SQRT(2.D0)*M1/(RDIST**4)
	    DR(1,3)=-SQRT(6.D0)*(M1/(2.D0*RDIST**4))
          DR(2,1)=DR(1,2)
          DR(2,2)=2.D0*M1/(RDIST**4)
	    DR(2,3)=-SQRT(3.D0)*M1/(2.D0*RDIST**4)
          DR(1,3)=DR(3,1)
          DR(2,3)=DR(3,2)
          DR(3,3)=0.D0

            CALL ZHEEVJ3(H,Q,W,RDIST,M1,DELTAE)
	     
	    L=1
	    DO J=2,3
	       IF (W(J) .LT. W(L)) THEN
		   L=J
	       ENDIF
            ENDDO	

   20	    ULR=-W(L)

	    DO I=1,3	    
		EIGVEC(I,1) = Q(I,L)
	    ENDDO	

   30       DEIGM1 =- MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM1,EIGVEC))           
   40       DEIGR =- MATMUL(TRANSPOSE(EIGVEC),MATMUL(DR,EIGVEC))

c            WRITE(25,600) RDIST ,ULR 
c  600 FORMAT(2D16.7)
c
	    WRITE(26,601) RDIST , DEIGM1, DEIGR 
  601 FORMAT(3D16.7)	

         Modulus = SQRABS(Z) 

      RETURN

      CONTAINS

* ----------------------------------------------------------------------------
      SUBROUTINE ZHEEVJ3(H, Q, W,RDIST,M1,DELTAE)
* ----------------------------------------------------------------------------
      REAL*8           RDIST,M1,DELTAE

      REAL*8           H(3,3),Q(3,3)
      DOUBLE PRECISION W(3)

      INTEGER          N
      PARAMETER        ( N = 3 )

      DOUBLE PRECISION SD, SO
      REAL*8           S, T
      DOUBLE PRECISION C
      DOUBLE PRECISION G, B, Z
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

      DOUBLE PRECISION FUNCTION SQRABS

* initialize matrix to 0.d0

        DO I=1,3
                    H(I,I)=0.0D0
                    ENDDO

* initialize matrix elements
 
            H(1,1)=-M1/(3.D0*(RDIST)**3)
            H(1,2)=(SQRT(2.D0)/3.D0)*M1/(RDIST)**3
            H(2,1)=H(1,2)
            H(1,3)=(SQRT(6.D0)/6.D0)*M1/(RDIST)**3
            H(3,1)=H(1,3)
            H(2,2)=-2.D0/3.D0*M1/RDIST**3 + DELTAE
            H(2,3)=(SQRT(3.D0)/6)*M1/RDIST**3
            H(3,2)=H(2,3)
            H(3,3)=DELTAE


*     Initialize Q to the identitity matrix
*     --- This loop can be omitted if only the eigenvalues are desired ---
      DO  X = 1, N
          Q(X,X) = 1.0D0
          DO  Y = 1, X-1
              Q(X, Y) = 0.0D0
              Q(Y, X) = 0.0D0
              ENDDO
          ENDDO

*     Initialize W to diag(A)
      DO  X = 1, N
          W(X) = DREAL(H(X, X))
          ENDDO

*     Calculate SQR(tr(A))
      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2

*     Main iteration loop
      DO 40 I = 1, 50
*       Test for convergence
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(DREAL(H(X, Y)))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

*       Do sweep
        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(DREAL(H(X, Y))) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X))
     $                    .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              H(X, Y) = 0.0D0
            ELSE IF (ABS(DREAL(H(X, Y)))
     $                 .GT. THRESH) THEN
*             Calculate Jacobi transformation
              B = W(Y) - W(X)
              IF ( ABS(B) + G .EQ. ABS(B) ) THEN
                T = H(X, Y) / B
              ELSE
                IF (B .LE. 0.0D0) THEN
                  T = -2.0D0 * H(X, Y)
     $                    /(SQRT(B**2 + 4.0D0 * SQRABS(H(X, Y))) - B)
                ELSE IF (B .EQ. 0.0D0) THEN
                  T = H(X, Y) * (1.0D0 / ABS(H(X, Y)))
                ELSE
                  T = 2.0D0 * H(X, Y)
     $                    /(SQRT(B**2 + 4.0D0 * SQRABS(H(X, Y))) + B)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + SQRABS(T) )
              S = T * C
              Z = DREAL(T * (H(X, Y)))

*             Apply Jacobi transformation
              H(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = H(R, X)
                H(R, X) = C * T - (S) * H(R, Y)
                H(R, Y) = S * T + C * H(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = H(X, R)
                H(X, R) = C * T - S * (H(R, Y))
                H(R, Y) = S * (T) + C * H(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = H(X, R)
                H(X, R) = C * T - S * H(Y, R)
                H(Y, R) = (S) * T + C * H(Y, R)
   90         CONTINUE

*             Update eigenvectors
*             --- This loop can be omitted if only the eigenvalues are desired ---
              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - (S) * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

      PRINT *, "ZHEEVJ3: No convergence."

      END SUBROUTINE ZHEEVJ3

* ----------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION SQRABS(Z)
* ----------------------------------------------------------------------------
* Calculates the squared absolute value of a complex number Z
* ----------------------------------------------------------------------------
*     .. Parameters ..
      REAL*8 Z

      SQRABS = DREAL(Z)**2
      RETURN

      END FUNCTION SQRABS

      END SUBROUTINE AF3X3 
