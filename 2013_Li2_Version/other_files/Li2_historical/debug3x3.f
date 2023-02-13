      PROGRAM debug3x3 

      REAL*8   H(3,3),DM1(3,3),DM3(3,3),DM5(3,3),DR(3,3),DDe(3,3)
c      REAL*8   Q(3,3),CmVAL(MXMLR)
      REAL*8   Q(3,3)
      REAL*8   DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),
     1                    DEIGDe(1,1)
      DOUBLE PRECISION W(3)
      DOUBLE PRECISION EIGVEC(3,1)
      REAL*8           RDIST,RE,De,RE3,RE6,C6adj,C9adj,DELTAE,ULR 
      REAL*8           RESID(3,1)

      REAL*8           M1,M2,M3,M5,DIFFER

      INTEGER          I,J,L,K,S
      INTEGER          ISTATE,MXMLR
      
      DOUBLE PRECISION MODULUS
      REAL*8           Z



c M2 = C6sigma, M3 = C6 adj

      De = 7091.6D0
      M1 = 357590.D0
      M2 = 14918100.D0 
      M3 = M2+(M1**2)/(4*De)
      M5 = 370200000.D0 
      DELTAE = 0.335338d0 



      DO 91 S = 1, 7997
           RDIST = 1.5D0+0.5D0*S

            CALL ZHEEVJ3(H, Q, W,RDIST,M1,M2,M3,M5,De,DELTAE) 
  
            DR(1,1) = dble((6 *M1*RDIST**5*De + 12 * RDIST **
     1        2 * M2 * De + 3 * RDIST ** 2 * M1 ** 2 + 16 * M5 * De) /
     2        RDIST ** 9 / De) / 0.6D1
            DR(1,2) = -dble(6 * M1 * RDIST ** 5 * De + 12 * RDIST ** 2
     1        * M2 * De + 3 * RDIST ** 2 * M1 ** 2 + 16 * M5 * De) *
     2         sqrt(0.2D1) / dble(RDIST ** 9) /
     3         dble(De) / 0.6D1
            DR(1,3) = -sqrt(0.6D1) * dble(M1) / dble(RDIST ** 4) / 0.2D1
            DR(2,1) = DR(1,2) 
            DR(2,2) = dble((6 * M1 * RDIST ** 5 * De + 12 * RDIST ** 2 *
     1          M2 * De + 3 * RDIST ** 2 * M1 ** 2 + 16 * M5 * De) /
     2          RDIST ** 9 / De) / 0.3D1
            DR(2,3) = -dble(M1) * sqrt(0.3D1) / dble(RDIST ** 4) / 0.2D1
            DR(3,1) = -sqrt(0.6D1) * dble(M1) / dble(RDIST ** 4) / 0.2D1
            DR(3,2) = -dble(M1) * sqrt(0.3D1) / dble(RDIST ** 4) / 0.2D1
            DR(3,3) = 0.d0


            DM1(1,1) = -dble((2 * RDIST ** 3 * De + M1) / RDIST ** 6 /
     1       De) / 0.6D1
            DM1(1,2) = dble(2 * RDIST ** 3 * De + M1) / dble(RDIST ** 6)
     1       / dble(De) * sqrt(0.2D1) / 0.6D1
            DM1(1,3) = sqrt(0.6D1) / dble(RDIST ** 3) / 0.6D1
            DM1(2,1) = dble(2 * RDIST ** 3 * De + M1) / dble(RDIST ** 6)
     1       / dble(De) * sqrt(0.2D1) / 0.6D1
            DM1(2,2) = -dble((2 * RDIST ** 3 *De + M1) / RDIST ** 6 / De
     1       ) / 0.3D1
            DM1(2,3) = 0.1D1 / dble(RDIST ** 3) * sqrt(0.3D1) / 0.6D1
            DM1(3,1) = sqrt(0.6D1) / dble(RDIST ** 3) / 0.6D1
            DM1(3,2) = 0.1D1 / dble(RDIST ** 3) * sqrt(0.3D1) / 0.6D1
            DM1(3,3) = 0.d0
 
          DM3(1,1)=-1.d0/(3.d0*RDIST**6)
          DM3(1,2)=-SQRT(2.d0)*DM3(1,1)
          DM3(1,3)=0.D0
          DM3(2,1)=DM3(1,2)
          DM3(2,2)=2.d0*DM3(1,1)
          DM3(2,3)=0.D0
          DM3(3,1)=DM3(1,3)
          DM3(3,2)=DM3(2,3)
          DM3(3,3)=0.D0

          DM5(1,1)=DM3(1,1)/(RDIST**2)
          DM5(1,2)=DM3(1,2)/(RDIST**2)
          DM5(1,3)=0.D0
          DM5(2,1)=DM3(1,2)
          DM5(2,2)=DM3(2,2)/(RDIST**2)
          DM5(2,3)=0.D0
          DM5(3,1)=DM5(1,3)
          DM5(3,2)=DM5(2,3)
          DM5(3,3)=0.D0

          DDe(1,1)=M1**2/(12.D0*RDIST**6*De**2)
          DDe(1,2)=-SQRT(2.D0)*DDe(1,1)
          DDe(1,3)=0.D0
          DDe(2,1)=DDe(1,2)
          DDe(2,2)=2.D0*DDe(1,1)
          DDe(2,3)=0.d0
          DDe(3,1)=DDe(1,3)
          DDe(3,2)=DDe(2,3)
          DDe(3,3)=0.D0



         L=1
           DO J=2,3
           IF (W(J) .LT. W(L)) THEN
         L=J
       ENDIF
            ENDDO

   20    ULR=-W(L)

      DO I=1,3    
          EIGVEC(I,1) = Q(I,L)
      ENDDO


      WRITE(6,601) RDIST ,ULR,EIGVEC 
  601 FORMAT(6D16.7)
   
      DIFFER = MAXVAL(MATMUL(H,EIGVEC) + EIGVEC*ULR)

      WRITE(7,602) RDIST,ULR,DIFFER
  602 FORMAT(3D16.7)

   30       DEIGM1 = - MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM1,EIGVEC))   
   40       DEIGR  = - MATMUL(TRANSPOSE(EIGVEC),MATMUL(DR,EIGVEC))
   50       DEIGDe = - MATMUL(TRANSPOSE(EIGVEC),MATMUL(DDe,EIGVEC))

      WRITE(9,605) RDIST,DEIGM1,DEIGR,DEIGDe
  605 FORMAT(4D26.16)

         Modulus = SQRABS(Z) 

   91 CONTINUE


      CONTAINS

* ----------------------------------------------------------------------------
      SUBROUTINE ZHEEVJ3(H, Q, W,RDIST,M1,M2,M3,M5,De,DELTAE)
* ----------------------------------------------------------------------------
      REAL*8           RDIST,M1,M2,M3,M5,De,DELTAE

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
 
       H(1,1)=-(1.D0/3.D0)*(M1/(RDIST**3)+M3/(RDIST**6)+M5/(RDIST**8))
       H(1,2)=-(SQRT(2.D0))*H(1,1)
       H(2,1)=H(1,2)
       H(1,3)=(SQRT(6.D0)/6.D0)*M1/(RDIST)**3
       H(3,1)=H(1,3)
       H(2,2)=2*H(1,1) + DELTAE
       H(2,3)=(SQRT(3.D0)/6)*M1/RDIST**3
       H(3,2)=H(2,3)
       H(3,3)=DELTAE

      WRITE(8,603) H 
  603 FORMAT(9D26.16)



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

      END PROGRAM  debug3x3
