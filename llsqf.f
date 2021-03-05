!> Program for performing linear least squares fits using orthogonal decomposition of the Design (partial derivative) matrix.
!! It is designed for the data sets of modest size where it is convenient to generate and store the complete
!! partial derivative matrix prior to calling LLSQF.  If this is not the case, subroutine version LLSQFVL,
!! which generates this partial derivative array one row at a time through calls to a user-supplied subroutine,
!! should be used.
!!
!! On Entry:
!!-----------------------------------------------------------------------
!! \f$ NDATA \f$ is the number of data to be fitted (\f$ \leq MXDATA \f$).
!!
!! \f$ NPARM \f$ is the number of parameters to be varied (\f$ \leq MXPARM \f$).
!!
!!      If NPARM is <= 0, assume YD(i) = YO(i) and calculate the (RMS dimensionless deviation) = DSE from them & YU(i).
!!
!! \f$ MXDATA \f$ & \f$ MXPARM \f$ are array dimension parameters (see below). Internal array sizes currently assume MXPARM \f$ \leq \f$ 60.
!!
!! \f$ YO(i) \f$ are the \f$ NDATA \f$ 'observed' data; for iterative non-linear fits these are: \f$ [Y(obs,i) - Y(trial,i)] \f$.
!!
!! \f$ YU(i) \f$ are the uncertainties in these \f$ YO(i) \f$ values.
!!
!! \f$ DYDP(i,j) \f$ is the partial derivative array \f$ \frac{dYO(i)}{dPV(j)} \f$.
!!
!! On Exit:
!!-----------------------------------------------------------------------
!! \f$ PV(j) \f$ are the fitted parameter values; for iterative non-linear fits these are the parameter changes.
!!
!! \f$ PU(j) \f$ are 95% confidence limit uncertainties in the \f$ PV(j) \f$'s.
!!
!! \f$ PS(j) \f$ are 'parameter sensitivities' for the \f$ PV(j) \f$'s, defined such that the RMS displacement of predicted data
!! due to rounding off parameter-j by \f$ PS(j) \f$ is \f$ \leq \frac{DSE}{10NPARM} \f$.
!!
!! \f$ DSE \f$ is predicted (dimensionless) standard error of the fit.
!!
!! \f$ YD(i) \f$ is the array of differences \f$ [YO(i) - Ycalc(i)] \f$.
!! \f$ CM(j,k) \f$ is the correlation matrix obtained by normalizing variance/covariance matrix:
!!  \f[
!!      CM(j,k) = \frac{CM(j,k)}{\sqrt{CM(j,j)*CM(k,k)}}
!!  \f]
!! The squared 95% confidence limit uncertainty in a property \f$ F({PV(j)}) \f$ defined in terms of the
!! fitted parameters \f$ {PV(j)} \f$ is (where the L.H.S. involves \f$ [row]*[matrix]*[column] \f$ multiplication):
!!  \f[
!!      [D(F)]^2 = [PU(1)*\frac{dF}{dPV}(1), PU(2)*\frac{dF}{dPV}(2), ...]*[CM(j,k)]*[PU(2)*\frac{dF}{dPV}(1), PU(2)*\frac{dF}{dPV}(2), ...]
!!  \f]
!! Externally dimension:  \f$ YO \f$, \f$ YU \f$ and \f$ YD \geq NDATA \f$ (say as \f$ MXDATA \f$), \f$ PV \f$, \f$ PU \f$ and \f$ PS \geq NPARM \f$ (say as \f$ MXPARM \f$),
!! \f$ DYDP \f$ with column length \f$ MXDATA \f$ and row length \f$ \geq NPARM CM \f$ as a square matrix with column & row length \f$ MXPARM \f$.
!!
      SUBROUTINE LLSQF(NDATA,NPARM,MXDATA,MXPARM,YO,YU,DYDP,YD,PV,PU,PS,
     1                 CM,DSE)
c**  Program for performing linear least squares fits using orthogonal
c  decomposition of the Design (partial derivative) matrix.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2000  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** This version of the program is designed for the data sets of modest
c  size where it is convenient to generate and store the complete
c  partial derivative matrix prior to calling LLSQF.  If this is not the
c  case, subroutine version LLSQFVL, which generates this partial
c  derivative array one row at a time through calls to a user-supplied
c  subroutine, should be used.
c
c** On entry: NDATA  is the number of data to be fitted (.le.MXDATA)
c             NPARM  the number of parameters to be varied (.le.MXPARM)
c                 If NPARM.le.0 , assume  YD(i)=YO(i)  and calculate the
c                 (RMS dimensionless deviation)=DSE  from them & YU(i)
c             MXDATA & MXPARM are array dimension parameters (see below)
c                 Internal array sizes currently assume  MXPARM .le. 60
c             YO(i)  are the NDATA 'observed' data;  for iterative
c                  non-linear fits these are:  [Y(obs,i) - Y(trial,i)]
c             YU(i)  are the uncertainties in these YO(i) values
c             DYDP(i,j)  is the partial derivative array  dYO(i)/dPV(j)
c
c** On Exit: PV(j)  are the fitted parameter values;  for iterative
c                  non-linear fits these are the parameter changes
c            PU(j) are 95% confidence limit uncertainties in the PV(j)'s
c            PS(j) are 'parameter sensitivities' for the PV(j)'s, defined
c               such that the RMS displacement of predicted data  due to
c               rounding off parameter-j by PS(j) is .le. DSE/10*NPARM
c            DSE  is predicted (dimensionless) standard error of the fit
c            YD(i) is the array of differences  [YO(i) - Ycalc(i)]
c            CM(j,k)  is the correlation matrix obtained by normalizing
c   variance/covariance matrix:  CM(j,k) = CM(j,k)/SQRT[CM(j,j)*CM(k,k)]
c** The squared 95% confidence limit uncertainty in a property F({PV(j)})
c  defined in terms of the fitted parameters {PV(j)} is (where the
c  L.H.S. involves  [row]*[matrix]*[column]  multiplication):
c  [D(F)]^2 = [PU(1)*dF/dPV(1), PU(2)*dF/dPV(2), ...]*[CM(j,k)]*
c                              [PU(2)*dF/dPV(1), PU(2)*dF/dPV(2), ...]
c** Externally dimension:  YO, YU and YD  .ge. NDATA (say as MXDATA),
c             PV, PU  and  PS  .ge.  NPARM (say as MXPARM),
c             DYDP  with column length MXDATA and row length .ge. NPARM
c             CM   as a square matrix with column & row length  MXPARM
c  Authors: Robert J. Le Roy  &  Michael Dulick, Department of Chemistry
c    U. of Waterloo, Waterloo, Ontario  N2L 3G1.    Version of: 07/10/00
c***********************************************************************
      INTEGER I,J,K,L,IDF,NDATA,MXDATA,NPARM,MXPARM
      REAL*8  YO(NDATA), YU(NDATA), YD(NDATA), PV(NPARM), PU(NPARM),
     1   PS(NPARM), DYDP(MXDATA,NPARM), CM(MXPARM,MXPARM), DSE,
     2   PX(60), F95(10), TFACT, S, U
      DATA F95/12.7062D0,4.3027D0,3.1824D0,2.7764D0,2.5706D0,2.4469D0,
     1  2.3646D0,2.3060D0,2.2622D0,2.2281D0/
c
      IF((NDATA.GT.MXDATA).OR.(NPARM.GT.MXPARM).OR.(NPARM.GT.60)
     1                    .OR.(NPARM.GT.NDATA)) THEN
c** If array dimensioning inadequate, print warning & then STOP
          WRITE(6,601) NDATA,MXDATA,NPARM,MXPARM
          STOP
          ENDIF
      IF(NPARM.LE.0) THEN
c** If no parameters varied - simply calculate RMS deviation = DSE
          DSE= 0.D0
          DO  I= 1,NDATA
              YD(I)= YO(I)
              DSE= DSE+ (YD(I)/YU(I))**2
              ENDDO
          DSE= DSQRT(DSE/DFLOAT(NDATA))
          RETURN
          ENDIF
c** TFACT  is 95% student t-value for (NDATA-NPARM) degrees of freedom.
c [Approximate expression for (NDATA-NPARM).GT.10 accurate to ca. 0.002]
      TFACT= 0.D0
      IF(NDATA.GT.NPARM) THEN
          IDF= NDATA-NPARM
          IF(IDF.GT.10) TFACT= 1.960D0*DEXP(1.265D0/DFLOAT(IDF))
          IF(IDF.LE.10) TFACT= F95(IDF)
        ELSE
          TFACT= 0.D0
        ENDIF
      DO  I = 1,NPARM
          PS(I) = 0.D0
          PU(I) = 0.D0
          PX(I) = 0.D0
          DO  J = 1,NPARM
              CM(I,J) = 0.D0
              ENDDO
          ENDDO
c
c** Begin by forming the Jacobian Matrix from the input partial
c  derivative matrix DYDP.  For VERY large data sets, these partial
c  derivatives may be generated inside this loop (see version LLSQFVL).
      DO  I = 1,NDATA
          S = 1.D0 / YU(I)
          U = YO(I) * S
          DO  J = 1,NPARM
              PV(J) = DYDP(I,J) * S
              ENDDO
          CALL QQROD(NPARM,MXPARM,MXPARM,CM,PV,PX,U,PS,PU)
          ENDDO
c
c** Compute the inverse of  CM
      CM(1,1) = 1.D0 / CM(1,1)
      DO  I = 2,NPARM
          L = I - 1
          DO  J = 1,L
              S = 0.D0
              DO  K = J,L
                  S = S + CM(K,I) * CM(J,K)
                  ENDDO
              CM(J,I) = -S / CM(I,I)
              ENDDO
          CM(I,I) = 1.D0 / CM(I,I)
          ENDDO
c
c** Solve for parameter values  PV(j)
      DO  I = 1,NPARM
          J = NPARM - I + 1
          PV(J) = 0.D0
          DO  K = J,NPARM
              PV(J) = PV(J) + CM(J,K) * PX(K)
              ENDDO
          ENDDO
c
c** Get (upper triangular) "dispersion Matrix" [variance-covarience
c  matrix  without the sigma^2 factor].
      DO  I = 1,NPARM
          DO  J = I,NPARM
              U = 0.D0
              DO  K = J,NPARM
                  U = U + CM(I,K) * CM(J,K)
                  ENDDO
              CM(I,J) = U
              ENDDO
          ENDDO
c** Generate core of Parameter Uncertainties  PU(j) and (symmetric)
c   correlation matrix  CM
      DO  J = 1,NPARM
          PU(J) = DSQRT(CM(J,J))
          DO  K= J,NPARM
              CM(J,K)= CM(J,K)/PU(J)
              ENDDO
          DO  K= 1,J
              CM(K,J)= CM(K,J)/PU(J)
              CM(J,K)= CM(K,J)
              ENDDO
          PX(J)= 0.d0
          ENDDO
c
c** Generate differences:   YD(i) = [YO(i) - Ycalc(i)] , standard error
c  DSE = sigma^2,  and prepare to calculate Parameter Sensitivities PS
      DSE= 0.D0
      DO  I = 1,NDATA
          S = 1.D0 / YU(I)
          U = 0.D0
          DO  J = 1,NPARM
              PX(J)= PX(J)+ (DYDP(I,J)*S)**2
              U = U + DYDP(I,J) * PV(J)
              ENDDO
          YD(I) = YO(I) - U
          DSE= DSE+ (S*YD(I))**2
          ENDDO
      IF(NDATA.GT.NPARM) THEN
          DSE= DSQRT(DSE/(NDATA-NPARM))
        ELSE
          DSE= 0.d0
        ENDIF
c** Use DSE to get final (95% confid. limit) parameter uncertainties PU
c** Calculate 'parameter sensitivities', changes in PV(j) which would
c  change predictions of input data by an RMS average of  DSE*0.1/NPARM
      U= DSE*0.1d0/DFLOAT(NPARM)
      S= DSE*TFACT
      DO  J = 1,NPARM
          PU(J)= S* PU(J)
          PS(J)= U*DSQRT(NDATA/PX(J))
          ENDDO
c
      RETURN
  601 FORMAT(/' *** Dimensioning problems in LLSQF *** (NDATA, MXDATA, N
     1PARM, MXPARM)  =  (',I5,4(' ,',I5),' )')
      END
c***********************************************************************
!> Performs orthogonal decomposition of the linear least-squares equation:
!!  \f[
!!      (J*X = F) \rightarrow (A*X = B(Transpose)*F)
!!  \f]
!! where \f$ J \f$ is the Jacobian in which the first \f$ N \f$ rows and columns are transformed
!! to the upper triangular matrix \f$ A (J = B*A) \f$, X is the independent variable vector, and F is
!! the dependent variable vector. The transformation is applied to one row of the Jacobian matrix at a time.
!!
!! Parameters:
!!-----------------------------------------------------------------------
!! \f$ N \f$ (Integer) is the dimension of \f$ A \f$ to be transformed.
!!
!! \f$ NR \f$ (Integer) is the row dimension of a declared in calling program.
!!
!! \f$ NC \f$ (Integer) is the column dimension of \f$ F \f$ declared in calling program.
!!
!! \f$ A \f$ (Real*8 array of dimensions \f$ \geq N*N \f$) upper triangular transformation matrix.
!!
!! \f$ R \f$ (Real*8 linear array of dimension \f$ \geq N \f$) row of Jacobian to be added.
!!
!! \f$ F \f$ (Real*8 linear array \f$ \geq \f$ to the row dimension of the Jacobian) transformed dependent variable matrix.
!!
!! \f$ B \f$ (Real*8) value of \f$ F \f$ that corresponds to the added Jacobian row.
!!
!! \f$ GC \f$ (Real*8 linear array \f$ \geq N \f$) givens cosine transformations.
!!
!! \f$ GS \f$ (Real*8 linear array \f$ \geq N \f$) givens sine transformations.
	SUBROUTINE QQROD(N,NR,NC,A,R,F,B,GC,GS)
C** Performs ORTHOGONAL DECOMPOSITION OF THE LINEAR LEAST-SQUARES    
C            EQUATION J * X = F TO A * X = B(TRANSPOSE) * F WHERE   
C            J IS THE JACOBIAN IN WHICH THE FIRST N ROWS AND COLUMNS
C            ARE TRANSFORMED TO THE UPPER TRIANGULAR MATRIX A      
C            (J = B * A), X IS THE INDEPENDENT VARIABLE VECTOR, AND
C            F IS THE DEPENDENT VARIABLE VECTOR. THE TRANSFORMATION
C            IS APPLIED TO ONE ROW OF THE JACOBIAN MATRIX AT A TIME.
C  PARAMETERS :                                                   
C      N   -  (INTEGER) DIMENSION OF A TO BE TRANSFORMED.        
C      NR  -  (INTEGER) ROW DIMENSION OF A DECLARED IN CALLING PROGRAM.
C      NC  -  (INTEGER) Column DIMENSION OF F DECLARED IN CALLING PROGRAM.
C      A   -  (REAL*8 ARRAY OF DIMENSIONS .GE. N*N) UPPER TRIANGULAR
C             TRANSFORMATION MATRIX.                               
C      R   -  (REAL*8 LINEAR ARRAY OF DIMENSION .GE. N) ROW OF    
C             JACOBIAN TO BE ADDED.                             
C      F   -  (REAL*8 LINEAR ARRAY .GE. TO THE ROW DIMENSION OF THE
C             JACOBIAN) TRANSFORMED DEPENDENT VARIABLE MATRIX.    
C      B   -  (REAL*8) VALUE OF F THAT CORRESPONDS TO THE ADDED  
C             JACOBIAN ROW.                                     
C     GC   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS COSINE TRANSFORMATIONS.
C     GS   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS SINE TRANSFORMATIONS. 
C--------------------------------------------------------------------
C  AUTHOR : MICHAEL DULICK, Department of Chemistry,
C           UNIVERSITY OF WATERLOO, WATERLOO, ONTARIO N2L 3G1
C--------------------------------------------------------------------
      INTEGER  I,J,K,N,NC,NR
      REAL*8 A(NR,NC), R(N), F(NR), GC(N), GS(N), B, Z(2)
      DO  I = 1,N
          Z(1) = R(I)
          J = I - 1
          DO  K = 1,J
              Z(2) = GC(K) * A(K,I) + GS(K) * Z(1)
              Z(1) = GC(K) * Z(1) - GS(K) * A(K,I)
              A(K,I) = Z(2)
              ENDDO
          GC(I) = 1.D0
          GS(I) = 0.D0
          IF((Z(1) .GT. 0.d0) .OR. (Z(1) .LT. 0.d0)) THEN
              IF(DABS(A(I,I)) .LT. DABS(Z(1))) THEN
                  Z(2) = A(I,I) / Z(1)
                  GS(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
                  GC(I) = Z(2) * GS(I)
                ELSE
                  Z(2) = Z(1) / A(I,I)
                  GC(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
                  GS(I) = Z(2) * GC(I)
                ENDIF
              A(I,I) = GC(I) * A(I,I) + GS(I) * Z(1)
              Z(2) = GC(I) * F(I) + GS(I) * B
              B = GC(I) * B - GS(I) * F(I)
              F(I) = Z(2)
              ENDIF
          ENDDO
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

