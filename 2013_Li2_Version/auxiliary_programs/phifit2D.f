c***********************************************************************
c************  Program  2D_PHIFIT_1.1  dated  27 February 2007  ********
c***********************************************************************
c* Program to fit read-in potential fx. values {THE(i),RTP(i),VTP(i)} to
c  an MLRp(NS,NL) potential function form in which the well depth D_e,
c  equilibrium distance R_e, and exponent expansion coefficients \phi_i
c  are all expressed as Legendre expansion functions of the Jacobi 
c  coordinate angle \theta.
c** See  'phiFIT' Manual in http://leroy.uwaterloo.ca/programs/  
c   for documentation on the 1-D version of this code, and the 
c   definitions of most input variables.
c***********************************************************************
      INTEGER MXDATA, MXPARM, NPMAX
      PARAMETER (MXDATA=15000, MXPARM=60, NPMAX=20)
      INTEGER IROUND,ROBUST,LPRINT,NPARM,NTP
      REAL*8 UNC,TSTPS,TSTPU,DSE
      REAL*8 CM(MXPARM,MXPARM),VTP(MXDATA),
     1  UVTP(MXDATA),YD(MXDATA),THTH,RX,VX,DR0,DD,DS,DXX,DYY 
c-----------------------------------------------------------------------
      INTEGER I,J,K,IP,IFXP(MXPARM),
     1  p,NDE,NRE,NCN,MCM,NS,NL,NPOW,NPHI(NPMAX),
     2  NPS
      REAL*8  PV(MXPARM),PU(MXPARM),PS(MXPARM)
      REAL*8  RREF,RCN2CN,RCMCN,RCM2CN,CN,PI
      REAL*8 RTP(MXDATA),THETA(MXDATA)
      COMMON /DATABLK/RTP,THETA,RREF,CN,RCN2CN,RCMCN,RCM2CN,
     1 p,NPHI,NDE,NRE,NCN,MCM,NS,NL
c-----------------------------------------------------------------------
      ROBUST= 0
      DR0= 0.5d0
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** UNC  is a reference value for the uncertainties assigned to 
c        the potential energies being fitted to.  See comments in 
c        Lines 118-129 below.
c** IROUND  integer set .ne.0  to have program round parameters 
c          (see phiFIT manual)
c** LPRINT  controls print in fitting subroutins: (see phiFIT manual)
c** NDE and NRE are the numbers of Legendre expansion coefficients used
c          to represent  De and Re
c** {p,NS,NL,Rref} define the exponent coefficient expansion and expansion
c  variable y_p(r)  (see phiFIT manual).
c=======================================================================
      READ(5,*) UNC, IROUND, LPRINT
      READ(5,*) NDE, NRE   
      READ(5,*) p, NS, NL, RREF
c=======================================================================
      NPOW=NL+1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** NPHI(i) is the number of Legendre terms associated with exponent 
c           expansion parameter \phi(i-1).
c** NCN and CN  are the power and coefficient of the leading long-range
c          inverse-power term in the potential
c** MCM  is the power of the second longest-range inverse-power term
c** RCN2CN  is the fixed ratio  C6(2)/C6(0)
c** RCMCN   is the fixed ratio  C8(0)/C6(0)
c** RCM2CN  is the fixed ratio  C8(2)/C6(0)
c=======================================================================
      READ(5,*) (NPHI(I),I=1,NPOW)
      READ(5,*) NCN, CN, MCM, RCN2CN,RCMCN,RCM2CN 
c=======================================================================
      PI=DACOS(-1.0D0)
      DO I=1,NDE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Read initial trial values of Legendre expansion coefficients for D_e
c=======================================================================
          READ(5,*) PV(I),IFXP(I)
c=======================================================================
          ENDDO
      IP=NDE
      DO I=1,NRE
          IP=IP+1  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Read initial trial values of Legendre expansion coefficients for R_e
c=======================================================================
          READ(5,*) PV(IP),IFXP(IP)
c=======================================================================
          ENDDO
      IP=NDE+NRE
      NPS=0
      DO J=1,NPOW
          DO I=1,NPHI(J)
              IP=IP+1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Read initial trial values of polynomial coefficients \phi(J,I) 
c   defining the exponent \phi(R) of the MLR potential.  Outer loop in J
c   over the polynomial order in y_p(R);  inner loop in I over the
c   Legendre order for that value of J.
c=======================================================================
              READ(5,*) PV(IP), IFXP(IP)
c=======================================================================
              ENDDO
          NPS=NPS+NPHI(J)
          ENDDO
      NPARM= NDE+ NRE+ NPS + 1
      IFXP(NPARM)= 1
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Read the (fixed) energy defining the potential asymptote
c========================================================================
      READ(5,*) PV(NPARM)
c========================================================================
c** Read the turning points to be fitted to
      K=0
      DO I=1, 9999
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now read all input potential energy function values, one per line,
c      to the end if the data file.   Use Jacobi coordinates:  
c       THTH= angle in degrees,
c       RX= centre of mass distance,  VX= potential function value
c========================================================================
          READ(5,*,end=66) THTH,RX,VX
c========================================================================
          IF(VX.LE.1000.)THEN
              K=K+1
              THETA(K)=DCOS(THTH*PI/180.)
              RTP(K)=RX
              VTP(K)=VX
              IF(VTP(K).LT.PV(NPARM))THEN
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c** For potential function values in the well region, assign a constant
c   uncertainty of 0.1*UNC
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  UVTP(K)=UNC*0.1
                ELSE
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c** For potential function values in the repulsive wall region, assign 
c   an uncertainty which is  0.1*UNC  at the crossing point, but 
c   increases with energy to a limiting ratio of 2%  of the energy above
c   the asymptote
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  UVTP(K)=(VTP(K) - PV(NPARM) + 5.)/50.
                ENDIF
              ENDIF
          ENDDO
66    NTP=K
      WRITE(6,602) NTP, NCN, CN, MCM, RCMCN
c=======================================================================
c** Now ... do direct non-linear fit to potential values ... first with
c      Re and/or VMIN and De fixed, and then freeing them up too ...
c=======================================================================
      CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,LPRINT,IFXP,
     1                      VTP,UVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
      WRITE(6,610)DSE,(J,PV(J),PU(J),PS(J),J=1,NPARM)
      WRITE(6,611)DSE
      WRITE(6,612)NDE
      IP=0
      DO I=1,NDE
          IP=IP+1
          WRITE(6,622)I-1,IP,PV(IP),PU(IP),PS(IP)
          ENDDO
      WRITE(6,613)NRE
      IP=0
      IP=NDE
      DO I=1,NRE
          IP=IP+1
          WRITE(6,623)I-1,IP,PV(IP),PU(IP),PS(IP)
          ENDDO
      WRITE(6,614) NPOW
      IP=0
      IP=NDE+NRE
      NPS=0
      DO I=1,NPOW
          DO J=1,NPHI(I)
              IP=IP+1
              WRITE(6,624)I-1,IP,PV(IP),PU(IP),PS(IP)
              ENDDO
          NPS=NPS+NPHI(I)
          ENDDO
      WRITE(6,615)1
c     for Vasy
      IP=NDE+NRE+NPS
      IP=IP+1
      WRITE(6,625)1,IP,PV(IP),PU(IP),PS(IP)
      DS=0.0D0
      DO I=1,NTP
          DO J=I+1,NTP
              DXX=RTP(I)*THETA(I)-RTP(J)*THETA(J) 
              DYY=RTP(I)*DSIN(DACOS(THETA(I)))-
     &                                 RTP(J)*DSIN(DACOS(THETA(J)))
              DD=SQRT(DXX**2+DYY**2)
cccc          IF(DD.LT.DR0)THEN
              DS=DS+YD(I)*YD(J)/(UVTP(I)*UVTP(J))*DEXP(-DD**2/DR0**2)
cccc          ENDIF
              ENDDO
          ENDDO 
      WRITE(6,666)DS
      WRITE(6,661) (DACOS(THETA(I))*180./PI,
     1 RTP(I),VTP(I),VTP(I)+YD(I),YD(I),YD(I)/UVTP(I),I=1,NTP)
      STOP
  602   FORMAT('Determine MLR using', i5, 'points, NCN=',i2,' CN=', 
     1 f11.4,'MCM=',i2,'RCMCN=',f11.4)
  611  FORMAT(' Direct non-linear fit DSE=',1Pd11.5)
  612  FORMAT(' Direct non-linear fit De with ',i2,' parameter')
  613  FORMAT(' Direct non-linear fit Re with ',i2,' parameter')
  614  FORMAT(' Direct non-linear fit Phi with ',i2,' parameter')
  615  FORMAT(' Direct non-linear fit Vasy with ',i2,' parameter')
  622  FORMAT('   De_'i2,' para_{',i2,'} =',1pd17.9,' (+/-',1pd8.1,')
     1    PS=',1pd8.1)
  623  FORMAT('   Re_'i2,' para_{',i2,'} =',1pd17.9,' (+/-',1pd8.1,')
     1    PS=',1pd8.1)
  624  FORMAT('  phi_'i2,' para_{',i2,'} =',1pd17.9,' (+/-',1pd8.1,')
     1    PS=',1pd8.1)
  625  FORMAT(' Vasy_'i2,' para_{',i2,'} =',1pd17.9,' (+/-',1pd8.1,')
     1    PS=',1pd8.1)
  610  FORMAT(/' Direct non-linear fit DSE=',1Pd11.5/
     1  ('   para_{',i2,'} =',d17.9,' (+/-',d8.1,')   PS=',d8.1))
  661 FORMAT(1x,2F8.4,4F15.5)
  666 FORMAT('appraise non-systematical change factor for fit=',1Pd11.5)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DYIDPJ(IDAT,NDATA,NPARM,IFXP,YC,PV,PD,PS,RMSR)
      INTEGER MXDATA, MXPARM, NPMAX
      PARAMETER (MXDATA=15000, MXPARM=60, NPMAX=20)
      INTEGER  I,J,K,IP,IDAT,NDATA,NPARM,IFXP(MXPARM),
     1  p,NDE,NRE,NCN,MCM,MMN,NS,NL,NPOW,NPHI(NPMAX),
     2  NPS
      REAL*8  YC,PV(MXPARM),PD(MXPARM),PS(MXPARM),RMSR,
     1  Re,De,Vasy,RREF,AREF,AREFp,Rep,CN,RCN2CN,RCMCN,RCM2CN,VLRe,
     2 dVLRedRe,phiINF,RTPp,yp,ype,dype,yPOW,XP,SUM,DSUM,VLR,XPW,
     3 DER,XDE,XRE
      REAL*8 Pn(0:NPMAX+1),PHI(NPMAX) 
      REAL*8 RTP(MXDATA),THETA(MXDATA)
      COMMON /DATABLK/RTP,THETA,RREF,CN,RCN2CN,RCMCN,RCM2CN,
     1 p,NPHI,NDE,NRE,NCN,MCM,NS,NL
c=======================================================================
c  For the case of an  MLR_{p}  potential ...
c-----------------------------------------------------------------------
      Pn(0)=1.
      Pn(1)=THETA(IDAT)
      DO I=1,NPMAX
          Pn(I+1)=((2*I+1)/(I+1))*THETA(IDAT)*Pn(I)-(I/(I+1))*Pn(I-1)        
          ENDDO 
c caculate the derivative of the parameters of De not including
c the coefficient before, only the Legendre expansion.
      De=0.0d0
      DO I=0,NDE-1
          De=De+Pn(I)*PV(I+1)
          PD(I+1)=Pn(I)
          ENDDO 
c caculate the derivative of the parameters of Re not including
c the coefficient before, only the Legendre expansion.
      Re=0.0d0
      IP=NDE
      DO I=0,NRE-1
          IP=IP+1
          Re=Re+Pn(I)*PV(IP)
          PD(IP)=Pn(I)
          ENDDO  
      AREF= RREF
      IF(RREF.LE.0.d0) AREF= Re
      AREFp= AREF**p
      Rep= Re**p
      VLRe= CN*(1.+Pn(2)*RCN2CN)/Re**NCN
      dVLRedRe= -NCN*VLRe/Re
      IF(MCM.GT.NCN) THEN
          MMN= MCM - NCN
          IF(p.LE.MMN) MMN= 0
          IF(MMN.GT.0) THEN
              VLRe= (CN/Re**NCN)*(1.+RCN2CN*Pn(2)+(RCMCN+RCM2CN*Pn(2))/
     1                                                        Re**MMN)
              dVLRedRe= dVLRedRe - MCM*CN*(RCMCN+RCM2CN*Pn(2))/
     1                                                     Re**(MCM+1)
              ENDIF
          phiINF= DLOG(2.d0*De/VLRe)
          ENDIF
      RTPp= RTP(IDAT)**p
      yp= (RTPp - AREFp)/(RTPp + AREFp)
      ype= (RTPp - Rep)/(RTPp + Rep)
c caculate the derivative of the parameters of PHI(N) not including
c the coefficient before, the Legendre expansion and exponent expansion.
      NPOW= NS+1
      IF(RTP(IDAT).GE.Re) NPOW= NL+1
      yPOW= 1.d0 - yp
      SUM=0.0 
      DSUM=0.0 
      NPS=0
      IP=NDE+NRE
      DO J=1,NPOW
          IP=IP+1
          PHI(J)= PV(IP)*Pn(0)
          PD(IP)=yp**(J-1)
          DO  K=2,NPHI(J)
              IP=IP+1
              PHI(J)=PHI(J)+ PV(IP)*Pn(K-1)
              PD(IP)=Pn(K-1)*yp**(J-1)
              ENDDO
          NPS=NPS+NPHI(J)
          SUM=SUM+PHI(J)*yp**(J-1)
          IF(RREF.LE.0.D0) DSUM=DSUM+yPOW*PHI(J)*(J-1)*yp**(J-1)
          ENDDO
calculate the derivative of the parameters of Vasy 
      IP=NDE+NRE+NPS
      PD(IP+1)=1.0D0
      Vasy=PD(IP+1)*PV(IP+1) 
      XP= SUM*yPOW+ phiINF*yp
      VLR= CN*(1.+RCN2CN*Pn(2))/RTP(IDAT)**NCN
      IF(MMN.GT.0) THEN
          VLR= (CN/RTP(IDAT)**NCN)*(1.+RCN2CN*Pn(2)+
     &                            (RCMCN+RCM2CN*Pn(2))/RTP(IDAT)**MMN)
          ENDIF
      XPW= DEXP(-XP*ype) * VLR/VLRe
      YC= De*(1.d0 - XPW)**2-De+Vasy 
c     write(11,99)IDAT, RTP(IDAT),THETA(IDAT),YC
c99   format(1x,I5,2F8.3,F15.7)
      DER= 2.d0*De*(1.d0- XPW)*XPW
      yPOW= DER*ype*(1.d0- yp)
c  calculate derivatives of all parameters include the coefficents before
c  for De
      XDE=((1.d0- XPW)**2 + DER/De*ype*yp-1.D0)
      DO I=1,NDE
          PD(I)=XDE*PD(I)
          ENDDO 


c for Re
      dype= -0.5d0*(p/Re)*(1.d0 - yp**2)
      IF(RREF.LE.0.d0) THEN
          DSUM= phiINF - SUM + DSUM 
        ELSE
          DSUM= 0.d0
        ENDIF
      XRE= DER*(dype*(XP + ype*DSUM) + (1.d0 - ype*yp)*dVLRedRe/VLRe)
      IP=NDE
      DO I=1,NRE
          IP=IP+1
          PD(IP)=XRE*PD(IP)  
          ENDDO 
c for phi(n)
      IP=NDE+NRE
      DO J=1,NPOW
          DO K=1,NPHI(J)
              IP=IP+1
              PD(IP)=PD(IP)*yPOW
              ENDDO
          ENDDO
c  for Vasy 
c  the same as before      
      RETURN
      END
c-----------------------------------------------------------------------
