c***********************************************************************
      SUBROUTINE DYIDPJ(IDAT,NDATA,NPARM,IFXP,YC,PV,PD,PS,RMSR)
c** Subroutine to calculate potential function value YC at distance
c  RDIST= RTP(IDAT), and its partial derivatives w.r.t. the various 
c  potential parameters.  If  IDAT.LE.1  generate a new set of 
c  internal potential variables, while if  IDAT > 1  use SAVED values
c... [Must ensure that calculations based on the current UPDATED PV(j)]
c------------------------------------------------------------------------
      INTEGER MXDATA, MXPARM, MXMLR
      PARAMETER (MXDATA=1501, MXPARM=43, MXMLR= 15)
      INTEGER  i,j,m,IDAT, NPOW,NPARM,NDATA, IFXP(MXPARM),JFXRe,JFXDe,
     1  JFXVMIN
      REAL*8  YC,PV(NPARM),PD(NPARM),PS(NPARM),DM(MXMLR),DMP(MXMLR),
     1 DMPP(MXMLR),RMSR,RTPp,RTPq,Rep,AREF,AREFp,AREFq,ype,dype,
     2 betaINF,yp,yq,yPOW,XP,XPW,DER,TCM,UM,TTMM,DERP,SUM,DSUM,AA,BB,
     3 FCT,FCT2,ULR,ULRe,dULRe,d2ULRe,DDER,T0,T0P,T1,RE3,RE6,RE8,RTP3,
     4 RTP6,RTP8,RDIST,C6adj,C9adj,dAAdRe,dBBdRe ,BETAN,beta0,C1tst
      REAL*8  DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),DEIGRe(1,1)
     1  ,DEIGDe(1,1)
c-----------------------------------------------------------------------
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,IDF,IDSTT,KDER,NCMM,MMLR(MXMLR),
     1                                              p,q,NS,NL,NPHI,SAP
      REAL*8 Re,De,VMIN,RREF,M2,as,bs,rhoAB,
     1  CmVAL(MXMLR),RTP(MXDATA),SAS(MXDATA,MXPARM)
      COMMON /DATABLK/Re,De,VMIN,RREF,M2,as,bs,rhoAB,CmVAL,RTP,SAS,PSEL,
     1 IFXRe,IFXDe,IFXVMIN,IDF,IDSTT,KDER,NCMM,MMLR,p,q,NS,NL,NPHI,SAP
c-----------------------------------------------------------------------
      SAVE JFXRe,JFXDe,JFXVMIN, AREF,AREFp,AREFq,Rep,C6adj,C9adj,
     1  betaINF,BETA0, AA,BB,ULRe,dULRe
c=======================================================================
      IF(ABS(IDAT).LE.1) THEN
          JFXRe= IFXP(NPHI+1)
          JFXDe= IFXP(NPHI+2)
          JFXVMIN= IFXP(NPHI+3)
          ENDIF
      RDIST= RTP(IDAT)
      DO  j=1,NPARM
          PD(j)= 0.d0
          ENDDO
c=======================================================================
c** For case of an  EMO_p  potential
c-----------------------------------------------------------------------
      IF(PSEL.EQ.1) THEN
          IF(ABS(IDAT).LE.1) THEN
              IF(JFXRe.LE.0) Re= PV(NPHI+1)
              IF(JFXDe.LE.0) De= PV(NPHI+2)
              IF(JFXVMIN.LE.0) VMIN= PV(NPHI+3)
              AREF= RREF
              IF(RREF.LE.0) AREF= Re
              AREFp= AREF**p
              ENDIF
          RTPp= RDIST**p
          yp= (RTPp - AREFp)/(RTPp + AREFp)
          NPOW= NS+1
          IF(RDIST.GE.Re) NPOW= NL+1
          yPOW= 1.d0
          SUM= PV(1)
          DSUM= 0.d0
          IF(NPOW.GE.2) THEN
              DO  j= 2,NPOW
                  IF(RREF.LE.0.d0) DSUM= DSUM+ (j-1)*PV(j)*yPOW
                  yPOW= yPOW*yp
                  SUM= SUM+ yPOW*PV(j)
                  ENDDO
              ENDIF
          XP= DEXP(-SUM*(RDIST- Re))
          YC= De*(1.d0 - XP)**2 + VMIN
          DER= 2.d0*De*(1.d0- XP)*XP
          DERP= DER*(RDIST- Re)
          DO  j= 1,NPOW
              PD(j)= DERP
              DERP= DERP*yp
              ENDDO
c** If appropriate, also get partial derivative w.r.t. Re
          IF(JFXRE.LE.0) THEN
              IF(RREF.LE.0.d0) SUM= SUM +(RDIST- Re)*DSUM*0.5d0*
     1                                             (p/Re)*(1.d0-yp**2)
              PD(NPHI+1)= -DER*SUM
              ENDIF
c** If appropriate, also get partial derivative w.r.t. De & VMIN
          IF(JFXDe.LE.0) PD(NPHI+2)= (1.d0- XP)**2
          IF(JFXVMIN.LE.0) PD(NPHI+3)= 1.d0
          ENDIF

c=======================================================================
c  For the case of an  MLR_{p}  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.2) THEN
          IF(ABS(IDAT).LE.1) THEN
              IF(JFXRe.LE.0) Re= PV(NPHI+1)
              IF(JFXDe.LE.0) De= PV(NPHI+2)
              IF(JFXVMIN.LE.0) VMIN= PV(NPHI+3)
              AREF= RREF
              IF(RREF.LE.0.d0) AREF= Re
              AREFp= AREF**p
              AREFq= AREF**q
              Rep= Re**p
              IF((NCMM.GE.4).AND.(MMLR(2).LE.0)) THEN
c** For Aubert-Frecon based  coupled-states ULR(r) for Li(2S)+Li(2P)
                  RE3= 1.d0/Re**3
                  RE6= RE3*RE3
                  RE8= RE6/Re**2
                  C6adj= CmVAL(3) + CmVAL(1)**2/(4.d0*DE)
                  C9adj= 0.5d0*CmVAL(1)*C6adj/De
                  IF(MMLR(2).EQ.0) THEN
c ... for Aubert-Frecon Li2(A) 2x2 {3,0,6,6,8,8} case ...
                      T1= (0.5d0*CmVAL(1)+ (C6adj- CmVAL(4))*RE3)
     1                                                       *RE3/3.d0
                      IF(NCMM.GT.4) THEN
                          T1= T1+ (CmVAL(5)- CmVAL(6))*RE8/3.d0
                          ENDIF
                      T0= DSQRT((T1- CmVAL(2))**2 + 8.d0*T1**2)
                      ULRe= 0.5d0*( - CmVAL(2) + (1.5d0*CmVAL(1)
     1         + (C6adj + CmVAL(4))*RE3)*RE3) + 0.5d0*T0 + C9adj*RE6*RE3
                      IF(NCMM.GT.4) 
     1                        ULRe= ULRe+0.5d0*(CmVAL(5)+CmVAL(6))*RE8
                      T0P= (9.d0*T1-CmVAL(2))/T0
                      dULRe= -RE3*(0.25d0*CmVAL(1)*(9.d0 + T0P)
     1               + RE3*(C6adj*(3.d0 + T0P) + CmVAL(4)*(3.d0 - T0P)
     2                                           + RE3*9.d0*C9adj))/Re
                      IF(NCMM.GT.4) THEN
                          dULRe= dULRe -RE8*4.d0*(CmVAL(5) 
     1                *(3.d0 + T0P) + CmVAL(6)*(3.d0 - T0P))/(3.d0*Re)
                          ENDIF
                      ENDIF 
                  IF(MMLR(2).EQ.-1) THEN
c ... for Aubert-Frecon Li2(c) 3x3 {3,0,6,6,8,8} case ...
                      CALL AF3x3potret(Re,CmVAL(2),CmVAL(1),C6adj,
     1             CmVAL(5),De,ULRe,DEIGM1,DEIGM3,DEIGM5,DEIGR,DEIGDe)

ccc                  c9adj= 0.0

                      ULRe= ULRe + C9adj*RE6*RE3
                      dULRe= DEIGR(1,1) - 9.d0*C9adj*RE8/Re**2
                      ENDIF
                ELSE
c** For normal inverse-power sum MLR/MLJ case, with or without damping
                  KDER=1
                  IF(rhoAB.GT.0.d0) CALL dampF(Re,rhoAB,NCMM,MMLR,IDF,
     1                                         IDSTT,KDER,DM,DMP,DMPP)
                  ULRe= 0.d0
                  dULRe= 0.d0
                  DO  m= 1,NCMM
                      T0= CmVAL(m)/Re**MMLR(m)
                      IF(rhoAB.GT.0.d0) T0= T0*DM(m)
                      ULRe= ULRe + T0
                      dULRe= dULRe - T0*MMLR(m)/Re
                      IF(rhoAB.GT.0.d0) dULRe= dULRe + T0*DMP(m)/DM(m)
                      ENDDO
                ENDIF
              betaINF= DLOG(2.d0*De/ULRe)
              IF(IDSTT.GT.1) THEN
c** For MLR with DS(s=-1/2) damping constrained to have C1= Z1*Z2= IDSTT
c   then fix the  beta(N) value ...
                  BETA0= 0.d0
                  DO  m= 1,NCMM
                      BETA0= BETA0 + CmVAL(m)*DSQRT(3.69d0*rhoAB/
     1                                         MMLR(m))**(2*MMLR(m)-1)
                      ENDDO
                  C1tst= BETA0/uLre
                  BETAN= DLOG(DSQRT(4.d0*De*IDSTT*116140.97d0)/BETA0)
                  BETA0= DLOG(DSQRT(IDSTT*116140.97d0/De)*uLRe/BETA0) 
                  C1tst= De*(DEXP(BETA0)*C1tst)**2
c*** For constrained C1/r case, now define  beta(N)
                  BETAN= -0.5D0*BETAN*(-1)**NPOW
                  T0= 1.d0
                  NPOW= NPHI- 1
                  DO  j= NPOW, 1, -1
                      BETAN= BETAN + T0*PV(j)
                      T0= -T0
                      ENDDO
                  PV(NPHI)= BETAN
                  IFXP(NPHI)= 1
                  WRITE(6,666) C1tst/116140.97d0, NPHI,PV(NPHI)
  666 FORMAT(/'  Test  C1(sr)=',1Pd10.3,'    PV(',i3,')=',1PD16.8/)
                  ENDIF
              ENDIF
          RTPp= RDIST**p
          RTPq= RDIST**q
          yp= (RTPp - AREFp)/(RTPp + AREFp)
          yq= (RTPq - AREFq)/(RTPq + AREFq)
          ype= (RTPp - Rep)/(RTPp + Rep)
          IF(SAP.GT.0) THEN
c*** Case of Pashov natural spline exponent ....
c... Now, use a spline through the exponent values defined by the input
c    points to generate values of that exponent at the desired
c    spline-definition points
              NPOW= NPHI
              XP= 0.d0
              DO  J= 1,NPOW
                  PD(J)= SAS(IDAT,J)
                  XP= XP + PV(J)*PD(J)
                  ENDDO
            ELSE
c... For conventional case of a constrained polynomial exponent function
              NPOW= NS+1
              IF(RDIST.GE.Re) NPOW= NL+1
              yPOW= 1.d0 - yp
              SUM= PV(1)*yPOW
              DSUM= 0.d0
              IF(NPOW.GE.2) THEN
                  DO  j= 2,NPOW
                      IF(RREF.LE.0.d0) DSUM= DSUM + PV(j)*(j-1)*yPOW
                      yPOW= yPOW*yq
                      SUM= SUM+ yPOW*PV(j)
                      ENDDO
                  ENDIF
              XP= SUM + betaINF*yp
            ENDIF
          IF((NCMM.GE.4).AND.(MMLR(2).LE.0)) THEN
c** For Aubert-Frecon based  uLR(r)
              RTP3= 1.d0/RDIST**3
              RTP6= RTP3*RTP3
              RTP8= RTP6/RDIST**2
              IF(MMLR(2).EQ.0) THEN
c ... extension for Aubert-Frecon Li2(A) {3,0,6,6,8,8} case ...
                  T1= (0.5d0*CmVAL(1) + (C6adj - CmVAL(4))*RTP3)*RTP3/3.d0
                  IF(NCMM.GT.4) THEN
                      T1= T1+ (CmVAL(5)- CmVAL(6))*RTP8/3.d0
                      ENDIF
                  T0= DSQRT((T1- CmVAL(2))**2 + 8.d0*T1**2)
                  ULR= 0.5d0*( - CmVAL(2) + (1.5d0*CmVAL(1) + (C6adj 
     1            + CmVAL(4))*RTP3)*RTP3) + 0.5d0*T0 + C9adj*RTP3*RTP6
                  IF(NCMM.GT.4) ULR= ULR+0.5d0*(CmVAL(5)+ CmVAL(6))*RTP8
c... SKIP Re derivative corrections for all?
                  T0P= (9.d0*T1-CmVAL(2))/T0
c             dULRdRe= -RTP3*(0.25d0*CmVAL(1)*(9.d0 + T0P)
c    1      + RTP3*(CmVAL(3)*(3.d0 + T0P) + CmVAL(4)*(3.d0 - T0P)))/Re
                  ENDIF
              IF(MMLR(2).EQ.-1) THEN
c ... for Aubert-Frecon Li2(c) 3x3 {3,0,6,6,8,8} case ...
                  CALL AF3x3potret(RDIST,CmVAL(2),CmVAL(1),C6adj,
     1             CmVAL(5),De,ULR,DEIGM1,DEIGM3,DEIGM5,DEIGR,DEIGDe)
                  ULR= ULR + C9adj*RTP6*RTP3
                  ENDIF
              ENDIF
          IF((NCMM.LE.1).OR.(MMLR(2).GT.0)) THEN
c** For normal inverse-power sum MLR/MLJ case, with or without damping
              KDER= 0
              IF(rhoAB.GT.0.d0) CALL dampF(RDIST,rhoAB,NCMM,MMLR,IDF,
     1                                         IDSTT,KDER,DM,DMP,DMPP)
              ULR= 0.d0
              DO  m= 1,NCMM
                  T0= CmVAL(m)/RDIST**MMLR(m)
                  IF(rhoAB.GT.0.d0) T0= T0*DM(m)
                  ULR= ULR + T0
                  ENDDO
              ENDIF
          XPW= DEXP(-XP*ype) * ULR/ULRe
          YC= De*(1.d0 - XPW)**2 + VMIN
          DER= 2.d0*De*(1.d0- XPW)*XPW
          yPOW= DER*ype*(1.d0- yp)
          IF(SAP.GT.0) THEN
c... finalize derivative w.r.t. exponent beta-function spline points ...
              DO  J= 1,NPOW
                  PD(J)= PD(J)*DER*ype
                  ENDDO 
            ELSE
c... finalize derivative w.r.t. exponent polynomial coefficient ....
              DO  j= 1,NPOW
                  PD(j)= yPOW
                  yPOW= yPOW*yq
                  ENDDO
            ENDIF
c** If appropriate, also get partial derivative w.r.t. De & VMIN
          IF(JFXDe.LE.0) PD(NPHI+2)= (1.d0- XPW)**2 + DER*ype*yp/De
          IF(JFXVMIN.LE.0) PD(NPHI+3)= 1.d0
c** If appropriate, also get partial derivative w.r.t. Re
          IF(JFXRe.LE.0) THEN
              dype= -0.5d0*(p/RE)*(1.d0 - yp**2)
              IF(RREF.LE.0.d0) THEN
                  DSUM= betaINF - SUM/(1.d0-yp) + DSUM 
                ELSE
                  DSUM= 0.d0
                ENDIF
              PD(NPHI+1)= DER*(dype*(XP + ype*DSUM) 
     1                               + (1.d0 - ype*yp)*dULRe/ULRe )
              ENDIF
          ENDIF

c=======================================================================
c** For the case of a  DELR_p  potential
c-----------------------------------------------------------------------
      IF(PSEL.EQ.3) THEN
          IF(ABS(IDAT).LE.1) THEN
              IF(JFXRe.LE.0) Re= PV(NPHI+1)
              IF(JFXDe.LE.0) De= PV(NPHI+2)
              IF(JFXVMIN.LE.0) VMIN= PV(NPHI+3)
              AREF= RREF
              IF(RREF.LE.0) AREF= Re
              AREFp= AREF**p
              ULRe= 0.d0
              dULRe= 0.d0
              d2ULRe= 0.d0
c-----------------------------------------------------------------------
c** Evaluate uLR & its first 2 deriv. at  Re ... 
              KDER=2
              IF(rhoAB.GT.0.d0) CALL dampF(Re,rhoAB,NCMM,MMLR,IDF,
     1                                         IDSTT,KDER,DM,DMP,DMPP)
              ULRe= 0.d0
              dULRe= 0.d0
              d2ULRe= 0.d0
              DO  m= 1,NCMM
                  T0= CmVAL(m)/Re**MMLR(m)
                  ULRe= ULRe + T0*DM(m)
                  dULRe= dULRe + T0*(DMP(m) - DM(m)*MMLR(m)/RE)
                  d2ULRe= d2ULRe + T0*(DMPP(m) - 2.d0*MMLR(m)*DMP(m)/RE
     1                           + MMLR(m)*(MMLR(m)+1.d0)*DM(m)/RE**2)
                  ENDDO
c-----------------------------------------------------------------------
              AA= De - ULRe - dULRe/PV(1)
              BB= 2.d0*(De - ULRe) - dULRe/PV(1)
              dAAdRe= -dULRe - d2ULRe/PV(1)
              dBBdRe= -2.d0*dULRe - d2ULRe/PV(1)
              ENDIF
          RTPp = RDIST**p
          yp= (RTPp - AREFp)/(RTPp + AREFp)
          ULR= 0.d0
c ... evaluate uLR at the actual distance ...
          KDER=0
          CALL dampF(RDIST,rhoAB,NCMM,MMLR,IDF,IDSTT,KDER,DM,DMP,DMPP)
          ULR= 0.d0
          DO  m= 1,NCMM
              ULR= ULR + DM(m)*CmVAL(m)/RDIST**MMLR(m)
              ENDDO
          NPOW= NS+1
          IF(RDIST.GE.Re) NPOW= NL+1
          yPOW= 1.d0
          SUM= PV(1)
          DSUM= 0.d0
          IF(NPOW.GE.2) THEN
              DO  j= 2,NPOW
                  IF(RREF.LE.0.d0) DSUM= DSUM+ (j-1)*PV(j)*yPOW
                  yPOW= yPOW*yp
                  SUM= SUM+ yPOW*PV(j)
                  ENDDO
              ENDIF
          XP= DEXP(-SUM*(RDIST- Re))
          YC= (AA*XP - BB)*XP + De - ULR + VMIN
          DER= XP*(BB - 2.d0*AA*XP)
          DERP= DER*(RDIST- Re)
          DO  j= 1,NPOW
              PD(j)= DERP
              DERP= DERP*yp
              ENDDO
          PD(1)= PD(1) + XP*(XP-1.d0)*dULRe/PV(1)**2
c** If appropriate, also get partial derivative w.r.t. De & VMIN
          IF(JFXDe.LE.0) PD(NPHI+2)= (XP - 2.d0)*XP + 1.d0
          IF(JFXVMIN.LE.0) PD(NPHI+3)= 1.d0
          IF(JFXRe.LE.0) THEN
c** If appropriate, also get partial derivative w.r.t. Re
              IF(RREF.GT.0.d0) THEN
                  DSUM= 0.d0
                ELSE
cc                DSUM= - DSUM*0.5d0*(p/Re)*(1.d0-yp**2)*(RDIST - Re)
                  DSUM= - DSUM*0.5d0*(p/Re)*(1.d0-yp**2)
                ENDIF
              PD(NPHI+1)= DER*(DSUM - SUM) + XP*(dAAdRe*XP- dBBdRe)
              ENDIF
          ENDIF
c=======================================================================

c=======================================================================
c** For the case of a  GPEF(p,as,bs)  potential
c-----------------------------------------------------------------------
      IF(PSEL.EQ.4) THEN
          IF(ABS(IDAT).LE.1) THEN
              JFXVMIN= IFXP(NPHI+1)
              IF(JFXVMIN.LE.0) VMIN= PV(NPHI+1)
              JFXRe= IFXP(NPHI+2)
              IF(JFXRe.LE.0) Re= PV(NPHI+2)
              Rep= Re**p
              ENDIF
          RTPp= RDIST**p
          yp= (RTPp - Rep)/(as*RTPp + bs*Rep)
          yPOW= yp**2
          SUM= PV(1)*yPOW
          DSUM= 2.d0*PV(1)*yp
          IF(IFXDe.LE.0) yPOW= PV(1)*yPOW
          PD(1)= yp**2
          IF(NPHI.GT.1) THEN
              DO  j= 2,NPHI
                  DSUM= DSUM + (j+1)*PV(j)*yPOW
                  yPOW= yPOW* yp
                  PD(j)= yPOW
                  SUM= SUM+ PV(j)*yPOW
                  ENDDO
              ENDIF
          YC= SUM + VMIN
          IF(JFXVMIN.LE.0) PD(NPHI+1)= 1.d0
          IF(JFXRe.LE.0) THEN
              PD(NPHI+2)= -DSUM* (p/Re)*Rep*RTPp*(as+bs)/
     1                                            (as*RTPp +bs*Rep)**2
              ENDIF 
          ENDIF
c=======================================================================
c%%%%%%%
cc    if(IDAT.eq.1) then
cc          write(7,700) nparm,(PV(i),i=1,nparm)
cc          write(7,702)  (i,i=1,min0(5,nparm))
cc          ENDIF
cc    write(7,704) i,yc,(pd(i),i=1,nparm)
cc700 FORMAT(///' Partial derivatives for',i3,' input parameters:'/
cc   1  (1P5D16.8))
cc702 FORMAT('  I    YC   ',5('     PD(',i1,')   ':))
cc704 format(i3,f9.2,1P5D13.5:/(12x,5d13.5:))
c%%%%%%%
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      double precision function Scalc(x,m,n,y,rKL,LMAX)
c** At the position 'x', evaluate the m'th Sm(x) function contributing 
c  the definition of the the natural cubic spline defined by 
c  function values at the  n  points  y(i) [i=1,n]
      INTEGER  LMAX,I,K,KK,M,N
      REAL*8  x,y1,y2,y(1:LMAX),rKL(1:LMAX,1:LMAX)
      k= 0
      kk= 0
      do i=2,n
c... select interval
          if ((x.gt.y(i-1)).and.(x.le.y(i)))  k=i
          end do
      if (x.lt.y(1)) then
          k=2
          kk=1
          end if
      if (x.gt.y(n)) then
          k=n
          kk=1
          end if
      if(x.eq.y(1)) k=2
      y1=y(k-1)
      y2=y(k)
      Scalc= 0.d0
      IF(kk.eq.0) 
     1    Scalc= rKL(m,k)*((y1-x)*(((y1-x)/(y1-y2))**2-1)/6)*(y1-y2)
     2         + rKL(m,k-1)*((x-y2)*(((x-y2)/(y1-y2))**2-1)/6)*(y1-y2)
      IF(k.EQ.m) Scalc= Scalc + (y1-x)/(y1-y2)
      IF(k-1.EQ.m) Scalc= Scalc + (x-y2)/(y1-y2)
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      subroutine Lkoef(n,x,A,LMAX)   
c*** Based on nespl subroutine          
      INTEGER LMAX
      INTEGER I,J,N,INDX(1:LMAX)
      REAL*8 X(1:LMAX),A(1:LMAX,1:LMAX),B(1:LMAX,1:LMAX), d
c
      A(1,1)=(x(3)-x(1))/3
      A(1,2)=(x(3)-x(2))/6
      do i=2,n-3
          A(i,i-1)=(x(i+1)-x(i))/6
          A(i,i)=(x(i+2)-x(i))/3
          A(i,i+1)=(x(i+2)-x(i+1))/6
          end do
      A(n-2,n-3)=(x(n-1)-x(n-2))/6
      A(n-2,n-2)=(x(n)-x(n-2))/3  
      do i=1,n-2
          B(i,i)=1/(x(i+1)-x(i))
          B(i,i+1)=-1/(x(i+2)-x(i+1))-1/(x(i+1)-x(i))
          B(i,i+2)=1/(x(i+2)-x(i+1))
          end do  
      call ludcmp(A,n-2,LMAX,indx,d)
      do i=1,n 
          call lubksb(A,n-2,LMAX,indx,B(1,i))
          end do 
      do i=1,n-2
          do j=1,n
              A(j,i+1)=B(i,j)
              end do
          end do 
      do i=1,n
          A(i,1)=0.0
          A(i,n)=0.0
          end do
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      double precision d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      double precision aamax,dum,sum,vv(NMAX)
      d=1.d0
      do  i=1,n
          aamax=0.
          do  j=1,n
              if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
              enddo
          if (aamax.eq.0.) pause 'singular matrix in ludcmp'
          vv(i)=1.d0/aamax
          enddo
      do  j=1,n
          do  i=1,j-1
              sum=a(i,j)
              do  k=1,i-1
                  sum=sum-a(i,k)*a(k,j)
                  enddo
              a(i,j)=sum
              enddo
          aamax=0.
          do  i=j,n
              sum=a(i,j)
              do  k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
                  enddo
              a(i,j)=sum
              dum=vv(i)*abs(sum)
              if (dum.ge.aamax) then
                  imax=i
                  aamax=dum
                  endif
              enddo
          if(j.ne.imax)then
              do  k=1,n
                  dum=a(imax,k)
                  a(imax,k)=a(j,k)
                  a(j,k)=dum
                  enddo
              d=-d
              vv(imax)=vv(j)
              endif
          indx(j)=imax
          if(a(j,j).eq.0.)a(j,j)=TINY
              if(j.ne.n)then
                  dum=1./a(j,j)
                  do  i=j+1,n
                      a(i,j)=a(i,j)*dum
                      enddo
                  endif
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER i,ii,j,ll, n,np,indx(n)
      double precision a(np,np),b(n), sum
      ii=0
      do  i=1,n
          ll=indx(i)
          sum=b(ll)
          b(ll)=b(i)
          if (ii.ne.0)then
              do  j=ii,i-1
                  sum=sum-a(i,j)*b(j)
                  enddo
            else if (sum.ne.0.) then
              ii=i
            endif
          b(i)=sum
          enddo
      do  i=n,1,-1
          sum=b(i)
          do  j=i+1,n
              sum=sum-a(i,j)*b(j)
              enddo
          b(i)=sum/a(i,i)
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
