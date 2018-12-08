c***********************************************************************
c***********  Program  betaFIT_2.0  dated  7 June 2009  ****************
c***********************************************************************
c* Program to fit NTP read-in potential fx. values {RTP(i),VTP(i)} 
c  (with or without individual weights) to a chosen analytic form.
c***********************************************************************
      INTEGER MXDATA, MXPARM, MXMLR
      PARAMETER (MXDATA=1501, MXPARM=43, MXMLR= 15)
      INTEGER I,J,INFL,ITER,IROUND,ROBUST,prFIT,prNLL,prDIFF,M,NPARM,
     1  NTP,NLIN,NDGF,IFXP(MXPARM)
      REAL*8 BETA(0:MXPARM),PV(MXPARM),PU(MXPARM),PS(MXPARM),
     1 CM(MXPARM,MXPARM),DYDP(MXDATA,MXPARM),VTP(MXDATA),
     2 uVTP(MXDATA),betay(MXDATA),Ubetay(MXDATA),YD(MXDATA),
     3 ypSAP(MXPARM),xSAP(MXDATA), rKL(1:MXDATA,1:MXDATA),
     4 DM(MXMLR),DMP(MXMLR),DMPP(MXMLR),
     3 betaINF,UNC,yPOW,DSE,TSTPS,TSTPU,DSEB, bDAMPR,TCM,UM,
     4 Rep,AREF,AREFp,AREFq,RTPp,RTPq, AA,BB,ULR,dULR,FCT,RAT,UMAX,
     5 XX,YY,YH,yp,yq,ypRE,ReDE, ReIN,DeIN,VMINin,ULRe,dULRe,RE3,RE6,
     6 RE8,T0,T1,C6adj,C9adj,RTP3,RTP6,RTP8,RH,RR,RB,RBB,VV,VB,VBB,
     7 SCALC,Rsap,DRMSD
      CHARACTER*4  NNAME,NAME(4)
      DATA NAME/' EMO',' MLR','DELR','GPEF'/
c-----------------------------------------------------------------------
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR(MXMLR),p,q,NS,NL,
     1                                                        NPHI,SAP
      REAL*8 Re,De,VMIN,RREF,M2,as,bs,bDAMP,
     1  CmVAL(MXMLR),RTP(MXDATA),SAS(MXDATA,MXPARM)
      COMMON /DATABLK/Re,De,VMIN,RREF,M2,as,bs,bDAMP,CmVAL,RTP,SAS,PSEL,
     1          IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR,p,q,NS,NL,NPHI,SAP
c-----------------------------------------------------------------------
      ROBUST= 0
      IDF= 1
c-----------------------------------------------------------------------
c** PSEL  specifies the type of potential being fitted to:
c     PSEL= 1 for EMO;   PSEL=2 for an MLR (or MLJ);   PSEL=3  for DELR ;
c     PSEL= 4  for GPEF
c* NPT  is the number of read-in potential points being fitted to.
c* UNC  is the energy uncertainty associated with the potential points
c     (plausibly ca. 0.1 cm-1 for RKR).  To weight each point separately
c     set  UNC < 0.0  and read in a separate uncertainty for each point
c  IROUND  specifies the level of rounding inside NLLSSRR if:
c          > 0 : requires that Sequential Rounding & Refitting be
c                performed, with each parameter being rounded at the
c                IROUND'th sig. digit of its local uncertainty.
c          <=0 : simply stops after full convergence (without rounding).
c* prFIT = 1  causes printout of results of preliminary linearized fit 
c       and (as appropriate) other non-final fits.  Normally set prFIT=0
c  prFIT > 1  specifies the level of printing inside fit subroutine NLLSSRR
c            >= 1  print converged parameters, PU & PS's
c            >= 2  also print parameter change each rounding step
c            >= 3  also indicate nature of convergence
c            >= 4  also print convergence tests on each cycle
c            >= 5  also parameters changes & uncertainties, each cycle
c  prDIFF (integer) > 0  causes fit residuals to be printed after each
c                        complete fit:  no print if  prDIFF.LE.0
c** Re & De  are assumed potential well minimum position & well depth.
c  ...  note that  De  is a dummy variable for a GPEF function
c** VMIN is the minimum of the potential function defined by the read-in 
c   points [typically =0 for RKR potential, but non-zero for ab initio]
c** To fix  Re, De or VMIN  unchanged at read-in values, set (integer)
c   IFXRe, IFXDe, and/or IFXVMIN > 0 ;  to fit them, set value =0
c   [Normally set  IFXVMIN= 0 !!]
c=======================================================================
      READ(5,*) PSEL, NTP, UNC, IROUND, prFIT, prDIFF
      READ(5,*) Re, De, VMIN
      READ(5,*) IFXRe, IFXDe, IFXVMIN
c=======================================================================
      prNLL= 0
      IF(prFIT.ge.2) prNLL= 5
      NNAME= NAME(PSEL)
      ReIN= Re
      DeIN= De
      VMINin= VMIN
      IF(PSEL.LE.3) WRITE(6,600) NNAME, VMIN, Re, De
      IF(PSEL.GT.3) WRITE(6,600) NNAME, VMIN, Re
c** For an MLR_p potential (PSEL=2) read number of long-range terms NCMM
c   to define long-range potential tail:  V(r)= De - \sum{CmVAL/r^MMLR}
c   with     beta(r)= y_{p} beta_inf + [1 - y_{p}] Sum{ beta_i y{_q}q^i }
c** bDAMP (real number) controls whether, and if 'yes' the strength of,
c      a damping function is included with each inverse-power long-range
c      term.  For  bDAMP.le.0  damping is ignored.  For  bDAMP > 0  it
c      is the coefficient in the damping function exponent. 
c** SAP is an integer: > 0  to use Pashov natural spline in MLR exponent
c                    : .le. 0  for normal constrained polynomial exponent
c*  For an MLR potential with  SAP > 0,  Rsap  is a negative{!} number
c      (-1 .le. Rsap < 0) which specifies the lower bound on the yp values 
c      selected to define the MLR exponent spline.  
c** For each long-range term read power  MMLR(i)  & coefficient CmVAL(i)
c** For special Aubert-Frecon case,  NCMM= 4,  MMLR= {3,0,6,6} or 
c  {3,0,6,6,8,8} and coefficients are: CmVAL(1)= C3(sig), CmVAL(2)= ASO,
c  CmVAL(3)= C6(Sig), CmVAL(4)= C6(pi)  for NCMM=4, and for NCMM= 6  add
c  coefficients CmVAL(5)= C8(sig) & CmVAL(6)= C8(pi)
c=======================================================================
      IF(PSEL.EQ.2) THEN
          READ(5,*) NCMM, bDAMP, SAP, Rsap
          READ(5,*) (MMLR(m), CmVAL(m), m= 1,NCMM)
          IF((NCMM.GE.4).AND.(MMLR(2).EQ.0)) THEN
c** For Lyon treatment of A-state alkali dimers ...
              WRITE(6,602) CmVAL(2),CmVAL(1),(CmVAL(i),i=3,NCMM)
            ELSE
              IF(bDAMP.LE.0) WRITE(6,698) (MMLR(m),CmVAL(m),m= 1,NCMM)
              IF(bDAMP.GT.0) WRITE(6,696) bDAMP,(MMLR(m),CmVAL(m),
     1                                                      m= 1,NCMM)
            ENDIF
          ENDIF
c=======================================================================
      IF(PSEL.EQ.3) THEN
c* For a DELR potential, NCMM is the number of long-range terms in the
c  sum  u_{LR}(r)= \sum_i{D_i(r) CmVAL(i)/r^i} , where i= MMLR(i) and
c  bDAMP  is the system-dependent scaling factor in the damping function.
c  Attractive terms have positive  CmVAL  values.
c* If  IDF= 1  use the Tang-Toennies damping function [bDAMP= 2.78*RHOd]
c   D_m(r)= 1 - exp{-bDAMP*r} \Sum_{k=0}^{m-1} (bDAMP*r)**{k}/k!
c Also ...  BETA(0) is the initial trial value of  \beta_0  used to
c    generate initial trial values of the A & B potential parameters.
c  * If the read-in  BETA(0) \leq 0.0 , the program uses a preliminary 
c    EMO_{p} fit to generate an estimate of  BETA(0); this usually works
c    if  u_{LR} is attractive, so that the potential has no barrier.
c  * If the potential has a barrier, one must determine an initial trial
c    value another way.  A one way would be to do an EMO_{p} fit to the
c    potential points, treating the barrier maximum as dissociation
c=======================================================================
          READ(5,*) NCMM, bDAMP
          READ(5,*) (MMLR(m), CmVAL(m), m= 1,NCMM)
c=======================================================================
          BETA(0)= 1.d0
          WRITE(6,696) bDAMP,(MMLR(m),CmVAL(m),m= 1,NCMM)
          ENDIF
c** For a GPEF potential, read coefficients to define expansion vble:
c       y = (r^p - Re^p)/(as*r^p + bs*Re^p)  where  p, as & bs all fixed
c=======================================================================
      IF(PSEL.EQ.4) THEN
          READ(5,*) as, bs 
          WRITE(6,605) as,bs
          ENDIF
c=======================================================================
c** Read the turning points to be fitted to
c=======================================================================
      IF(UNC.GT.0.d0) READ(5,*) (RTP(i), VTP(i),i= 1,NTP)
      IF(UNC.LE.0.d0) READ(5,*) (RTP(i), VTP(i),uVTP(i),i= 1,NTP)
c=======================================================================
      IF(UNC.GT.0.d0) THEN
          WRITE(6,606) NTP,UNC,(RTP(i),VTP(i),i= 1,NTP)
          DO  i= 1,NTP
              uVTP(i)= UNC
              ENDDO
        ELSE
          WRITE(6,609) NTP,VMIN,(RTP(i),VTP(i),uVTP(i),i= 1,NTP)
        ENDIF
      WRITE(6,608)
      IF((PSEL.EQ.2).AND.(NCMM.GE.4).AND.(MMLR(2).EQ.0)) THEN
          DO  I= 1,NTP
              VTP(I)= VTP(I) - 0.5d0*CmVAL(2)
              ENDDO
          ENDIF
c
  696 FORMAT(' uLR(r)  is a damped inverse-power sum with',6x,'bDAMP=',
     1  F7.4/(49x,'C',i2,'=',1PD15.8:))
  698 FORMAT(' uLR(r)  is a simple inverse-power sum with',6x,'C',I2,
     1  '=',1PD15.8:/(49x,'C',i2,'=',D15.8:))
  600 FORMAT(' Fit an ',A4,'p  potential function to the input points'/
     1  1x,26('==')/5x,'with initial   VMIN=',f13.4,'   Re=',f11.8:,
     2  '   De=',f11.4)
  601 FORMAT('     with expansion variable   y_',i1,'(r)= [r^',i1,' -',
     1 f7.4,'^',i1,']/[r^',i1,' +',f7.4,'^',i1,']' )
  602 FORMAT(' Use Lyon 2x2  uLR(r)  with   Aso=',F10.6,'   C_3(Sigma)='
     1  ,1PD15.7:/47x,'C_6(Sigma)=',D15.7:/47x,'C_6(Pi)   =',D15.7:/
     2  47x,'C_8(Sigma)=',1PD15.7:/47x,'C_8(Pi)   =',D15.7)
  603 FORMAT('     with expansion variable   y_',i1,'(r)= [r^',i1,
     1  ' - Re^',i1,']/[r^',i1,' + Re^',i1,']' )
  652 FORMAT( '   & define beta(y(r)) as a natural spline through points
     1 at the',i4,'  yp values:'/(2x,7F11.7))
  605 FORMAT(' GPEF expansion variable is:   y= (R^p - Re^p)/(',F9.5,
     1  '*R^p',1x,SP,F9.5,'*Re^p)')
  606 FORMAT(/' Fit to',I5,' input turning points assuming common uncert
     1ainty  u(VTP)=',1PD9.2/1x,39('--')/
     2  4('    RTP',7x,'VTP   ')/1x,39('--')/(4(0pF9.5,f11.3)))
  607 FORMAT('      which yields initial values of   AA=',1PD14.7,
     1   '   BB=',D14.7)
  608 FORMAT(1x,39('--')/)
  609 FORMAT(/' Fit to',I5,' input turning points with initial energy mi
     1nimum   VMIN=',f11.4/1x,32('--')/2(5x,'RTP',8x,'VTP',6x,'unc',4x)/
     2  1x,32('--')/(2(0PF10.5,F12.4,1PD10.2)))
  610 FORMAT(' Shift read-in VMIN value below lowest input potential val
     1ue',f11.4)
  614 FORMAT('  uLR multiplies  Cm/r^m  by incomplete gamma function  P(
     1m-1)')
c=======================================================================
c** Now ... loop over different {p,q,NS,NL} combinations till end of data
c**  p and q are powers used to define radial variable in the exponent
c       beta(r)= yp*betaINF + [1 - yp]*sum{beta_i*yq^i}  where  
c       ya=(R^a - AREF^a)/(R^a + AREF^a)   for  a= p  or  q
c** Read powers NS used in  \beta(y)  expansion for  r > Re  and 
c               NL used in  \beta(y)  expansion for  r = Re  or  > Re 
c** For an MLR potential with  SAP > 0,  (NS+NL+1) is the number of yp
c      values to be used to define the exponent spline used for beta(y):
c   NS specifies the number of (equally spaced) yp values for  yp < 0
c   NL specifies the number of yp values from  yp= 0  to  yp= 1
c** For a GPEF potential, consider powers ranging from  NS(.ge.0) to NL
c* RREF   defines the reference distance in the expansion variable
c      - for  RREF.le.0 , define parameter  RREF = Re
c      - for  RREF.gt.0 , fix parameter  RREF   at its read-in value
c-----------------------------------------------------------------------
   10 READ(5,*, END= 999) p, q, NS, NL, RREF
c-----------------------------------------------------------------------
      Re= ReIN
      De= DeIN
      VMIN= VMINin
      IF((p.LE.0).OR.(NS.LT.0).OR.(NL.LT.0)) GOTO 999
      NPHI= MAX(NS,NL)+ 1
      IF(PSEL.LE.3) WRITE(6,600) NNAME, VMIN, Re, De
      IF(PSEL.GE.4) WRITE(6,600) NNAME, VMIN, Re
      IF(PSEL.EQ.2) THEN
          IF(SAP.GT.0) WRITE(6,650) NS,NL
          IF(SAP.LE.0) THEN
c** For 'conventional' exponent polynomial \beta(r)
              IF((NCMM.GE.4).AND.(MMLR(2).EQ.0)) THEN
c... For Lyon treatment of A-state alkali dimers ...
                  WRITE(6,602) CmVAL(2),CmVAL(1),(CmVAL(m),m=3,NCMM)
                ELSE
c... For 'conventional' inverse-power sum uLR(r)
                  IF(bDAMP.LE.0) WRITE(6,698) (MMLR(m),CmVAL(m),
     1                                                      m= 1,NCMM)
                  IF(bDAMP.GT.0) WRITE(6,696) bDAMP,(MMLR(m),CmVAL(m),
     1                                                      m= 1,NCMM)
                ENDIF
              WRITE(6,634) p,p,q
              ENDIF
          ENDIF
      IF(PSEL.EQ.3) THEN
          WRITE(6,696) bDAMP, (MMLR(m),CmVAL(m),m=1, NCMM)
          IF(IDF.EQ.1) WRITE(6,614)
          ENDIF
      IF(PSEL.LE.3) THEN
          IF(RREF.gt.0.d0) THEN
              AREF= RREF
              WRITE(6,601) q,q,AREF,q,q,AREF,q
            ELSE
              AREF= Re
              WRITE(6,603) q,q,q,q,q
            ENDIF
          IF((PSEL.EQ.2).AND.(SAP.GT.0)) THEN
              YH= -Rsap/NS
              ypSAP(1)= Rsap
              DO  I= 2,NS
                  ypSAP(I)= ypSAP(I-1) + YH
                  ENDDO
              NPHI= NS+ 1
              ypSAP(NPHI)= 0.d0
              YH= 1.d0/NL
              DO  I= 1,NL-1
                  ypSAP(NPHI+I)= ypSAP(NPHI+I-1) + YH
                  ENDDO
              NPHI= NPHI+ NL
              ypSAP(NPHI)= 1.d0
              WRITE(6,652) NPHI,(ypSAP(i), i= 1,NPHI)
              ENDIF
          AREFp= AREF**p
          AREFq= AREF**q
          ENDIF
      DSE= VMIN
c** Scan input potential array to ensure  VMIN .le. {lowest input point}
      DO  i= 1,NTP
          DSE= DMIN1(DSE,VTP(i))
          ENDDO
      IF(DSE.LT.VMIN) THEN
          VMIN= DSE - 0.1d0
          WRITE(6,610) VMIN
          ENDIF

c=======================================================================
c** Preliminary linearized fit for an  EMOp  potential ...
c-----------------------------------------------------------------------
      IF((PSEL.EQ.1).OR.(PSEL.EQ.3).AND.(beta(0).le.0)) THEN
c ... first define ordinate array
          DO  i= 1,NTP
              RTPp= RTP(i)**p
              RTPq= RTP(i)**q
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              yq= (RTPq - AREFq)/(RTPq + AREFq)
              IF(RTP(i).GT.Re) THEN
                  betay(i)= - DLOG(1.d0 - DSQRT((VTP(i)-VMIN)/De))
                  IF(VTP(i).GT.uVTP(i)) THEN
                      Ubetay(i)= 0.5d0*uVTP(i)/(DSQRT((VTP(i)-VMIN)*De)
     1                                                - (VTP(i)-VMIN))
                    ELSE
                      Ubetay(i)= DSQRT(uVTP(i)/De)
                    ENDIF
                ELSE
                  betay(i)= - DLOG(1.d0 + DSQRT((VTP(i)-VMIN)/De))
                  IF(VTP(i).GT.uVTP(i)) THEN
                      Ubetay(i)= 0.5d0*uVTP(i)/(DSQRT((VTP(i)-VMIN)*De)
     1                                                + (VTP(i)-VMIN))
                    ELSE
                      Ubetay(i)= DSQRT(uVTP(i)/De)
                    ENDIF
                ENDIF
c ... next create partial derivative array for linearized fit ...
              yPOW= (RTP(i)- Re)
              DO  j= 1, NPHI
                  DYDP(i,j)= yPOW
                  IF((RTP(i).GT.Re).AND.(j.GT.NL+1)) DYDP(i,j)= 0.d0
                  IF((RTP(i).LE.Re).AND.(j.GT.NS+1)) DYDP(i,j)= 0.d0
                  yPOW= yPOW*yp
                  ENDDO
c%%
cc        if(i.eq.1) write(8,700) 
cc        write(8,702) rtp(i),yp,vtp(i),betay(i),Ubetay(i),
cc   1                                           (dydp(i,j),j=1,nphi)
cc700 format('  RTP     yp      VTP      beta*y      unc(beta*y) :',
cc   1  ' {dY/dp}')
cc702 format(f6.3,f8.4,f9.2,1P2d13.5:/(14x,5d13.5))
c%%
              ENDDO
          NPARM= NPHI
          CALL LLSQF(NTP,NPARM,MXDATA,MXPARM,betay,Ubetay,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
          IF(prFIT.GT.0) WRITE(6,620) NNAME,p,NS,NL,DSE,
     1                        ('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPARM)
          ENDIF

c=======================================================================
      IF(PSEL.EQ.2) THEN
c-----------------------------------------------------------------------
c*** Preliminary linearized fit for an  MLRp  potential ...
          IF(p.LE.(MMLR(NCMM)-MMLR(1))) WRITE(6,616) p, NCMM,
     1                                              MMLR(NCMM)-MMLR(1)
c** First define array of exponent values with uncertainties defined by 
c  the assumption that all potential values have equal uncertainties UNC
          IF((NCMM.GE.4).AND.(MMLR(2).EQ.0)) THEN
c** Aubert-Frecon based ULR(r)
              RE3= 1.d0/Re**3
              RE6= RE3*RE3
              C6adj= CmVAL(3) + 0.25D0*CmVAL(1)**2/De
              C9adj= 0.5d0*CmVAL(1)*C6adj/De
              T1= (0.5d0*CmVAL(1)+ (C6adj-CmVAL(4))*RE3)*RE3/3.d0
              IF(NCMM.GT.4) THEN
c ... extension for Aubert-Frecon Li2(A) {3,0,6,6,8,8} case ...
                  RE8= RE6/Re**2
                  T1= T1+ (CmVAL(5)- CmVAL(6))*RE8/3.d0
                  ENDIF
              T0= DSQRT((T1- CmVAL(2))**2 + 8.d0*T1**2)
              ULRe= 0.5d0*( -CmVAL(2)+ RE3*(1.5d0*CmVAL(1) 
     1           + RE3*(C6adj + CmVAL(4)))) + 0.5d0*T0 + C9adj*RE6*RE3
              IF(NCMM.GT.4) THEN
                  ULRe= ULRe+ 0.5d0*(CmVAL(5)+ CmVAL(6))*RE8
                  ENDIF
            ELSE
c** For normal inverse-power sum MLR/MLJ case, with or without damping
              IF(bDAMP.GT.0.d0) CALL DAMPIG(Re,bDAMP,MMLR(NCMM),DM,DMP,
     1                                                           DMPP)
              ULRe= 0.d0
              DO  m= 1,NCMM
                  T0= CmVAL(m)/Re**MMLR(m)
                  IF(bDAMP.GT.0.d0) T0= T0*DM(MMLR(m))
                  ULRe= ULRe + T0
                  ENDDO
            ENDIF
          betaINF= DLOG(2.d0*De/ULRe)
          WRITE(6,619) betaINF
          Rep= RE**p
          DO  i= 1, NTP
              RTPp= RTP(i)**p
              RTPq= RTP(i)**q
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              yq= (RTPq - AREFq)/(RTPq + AREFq)
              ypRE= (RTPp - Rep)/(RTPp + Rep)
              xSAP(i)= ypRE
              IF((NCMM.GE.4).AND.(MMLR(2).EQ.0)) THEN
c... for Aubert-Frecon {3,0,6,6} case ...
                  RTP3= 1.d0/RTP(i)**3
                  RTP6= RTP3*RTP3
                  T1= (0.5d0*CmVAL(1)+ (C6adj- CmVAL(4))*RTP3)*RTP3/3.d0
                  IF(NCMM.GT.4) THEN
c ... extension for Aubert-Frecon Li2(A) {3,0,6,6,8,8} case ...
                      RTP8= RTP6/RTP(i)**2
                      T1= T1+ (CmVAL(5)- CmVAL(6))*RTP8/3.d0
                      ENDIF
                  T0= DSQRT((T1- CmVAL(2))**2 + 8.d0*T1**2)
                  ULR= 0.5d0*( -CmVAL(2) + RTP3*(1.5d0*CmVAL(1)
     1                       + RTP3*(CmVAL(3) + CmVAL(4)))) + 0.5d0*T0
     2                                               + C9adj*RTP3*RTP6
                   IF(NCMM.GT.4) ULR=ULR+ 0.5d0*(CmVAL(5)+CmVAL(6))*RTP8
                ELSE
c... for normal MLR/MLJ case ... with or without damping
                  IF(bDAMP.GT.0.d0) 
     1                CALL DAMPIG(RTP(i),bDAMP,MMLR(NCMM),DM,DMP,DMPP)
                  ULR= 0.d0
                  DO  m= 1,NCMM
                      T0= CmVAL(m)/RTP(i)**MMLR(m)
                      IF(bDAMP.GT.0.d0) T0= T0*DM(MMLR(m))
                      ULR= ULR + T0
                      ENDDO
                ENDIF
              IF(RTP(i).GT.Re) THEN
                  betay(i)= - DLOG((1.d0 - DSQRT((VTP(i)-VMIN)/De))
     1                                                      *ULRe/ULR)
                  IF((VTP(i)-VMIN).GT.uVTP(i)) THEN
                      Ubetay(i)= 0.5d0*uVTP(i)
     1                      /(DSQRT((VTP(i)-VMIN)*De) - (VTP(i)-VMIN))
                    ELSE
                      Ubetay(i)= DSQRT(uVTP(i)/De)
                    ENDIF
                ELSE
                  betay(i)= - DLOG((1.d0 + DSQRT((VTP(i)-VMIN)/De))
     1                                                      *ULRe/ULR)
                  IF((VTP(i)-VMIN).GT.uVTP(i)) THEN
                      Ubetay(i)= 0.5d0*uVTP(i)
     1                      /(DSQRT((VTP(i)-VMIN)*De) + (VTP(i)-VMIN))
                    ELSE
                      Ubetay(i)= DSQRT(uVTP(i)/De)
                    ENDIF
                ENDIF
              IF(SAP.LE.0) THEN
c** Subtract the \beta_\infty term to yield polynomial for fitting
c... For Huang's MLR exponent function
                  betay(i)= betay(i)- betaINF*yp*ypRE
                  yPOW= ypRE*(1.d0- yp)
c... then create partial derivative array for linearized fit ...
                  DO  j= 1, NPHI
                      DYDP(i,j)= yPOW
                      IF((RTP(i).GT.Re).AND.(j.GT.NL+1)) DYDP(i,j)= 0.d0
                      IF((RTP(i).LE.Re).AND.(j.GT.NS+1)) DYDP(i,j)= 0.d0
                      yPOW= yPOW*yq
                      ENDDO
c%%%
cc                if(i.eq.1) write(8,700) 
cc                        write(8,702) rtp(i),yp,vtp(i),betay(i),
cc   1                                  Ubetay(i),(dydp(i,j),j=1,nphi)
cc                write(8,800) rtp(i),yp,ypRE,ULR,betay(i)
cc800 Format( f7.4,2f12.8,4(1Pd15.7))
c%%%
                  ENDIF
cc            write(6,800) rtp(i),yp,ypRE,vlr,betay(i), betay(i)/ypRE
cc800 Format( f7.4,2f12.8,5(1Pd15.7))
              ENDDO
          IF(SAP.LE.0) THEN
              CALL LLSQF(NTP,NPHI,MXDATA,MXPARM,betay,Ubetay,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
              IF(prFIT.GT.0) THEN
                  IF(RREF.GE.0.d0) WRITE(6,621) NNAME,p,q,RREF,NS,NL,
     1                    DSE,('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPHI)
                  IF(RREF.LT.0.d0) WRITE(6,623) NNAME,p,q,NS,NL,
     1                    DSE,('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPHI)
                  ENDIF
              ENDIF
          IF(SAP.GT.0) THEN
c** For Pashov exponent spline, use spline through linearized input 
c  exponent values to define spline point \beta(i) values !
              xSAP(NTP+1)= 1.d0
              betay(NTP+1)= betaINF
              CALL Lkoef(NTP+1,xSAP,rKL,MXDATA)
              DO  I= 1,NPHI
                  XX= ypSAP(I)
c... Now, use a spline through the exponent values defined by the input 
c    points to generate values of that exponent at the desired 
c    spline-definition points
                  PV(I)= 0.d0
                  DO  m= 1,NTP+1
                      PV(I)= PV(I) + 
     1                Scalc(XX,m,NTP+1,xSAP,rKL,MXDATA)*betay(m)/xSAP(m)
                      ENDDO
                  ENDDO
              NLIN= (NPHI+1)/2
              IF(prFIT.GT.0) WRITE(6,653) NNAME,p,
     1                      ((ypSAP(I),PV(I),I= J,NPHI,NLIN),J=1,NLIN)
  653 FORMAT(/' Linearized ',A4,'{p=',i1,'}-SAP treatment yields:'/
     1  (2('    ypSAP=',f10.6,'  beta(y)='f10.5) ))
c
c... Finally, create the fixed array of S(m,x) to deine the exponent
c   and its partial derivatives in subsequent fits ...
              CALL Lkoef(NPHI,ypSAP,rKL,MXDATA)
              DO  I= 1,NTP
                  DO  m= 1,NPHI
                      SAS(I,m)= Scalc(xSAP(I),m,NPHI,ypSAP,rKL,MXDATA)
                      ENDDO
                  ENDDO
              ENDIF
          ENDIF

c=======================================================================
c** Preliminary linearized fit for a  DELR  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.3) THEN
c ... NOTE need iteration to determine self-consistent  beta(0)  value
c  First generate  A & B  from input Re, De and trial beta(0)
          ITER= 0
   40     ULRe= 0.d0
          dULRe= 0.d0
          CALL DAMPIG(Re,bDAMP,MMLR(NCMM),DM,DMP,DMPP)
          DO  m= 1,NCMM
              FCT= dexp(-bDAMP*Re/MMLR(m)
     1               - 0.025*(bDAMP*Re)**2/DSQRT(DFLOAT(MMLR(m))))
              AA= CmVAL(m)/Re**MMLR(m)
              ULRe= ULRe+ AA*DM(MMLR(m))
              dULRe= dULRe+ AA*(DMP(MMLR(m)) - MMLR(m)*DM(MMLR(m))/Re)
              ENDDO
          AA= De - ULRe - dULRe/beta(0)
          BB= 2.d0*(De - ULRe) - dULRe/beta(0)
          WRITE(6,607) AA,BB
          RAT= 0.5d0*BB/AA
          UMAX= DSQRT(RAT**2 + (uVTP(i) + ULRe - DE)/AA)
          ReDE= Re- dlog(RAT)/beta(0)
          DO  i= 1,NTP
              ULR= 0.d0
              RTPp= RTP(i)**p
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              CALL DAMPIG(RTP(i),bDAMP,MMLR(NCMM),DM,DMP,DMPP)
              DO  m= 1,NCMM
                  ULR= ULR + DM(MMLR(m))*CmVAL(m)/RTP(i)**MMLR(m)
                  ENDDO
              FCT= (VTP(i) - VMIN + ULR - De)/AA + RAT**2
              IF(FCT.LT.0.d0) THEN
c** If estimate of ReDE off a bit and  FCT < 0 , ignore & deweight point
                    betay(i)= 0.d0
                    Ubetay(i)= 9.d99
                    GOTO 44
                    ENDIF
              FCT= DSQRT(FCT)
              IF(RTP(i).GT.ReDE) THEN
                  IF(RAT.GT.FCT) THEN
                      betay(i)= - DLOG(RAT - FCT)
                      IF((VTP(i)-VMIN).GT.uVTP(i)) THEN
                          Ubetay(i)= 0.5d0*uVTP(i)/(AA*(RAT- FCT)*FCT)
                        ELSE
                          Ubetay(i)= uVTP(i)/(AA*UMAX*(RAT- UMAX))
                        ENDIF
                    ELSE
c ... deweight away points for which \ln argument would be negative
                      betay(i)= 0.d0
                      Ubetay(i)= 9.d9
                    ENDIF
                ELSE
                  betay(i)= - DLOG(RAT + FCT)
                  IF((VTP(i)-VMIN).GT.uVTP(i)) THEN
                      Ubetay(i)= 0.5d0*uVTP(i)/(AA*(RAT+ FCT)*FCT)
                    ELSE
                      Ubetay(i)= uVTP(i)/(AA*UMAX*(RAT+ UMAX))
                    ENDIF
                ENDIF
c ... now create partial derivative array for linearized fit ...
   44         yPOW= (RTP(i)- Re)
              DO  j= 1, NPHI
                  DYDP(i,j)= yPOW
                  IF((RTP(i).GT.Re).AND.(j.GT.NL+1)) DYDP(i,j)= 0.d0
                  IF((RTP(i).LE.Re).AND.(j.GT.NS+1)) DYDP(i,j)= 0.d0
                  yPOW= yPOW*yp
                  ENDDO
              ENDDO
          NPARM= MAX(NS,NL)+ 1
          CALL LLSQF(NTP,NPARM,MXDATA,MXPARM,betay,Ubetay,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
          IF(DABS(PV(1)-BETA(0)).GT.PS(1)) THEN
              WRITE(6,644) beta(0),PV(1),PV(1)-beta(0),DSE
              beta(0)= PV(1)
              ITER= ITER+ 1
              IF(ITER.LE.10) GOTO 40
              WRITE(6,646) ITER
            ELSE
              WRITE(6,648) PV(1),PV(1)-beta(0)
            ENDIF
          IF(prFIT.GT.0) WRITE(6,620) NNAME,p,NS,NL,DSE,
     1                        ('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPARM)
          ENDIF 

      IF(PSEL.LE.3)  THEN
c======================================================================
c** Now ... do direct non-linear fit to potential values ... first with
c      Re and/or VMIN and De fixed, and then freeing them up too ...
c=======================================================================
c* FIRST optimize  BETA(j)'s (and VMIN) with  Re and De held fixed!
          DO  j= 1,NPHI
              BETA(j-1)= PV(j)
              IFXP(j)= 0
              ENDDO
          IF(SAP.GT.0) THEN
              IFXP(NPHI)= 1
              ENDIF
          NPARM= NPHI+ 3
          DO  j= NPHI+1, NPARM
              IFXP(j)= 1
              ENDDO
          NDGF= NTP- NPHI
          IF(IFXVMIN.LE.0) THEN
              IFXP(NPARM)= 0
              NDGF= NDGF-1
              ENDIF
          PV(NPHI+1)= Re
          PV(NPHI+2)= De
          PV(NPHI+3)= VMIN
          CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
          IF(prFIT.GT.0) THEN
              DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
              IF(SAP.LE.0) THEN
                  IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,q,AREF,NS,
     1                               NL,DRMSD,0,PV(1),PU(1),PS(1),DSE,
     2                                (j-1,PV(j),PU(j),PS(j),j=2,NPHI)
                  IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,q,NS,NL,
     1                                  DRMSD,0,PV(1),PU(1),PS(1),DSE,
     2                                (j-1,PV(j),PU(j),PS(j),j=2,NPHI)
                ELSE
                  IF(RREF.GT.0.d0) WRITE(6,626) NNAME,p,AREF,NS,
     1              NL,DRMSD,(j,ypSAP(j),j,PV(j),PU(j),PS(j),j=1,NPHI)
                  IF(RREF.LE.0.d0) WRITE(6,628) NNAME,p,NS,NL,
     1                 DRMSD,(j,ypSAP(j),j,PV(j),PU(j),PS(j),j=1,NPHI)
                 ENDIF
              IF(IFXVMIN.LE.0) THEN
                  WRITE(6,660) PV(NPHI+3),PU(NPHI+3),PS(NPHI+3)
                  VMIN= PV(NPHI+3)
                  ENDIF
              ENDIF
c ... the, if appropriate, set  Re  free too ...
          IF(IFXRe.LE.0) THEN
              NDGF= NDGF- 1
              IFXP(NPHI+1)= 0
              CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              Re= PV(NPHI+1)
              DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
              IF(prFIT.GT.0) THEN
                  IF(SAP.LE.0) THEN
                      IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,q,AREF,
     1                            NS,NL,DRMSD,0,PV(1),PU(1),PS(1),DSE,
     2                                (j-1,PV(j),PU(j),PS(j),j=2,NPHI)
                      IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,q,NS,
     1                               NL,DRMSD,0,PV(1),PU(1),PS(1),DSE,
     2                                (j-1,PV(j),PU(j),PS(j),j=2,NPHI)
                    ELSE
                      IF(RREF.GT.0.d0) WRITE(6,626) NNAME,p,AREF,NS,NL,
     1                 DRMSD,(j,ypSAP(j),j,PV(j),PU(j),PS(j),j=1,NPHI)
                      IF(RREF.LE.0.d0) WRITE(6,628) NNAME,p,NS,NL,DRMSD,
     1                       (j,ypSAP(j),j,PV(j),PU(j),PS(j),j=1,NPHI)
                    ENDIF
                  WRITE(6,662) PV(NPHI+1),PU(NPHI+1),PS(NPHI+1)
                  IF(IFXVMIN.LE.0)
     1                   WRITE(6,660) PV(NPHI+3),PU(NPHI+3),PS(NPHI+3)
                  ENDIF
              ENDIF
c ... then with Re fixed again, free De & VMIN (as well as the beta's)
          IF(IFXDe.LE.0) THEN
              DSEB= DSE
              IFXP(NPHI+1)= 1
              NDGF= NTP - NPHI
              IF(IFXDe.LE.0) THEN
                  IFXP(NPHI+2)= 0
                  NDGF= NDGF-1
                  ENDIF
              IF(IFXVMIN.LE.0) THEN
                  IFXP(NPHI+3)= 0
                  NDGF= NDGF-1
                  ENDIF
              CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              IF(IFXDe.LE.0) De= PV(NPHI+2)
              IF(IFXVMIN.LE.0) VMIN= PV(NPHI+3)
              IF((prFIT.GT.0).OR.(DSE.GT.DSEB*1.01)) THEN
                  IF(DSE.GT.DSEB*1.01) WRITE(6,654) DSEB,DSE
 
                  DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
                  IF(SAP.LE.0) THEN
                      IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,q,AREF,NS,
     1                               NL,DRMSD,0,PV(1),PU(1),PS(1),DSE,
     2                                (j-1,PV(j),PU(j),PS(j),j=2,NPHI)
                      IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,q,NS,NL,
     1                                  DRMSD,0,PV(1),PU(1),PS(1),DSE,
     2                                (j-1,PV(j),PU(j),PS(j),j=2,NPHI)
                    ELSE
                      IF(RREF.GT.0.d0) WRITE(6,626) NNAME,p,AREF,NS,NL,
     1                 DRMSD,(j,ypSAP(j),j,PV(j),PU(j),PS(j),j=1,NPHI)
                      IF(RREF.LE.0.d0) WRITE(6,628) NNAME,p,NS,NL,DRMSD,
     1                       (j,ypSAP(j),j,PV(j),PU(j),PS(j),j=1,NPHI)
                    ENDIF
                  IF(IFXRe.LE.0) 
     1                WRITE(6,662) PV(NPHI+1),PU(NPHI+1),PS(NPHI+1)
                  IF(IFXDe.LE.0)
     1                WRITE(6,630) PV(NPHI+2),PU(NPHI+2),PS(NPHI+2)
                  IF(IFXVMIN.LE.0)
     1                WRITE(6,660) PV(NPHI+3),PU(NPHI+3),PS(NPHI+3)
                  ENDIF
              ENDIF
c ... and finally ... fit to all three of  VMIN, De and Re
          IFXP(NPHI+1)= IFXRe
          IFXP(NPHI+2)= IFXDe
          IFXP(NPHI+3)= IFXVMIN
          PV(NPHI+1)= Re
          PV(NPHI+2)= De
          PV(NPHI+3)= VMIN
          NDGF= NTP- NPHI
          IF(IFXVMIN.LE.0) NDGF= NDGF-1
          IF(IFXDE.LE.0) NDGF= NDGF-1
          IF(IFXRE.LE.0) NDGF= NDGF-1
          CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
          DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
          IF(SAP.LE.0) THEN
              IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,q,AREF,NS,NL,DRMSD,
     1                                        0,PV(1),PU(1),PS(1),DSE,
     2                                (j-1,PV(j),PU(j),PS(j),j=2,NPHI)
              IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,q,NS,NL,DRMSD,
     1                                        0,PV(1),PU(1),PS(1),DSE,
     2                                (j-1,PV(j),PU(j),PS(j),j=2,NPHI)
            ELSE
              IF(RREF.GT.0.d0) WRITE(6,626) NNAME,p,AREF,NS,NL,DRMSD,
     1                       (j,ypSAP(j),j,PV(j),PU(j),PS(j),j=1,NPHI)
              IF(RREF.LE.0.d0) WRITE(6,628) NNAME,p,NS,NL,DRMSD,
     1                       (j,ypSAP(j),j,PV(j),PU(j),PS(j),j=1,NPHI)
              IF(RREF.GT.0.d0) WRITE(6,629) RREF
            ENDIF
          IF(PSEL.EQ.3) BETA(0)= PV(1)
          WRITE(6,662) PV(NPHI+1),PU(NPHI+1),PS(NPHI+1)
          WRITE(6,630) PV(NPHI+2),PU(NPHI+2),PS(NPHI+2)
          WRITE(6,660) PV(NPHI+3),PU(NPHI+3),PS(NPHI+3)
ccc Print [calc.-obs.]
          IF(prDIFF.gt.0) WRITE(6,730) (RTP(I),YD(I),YD(I)/uVTP(I),
     1                                                       I= 1,NTP)
  730 FORMAT(1x,39('==')/3x,3(3x,'RTP',4x,'[c-o] [c-o]/unc')/
     1  1x,39('--')/(1x,3(f10.5,f8.4,f7.2)))
ccc
          IF(IFXRe.LE.0) Re= PV(NPHI+1)
          IF(IFXDe.LE.0) De= PV(NPHI+2)
          IF(IFXVMIN.LE.0) VMIN= PV(NPHI+3)
          IF(PSEL.EQ.3) beta(0)= PV(1)
          ENDIF

c=======================================================================
c*** For case of a  GPEF  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.4) THEN
          DO  i= 1,NTP
              uVTP(i)= uVTP(i)
              betay(i)= VTP(i)
              ENDDO
          DO  NPHI= NS+1,NL+1
c*** Loop over expansion orders from NS to NL
              IFXP(NPHI+1)= IFXVMIN
              IFXP(NPHI+2)= IFXRe
              IFXP(NPHI+3)= 1
c.... first, do fully linearized fit to get trial expansion coefficients
              Rep= Re**p
              NPARM= NPHI
              IF(IFXVMIN.LE.0) NPARM= NPARM+ 1
              DO  i= 1, NTP
                  RTPp= RTP(i)**p
                  yp= (RTPp - Rep)/(as*RTPp + bs*Rep)
                  yPOW= yp
                  DO  j=1,NPHI
                      yPOW= yPOW*yp
                      DYDP(i,j)= yPOW
                      ENDDO
                  IF(IFXVMIN.LE.0) DYDP(i,NPARM)= 1.d0
                  ENDDO
              CALL LLSQF(NTP,NPARM,MXDATA,MXPARM,betay,uVTP,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
              IF(IFXVMIN.LE.0) VMIN= PV(NPHI+1)
              IF(prFIT.GT.0) THEN
                  WRITE(6,620) NNAME,p,NS,NL,DSE,
     1                         ('  c',j-1,PV(j),PU(j),PS(j),j= 1,NPHI)
                  IF(IFXVMIN.LE.0) WRITE(6,660) PV(NPHI+1),PU(NPHI+1),
     1                                                      PS(NPHI+1)
                  ENDIF
c.... then, proceed with fit to non-linear form
              NPARM= NPHI+2
              PV(NPHI+2)= Re
              IFXDe= 1
              IF(IFXRe.LE.0) THEN
c... If Re is to be free, first optimize it in fit to initial form
                  CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,prNLL,
     1                   IFXP,VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
                  IF(prFIT.GT.0) THEN
                      WRITE(6,636) 
                      WRITE(6,658) p,as,bs,DSE,(i-1,PV(i),PU(i),PS(i),
     1                                                       i=1,NPHI)
                      IF(IFXVMIN.LE.0) WRITE(6,660) PV(NPHI+1),
     1                                           PU(NPHI+1),PS(NPHI+1)
                      WRITE(6,662) PV(NPHI+2),PU(NPHI+2),PS(NPHI+2)
                      ENDIF
                  ENDIF
              IF(NPHI.GE.2) THEN
                  DO  j=2, NPHI
                      PV(j)= PV(j+1)/PV(1)
                      ENDDO
                  ENDIF
c ... IFXDe is a flag indicating fit to final  c0*y**2(1 + c1*y + ... )
              IFXDe= 0
              CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              WRITE(6,658) p,as,bs,DSE,(i-1,PV(i),PU(i),PS(i),i=1,NPHI)
              IF(IFXVMIN.LE.0) WRITE(6,660) PV(NPHI+1),PU(NPHI+1),
     1                                                     PS(NPHI+1)
              IF(IFXRE.LE.0) WRITE(6,662) PV(NPHI+2),PU(NPHI+2),
     1                                                     PS(NPHI+2)
              ENDDO
          ENDIF
c*** Step inward and check whether the repulsive inner potential wall
c   has inflection point or turnover.
      IF(SAP.LE.0) THEN
          RH= RTP(1)*1.0d-2
          RR= RTP(1)
          VV= VTP(1)
          RB= RR + RH
          J= -1
          CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VB,PV,PU,PS,RB)
          INFL= 0
          DO  I= 1,99
              J= J- 1
              RBB= RB
              VBB= VB
              RB= RR
              VB= VV
              RR= RR - RH
              IF(RR.LT.0.2d0) EXIT
              CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VV,PV,PU,PS,RR)
c
c     write(6,666) rr,vv, (vv-vb)/RH
c 666 format(f8.4,3f16.4)
c
              IF(INFL.GT.0) THEN
                  IF(VV.LT.VB) THEN
c... warning of potential inner-wall turnover
                      WRITE(6,672) RB,VB
                      WRITE(6,608)
                      GOTO 10
                      ENDIF
                ELSE
                  IF((VV-VB).LE.(VB-VBB)) THEN
c... warning of potential inner-wall inflection 
                      INFL= 1
                      WRITE(6,670) RB,VB
                      ENDIF
                ENDIF
              ENDDO
          ENDIF
      WRITE(6,608)
      GOTO 10
c-----------------------------------------------------------------------
  616 FORMAT(//' *** WARNING:   p=',i2,' .LE. [MMLR(',i1,')-MMLR(1)]=',
     1  i2,'  CAUTION !!!!!! ******')
  619 FORMAT(' Linearized fit uses    beta(INF)=',f10.6)
  620 FORMAT(/' Linearized ',A4,'{p=',i1,'; NS=',i2,', NL=',i2,'} fit yi
     1elds   DSE=',1Pd9.2/(3x,a4,'_{',i2,'}=',d19.11,' (+/-',d8.1,')   P
     2S=',d8.1))
  621 FORMAT(/' Linearized ',A4,'{p=',i1,', q=',i1,'; Rref=',f5.2,
     1 '; NS=',i2,', NL=',i2,'} fit yields   DSE=',1Pd9.2/(3x,a4,'_{',
     2 i2,'}=',d19.11,' (+/-',d8.1,')   PS=',d8.1))
  623 FORMAT(/' Linearized ',A4,'{p=',i1,', q=',i1,';  Rref=Re;  NS=',
     1 i2,', NL=',i2,'} fit yields   DSE=',1Pd9.2/(3x,a4,'_{',i2,'}=',
     2 d19.11,' (+/-',d8.1,')   PS=',d8.1))
  622 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',i1,'; Rref=',f5.2,
     1 ' ; NS=',i2,', NL=',I2,'}  potl:   dd=',1Pd9.2/
     2 '   beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,4x,'DSE=',
     3   D9.2/('   beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  624 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',I1,
     1 '; Rref= Re ; NS=',i2,', NL=',I2,'}  potl:     dd=',1Pd9.2/
     2 '   beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,4x'DSE=',
     3  D9.2/('   beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  626 FORMAT(/' Direct fit to ',A4,'{p=',i1,'; Rref=',f5.2,' ; NS=',i2,
     1  ', NL=',I2,'}  potential:    dd=',1Pd9.2/(' ypSAP{',i2,'}=',
     20PF11.7,'   beta_{',i2,'}=',1Pd20.12,' (+/-',d8.1,')   PS=',d8.1))
  628 FORMAT(/' Direct fit to ',A4,'{p=',i1,'; Rref= Re ; NS=',i2,
     1  ', NL=',I2,'}  potential:   DSE=',1Pd9.2/(' ypSAP{',i2,'}=',
     20PF11.7,'   beta_{',i2,'}=',1Pd20.12,' (+/-',d8.1,')   PS=',d8.1))
  629 FORMAT(8x,'Rref ='f15.12)
  630 FORMAT(9x,'De =',f13.6,' (+/-',f12.6,')      PS=',1pd8.1)
  634 FORMAT(' Use polynomial exponent fx:  beta(R)= betaINF*y_',i1,
     1 ' + (1-y_',I1,')*Sum{beta_i*[y_',i1,']^i}')
  650 FORMAT(' Use Pashov natural spline exponent based on', i4,'  yp va
     1lues for  y < 0'/41x,'and',i4,'  yp values for  y > 0')
  636 FORMAT(/' First perform full non-linear GPEF fit without taking ou
     1t common factor of c_0')
  644 FORMAT('   Update  beta_0  from',f9.6,'   to',f9.6,'   by',
     1  1Pd9.1,' :   DSE=',1PD8.1)
  646 FORMAT(' !!! CAUTION !!! Iteration to optimize  beta(0)  not conve
     1rged after',i3,' tries')
  648 FORMAT('   Converge on   beta_0=',f9.6,'   Next change=',1Pd9.1)
  654 FORMAT(/' *** PROBLEM *** freeing De makes DSE increase from',
     1  1PD9.2,' to',D9.2)
  658 FORMAT(/' Fit to GPEF{p=',i1,'} potential with   As=',F5.2,
     1  '   Bs=',F5.2,'   yields   DSE=',1Pd9.2/
     2  (5x,'c_{',i2,'} =',d18.10,' (+/-',d8.1,')     PS=',d8.1))
  660 FORMAT(7x,'VMIN =',f13.5,' (+/-',f12.6,')      PS=',1pd8.1)
  662 FORMAT(9x,'Re =',f13.9,' (+/-',f12.9,')      PS=',1pd8.1)
  670 FORMAT(1x,39('--')/ ' *** CAUTION *** inner wall has inflection at
     1   R=',F6.3,'   V=',1PD12.4)
  672 FORMAT(28x'and turns over at   R=',F6.3,'   V=',1PD12.4)
  999 STOP
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DYIDPJ(IDAT,NDATA,NPARM,IFXP,YC,PV,PD,PS,RMSR)
c** Subroutine to calculate potential function value YC at distance
c  RDIST= RTP(IDAT), and its partial derivatives w.r.t. the various 
c  potential parameters.  If  IDAT.LE.0, evaluate them at the distance
c  input as RDIST= RMSR.  If  ABS{IDAT}.LE.1  generate a new set of 
c  internal potential variables, while if  ABS{IDAT} > 1  use SAVEd values
c... [Must ensure that calculations based on the current UPDATED PV(j)]
c------------------------------------------------------------------------
      INTEGER MXDATA, MXPARM, MXMLR
      PARAMETER (MXDATA=1501, MXPARM=43, MXMLR= 15)
      INTEGER  i,j,m,IDAT, NPOW,NPARM,NDATA, IFXP(MXPARM),JFXRe,JFXDe,
     1  JFXVMIN
      REAL*8  YC,PV(NPARM),PD(NPARM),PS(NPARM),DM(MXMLR),DMP(MXMLR),
     1 DMPP(MXMLR),bDAMPR,RMSR,RTPp,RTPq,Rep,AREF,AREFp,AREFq,ype,dype,
     2 betaINF,yp,yq,yPOW,XP,XPW,DER,TCM,UM,TTMM,DERP,SUM,DSUM,AA,BB,
     3 FCT,FCT2,ULR,ULRe,dULRe,d2ULRe,DDER,T0,T0P,T1,RE3,RE6,RE8,RTP3,
     4 RTP6,RTP8,RDIST,C6adj,C9adj,dAAdRe,dBBdRe
c-----------------------------------------------------------------------
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR(MXMLR),p,q,NS,NL,
     1                                                        NPHI,SAP
      REAL*8 Re,De,VMIN,RREF,M2,as,bs,bDAMP,
     1  CmVAL(MXMLR),RTP(MXDATA),SAS(MXDATA,MXPARM)
      COMMON /DATABLK/Re,De,VMIN,RREF,M2,as,bs,bDAMP,CmVAL,RTP,SAS,PSEL,
     1          IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR,p,q,NS,NL,NPHI,SAP
c-----------------------------------------------------------------------
      SAVE JFXRe,JFXDe,JFXVMIN, AREF,AREFp,Rep,betaINF,AA,BB, ULRe,
     1  dULRe
c=======================================================================
      IF(ABS(IDAT).LE.1) THEN
          JFXRe= IFXP(NPHI+1)
          JFXDe= IFXP(NPHI+2)
          JFXVMIN= IFXP(NPHI+3)
          ENDIF
      RDIST= RTP(IDAT)
      IF(IDAT.LE.0) RDIST= RMSR
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
              IF((NCMM.GE.4).AND.(MMLR(2).EQ.0)) THEN
c** For Aubert-Frecon based  ULR(r)
                  RE3= 1.d0/Re**3
                  RE6= RE3*RE3
                  C6adj= CmVAL(3) + CmVAL(1)**2/(4.d0*DE)
                  C9adj= 0.5d0*CmVAL(1)*C6adj/De
                  T1= (0.5d0*CmVAL(1)+ (C6adj- CmVAL(4))*RE3)*RE3/3.d0
                  IF(NCMM.GT.4) THEN
c ... extension for Aubert-Frecon Li2(A) {3,0,6,6,8,8} case ...
                      RE8= RE6/Re**2
                      T1= T1+ (CmVAL(5)- CmVAL(6))*RE8/3.d0
                      ENDIF
                  T0= DSQRT((T1- CmVAL(2))**2 + 8.d0*T1**2)
                  ULRe= 0.5d0*( - CmVAL(2) + (1.5d0*CmVAL(1)
     1         + (C6adj + CmVAL(4))*RE3)*RE3) + 0.5d0*T0 + C9adj*RE6*RE3
                  IF(NCMM.GT.4) ULRe= ULRe+0.5d0*(CmVAL(5)+CmVAL(6))*RE8
                  T0P= (9.d0*T1-CmVAL(2))/T0
                  dULRe= -RE3*(0.25d0*CmVAL(1)*(9.d0 + T0P)
     1               + RE3*(C6adj*(3.d0 + T0P) + CmVAL(4)*(3.d0 - T0P)
     2                                           + RE3*9.d0*C9adj))/Re
                  IF(NCMM.GT.4) THEN
                      dULRe= dULRe -RE8*4.d0*(CmVAL(5) 
     1                *(3.d0 + T0P) + CmVAL(6)*(3.d0 - T0P))/(3.d0*Re)
                      ENDIF
                ELSE
c** For normal inverse-power sum MLR/MLJ case, with or without damping
                  IF(bDAMP.GT.0.d0) CALL DAMPIG(Re,bDAMP,MMLR(NCMM),DM,
     1                                                       DMP,DMPP)
                  ULRe= 0.d0
                  dULRe= 0.d0
                  DO  m= 1,NCMM
                      T0= CmVAL(m)/Re**MMLR(m)
                      IF(bDAMP.GT.0.d0) T0= T0*DM(MMLR(m))
                      ULRe= ULRe + T0
                      dULRe= dULRe - T0*MMLR(m)/Re
                      IF(bDAMP.GT.0.d0)
     1                dULRe= dULRe + T0*DMP(MMLR(m))/DM(MMLR(m))
                      ENDDO
                ENDIF
              betaINF= DLOG(2.d0*De/ULRe)
              ENDIF
          RTPp= RDIST**p
          RTPq= RDIST**q
          yp= (RTPp - AREFp)/(RTPp + AREFp)
          yq= (RTPq - AREFq)/(RTPq + AREFq)
          ype= (RTPp - Rep)/(RTPp + Rep)
          IF(SAP.GT.0) THEN
c*** Case of Pashov natural spline exponent ....
              NPOW= NPHI
c... Now, use a spline through the exponent values defined by the input
c    points to generate values of that exponent at the desired
c    spline-definition points
              XP= 0.d0
              DO  J= 1,NPOW
                  PD(J)= SAS(IDAT,J)
                  XP= XP + PV(J)*PD(J)
                  ENDDO
              ENDIF
          IF(SAP.LE.0) THEN
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
          IF((NCMM.GE.4).AND.(MMLR(2).EQ.0)) THEN
c** For Aubert-Frecon based  ULR(r)
              RTP3= 1.d0/RDIST**3
              RTP6= RTP3*RTP3
              T1= (0.5d0*CmVAL(1) + (C6adj - CmVAL(4))*RTP3)*RTP3/3.d0
              IF(NCMM.GT.4) THEN
c ... extension for Aubert-Frecon Li2(A) {3,0,6,6,8,8} case ...
                  RTP8= RTP6/RDIST**2
                  T1= T1+ (CmVAL(5)- CmVAL(6))*RTP8/3.d0
                  ENDIF
              T0= DSQRT((T1- CmVAL(2))**2 + 8.d0*T1**2)
              ULR= 0.5d0*( - CmVAL(2) + (1.5d0*CmVAL(1) + (C6adj 
     1            + CmVAL(4))*RTP3)*RTP3) + 0.5d0*T0 + C9adj*RTP3*RTP6
              IF(NCMM.GT.4) ULR= ULR + 0.5d0*(CmVAL(5)+ CmVAL(6))*RTP8
c... SKIP Re derivative corrections for all?
              T0P= (9.d0*T1-CmVAL(2))/T0
c             dULRdRe= -RTP3*(0.25d0*CmVAL(1)*(9.d0 + T0P)
c    1      + RTP3*(CmVAL(3)*(3.d0 + T0P) + CmVAL(4)*(3.d0 - T0P)))/Re
              ENDIF
          IF((NCMM.LE.1).OR.(MMLR(2).GT.0)) THEN
c** For normal inverse-power sum MLR/MLJ case, with or without damping
              IF(bDAMP.GT.0.d0) CALL DAMPIG(RDIST,bDAMP,MMLR(NCMM),
     1                                                    DM,DMP,DMPP)
              ULR= 0.d0
              DO  m= 1,NCMM
                  T0= CmVAL(m)/RDIST**MMLR(m)
                  IF(bDAMP.GT.0.d0) T0= T0*DM(MMLR(m))
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
              CALL DAMPIG(Re,bDAMP,MMLR(NCMM),DM,DMP,DMPP)
              ULRe= 0.d0
              dULRe= 0.d0
              d2ULRe= 0.d0
              DO  m= 1,NCMM
                  T0= CmVAL(m)/Re**MMLR(m)
                  ULRe= ULRe + T0*DM(MMLR(m))
                  dULRe= dULRe + T0*(DMP(MMLR(m))
     1                                   - DM(MMLR(m))*MMLR(m)/RE)
                  d2ULRe= d2ULRe + T0*(DMPP(MMLR(m))
     1                 - 2.d0*MMLR(m)*DMP(MMLR(m))/RE 
     2                 + MMLR(m)*(MMLR(m)+1.d0)*DM(MMLR(m))/RE**2)
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
          CALL DAMPIG(RDIST,bDAMP,MMLR(NCMM),DM,DMP,DMPP)
          ULR= 0.d0
          DO  m= 1,NCMM
              ULR= ULR + DM(MMLR(m))*CmVAL(m)/RDIST**MMLR(m)
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

c***********************************************************************
      SUBROUTINE DAMPIG(r,b,MMAX,DM,DMP,DMPP)
c** Subroutine to generate values DM(m) and the first and second radial
c  derivatives DMP(m) and DMPP(m) of normalized incomplete gamma 
c  functions of orders m=0 to MMAX, at radial distance 'r', for damping
c  parameter 'b'.  NOTE that  DM(m)= {Tang-Toennies function}(m+1).
** NOTE ... this algorithm becomes unstable for small 'b*r' (say b*r < 0.2)
c***********************************************************************
      INTEGER MMAX,I,m
      REAL*8 b,r,br,XP,SSm,SSm1,SSm2,TK,DM(MMAX),DMP(MMAX),DMPP(MMAX)
      br= b*r
      XP= DEXP(-br)
      SSm= 0.d0
      SSm1= 0.d0
      SSm2= 0.d0
      TK= 1.d0
      DO  m=1, MMAX
          SSm2= SSm1
          SSm1= SSm
          SSm= SSm+ TK
          DM(m)= 1.d0 - XP*SSM
          DMP(m)= b*XP*(SSm - SSm1)
          DMPP(m)= b**2 *XP*(2.d0*SSm1 - SSm - SSm2)
          TK= TK*br/DFLOAT(m)
          ENDDO
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE NLLSSRR(NDATA,NPTOT,NPMAX,IROUND,ROBUST,LPRINT,IFXP,
     1                           YO,YU,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c**  Program for performing linear or non-linear least-squares fits and
c  (if desired) automatically using sequential rounding and refitting 
c  to minimize the numbers of parameter digits which must be quoted [see
c  R.J. Le Roy, J.Mol.Spectrosc. 191, 223-231 (1998)].         21/08/04
c
c  21/08/04 test version ... attempting to stablize non-linear fits by
c               scaling back parameter changes by factor of 4 or 16
c         [ & corrects  CM  for constrained-parameter cases!]
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c             COPYRIGHT 1998-2009  by  Robert J. Le Roy                +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Program uses orthogonal decomposition of the "design" (partial 
c  derivative) matrix for the core locally linear (steepest descent) 
c  step, following a method introduced (to me) by Dr. Michael Dulick. 
c** If no parameters are free (NPTOT=0), simply return RMS(residuals) as
c  calculated from the input parameter values {PV(j)}.
c** A user MUST SUPPLY subroutine  DYIDPJ  to generate the predicted
c  value of each datum and the partial derivatives of each datum w.r.t.
c  each parameter (see below) from the current trial parameters.
c
c** On entry: 
c    NDATA  is the number of data to be fitted 
c    NPTOT  the total number of parameters in the model (.le.NPMAX).
c           If NPTOT.le.0 , assume  YD(i)=YO(i)  and calculate the (RMS 
c           dimensionless deviation)=DSE  from them & YU(i) 
c    NPMAX is the maximum number of model parameters allowed by current
c          external array sizes.  Should set internal NPINTMX = NPMAX 
c          (may be freely changed by the user).
c    IROUND .ne. 0  causes Sequential Rounding & Refitting to be 
c             performed, with each parameter being rounded at the 
c            |IROUND|'th sig. digit of its local incertainty.
c        > 0  rounding selects in turn remaining parameter with largest
c             relative uncertainy
c        < 0  round parameters sequentially from last to first
c        = 0  simply stops after full convergence (without rounding).
c    ROBUST > 0  causes fits to use Watson's ``robust'' weighting  
c        1/[u^2 +{(c-o)^2}/3].  ROBUST > 1 uses normal 1/u^2 on first
c        fit cycle and  'robust' on later cycles.
c    LPRINT  specifies the level of printing inside NLLSSRR
c          if: =  0, no print except for failed convergence.
c               < 0  only converged, unrounded parameters, PU & PS's
c              >= 1  print converged parameters, PU & PS's
c              >= 2  also print parameter change each rounding step
c              >= 3  also indicate nature of convergence
c              >= 4  also print convergence tests on each cycle
c              >= 5  also parameters changes & uncertainties, each cycle
c    IFXP(j)  specifies whether parameter  j  is to be held fixed
c           [IFXP > 0] or to be freely varied in the fit [IFXP= 0]
c    YO(i)  are the NDATA 'observed' data to be fitted  
c    YU(i)  are the uncertainties in these YO(i) values
c    PV(j)  are initial trial parameter values (for non-linear fits);  
c           should be set at zero for initially undefined parameters.
c
c** On Exit:   
c    YD(i)  is the array of differences  [Ycalc(i) - YO(i)]
c    PV(j)  are the final converged parameter values
c    PU(j)  are 95% confidence limit uncertainties in the PV(j)'s
c    PS(j)  are 'parameter sensitivities' for the PV(j)'s, defined such 
c           that the RMS displacement of predicted data  due to rounding
c           off parameter-j by PS(j) is .le. DSE/10*NPTOT
c    CM(j,k)  is the correlation matrix obtained by normalizing variance
c           /covariance matrix:  CM(j,k) = CM(j,k)/SQRT[CM(j,j)*CM(k,k)]
c    TSTPS = max{|delta[PV(j)]/PS(j)|}  is the parameter sensitivity 
c          convergence test:  delta[PV(j)] is last change in parameter-j
c    TSTPU = max{|delta[PV(j)]/PU(j)|}  is the parameter uncertainty 
c          convergence test:  delta[PV(j)] is last change in parameter-j
c    DSE    is the predicted (dimensionless) standard error of the fit
c
c  NOTE that the squared 95% confidence limit uncertainty in a property 
c  F({PV(j)}) defined in terms of the fitted parameters {PV(j)} (where
c  the L.H.S. involves  [row]*[matrix]*[column]  multiplication) is:
c  [D(F)]^2 = [PU(1)*dF/dPV(1), PU(2)*dF/dPV(2), ...]*[CM(j,k)]*
c                              [PU(2)*dF/dPV(1), PU(2)*dF/dPV(2), ...]
c
c** Externally dimension:  YO, YU and YD  .ge. NDATA 
c             PV, PU  and  PS  .ge.  NPTOT (say as NPMAX), 
c             CM   as a square matrix with column & row length  NPMAX
c***********************************************************************
      INTEGER NPINTMX
      PARAMETER (NPINTMX=2000)
      INTEGER I,J,K,L,IDF,ITER,NITER,IROUND,ISCAL,JROUND,LPRINT,NDATA,
     1 NPTOT,NPMAX,NPARM,NPFIT,JFIX,QUIT,ROBUST,
     2 IFXP(NPMAX),JFXP(NPINTMX)
      REAL*8  YO(NDATA), YU(NDATA), YD(NDATA), PV(NPTOT), PU(NPTOT), 
     1 PS(NPTOT),PSS(NPINTMX),PC(NPINTMX),PCS(NPINTMX),PX(NPINTMX),
     2 PY(NPINTMX),CM(NPMAX,NPMAX), F95(10),
     3 RMSR, RMSRB, DSE, TSTPS, TSTPSB, TSTPU, TFACT, S, UU, Zthrd
      DATA F95/12.7062D0,4.3027D0,3.1824D0,2.7764D0,2.5706D0,2.4469D0,
     1  2.3646D0,2.3060D0,2.2622D0,2.2281D0/
c
      IF((NPTOT.GT.NPMAX).OR.(NPTOT.GT.NPINTMX)
     1                    .OR.(NPTOT.GT.NDATA)) THEN
c** If array dimensioning inadequate, print warning & then STOP
          WRITE(6,602) NPTOT,NPINTMX,NPMAX,NDATA
          STOP
          ENDIF
      Zthrd= 0.d0
      IF(ROBUST.GE.2) Zthrd= 1.d0/3.d0
      TSTPS= 0.d0
      RMSR= 0.d0
      NITER= 0
      QUIT= 0
      NPARM= NPTOT
      DO J= 1, NPTOT
          PS(J)= 0.d0
          JFXP(J)= IFXP(J)
          IF(IFXP(J).GT.0) NPARM= NPARM- 1
          ENDDO
      NPFIT= NPARM
      JROUND= IABS(IROUND)
c=======================================================================
c** Beginning of loop to perform rounding (if desired).  NOTE that in 
c  sequential rounding, NPARM is the current (iteratively shrinking) 
c  number of free parameters. 
    6 IF(NPARM.GT.0) TSTPS= 9.d99
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
c======================================================================
c** Begin iterative convergence loop:  try for up to 30 cycles
      DO 50 ITER= 1, 30
          ISCAL= 0
          NITER= NITER+ 1
          DSE= 0.d0 
          TSTPSB= TSTPS
          RMSRB= RMSR
c** Zero out various arrays
   10     IF(NPARM.GT.0) THEN
c%%           DO  I = 1,NPARM
              DO  I = 1,NPTOT
c** PSS is the array of Saved Parameter Sensitivities from previous 
c   iteration to be carried into dyidpj subroutine - used in predicting
c   increment for derivatives by differences.
                  PSS(I)= PS(I)
c** PCS is the saved array of parameter changes from previous iteration
c   to be used (if necessary) to attempt to stablize fit
                  PCS(I)= PC(I)
                  PS(I) = 0.D0
                  PU(I) = 0.D0
                  PX(I) = 0.D0
                  PY(I) = 0.D0
                  DO  J = 1,NPARM
                      CM(I,J) = 0.D0
                      ENDDO
                  ENDDO
              ENDIF
c
c========Beginning of core linear least-squares step====================
c
c** Begin by forming the Jacobian Matrix from partial derivative matrix
          DO  I = 1,NDATA
c** User-supplied subroutine DYIDPJ uses current (trial) parameter 
c  values {PV} to generate predicted datum # I [y(calc;I)=UU] and its
c  partial derivatives w.r.t. each of the parameters, returning the 
c  latter in 1-D array PC.  See dummy sample version at end of listing.
c* NOTE 1: if more convenient, DYIDPJ could prepare the y(calc) values 
c     and derivatives for all data at the same time (when I=1), but only
c     returned the values here one datum at a time (for I > 1).]
c* NOTE 2: the partial derivative array PC returned by DYIDPJ must have
c     an entry for every parameter in the model, though for parameters 
cyPOW     which are held fixed [JFXP(j)=1], those PC(j) values are ignored.
              CALL DYIDPJ(I,NDATA,NPTOT,JFXP,UU,PV,PC,PSS,RMSR)
              IF(NPARM.LT.NPTOT) THEN
c** For constrained parameter or sequential rounding, collapse partial 
c   derivative array here
                  DO  J= NPTOT,1,-1
                      IF(JFXP(J).GT.0) THEN
                          IF(J.LT.NPTOT) THEN
                              DO  K= J,NPTOT-1
                                  PC(K)= PC(K+1)
                                  ENDDO
                              ENDIF
                          PC(NPTOT)= 0.d0
                          ENDIF
                      ENDDO
                  ENDIF
              YD(I)= UU - YO(I)
              S = 1.D0 / YU(I)
              IF(Zthrd.GT.0.d0) S= 1.d0/DSQRT(YU(I)**2 + Zthrd*YD(I)**2)
              UU = - YD(I) * S
              DSE= DSE+ UU*UU
              IF(NPARM.GT.0) THEN
                  DO  J = 1,NPARM
                      PC(J) = PC(J)*S
                      PS(J) = PS(J)+ PC(J)**2
                      ENDDO
                  CALL QROD(NPARM,NPMAX,NPMAX,CM,PC,PU,UU,PX,PY)
                  ENDIF
              ENDDO
          RMSR= DSQRT(DSE/NDATA)
          IF(NPARM.LE.0) GO TO 60
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
c** Solve for parameter changes  PC(j)
          DO  I = 1,NPARM
              J = NPARM - I + 1
              PC(J) = 0.D0
              DO  K = J,NPARM
                  PC(J) = PC(J) + CM(J,K) * PU(K)
                  ENDDO
              ENDDO
c
c** Get (upper triangular) "dispersion Matrix" [variance-covarience 
c  matrix  without the sigma^2 factor].
          DO  I = 1,NPARM
              DO  J = I,NPARM
                  UU = 0.D0
                  DO  K = J,NPARM
                      UU = UU + CM(I,K) * CM(J,K)
                      ENDDO
                  CM(I,J) = UU
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
              ENDDO
c
c** Generate standard error  DSE = sigma^2,  and prepare to calculate 
c  Parameter Sensitivities PS
          IF(NDATA.GT.NPARM) THEN
              DSE= DSQRT(DSE/(NDATA-NPARM))
            ELSE
              DSE= 0.d0
            ENDIF
c** Use DSE to get final (95% confid. limit) parameter uncertainties PU
c** Calculate 'parameter sensitivities', changes in PV(j) which would 
c  change predictions of input data by an RMS average of  DSE*0.1/NPARM
          UU= DSE*0.1d0/DFLOAT(NPARM)
          S= DSE*TFACT
          DO  J = 1,NPARM
              PU(J)= S* PU(J)
              PS(J)= UU*DSQRT(NDATA/PS(J))
              ENDDO

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          IF((ITER.GT.1).AND.(RMSR.GT.2.0d0*RMSRB).AND.(ISCAL.LE. 3))
     1                                                            THEN
c** LeRoy's Marquardt-like attempt to damp changes if RMSR increases ...
              ISCAL= ISCAL+ 1
              IF(LPRINT.GE.0) THEN
                  WRITE(6,620) ITER,RMSR,RMSR/RMSRB,ISCAL
  620 FORMAT(' At Iteration',i3,'  RMSD=',1PD8.1,'  RMSD/RMSDB=',D8.1,
     1 "  Scale PC by  (1/4)**",i1)
ccc               WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),J=1,NPTOT)
                  ENDIF
              DO  J= 1,NPTOT
                  PC(J)= 0.25d0*PCS(J)
                  PV(J)= PV(J)- 3.d0*PC(J)
                  ENDDO
              GOTO 10
              ENDIF
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c========End of core linear least-squares step==========================
c ... early exit if Rounding cycle finished ... 
          IF(QUIT.GT.0) GO TO 54
c
c** Next test for convergence 
          TSTPS= 0.D0
          TSTPU= 0.D0
          DO  J= 1, NPARM
              TSTPS= MAX(TSTPS,DABS(PC(J)/PS(J)))
              TSTPU= MAX(TSTPU,DABS(PC(J)/PU(J)))
              ENDDO
          IF(LPRINT.GE.4) WRITE(6,604) ITER,RMSR,TSTPS,TSTPU
c** Now ... update parameters (careful about rounding)
          DO  J= 1,NPTOT
              IF(JFXP(J).GT.0) THEN
c** If parameter held fixed (by input or rounding process), shift values
c   of change, sensitivity & uncertainty to correct label.
                  IF(J.LT.NPTOT) THEN
                      DO  I= NPTOT,J+1,-1
                          PC(I)= PC(I-1)
                          PS(I)= PS(I-1)
                          PU(I)= PU(I-1)
                          ENDDO
                      ENDIF
                  PC(J)= 0.d0
                  PS(J)= 0.d0
                  PU(J)= 0.d0
                ELSE
                  PV(J)= PV(J)+ PC(J)
                ENDIF
              ENDDO
          IF(LPRINT.GE.5) WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),
     1                                                      J=1,NPTOT)
          IF(NITER.GT.1) THEN
c** Test for convergence:  for every parameter desire:
c  |parameter change| < |parameter sensitivity|,  but after iteration #5
c  STOP iterating if  Max{|change/sens.|} increases AND 
c  Max{|change/unc.|} < 0.01
              IF(TSTPS.GT.1.d0) THEN
                  IF((RMSR.GT.RMSRB).AND.(ITER.GT.5)) THEN
                      IF((TSTPU.LT.1.d-2).OR.((TSTPU.LT.0.5d0).AND.
     1                                             (ITER.GT.10))) THEN
                          IF(LPRINT.GE.3) WRITE(6,606) ITER,TSTPU,RMSR
                          GO TO 54
                          ENDIF
                      ENDIF
                ELSE
                  IF(LPRINT.GE.3) WRITE(6,608) ITER,TSTPS,RMSR
                  GO TO 54
                ENDIF
              ENDIF
cc        CALL FLUSH(6)
          IF(ROBUST.GT.0) Zthrd= 1.d0/3.d0
   50     CONTINUE
      WRITE(6,610) NPARM,NDATA,ITER,RMSR,TSTPS,TSTPU
c** End of iterative convergence loop for (in general) non-linear case.
c======================================================================
c
   54 IF(NPARM.LT.NPTOT) THEN
c** If necessary, redistribute correlation matrix elements to full 
c  NPTOT-element correlation matrix
          DO  J= 1,NPTOT
              IF(JFXP(J).GT.0) THEN  
c* If parameter J was held fixed
                  IF(J.LT.NPTOT) THEN
c ... then move every lower CM element down one row:
                      DO  I= NPTOT,J+1,-1
c ... For  K < J, just shift down or over to the right
                          IF(J.GT.1) THEN
                              DO  K= 1,J-1
                                  CM(I,K)= CM(I-1,K) 
                                  CM(K,I)= CM(I,K)
                                  ENDDO
                              ENDIF
c ... while for  K > J  also shift elements one column to the right
                          DO  K= NPTOT,J+1,-1
                              CM(I,K)= CM(I-1,K-1)
                              ENDDO
                          ENDDO
                      ENDIF
c ... and finally, insert appropriate row/column of zeros ....
                  DO  I= 1,NPTOT
                      CM(I,J)= 0.d0
                      CM(J,I)= 0.d0
                      ENDDO 
                  CM(J,J)= 1.d0
                  ENDIF
              ENDDO
          ENDIF
      IF(QUIT.GT.0) GOTO 60
      IF(NPARM.EQ.NPFIT) THEN
c** If desired, print unrounded parameters and fit properties
          IF(LPRINT.NE.0) THEN
              WRITE(6,616) NDATA,NPARM,RMSR,TSTPS
              WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),J=1,NPTOT)
              ENDIF
          ENDIF
      IF(IROUND.EQ.0) RETURN
c** Automated 'Sequential Rounding and Refitting' section:  round 
c  selected parameter, fix it, and return (above) to repeat fit.
      IF(IROUND.LT.0) THEN
c ... if IROUND < 0, sequentially round off 'last' remaining parameter
          DO  J= 1, NPTOT
              IF(JFXP(J).LE.0) THEN
                  JFIX= J
                  ENDIF
              ENDDO
        ELSE
c ... if IROUND > 0, sequentially round off remaining parameter with
c                    largest relative uncertainty.
c ... First, select parameter JFIX with the largest relative uncertainty
          K= 0
          TSTPS= 0.d0
          DO  J= 1,NPTOT
              IF(JFXP(J).LE.0) THEN
                  K= K+1
                  TSTPSB= DABS(PU(J)/PV(J))
                  IF(TSTPSB.GT.TSTPS) THEN
                      JFIX= J
                      TSTPS= TSTPSB
                      ENDIF
                  ENDIF
              ENDDO 
        ENDIF
      UU= PV(JFIX)
      CALL ROUND(JROUND,NPMAX,NPTOT,NPTOT,JFIX,PV,PU,PS,CM)
      JFXP(JFIX)= 1
      IF(LPRINT.GE.2)
     1            WRITE(6,614) JFIX,UU,PU(JFIX),PS(JFIX),PV(JFIX),RMSR
      NPARM= NPARM-1
      IF(NPARM.EQ.0) THEN
c** After rounding complete, make one more pass with all non-fixed 
c  parameters set free to get full correct final correlation matrix, 
c  uncertainties & sensitivities
          NPARM= NPFIT
          QUIT= 1
          DO  J= 1,NPTOT
              JFXP(J)= IFXP(J)
              ENDDO
c ... reinitialize for derivative-by-differences calculation
          RMSR= 0.d0
          ENDIF
      GO TO 6
c
c** If no parameters varied or sequential rounding completed - simply 
c   calculate DSE from RMS residuals and return.
   60 DSE= 0.d0
      IF(NDATA.GT.NPFIT) THEN
          DSE= RMSR*DSQRT(DFLOAT(NDATA)/DFLOAT(NDATA-NPFIT))
        ELSE
          DSE= 0.d0
        ENDIF
      IF(NPFIT.GT.0) THEN
          IF(LPRINT.GT.0) THEN
c** Print final rounded parameters with original Uncert. & Sensitivities
              IF(QUIT.LT.1) WRITE(6,616) NDATA, NPFIT, RMSR, TSTPS
              IF(QUIT.EQ.1) WRITE(6,616) NDATA, NPFIT, RMSR
              DO  J= 1, NPTOT
                  IF(JFXP(J).GT.0) THEN
c** If parameter held fixed (by rounding process), shift values of
c   change, sensitivity & uncertainty to correct absolute number label.
                      DO  I= NPTOT,J+1,-1
                          PC(I)= PC(I-1)
                          PS(I)= PS(I-1)
                          PU(I)= PU(I-1)
                          ENDDO
                      PC(J)= 0.d0
                      PS(J)= 0.d0
                      PU(J)= 0.d0
                      ENDIF
                  ENDDO
              WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),J=1,NPTOT)
              ENDIF
          ENDIF
      RETURN
c
  602 FORMAT(/' *** NLLSSRR problem:  [NPTOT=',i4,'] > min{NPINTMX=',
     1  i4,' NPMAX=',i4,', NDATA=',i6,'}')
  604 FORMAT(' After Cycle #',i2,':   DRMSD=',1PD10.3,'    test(PS)=',
     1  1PD8.1,'    test(PU)=',D8.1)
  606 FORMAT(' Effective',i3,'-cycle Cgce:  MAX{|change/unc.|}=',1PD8.1,
     1  ' < 0.01   DRMSD=',D10.3)
  608 FORMAT(' Full',i3,'-cycle convergence:  Max{|change/sens.|}=',
     1  1PD8.1,' < 1   DRMSD=',D10.2)
  610 FORMAT(' !! CAUTION !! fit of',i4,' parameters to',I6,' data not c
     1onverged after',i3,' Cycles'/5x,'DRMS(deviations)=',1PD10.3,
     2 '    test(PS) =',D9.2,'    test(PU) =',D9.2/1x,31('**'))
  612 FORMAT((4x,'PV(',i4,') =',1PD22.14,' (+/-',D8.1,')    PS=',d8.1,
     1  '   PC=',d8.1))
  614 FORMAT(' =',39('==')/' Round Off  PV(',i4,')=',1PD21.13,' (+/-',
     1 D9.2,')    PS=',d9.2/10x,'fix it as ',D21.13,'  & refit:  DRMS(de
     2viations)=',D10.3)
  616 FORMAT(/i6,' data fit to',i5,' param. yields   DRMS(devn)=',G11.4:
     1 '  tst(PS)=',1Pd8.1)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
	SUBROUTINE QROD(N,NR,NC,A,R,F,B,GC,GS)
C** Performs ORTHOGONAL DECOMPOSITION OF THE LINEAR LEAST-SQUARES    
C            EQUATION J * X = F TO A * X = B(TRANSPOSE) * F WHERE   
C            J IS THE JACOBIAN IN WHICH THE FIRST N ROWS AND COLUMNS
C            ARE TRANSFORMED TO THE UPPER TRIANGULAR MATRIX A      
C            (J = B * A), X IS THE INDEPENDENT VARIABLE VECTOR, AND
C            F IS THE DEPENDENT VARIABLE VECTOR. THE TRANSFORMATION
C            IS APPLIED TO ONE ROW OF THE JACOBIAN MATRIX AT A TIME.
C  Parameters:                                                   
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
      DO 10 I = 1,N
          Z(1) = R(I)
          J = I - 1
          DO  K = 1,J
              Z(2) = GC(K) * A(K,I) + GS(K) * Z(1)
              Z(1) = GC(K) * Z(1) - GS(K) * A(K,I)
              A(K,I) = Z(2)
              ENDDO
          GC(I) = 1.D0
          GS(I) = 0.D0
          IF(DABS(Z(1)) .LE. 0.D0) GOTO 10
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
   10     CONTINUE
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ROUND(IROUND,NPMAX,NPARM,NPTOT,IPAR,PV,PU,PS,CM)
c** Subroutine to round off parameter # IPAR with value PV(IPAR) at the
c  |IROUND|'th significant digit of:  [its uncertainty  PU(IPAR)] . 
c** On return, the rounded value replaced the initial value  PV(IPAR).
c** Then ... use the correlation matrix CM and the uncertainties PU(I)
c  in the other (NPTOT-1) [or (NPARM-1) free] parameters to calculate 
c  the optimum compensating changes PV(I) in their values.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 1998-2009  by  Robert J. Le Roy             +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER    IROUND,NPMAX,NPARM,NPTOT,IPAR,I,IRND,KRND
      REAL*8  PU(NPMAX),PS(NPMAX),PV(NPMAX),CM(NPMAX,NPMAX),CNST,
     1        CRND,XRND,FCT,Z0
      DATA Z0/0.d0/
      CNST= PV(IPAR)
      XRND= DLOG10(PU(IPAR))
c** If appropriate, base last rounding step on sensitivity (not uncert.)
      IF((NPARM.EQ.1).AND.(PS(IPAR).LT.PU(IPAR))) XRND= DLOG10(PS(IPAR))
c** First ... fiddle with log's to perform the rounding
      IRND= INT(XRND)
      IF(XRND.GT.0) IRND=IRND+1
      IRND= IRND- IROUND
      FCT= 10.D0**IRND
      CRND= PV(IPAR)/FCT
      XRND= Z0
c ... if rounding goes past REAL*8 precision, retain unrounded constant
      IF(DABS(CRND).GE.1.D+16) THEN
          WRITE(6,601) IROUND,IPAR
           RETURN
           ENDIF
      IF(DABS(CRND).GE.1.D+8) THEN
c ... to avoid problems from overflow of I*4 integers ...
          KRND= NINT(CRND/1.D+8)
          XRND= KRND*1.D+8
          CRND= CRND-XRND
          XRND= XRND*FCT
          END IF
      IRND= NINT(CRND)
      CNST= IRND*FCT+ XRND
c????????????????
c** Zero parameters more aggressively ... if unc. > 2* value
        if(dabs(PU(IPAR)/PV(IPAR)).GT.2.d0) then
            cnst= 0.d0
            endif
c????????????????
c** Now ... combine rounding change in parameter # IPAR, together with
c  correlation matrix CM and parameter uncertainties PU to predict
c  changes in other parameters to optimally compensate for rounding off
c  of parameter-IPAR.  Method pointed out by Mary Thompson (Dept. of
c  Statistics, UW),
      IF(IPAR.GT.1) THEN
          XRND= (CNST-PV(IPAR))/PU(IPAR)
          DO  I= 1,NPTOT
              IF(I.NE.IPAR) THEN
                  PV(I)= PV(I)+ CM(IPAR,I)*PU(I)*XRND
                  ENDIF
              ENDDO
          ENDIF
      PV(IPAR)= CNST
      RETURN
  601 FORMAT(' =',39('==')/' Caution:',i3,'-digit rounding of parameter-
     1',i2,' would exceed (assumed) REAL*8'/' ********   precision overf
     2low at 1.D+16, so keep unrounded constant')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
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
	SUBROUTINE QQROD(N,NR,NC,A,R,F,B,GC,GS)
C** Performs ORTHOGONAL DECOMPOSITION OF THE LINEAR LEAST-SQUARES    
C            EQUATION J * X = F TO A * X = B(TRANSPOSE) * F WHERE   
C            J IS THE JACOBIAN IN WHICH THE FIRST N ROWS AND COLUMNS
C            ARE TRANSFORMED TO THE UPPER TRIANGULAR MATRIX A      
C            (J = B * A), X IS THE INDEPENDENT VARIABLE VECTOR, AND
C            F IS THE DEPENDENT VARIABLE VECTOR. THE TRANSFORMATION
C            IS APPLIED TO ONE ROW OF THE JACOBIAN MATRIX AT A TIME.
C  Parameters:                                                   
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

