c***********************************************************************
c***********  Program  betaFIT_2.2  dated  04 April 2010  **************
c***********************************************************************
c* Program to fit NTP read-in potential fx. values {RTP(i),VTP(i)} 
c  (with or without individual weights) to a chosen analytic form.
c***********************************************************************
      INTEGER MXDATA, MXPARM, MXMLR
      PARAMETER (MXDATA=1501, MXPARM=43, MXMLR= 15)
      INTEGER I,J,INFL,ITER,IROUND,JROUND,ROBUST,prFIT,prNLL,prDIFF,M,
     1  NPARM,NTP,NLIN,NDGF,IFXP(MXPARM)
      REAL*8 BETA(0:MXPARM),PV(MXPARM),PU(MXPARM),PS(MXPARM),
     1 CM(MXPARM,MXPARM),DYDP(MXDATA,MXPARM),VTP(MXDATA),
     2 uVTP(MXDATA),betay(MXDATA),Ubetay(MXDATA),YD(MXDATA),
     3 ypSAP(MXPARM),xSAP(MXDATA), rKL(1:MXDATA,1:MXDATA),
     4 DM(MXMLR),DMP(MXMLR),DMPP(MXMLR),
     3 betaINF,UNC,yPOW,DSE,TSTPS,TSTPU,DSEB, TCM,UM,
     4 Rep,AREF,AREFp,AREFq,RTPp,RTPq, AA,BB,ULR,dULR,FCT,RAT,UMAX,
     5 XX,YY,YH,yp,yq,ypRE,ReDE, ReIN,DeIN,VMINin,ULRe,dULRe,RE3,RE6,
     6 RE8,T0,T1,C6adj,C9adj,RTP3,RTP6,RTP8,RH,RR,RB,RBB,VV,VB,VBB,
     7 SCALC,Rsap,DRMSD ,BETA0,BETAN
      REAL*8  DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),DEIGRe(1,1)
     1  ,DEIGDe(1,1)
      CHARACTER*4  NNAME,NAME(4)
      DATA NAME/' EMO',' MLR','DELR','GPEF'/
c-----------------------------------------------------------------------
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,IDF,IDSTT,KDER,NCMM,MMLR(MXMLR),
     1                                              p,q,NS,NL,NPHI,SAP
      REAL*8 Re,De,VMIN,RREF,M2,as,bs,rhoAB,
     1  CmVAL(MXMLR),RTP(MXDATA),SAS(MXDATA,MXPARM)
      COMMON /DATABLK/Re,De,VMIN,RREF,M2,as,bs,rhoAB,CmVAL,RTP,SAS,PSEL,
     1 IFXRe,IFXDe,IFXVMIN,IDF,IDSTT,KDER,NCMM,MMLR,p,q,NS,NL,NPHI,SAP
c-----------------------------------------------------------------------
      ROBUST= 0
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
      JROUND= 0
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
c** If rhoAB .LE. 0.0  have NO damping functions: all  Dm(R)= 1.0
c   If rhoAB > 0.0  it is the molecule-dependent radial scaling factor
c                    of Douketis et al. [JCP 76, 3057 (1982)]
c     rhoAB =  2*rhoA*rhoB/(rhoA+rhoB)   where  rhoA  is the ionization
c           potential ratio  (I_p^A/I_p^H)^0.66  for atom A vs. hydrogen
c  For rhoAB > 0.0,  IDF specifies damping s.th.  Dm(r)/r^m --> r^{IDF/2}
c                    IDSTT > 0  use Douketis et al. damping functions
c                    IDSTT > 1  use DS damping for s=-1/2 AND impose very
c                       short-range 1/r behaviour with  C1= IDSTT = Z1*Z2
c                    IDSTT .LE. 0  use Tang-Toennies damping functions
c   {NOTE:  IDF & IDSTT required only for testing - remove em later?}
c
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
          READ(5,*) NCMM, rhoAB, IDSTT, IDF, SAP, Rsap
          DO m= 1, NCMM
              READ(5,*) MMLR(m), CmVAL(m)
              ENDDO
c=======================================================================
          IF((NCMM.GE.4).AND.(MMLR(2).EQ.0)) THEN
c** For Lyon treatment of A-state alkali dimers ...
              WRITE(6,602) CmVAL(2),CmVAL(1),(CmVAL(i),i=3,NCMM)
            ELSE
c** For 'normal' inverse-power sum uL, with/without damping
              IF(rhoAB.LE.0) THEN
                  WRITE(6,698) (MMLR(m),CmVAL(m),m= 1,NCMM)
                ELSE
                  IF(IDSTT.GT.0) WRITE(6,696) rhoAB,IDF,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
                  IF(IDSTT.LE.0) WRITE(6,694) rhoAB,IDF,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
                ENDIF
            ENDIF
          ENDIF
      IF(PSEL.EQ.3) THEN
c* For a DELR potential, NCMM is the number of long-range terms in the
c  sum  u_{LR}(r)= \sum_i{D_i(r) CmVAL(i)/r^i} , where i= MMLR(i) and
c* rhoAB   is the system-dependent scaling factor in the damping function.
c          and it and  IDF  and  IDSTT  are defined as above.
c  Attractive terms have positive  CmVAL  values.
c=======================================================================
          READ(5,*) NCMM, rhoAB, IDF, IDSTT
          READ(5,*) (MMLR(m), CmVAL(m), m= 1,NCMM)
c=======================================================================
c* IDF  specifies the version of the incomplete-gamma damping function 
c       to be used where limiting behaviour of  Dm(r)/r^m  is  r^{IDF}
c    Our Li2(B) study used  IDF=+1, but  IDF=0  would also be accptable
          IDF= +1
c-----------------------------------------------------------------------
          BETA(0)= 1.d0
          IF(IDSTT.GT.0) WRITE(6,696) rhoAB,IDF,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
          IF(IDSTT.LE.0) WRITE(6,694) rhoAB,IDF,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
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
  690 FORMAT(' uLR(r)  is a damped TT inverse-power sum with',6x,
     1'   rhoAB=',  F7.4,'   IDF=',i3/(49x,'C',i2,'=',1PD15.8:))
  694 FORMAT(' uLR inverse-power terms incorporate TT damping with   rho
     1AB=',f10.7/8x,'defined to give very short-range damped LR-term beh
     2aviour   r^{',i2,'/2}'/(49x,'C',i2,'=',1PD15.8:))
  696 FORMAT(' uLR inverse-power terms incorporate DS damping with   rho
     1AB=',f10.7/8x,'defined to give very short-range damped LR-term beh
     2aviour  r^{',i2,'/2}'/(49x,'C',i2,'=',1PD15.8:)) 
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
 6022 FORMAT(' Use Lyon 3x3  uLR(r)  with   Aso=',F10.6,'   C_3(Sigma)='
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
     1ue',f12.4)
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
   10 READ(5,*,END= 999) p, q, NS, NL, RREF
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
              IF((NCMM.GE.4).AND.(MMLR(2).LE.0)) THEN
c... For Lyon treatment of A-state alkali dimers ...
                  IF(MMLR(2).EQ.0)
     1              WRITE(6,602) CmVAL(2),CmVAL(1),(CmVAL(m),m=3,NCMM)
                  IF(MMLR(2).EQ.-1)
     1              WRITE(6,6022) CmVAL(2),CmVAL(1),(CmVAL(m),m=3,NCMM)
                ELSE
c... For 'conventional' (damped or non-damped) inverse-power sum uLR(r)
                  IF(rhoAB.LE.0) THEN
                      WRITE(6,698) (MMLR(m),CmVAL(m),m= 1,NCMM)
                    ELSE
                      IF(IDSTT.GT.0) WRITE(6,696) rhoAB,IDF,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
                      IF(IDSTT.LE.0) WRITE(6,694) rhoAB,IDF,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
                    ENDIF
                ENDIF
              WRITE(6,634) p,p,q
              ENDIF
          ENDIF
      IF(PSEL.EQ.3) THEN
          IF(IDSTT.GT.0) WRITE(6,696) rhoAB,IDF,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
          IF(IDSTT.LE.0) WRITE(6,694) rhoAB,IDF,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
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
          IF((NCMM.GE.4).AND.(MMLR(2).LE.0)) THEN
c** Extension for Aubert-Frecon Li2(A) {3,0,6,6,8,8} case ...
              RE3= 1.d0/Re**3
              RE6= RE3*RE3
              RE8= RE6/Re**2
              C6adj= CmVAL(3) + 0.25D0*CmVAL(1)**2/De
              C9adj= 0.5d0*CmVAL(1)*C6adj/De
              IF(MMLR(2).EQ.0) THEN
c ... for Aubert-Frecon 2x2 case ...
                  T1= (0.5d0*CmVAL(1)+ (C6adj-CmVAL(4))*RE3)*RE3/3.d0
                  IF(NCMM.GT.4) THEN
                      T1= T1+ (CmVAL(5)- CmVAL(6))*RE8/3.d0
                      ENDIF
                  T0= DSQRT((T1- CmVAL(2))**2 + 8.d0*T1**2)
                  ULRe= 0.5d0*( -CmVAL(2)+ RE3*(1.5d0*CmVAL(1) 
     1           + RE3*(C6adj + CmVAL(4)))) + 0.5d0*T0 + C9adj*RE6*RE3
                  IF(NCMM.GT.4) THEN
                      ULRe= ULRe+ 0.5d0*(CmVAL(5)+ CmVAL(6))*RE8
                      ENDIF
                  ENDIF
              IF(MMLR(2).EQ.-1) THEN
c ... for Aubert-Frecon 3x3 Li2(c) {3,0,6,6,8,8} case ...
                  CALL AF3x3potret(Re,CmVAL(2),CmVAL(1),C6adj,CmVAL(5),
     1                      De,ULRe,DEIGM1,DEIGM3,DEIGM5,DEIGR,DEIGDe)

ccc             c9adj= 0.0

                  ULRe= ULRe + C9adj*RE6*RE3
                  ENDIF
            ELSE
c** For normal inverse-power sum MLR/MLJ case, with or without damping
              KDER=0
              IF(rhoAB.GT.0.d0) CALL dampF(Re,rhoAB,NCMM,MMLR,IDF,
     1                                         IDSTT,KDER,DM,DMP,DMPP)
              ULRe= 0.d0
              DO  m= 1,NCMM
                  T0= CmVAL(m)/Re**MMLR(m)
                  IF(rhoAB.GT.0.d0) T0= T0*DM(m)
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
              IF((NCMM.GE.4).AND.(MMLR(2).LE.0)) THEN
c ... extension for Aubert-Frecon {3,0,6,6,8,8} cases ...
                  RTP3= 1.d0/RTP(i)**3
                  RTP6= RTP3*RTP3
                  RTP8= RTP6/RTP(i)**2
                  IF(MMLR(2).EQ.0) THEN
c ... for Aubert-Frecon Li2(A) 2x2 {3,0,6,6,8,8} case ...
                      T1= (0.5d0*CmVAL(1)+ (C6adj- CmVAL(4))*RTP3)
     1                                                      *RTP3/3.d0
                      IF(NCMM.GT.4) THEN
                          T1= T1+ (CmVAL(5)- CmVAL(6))*RTP8/3.d0
                          ENDIF
                      T0= DSQRT((T1- CmVAL(2))**2 + 8.d0*T1**2)
                      ULR= 0.5d0*( -CmVAL(2) + RTP3*(1.5d0*CmVAL(1)
     1                       + RTP3*(CmVAL(3) + CmVAL(4)))) + 0.5d0*T0
     2                                               + C9adj*RTP3*RTP6
                      IF(NCMM.GT.4) 
     1                        ULR= ULR+ 0.5d0*(CmVAL(5)+CmVAL(6))*RTP8
                      ENDIF
                  IF(MMLR(2).EQ.-1) THEN
c ... for Aubert-Frecon Li2(A) 3x3 {3,0,6,6,8,8} case ...
                      CALL AF3x3potret(RTP(i),CmVAL(2),CmVAL(1),C6adj,
     1              CmVAL(5),De,ULR,DEIGM1,DEIGM3,DEIGM5,DEIGR,DEIGDe)
                      ULR= ULR + C9adj*RTP6*RTP3
                      ENDIF
                ELSE
c... for normal MLR/MLJ case ... with or without damping
                  IF(rhoAB.GT.0.d0) CALL dampF(RTP(i),rhoAB,NCMM,MMLR,
     1                                     IDF,IDSTT,KDER,DM,DMP,DMPP)
                  ULR= 0.d0
                  DO  m= 1,NCMM
                      T0= CmVAL(m)/RTP(i)**MMLR(m)
                      IF(rhoAB.GT.0.d0) T0= T0*DM(m)
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
c!!!!!!!!!
c!!!  RJL  DELR damping needs fixing up !!  [This assumes  DS L_{sr}=0]
c!!!!!!!!!
          KDER= 1
          CALL dampF(Re,rhoAB,NCMM,MMLR,IDF,IDSTT,KDER,DM,DMP,DMPP)
          DO  m= 1,NCMM
              FCT= dexp(-3.95d0*rhoAB*Re/MMLR(m)
     1               - 0.39*(rhoAB*Re)**2/DSQRT(DFLOAT(MMLR(m))))
              AA= CmVAL(m)/Re**MMLR(m)
              ULRe= ULRe+ AA*DM(MMLR(m))
              dULRe= dULRe+ AA*(DMP(m) - MMLR(m)*DM(m)/Re)
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
              CALL dampF(RTP(i),rhoAB,NCMM,MMLR,IDF,IDSTT,KDER,DM,DMP,
     1                                                           DMPP)
              DO  m= 1,NCMM
                  ULR= ULR + DM(m)*CmVAL(m)/RTP(i)**MMLR(m)
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
          IF(IDSTT.GT.1) IFXP(NPHI)= IDSTT
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
          CALL NLLSSRR(NTP,NPARM,MXPARM,JROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
          IF(prFIT.GT.0) THEN
              DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
              IF(SAP.LE.0) THEN
                  IF(RREF.GT.0.d0) THEN
                      WRITE(6,622) NNAME,p,q,AREF,NS,NL,DRMSD,0,PV(1),
     1                PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPHI)
cc                    WRITE(7,722) NNAME,p,q,AREF,NS,NL,DRMSD,
cc   1                                                (PV(j),j=1,NPHI)
                      ENDIF
                  IF(RREF.LE.0.d0) THEN
                      WRITE(6,624) NNAME,p,q,NS,NL,DRMSD,0,PV(1),PU(1),
     1                      PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPHI)
cc                    WRITE(6,724) NNAME,p,q,NS,NL,DRMSD,
cc   1                                                (PV(j),j=1,NPHI)
                      ENDIF
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
              CALL NLLSSRR(NTP,NPARM,MXPARM,JROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              Re= PV(NPHI+1)
              DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
              IF(prFIT.GT.0) THEN
                  IF(SAP.LE.0) THEN
                      IF(RREF.GT.0.d0) THEN
                          WRITE(6,622) NNAME,p,q,AREF,NS,NL,DRMSD,0,
     1          PV(1),PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPHI)
cc                        WRITE(7,722) NNAME,p,q,NS,NL,DRMSD,
cc   1                                                (PV(j),j=1,NPHI)
                          ENDIF
                      IF(RREF.LE.0.d0) THEN
                          WRITE(6,624) NNAME,p,q,NS,NL,DRMSD,0,PV(1),
     1                PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPHI)
cc                        WRITE(7,724) NNAME,p,q,NS,NL,DRMSD,
cc   1                                                (PV(j),j=1,NPHI)
                          ENDIF
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
              CALL NLLSSRR(NTP,NPARM,MXPARM,JROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              IF(IFXDe.LE.0) De= PV(NPHI+2)
              IF(IFXVMIN.LE.0) VMIN= PV(NPHI+3)
              IF((prFIT.GT.0).OR.(DSE.GT.DSEB*1.01)) THEN
                  IF(DSE.GT.DSEB*1.01) WRITE(6,654) DSEB,DSE
 
                  DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
                  IF(SAP.LE.0) THEN
                      IF(RREF.GT.0.d0) THEN
                          WRITE(6,622) NNAME,p,q,AREF,NS,NL,DRMSD,0,
     1          PV(1),PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPHI)
cc                        WRITE(7,722) NNAME,p,q,AREF,NS,NL,DRMSD,
cc   1                                                (PV(j),j=1,NPHI)
                          ENDIF
                      IF(RREF.LE.0.d0) THEN
                          WRITE(6,624) NNAME,p,q,NS,NL,DRMSD,0,PV(1),
     1                PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPHI)
cc                        WRITE(7,724) NNAME,p,q,NS,NL,DRMSD,
cc   1                                                (PV(j),j=1,NPHI)
                          ENDIF
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
              IF(RREF.GT.0.d0) THEN
                  WRITE(6,622) NNAME,p,q,AREF,NS,NL,DRMSD,0,
     1          PV(1),PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPHI)
                  WRITE(7,722) NNAME,p,q,AREF,NS,NL,DRMSD,
     1                                                (PV(j),j=1,NPHI)
                  WRITE(7,726)  De,Re
                  ENDIF
              IF(RREF.LE.0.d0) THEN
                  WRITE(6,624) NNAME,p,q,NS,NL,DRMSD,0,PV(1),
     1                PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPHI)
                  WRITE(7,724) NNAME,p,q,NS,NL,DRMSD,(PV(j),j=1,NPHI)
                  ENDIF
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
                  CALL NLLSSRR(NTP,NPARM,MXPARM,JROUND,ROBUST,prNLL,
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
c======================================================================
c*** Step inward and check whether the repulsive inner potential wall
c   has an inflection point or turnover.
      IF(SAP.LE.0) THEN
          IF(NTP.GE.MXDATA) THEN
              WRITE(6,674)  NTP,MXDATA
              GOTO 90
              ENDIF
          RH= RTP(1)*1.0d-2
          J= NTP+ 1
          RR= RTP(1)
          RTP(J)= RR
          CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VV,PV,PU,PS,RR)
          RB= RR + RH
          RTP(J)= RB
          CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VB,PV,PU,PS,RR)
          INFL= 0
          DO  I= 1,99
              RBB= RB
              VBB= VB
              RB= RR
              VB= VV
              RR= RR - RH
              RTP(J)= RR
              IF(RTP(J).LT.0.1d0) EXIT
              CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VV,PV,PU,PS,RR)
c
c     write(6,666) rr,vv, (vv-vb)/RH
c 666 format(f8.4,3f16.4)
c
              IF(INFL.GT.0) THEN
                  IF(VV.LT.VB) THEN
c... warning of potential inner-wall turnover
                      WRITE(6,672) RB,VB
                      GOTO 90
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
   90 WRITE(6,608)
      GOTO 10
c-----------------------------------------------------------------------
  616 FORMAT(//' *** WARNING:   p=',i2,' .LE. [MMLR(',i1,')-MMLR(1)]=',
     1  i2,'  CAUTION !!!!!! ******')
  619 FORMAT(' Linearized fit uses    beta(INF)=',f10.6)
  620 FORMAT(/' Linearized ',A4,'{p=',i1,'; NS=',i2,', NL=',i2,'} fit yi
     1elds   DSE=',1Pd9.2/(3x,a4,'_{',i2,'}=',d19.11,' (+/-',d8.1,')   P
     2S=',d8.1))
  621 FORMAT(/' Linearized ',A4,'{p=',i1,', q=',i2,'; Rref=',f5.2,
     1 '; NS=',i2,', NL=',i2,'} fit yields   DSE=',1Pd9.2/(3x,a4,'_{',
     2 i2,'}=',d19.11,' (+/-',d8.1,')   PS=',d8.1))
  623 FORMAT(/' Linearized ',A4,'{p=',i1,', q=',i2,';  Rref=Re;  NS=',
     1 i2,', NL=',i2,'} fit yields   DSE=',1Pd9.2/(3x,a4,'_{',i2,'}=',
     2 d19.11,' (+/-',d8.1,')   PS=',d8.1))
  622 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',i2,'; Rref=',f5.2,
     1 ' ; NS=',i2,', NL=',I2,'}  potl:  dd=',1Pd12.5/
     2 '   beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,3x,'DSE=',
     3   D12.5/('   beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  722 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',i2,'; Rref=',f5.2,
     1 ' ; NS=',i2,', NL=',I2,'}  potl:  dd=',1Pd12.5/(d20.12,'  0'))
  624 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',I2,
     1 '; Rref= Re ; NS=',i2,', NL=',I2,'}  potl:     dd=',1Pd9.2/
     2 '   beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,4x'DSE=',
     3  D9.2/('   beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  724 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',I2,
     1 '; Rref= Re ; NS=',i2,', NL=',I2,'}  potl:     dd=',1Pd9.2/
     2  (d20.12,'  0'))
  726 FORMAT('   De=',f14.6 /'   Re=',f11.7)
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
  674 FORMAT(/' *** {NTP=',I4,'}.GE.{MXDATA=',i4,'}: array size prevents
     1 inner-wall curvature check')
  999 STOP
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

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

c***********************************************************************
      SUBROUTINE dampF(r,rhoAB,NCMM,MMLR,IDF,IDSTT,KDER,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. R of the chosen version of the incomplete
c    gamma function damping function, for  m= 1 to MMAX.
c---------------------- RJL Version of 06 July 2010 --------------------
c-----------------------------------------------------------------------
c                 Upon Input
c* r - the radial distance in Angsroms (!) 
c* RHOab  'universal' scaling coefficient used for systems other than H_2
c       RHOab= 2*(RHOa*RHOb)/(RHOa+RHOb) where RHOa = (I_p^A/I_p^H)^0.66
c              where I_p^A is the ionization potential of atom A
c              and I_p^H is the ionization potential of atomic hydrogen
c* NCMM  the number of inverse-power terms to be considered
c* MMLR  are the powers of the NCMM inverse-power terms
c* IDF requires damping to be defined s.th.  Dm(r)/r^m --> r^{IDF/2}
c* IDSTT specifies damping function type:  > 0  use Douketis et al. form 
c                               if  IDSTT .LE. 0  use Tang-Toennies form
c* KDER:  if KDER.GT.0  the first derivative is also calculated 
c*        if KDER.GT.1  the second derivative is also calculated 
c-----------------------------------------------------------------------
c                 Upon Output
c  DM(m) - The value of the damping function for the long range term 
c          C_MMLR(m)/r^MMLR(m)    {m= 1, NCMM}
c  DMP(m) - The first derivative of the damping function  DM(m)
c  DMPP(m) - The second derivative of the damping function  DM(m)
c-----------------------------------------------------------------------
      INTEGER NCMM,NCMMax,MMLR(NCMM),IDF,IDSTT,KDER,IDFF,FIRST,
     1  Lsr,m,MM,MMAX
      REAL*8 r,rhoAB,bTT(-2:2),cDS(-4:0),bDS(-4:0),aTT,br,XP,YP,
     1  TK, DM(NCMM),DMP(NCMM),DMPP(NCMM),SM(-3:25),
     2  bpm(20,-2:0), cpm(20,-2:0),ZK
c------------------------------------------------------------------------
c  The following values for the numerical factors used in both TT and DS
c  were  normalized to the Hydrogen data presented
c  by Kreek and Meath in J.Chem.Phys. 50, 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c------------------------------------------------------------------------
       DATA bTT/2.10d0,2.44d0,2.78d0,3.126d0,3.471d0/
       DATA bDS/2.50d0,2.90d0,3.3d0,3.69d0,3.95d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.40d0,0.39d0/
       DATA FIRST/ 1/
       SAVE FIRST, bpm, cpm
c------------------------------------------------------------------------
      IF(RHOab.LE.0) THEN
          WRITE(6,602) RHOab
          STOP
          ENDIF
      IF(IDSTT.LE.0) THEN
c===========================================
c** For Tang-Toennies type damping functions
c===========================================
          IF((IDF.LT.-4).OR.(IDF.GT.4)) THEN
                WRITE(6,600) IDSTT,IDF
                STOP
                ENDIF
          Lsr= IDF/2
          MMAX= MMLR(NCMM) + Lsr - 1
          aTT= RHOab*bTT(Lsr)
          br= aTT*r
          XP= DEXP(-br)
          SM(-3)= 0.d0
          SM(-2)= 0.d0
          SM(-1)= 0.d0
          SM(0)=  1.d0
          TK= 1.d0
          IF(br.GT.0.5d0) THEN
              DO  m= 1,MMAX
                  TK= TK*br/DFLOAT(m)
                  SM(m)= SM(m-1)+ TK
                  ENDDO
              DO m= 1, NCMM
                  MM= MMLR(m) - 1 + Lsr
                  DM(m)= 1.d0 - XP*SM(MM)
                  IF(KDER.GT.0) THEN
                      DMP(m)= aTT*XP*(SM(MM) - SM(MM-1))
                      IF(KDER.GT.1) DMPP(m)= -aTT*aTT*XP*(SM(MM) 
     1                                     - 2.d0*SM(MM-1) + SM(MM-2))
                      ENDIF
                  ENDDO
c-----------------------------------------------------------------------
c  The above section handles the calculation of the value of the damping
c  function for most values of r.  However, at very small r that algorithm
c  becomes unstable due to numerical noise.  To avoid this, if the 
c  argument is very small it is re-evaluated as a finite sum ...
c-----------------------------------------------------------------------
            ELSE
              MMAX= MMAX+5
              DO  m= 1, MMAX
c... NOTE that here SM(m) is the m'th term  (b*r)^m/m!  [not a sum]
                  SM(m)= SM(m-1)*br/DFLOAT(m)
                  ENDDO
              DO  m= 1, NCMM
                  MM= MMLR(m) + Lsr
                  DM(m)= XP*(SM(MM)+ SM(MM+1)+ SM(MM+2)+ SM(MM+3) 
     1                                                     + SM(MM+4))
                  IF(KDER.GT.0) THEN
                      DMP(m)= aTT*XP*SM(m-1)
                      IF(KDER.GT.1)DMPP(m)= aTT*aTT*XP*(SM(m-2)-SM(m-1))
                      ENDIF
                  ENDDO
            ENDIF
          ENDIF
c
      IF(IDSTT.GT.0) THEN
c=======================================================================
c** For Douketis-Scoles-Marchetti-Zen-Thakkar type damping function ...
c=======================================================================
          IF((IDF.LT.-4).OR.(IDF.GT.0)) THEN
              WRITE(6,600) IDSTT,IDF
              STOP
              ENDIF
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  IDFF= -2,0
                      bpm(m,IDFF)= bDS(IDFF)/DFLOAT(m)
                      cpm(m,IDFF)= cDS(IDFF)/DSQRT(DFLOAT(m))
                      ENDDO
                  ENDDO
              FIRST= 0 
              ENDIF
          br= rhoAB*r
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,IDF) + cpm(MM,IDF)*br)*br)
              YP= 1.d0 - XP
              ZK= MM-1.d0
              DM(m)= YP**(MM-1)
c... Actually ...  DM(m)= YP**(MM + IDF/2)  :  set it up this way to 
c   avoid taking exponential of a logarithm for fractional powers (slow)
              IF(IDF.EQ.-4) THEN
                  ZK= ZK- 1.d0
                  DM(m)= DM(m)/YP
                  ENDIF
              IF(IDF.EQ.-3) THEN
                  ZK= ZK- 0.5d0
                  DM(m)= DM(m)/DSQRT(YP)
                  ENDIF
              IF(IDF.EQ.-1) THEN
                  ZK= ZK+ 0.5d0
                  DM(m)= DM(m)*DSQRT(YP)
                  ENDIF
              IF(IDF.EQ.0) THEN
                  ZK= MM
                  DM(m)= DM(m)*YP
                  ENDIF
              IF(KDER.GT.0) THEN
                  TK= bpm(MM,IDF) + 2.d0*cpm(MM,IDF)*br
                  DMP(m) = ZK*XP*rhoAB*TK*DM(m)/YP
                  IF(KDER.GT.1) THEN
c ... if desired ... calculate second derivative [for DELR case] {check this!}
                      DMPP(m)= (ZK-1.d0)*XP*TK*DMP(m)/YP
     1               - DMP(m)*TK + DMP(m)*2.d0*cpm(MM,IDF)*rhoAB**2/TK
                      ENDIF
                  ENDIF
              ENDDO   
          ENDIF  
      RETURN
  600 FORMAT(/,' *** ERROR ***  For  IDSTT=',i3,'   IDF=',i3,'  no dampi
     1ng function is defined')
  602 FORMAT( /,' ***ERROR ***  rhoAB=', F7.4,'  yields an invalid Dampi
     1ng Function definition')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
