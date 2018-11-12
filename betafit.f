c***********************************************************************
c***********  Program  betaFIT_2.0  dated  22 July 2008  ***************
c***********************************************************************
c* Program to fit NTP read-in potential fx. values {RTP(i),VTP(i)} to
c  a chosen analytic form.
c***********************************************************************
      INTEGER MXDATA, MXPARM, MXMLR
      PARAMETER (MXDATA=1501, MXPARM=40, MXMLR= 8)
      INTEGER I,J,INFL,ITER,IROUND,ROBUST,LPRINT,IWR,M,NPARM,NTP,
     1  NLIN,IFXP(MXPARM)
      REAL*8 BETA(0:MXPARM),PV(MXPARM),PU(MXPARM),PS(MXPARM),
     1 CM(MXPARM,MXPARM),DYDP(MXDATA,MXPARM),VTP(MXDATA),
     2 uVTP(MXDATA),betay(MXDATA),Ubetay(MXDATA),YD(MXDATA),
     3 ypSAP(MXPARM),xSAP(MXDATA), rKL(1:MXDATA,1:MXDATA),
     3 betaINF,UNC,yPOW,DSE,TSTPS,TSTPU,DSEB,TT(0:20),RHOdR,RHOp,TTM,
     4 Rep,AREF,AREFp,AREFq,RTPp,RTPq, AA,BB,ULR,dULR,FCT,RAT,UMAX,
     5 XX,YY,YH,yp,yq,fsw,ypRE,ReDE, ReIN,DeIN,VMINin ,ULRe,RE3,RE6,
     6 RE8,T0,T1,C6adj,C9adj,RTP3,RTP6,RTP8,RH,RR,RB,RBB,VV,VB,VBB,
     7 SCALC
      CHARACTER*4  NNAME,NAME(5)
      DATA NAME/' EMO',' MLJ',' MLR','DELR','GPEF'/
c-----------------------------------------------------------------------
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR(MXMLR),p,q,NS,NL,
     1                                                        NPHI,SAP
      REAL*8 Re,De,VMIN,RREF,Asw,Rsw,M2,as,bs,RHOd,
     1  CmVAL(MXMLR),RTP(MXDATA),SAS(MXDATA,MXPARM)
      COMMON /DATABLK/Re,De,VMIN,RREF,Asw,Rsw,M2,as,bs,
     1  RHOd,CmVAL,RTP,SAS,PSEL,IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR,p,q,
     2                                                  NS,NL,NPHI,SAP
c-----------------------------------------------------------------------
      ROBUST= 0
c-----------------------------------------------------------------------
c** PSEL  specifies the type of potential being fitted to:
c     PSEL=1 for EMO;   PSEL=2 for an MLR (or MLJ);   PSEL=3  for DELR ;
c     PSEL= 4  for GPEF
c* NPT  is the number of read-in potential points being fitted to.
c* UNC  is the energy uncertainty associated with the potential points
c     (plausibly ca. 0.1 cm-1 for RKR).
c  IROUND  specifies the level of rounding inside NLLSSRR if:
c          > 0 : requires that Sequential Rounding & Refitting be
c                performed, with each parameter being rounded at the
c                IROUND'th sig. digit of its local uncertainty.
c          <=0 : simply stops after full convergence (without rounding).
c* IWR > 0  causes printout of results of preliminary linearized fit and
c       (as appropriate) other non-final fits.  Notmally set  IWR= 0.
c  LPRINT  specifies the level of printing inside NLLSSRR
c        if: =  0, no print except for failed convergence [normal value]
c             < 0  only converged, unrounded parameters, PU & PS's
c            >= 1  print converged parameters, PU & PS's
c            >= 2  also print parameter change each rounding step
c            >= 3  also indicate nature of convergence
c            >= 4  also print convergence tests on each cycle
c            >= 5  also parameters changes & uncertainties, each cycle
c
c** Re & De  are assumed potential well minimum position & well depth.
c  ...  note that  De  is a dummy variable for a GPEF function
c** VMIN is the minimum of the potential function defined by the read-in 
c   points [typically =0 for RKR potential, but non-zero for ab initio]
c** To fix  Re, De or VMIN  unchanged at read-in values, set (integer)
c   IFXRe, IFXDe, and/or IFXVMIN > 0 ;  to fit them, set value =0
c   [Normally set  IFXVMIN= 0 !!]
c=======================================================================
      READ(5,*) PSEL, NTP, UNC, IROUND, IWR, LPRINT
      READ(5,*) Re, De, VMIN
      READ(5,*) IFXRe, IFXDe, IFXVMIN
c=======================================================================
      ReIN= Re
      DeIN= De
      VMINin= VMIN
c** For an MLR_p potential (PSEL=2) read number of long-range terms NCMM
c   to define long-range potential tail:  V(r)= De - \sum{CmVAL/r^MMLR}
c** To use the switching function  fsw(r)= 1/[exp{Asw*(r-Rsw)} + 1]  form
c   of exponent coefft  beta(y)= [1-fsw] beta_inf + fsw \sum {beta_i y^i} 
c   read positive values of Asw and Rsw; otherwise, set them .LE. 0 and
c        beta(r)= yp beta_inf + [1 - yp] Sum{ beta_i yp^i }
c** SAP is an integer: > 0  to use Pashov natural spline in MLR exponent
c                    : .le. 0  to use constrained polynomial exponent
c*  For an MLR potential with  SAP > 0,  Rsw  is a negative{!} number
c      (0 .le. Rsw < 1) which specifies the lower bound on the yp values 
c      selected to define the MLR exponent spline.  
c** For each long-range term read power  MMLR(i)  & coefficient CmVAL(i)
c** For special Aubert-Frecon case,  NCMM= 4,  MMLR= {3,0,6,6} or 
c  {3,0,6,6,8,8} and coefficients are: CmVAL(1)= C3(sig), CmVAL(2)= ASO,
c  CmVAL(3)= C6(Sig) CmVAL(4)= C6(pi) // oa, for NCMM= 6  with added 
c  coefficients CmVAL(5)= C8(sig) & CmVAL(6)= C8(pi)
c=======================================================================
      IF(PSEL.EQ.2) THEN
          READ(5,*) NCMM, Asw, Rsw,SAP
          READ(5,*) (MMLR(i), CmVAL(i), i= 1,NCMM)
          ENDIF
c=======================================================================
      IF(PSEL.EQ.3) THEN
c* For a DELR potential, NCMM is the number of long-range terms in the
c  sum  u_{LR}(r)= \sum_i{D_i(r) CmVAL(i)/r^i} , where i= MMLR(i)
c  and  RHOd  is the scaling factor \rho_d in the damping function.
c  Attractive terms have positive  CmVAL  values.
c* If  IDF= 1  use the Tang-Toennies damping function
c   D_m(r)= 1 - exp{-3.16*RHOd*r} \Sum_{k=0}^{NCMM} (3.16*RHOd*r)**{k}/k!
c* If  IDF= 2  use the Scoles damping function
c      D_m(r)= [1 - exp{-3.97(RHOd*r)/m - 0.39(RHOd*r)^2/sqrt(m)]^m
c Also ...  BETA(0) is the initial trial value of  \beta_0  used to
c    generate initial trial values of the A & B potential parameters.
c  * If the read-in  BETA(0) \leq 0.0 , the program uses a preliminary 
c    EMO_{p} fit to generate an estimate of  BETA(0); this usually works
c    if  u_{LR} is attractive, so that the potential has no barrier.
c  * If the potential has a barrier, one must determine an initial trial
c    value another way.  A one way would be to do an EMO_{p} fit to the
c    potential points, treating the barrier maximum as dissociation
c=======================================================================
          READ(5,*) NCMM, RHOd, IDF, BETA(0)
          READ(5,*) (MMLR(i), CmVAL(i), i= 1, NCMM)
c=======================================================================
          BETA(0)= 1.d0
          ENDIF
c** For a GPEF potential, read coefficients to define expansion vble:
c       y = (r^p - Re^p)/(as*r^p + bs*Re^p)  where  p, as & bs all fixed
c=======================================================================
      IF(PSEL.EQ.4) READ(5,*) as, bs 
c=======================================================================
c** Read the turning points to be fitted to
c=======================================================================
      READ(5,*) (RTP(i), VTP(i),i= 1,NTP)
c=======================================================================
      IF(PSEL.EQ.1) WRITE(6,600) Re, De
      IF(PSEL.EQ.2) THEN
          NNAME= NAME(2)
          IF(NCMM.GT.1) NNAME= NAME(3)
          WRITE(6,602) NNAME, Re, De, (MMLR(i),CmVAL(i),i= 1,NCMM)
          IF((NCMM.GE.4).AND.(MMLR(2).EQ.0)) THEN
c** For Lyon treatment of A-state alkali dimers ...
              WRITE(6,618) CmVAL(2),CmVAL(3),CmVAL(4)
              IF(NCMM.GT.4) WRITE(6,617) CmVAL(5),CmVAL(6)
              ENDIF
          ENDIF
      IF(PSEL.EQ.3) THEN
          RHOp= RHOd
          IF(IDF.EQ.1) RHOd= 3.16d0*RHOd
          WRITE(6,604) Re, De, RHOp, NCMM,(MMLR(j),CmVAL(j),j=1, NCMM)
          ENDIF
      IF(PSEL.EQ.4) WRITE(6,605) as,bs,Re
      WRITE(6,606) NTP,VMIN,UNC,(RTP(i),VTP(i),i= 1,NTP)
      WRITE(6,608)
      IF((PSEL.EQ.2).AND.(NCMM.GE.4).AND.(MMLR(2).EQ.0)) THEN
          DO  I= 1,NTP
              VTP(I)= VTP(I) - 0.5d0*CmVAL(2)
              ENDDO
          ENDIF
c
  600 FORMAT(' Determine  EMOp  exponent expansion coefficients'/
     1  1x,24('==')/' Start with   Re=',f11.8,'   De=',f11.4)
  601 FORMAT(' using exponent expansion variable  y_',i1,'(r)= [r^',i1,
     1 ' -',f8.4,'^',i1,']/[r^',i1,' +',f8.4,'^',i1,']' )
  602 FORMAT(' Fit an ',A4,'p  potential function to the input points'/
     1  1x,24('==')/' Start with   Re=',f11.8,'   De=',f11.4,'    C',
     2  i2,'=',1PD15.8:/(49x,'C',i2,'=',D15.8:))
  618 FORMAT(' Use Lyon 2x2 ULR(r) with   C_0=',F10.6,'   C_6(Sigma)=',
     1  1PD15.7/45x,'C_6(Pi)   =',D15.7)
  617 FORMAT(42x,'C_8(Sigma)=',1PD15.7/42x,'C_8(Pi)   =',D15.7)
  603 FORMAT(' using exponent expansion variable  y_',i1,'(r)= [r^',i1,
     1 ' - Re^',i1,']/[r^',i1,' + Re^',i1,']' )
  652 FORMAT( '   & define beta(y(r)) as a natural spline through points
     1 at the',i4,'  yp values:'/(2x,7F11.7))
  604 FORMAT(' Determine  DELRp potential exponent expansion coefficient
     1s'/1x,29('==')/' Start with   Re=',f11.8,'   De=',f11.4,
     2  '   where   RHOd=',f10.6/5x,'and  u_{LR}  has',i2,
     3 ' inverse-power terms with coefficients (+ve attractive):'/
     4  (1x,3('   C_{',i2,'}=',1Pd15.7:)))
  605 FORMAT(' Determine coefficient for a  GPEF{p}  polynomial potentia
     1l using the'/1x,29('==')/' expansion variable:   y= (R^p - Re^p)/(
     2',1Pd11.3,'*R^p',1x,SP,d11.3,'*Re^p)'/' with initial   Re=',
     3 0pf12.8 )
  606 FORMAT(/' Fit to',I5,' input turning points with initial energy mi
     1nimum   VMIN=',f11.4/'    assuming uncertainties of  u(VTP)=',
     2  1PD9.2,'  for all input potential values.'/1x,39('--')/
     3  4('    RTP',7x,'VTP   ')/1x,39('--')/(4(0pF9.5,f11.3)))
  607 FORMAT('      which yields initial values of   AA=',1PD14.7,
     1   '   BB=',D14.7)
  608 FORMAT(1x,39('--')/)
  610 FORMAT(' Shift read-in VMIN value to lowest input potential value'
     1   ,f11.4)
  612 FORMAT('  VLR component of potential uses Scoles-type damping func
     1tion')
  614 FORMAT('  VLR component of potential uses Tang-Toennies damping fu
     1nction')
c=======================================================================
c** Now ... loop over different {p,q,NS,NL} combinations till end of data
c**  p and q are powers used to define radial variable in the exponent
c       beta(r)= yp*betaINF + sum{beta_i*yq^i} where  
c       ya=(R^a - AREF^a)/(R^a + AREF^a)
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
      IF((p.LE.0).OR.(NS.LE.0).OR.(NL.LE.0)) GOTO 999
      NPHI= MAX(NS,NL)+ 1
      IF(PSEL.EQ.1) WRITE(6,600) Re, De
      IF(PSEL.EQ.2) THEN
          WRITE(6,602) NNAME, Re, De, (MMLR(i),CmVAL(i),i= 1,NCMM)
          IF(SAP.GT.0) WRITE(6,650) NS,NL
          IF(SAP.LE.0) THEN
              IF(Asw.GT.0.d0) WRITE(6,632) Asw, Rsw
              IF(Asw.LE.0.d0) WRITE(6,634) p,p,q
              ENDIF
          ENDIF
      IF(PSEL.EQ.3) THEN
          WRITE(6,604) Re, De, RHOp, NCMM,(MMLR(j),CmVAL(j),j=1, NCMM)
          IF(IDF.EQ.1) WRITE(6,614)
          IF(IDF.EQ.2) WRITE(6,612)
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
              YH= -Rsw/NS
              ypSAP(1)= Rsw
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
          WRITE(6,610) DSE
          VMIN= DSE
          ENDIF
c=======================================================================
c** Preliminary linearized fit for an  EMOp  potential ...
c-----------------------------------------------------------------------
      IF((PSEL.EQ.1).OR.(PSEL.EQ.3).AND.(beta(0).le.0)) THEN
c ... first define ordinate array
          NNAME= NAME(1)
          DO  i= 1,NTP
              RTPp= RTP(i)**p
              RTPq= RTP(i)**q
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              yq= (RTPq - AREFq)/(RTPq + AREFq)
              IF(RTP(i).GT.Re) THEN
                  betay(i)= - DLOG(1.d0 - DSQRT((VTP(i)-VMIN)/De))
                  IF(VTP(i).GT.UNC) THEN
                      Ubetay(i)= 0.5d0*UNC/(DSQRT((VTP(i)-VMIN)*De) - 
     1                                                  (VTP(i)-VMIN))
                    ELSE
                      Ubetay(i)= DSQRT(UNC/De)
                    ENDIF
                ELSE
                  betay(i)= - DLOG(1.d0 + DSQRT((VTP(i)-VMIN)/De))
                  IF(VTP(i).GT.UNC) THEN
                      Ubetay(i)= 0.5d0*UNC/(DSQRT((VTP(i)-VMIN)*De) +
     1                                                  (VTP(i)-VMIN))
                    ELSE
                      Ubetay(i)= DSQRT(UNC/De)
                    ENDIF
                ENDIF
              uVTP(i)= UNC
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
          IF(IWR.GT.0) WRITE(6,620) NNAME,p,NS,NL,DSE,
     1                        ('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPARM)
          ENDIF

c=======================================================================
      IF(PSEL.EQ.2) THEN
c-----------------------------------------------------------------------
c*** Preliminary linearized fit for an  MLRp  potential ...
c** First define array of exponent values with uncertainties defined by 
c  the assumption that all potential values have equal uncertainties UNC
          NNAME= NAME(3)
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
              WRITE(6,618) (CmVAL(i), i=2,4)
              IF(NCMM.GT.4) THEN
                  ULRe= ULRe+ 0.5d0*(CmVAL(5)+ CmVAL(6))*RE8
                  WRITE(6,617) CmVAL(5),CmVAL(6)
                  ENDIF
            ELSE
c** For normal inverse-power sum MLR/MLJ case
              IF(NCMM.EQ.1) NNAME= NAME(2)
              IF(p.LE.(MMLR(NCMM)-MMLR(1))) THEN
                  WRITE(6,616) p, NCMM,MMLR(NCMM)-MMLR(1)
                  ENDIF
              ULRe= 0.d0
              DO  i= 1,NCMM
                  ULRe= ULRe + CmVAL(i)/Re**MMLR(i)
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
                  T1= (0.5d0*CmVAL(1)+ RTP3*(C6adj- CmVAL(4)))*RTP3/3.d0
                  IF(NCMM.GT.4) THEN
c ... extension for Aubert-Frecon Li2(A) {3,0,6,6,8,8} case ...
                      RTP8= RTP6/RTP(i)**2
                      T1= T1+ (CmVAL(5)- CmVAL(6))*RE8/3.d0
                      ENDIF
                  T0= DSQRT((T1- CmVAL(2))**2 + 8.d0*T1**2)
                  ULR= 0.5d0*( -CmVAL(2) + RTP3*(1.5d0*CmVAL(1)
     1                       + RTP3*(CmVAL(3) + CmVAL(4)))) + 0.5d0*T0
     2                                               + C9adj*RTP3*RTP6
                  IF(NCMM.GT.4) ULR= ULR+ 0.5d0*(CmVAL(5)+CmVAL(6))*RTP8
                ELSE
c... for normal MLR/MLJ case ...
                  ULR= 0.d0
                  DO  j= 1, NCMM
                      ULR= ULR+ CmVAL(j)/RTP(i)**MMLR(j)
                      ENDDO
                ENDIF
              IF(RTP(i).GT.Re) THEN
                  betay(i)= - DLOG((1.d0 - DSQRT((VTP(i)-VMIN)/De))
     1                                                      *ULRe/ULR)
                  IF((VTP(i)-VMIN).GT.UNC) THEN
                      Ubetay(i)= 0.5d0*UNC
     1                      /(DSQRT((VTP(i)-VMIN)*De) - (VTP(i)-VMIN))
                    ELSE
                      Ubetay(i)= DSQRT(UNC/De)
                    ENDIF
                ELSE
                  betay(i)= - DLOG((1.d0 + DSQRT((VTP(i)-VMIN)/De))
     1                                                      *ULRe/ULR)
                  IF((VTP(i)-VMIN).GT.UNC) THEN
                      Ubetay(i)= 0.5d0*UNC
     1                      /(DSQRT((VTP(i)-VMIN)*De) + (VTP(i)-VMIN))
                    ELSE
                      Ubetay(i)= DSQRT(UNC/De)
                    ENDIF
                ENDIF
              uVTP(i)= UNC
              IF(SAP.LE.0) THEN
c** Subtract the \beta_\infty term to yield polynomial for fitting
                  IF(Asw.LE.0.d0) THEN
c... For Huang's MLR exponent function
                      betay(i)= betay(i)- betaINF*yp*ypRE
                      yPOW= ypRE*(1.d0- yp)
                    ELSE
c... For Photos' origonal MLJ exponent switching function
                      fsw= 1.d0/(DEXP(Asw*(RTP(i)- Rsw)) + 1.d0)
                      betay(i)= betay(i)- betaINF*ypRE*(1.d0 - fsw)
                      yPOW= ypRE*fsw
                    ENDIF
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
              IF(IWR.GT.0) WRITE(6,620) NNAME,p,NS,NL,DSE,
     1                        ('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPHI)
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
     1              Scalc(XX,m,NTP+1,xSAP,rKL,MXDATA)*betay(m)/xSAP(m)
                       ENDDO
                   ENDDO
              NLIN= (NPHI+1)/2
              IF(IWR.GT.0) WRITE(6,653) NNAME,p,
     1                        ((ypSAP(I),PV(I),I= J,NPHI,NLIN),J=1,NLIN)
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
          NNAME= NAME(4)
          ITER= 0
   40     ULR= 0.d0
          dULR= 0.d0
          IF(IDF.EQ.2) THEN
c... Using Scoles-type damping function ...
              DO  j= 1,NCMM
                  FCT= dexp(-3.97d0*RHOd*Re/MMLR(j)
     1                    - 0.39d0*(RHOd*Re)**2/DSQRT(DFLOAT(MMLR(j))))
                  AA= CmVAL(j)*((1.d0- FCT)/Re)**MMLR(j)
                  ULR= ULR+ AA
                  dULR= dULR- MMLR(j)*AA/Re
                  BB= FCT*(3.97d0*RHOd/MMLR(j) 
     1                      + 0.78d0*Re*RHOd**2/DSQRT(DFLOAT(MMLR(j))))
                  dULR= dULR+  MMLR(j)*BB*AA/(1.d0- FCT)
                  ENDDO
              ENDIF
          IF(IDF.EQ.1) THEN
c... Using Tang-Toennies damping function ...
              TT(0)= 1.d0
              yPOW= 1.d0
              RHOdR= RHOd*Re
              DO  j= 1,MMLR(NCMM)
                  yPOW= yPOW*RHOdR/DFLOAT(J)
                  TT(J)= TT(J-1)+ yPOW
                  ENDDO
              yPOW= DEXP(-RHOdR)
              ULR= 0.d0
              dULR= 0.d0
              DO  j= 1,NCMM
                  TTM= (1.d0- yPOW*TT(MMLR(j)))*CmVAL(j)/Re**MMLR(j)
                  ULR= ULR+ TTM
		      dULR= dULR+ yPOW*RHOd*(TT(MMLR(j)) - TT(MMLR(j)-1))
     1                          *CmVAL(j)/Re**MMLR(j) - MMLR(j)*TTM/Re
                  ENDDO
              ENDIF
          AA= De - ULR - dULR/beta(0)
          BB= 2.d0*(De - ULR) - dULR/beta(0)
          WRITE(6,607) AA,BB
          RAT= 0.5d0*BB/AA
          UMAX= DSQRT(RAT**2 + (UNC + ULR - DE)/AA)
          ReDE= Re- dlog(RAT)/beta(0)
          DO  i= 1,NTP
              ULR= 0.d0
              RTPp= RTP(i)**p
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              IF(IDF.EQ.2) THEN
c... if using Scoles damping function ...
                  DO  j= 1,NCMM
                      FCT= dexp(-3.97d0*RHOd*RTP(i)/MMLR(j)
     1               - 0.39d0*(RHOd*RTP(i))**2/DSQRT(DFLOAT(MMLR(j))))
                      ULR= ULR+ CmVAL(j)*((1.d0-FCT)/RTP(i))**MMLR(j)
                      ENDDO
                  ENDIF
              IF(IDF.EQ.1) THEN
c... if using Tang-Toennies damping function ...
                  RHOdR= RHOd*RTP(i)
                  yPOW= DEXP(-RHOdR)
                  TT(0)= yPOW
                  DO  j= 1,MMLR(NCMM)
                      yPOW= yPOW*RHOdR/DFLOAT(J)
                      TT(J)= TT(J-1)+ yPOW
                      ENDDO
                  DO  j=1,NCMM
                      ULR= ULR+CmVAL(j)*(1.d0-TT(MMLR(j)))/
     1                                                 RTP(i)**MMLR(j)
                      ENDDO
                  ENDIF
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
                      IF((VTP(i)-VMIN).GT.UNC) THEN
                          Ubetay(i)= 0.5d0*UNC/(AA*(RAT- FCT)*FCT)
                        ELSE
                          Ubetay(i)= UNC/(AA*UMAX*(RAT- UMAX))
                        ENDIF
                    ELSE
c ... deweight away points for which \ln argument would be negative
                      betay(i)= 0.d0
                      Ubetay(i)= 9.d9
                    ENDIF
                ELSE
                  betay(i)= - DLOG(RAT + FCT)
                  IF((VTP(i)-VMIN).GT.UNC) THEN
                      Ubetay(i)= 0.5d0*UNC/(AA*(RAT+ FCT)*FCT)
                    ELSE
                      Ubetay(i)= UNC/(AA*UMAX*(RAT+ UMAX))
                    ENDIF
                ENDIF
   44         uVTP(i)= UNC
c ... now create partial derivative array for linearized fit ...
              yPOW= (RTP(i)- Re)
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
          IF(IWR.GT.0) WRITE(6,620) NNAME,p,NS,NL,DSE,
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
          IF(IFXVMIN.LE.0) IFXP(NPARM)= 0
          PV(NPHI+1)= Re
          PV(NPHI+2)= De
          PV(NPHI+3)= VMIN
          CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,LPRINT,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
          IF(IWR.GT.0) THEN
              IF(SAP.LE.0) THEN
                  IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,q,AREF,NS,NL,
     1                            DSE,(j-1,PV(j),PU(j),PS(j),j=1,NPHI)
                  IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,q,NS,NL,DSE,
     1                                (j-1,PV(j),PU(j),PS(j),j=1,NPHI)
                ELSE
                  IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,p,AREF,NS,NL,
     1                              DSE,(j,PV(j),PU(j),PS(j),j=1,NPHI)
                  IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,p,NS,NL,DSE,
     1                                  (j,PV(j),PU(j),PS(j),j=1,NPHI)
                 ENDIF
              IF(IFXVMIN.LE.0)
     1                   WRITE(6,660) PV(NPHI+3),PU(NPHI+3),PS(NPHI+3)
              ENDIF
c ... the, if appropriate, set  Re  free too ...
          IF(IFXRe.LE.0) THEN
              IFXP(NPHI+1)= 0
              CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,LPRINT,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              Re= PV(NPHI+1)
              IF(IWR.GT.0) THEN
                  IF(SAP.LE.0) THEN
                      IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,q,AREF,NS,
     1                         NL,DSE,(j-1,PV(j),PU(j),PS(j),j=1,NPHI)
                      IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,q,NS,NL,
     1                            DSE,(j-1,PV(j),PU(j),PS(j),j=1,NPHI)
                    ELSE
                      IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,q,AREF,NS,
     1                         NL,DSE,(j,PV(j),PU(j),PS(j),j=1,NPHI)
                      IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,q,NS,NL,
     1                              DSE,(j,PV(j),PU(j),PS(j),j=1,NPHI)
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
              IF(IFXDe.LE.0) IFXP(NPHI+2)= 0
              IF(IFXVMIN.LE.0) IFXP(NPHI+3)= 1
              CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,LPRINT,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              IF(IFXDe.LE.0) De= PV(NPHI+2)
              IF(IFXVMIN.LE.0) VMIN= PV(NPHI+3)
              IF((IWR.GT.0).OR.(DSE.GT.DSEB*1.01)) THEN
                  IF(DSE.GT.DSEB*1.01) WRITE(6,654) DSEB,DSE
                  IF(SAP.LE.0) THEN
                      IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,q,AREF,NS,
     1                         NL,DSE,(j-1,PV(j),PU(j),PS(j),j=1,NPHI)
                      IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,q,NS,NL,
     1                            DSE,(j-1,PV(j),PU(j),PS(j),j=1,NPHI)
                    ELSE
                      IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,q,AREF,NS,
     1                           NL,DSE,(j,PV(j),PU(j),PS(j),j=1,NPHI)
                      IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,q,NS,NL,
     1                              DSE,(j,PV(j),PU(j),PS(j),j=1,NPHI)
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
          CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,LPRINT,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
          IF(SAP.LE.0) THEN
              IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,q,AREF,NS,NL,DSE,
     1                                (j-1,PV(j),PU(j),PS(j),j=1,NPHI)
              IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,q,NS,NL,DSE,
     1                                (j-1,PV(j),PU(j),PS(j),j=1,NPHI)
            ELSE
              IF(RREF.GT.0.d0) WRITE(6,626) NNAME,p,AREF,NS,NL,DSE,
     1                       (j,ypSAP(j),j,PV(j),PU(j),PS(j),j=1,NPHI)
              IF(RREF.LE.0.d0) WRITE(6,628) NNAME,p,NS,NL,DSE,
     1                       (j,ypSAP(j),j,PV(j),PU(j),PS(j),j=1,NPHI)
              IF(RREF.GT.0.d0) WRITE(6,629) RREF
            ENDIF
          IF(PSEL.EQ.3) BETA(0)= PV(1)
          WRITE(6,662) PV(NPHI+1),PU(NPHI+1),PS(NPHI+1)
          WRITE(6,630) PV(NPHI+2),PU(NPHI+2),PS(NPHI+2)
          WRITE(6,660) PV(NPHI+3),PU(NPHI+3),PS(NPHI+3)
ccc 
c         WRITE(7,730) (RTP(I),YD(I),I= 1,NTP)
c 730 FORMAT((4(f10.5,f8.3)))
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
          NNAME= NAME(5)
          DO  i= 1,NTP
              uVTP(i)= UNC
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
              IF(IWR.GT.0) THEN
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
                  CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,LPRINT,
     1                   IFXP,VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
                  IF(IWR.GT.0) THEN
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
              CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,LPRINT,IFXP,
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
          RB= RR+ RH
          J= -1
          CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VB,PV,PU,PS,RB)
          INFL= 0
          DO  I= 1,99
              J= J- 1
              RBB= RB
              VBB= VB
              RB= RR
              VB= VV
              RR= RR- RH
              CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VV,PV,PU,PS,RR)
              IF(INFL.GT.0) THEN
                  IF(VV.LT.VB) THEN
c... test for potential inner-wall turnover
                      WRITE(6,672) RB,VB
                      WRITE(6,608)
                      GOTO 10
                      ENDIF
                ELSE
                  IF((VV-VB).LE.(VB-VBB)) THEN
c... test for potential inner-wall inflection 
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
  620 FORMAT(/' Linearized ',A4,'{p=',i1,'} fit with   NS=',i2,'   NL=',
     1 i2,'   yields   DSE=',1Pd9.2/(3x,a4,'_{',i2,'} =',d19.11,
     2 ' (+/-',d8.1,')   PS=',d8.1))
  622 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',i1,'; Rref=',f5.2,
     1  ' ; NS=',i2,', NL=',I2,'}  potl:   DSE=',1Pd9.2/
     2  ('   beta_{',i2,'}=',d19.11,' (+/-',d8.1,')   PS=',d8.1))
  624 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',I1,
     1 '; Rref= Re ; NS=',i2,', NL=',I2,'}  potl:   DSE=',1Pd9.2/
     2 ('   beta_{',i2,'}=',d19.11,' (+/-',d8.1,')   PS=',d8.1))
  626 FORMAT(/' Direct fit to ',A4,'{p=',i1,'; Rref=',f5.2,' ; NS=',i2,
     1  ', NL=',I2,'}  potential:   DSE=',1Pd9.2/(' ypSAP{',i2,'}=',
     20PF11.7,'   beta_{',i2,'}=',1Pd19.11,' (+/-',d8.1,')   PS=',d8.1))
  628 FORMAT(/' Direct fit to ',A4,'{p=',i1,'; Rref= Re ; NS=',i2,
     1  ', NL=',I2,'}  potential:   DSE=',1Pd9.2/(' ypSAP{',i2,'}=',
     20PF11.7,'   beta_{',i2,'}=',1Pd19.11,' (+/-',d8.1,')   PS=',d8.1))
  629 FORMAT(8x,'Rref ='f15.12)
  630 FORMAT(10x,'De =',f13.6,' (+/-',f12.6,')   PS=',1pd8.1)
  632 FORMAT(' Use exponent switching function with   Asw=',F9.6,
     1  '   Rsw=',F9.6)
  634 FORMAT(' Use Huang exponent fx:  beta(R)= betaINF*y_',i2,
     1 ' + (1-y_',I1,')* Sum{beta_i*[y_',i1,']^i}')
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
     2  (6x,'c_{',i2,'} =',d17.9,' (+/-',d8.1,')   PS=',d8.1))
  660 FORMAT(8x,'VMIN =',f13.5,' (+/-',f12.6,')   PS=',1pd8.1)
  662 FORMAT(10x,'Re =',f13.9,' (+/-',f12.9,')   PS=',1pd8.1)
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
      PARAMETER (MXDATA=1501, MXPARM=40, MXMLR= 8)
      INTEGER  i,j,IDAT, NPOW,NPARM,NDATA, IFXP(MXPARM),JFXRe,JFXDe,
     1  JFXVMIN
      REAL*8  YC,PV(NPARM),PD(NPARM),PS(NPARM),TT(0:20),RHOdR,RMSR,RTPp,
     1  RTPq,Rep,AREF,AREFp,AREFq,ype,dype,betaINF,yp,yq,fsw,yPOW,XP,
     2  XPW,DER,TTM,TTMM,DERP,SUM,DSUM,AA,BB,FCT,ULR,ULRe,dULRe,d2ULRe,
     3  VCN,DDER,T0,T0P,T1,RE3,RE6,RE8,RTP3,RTP6,RTP8,dULRedRe,RDIST,
     4  C3VAL,C6adj,C9adj
c-----------------------------------------------------------------------
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR(MXMLR),p,q,NS,NL,
     1                                                        NPHI,SAP
      REAL*8 Re,De,VMIN,RREF,Asw,Rsw,M2,as,bs,RHOd,
     1  CmVAL(MXMLR),RTP(MXDATA),SAS(MXDATA,MXPARM)
      COMMON /DATABLK/Re,De,VMIN,RREF,Asw,Rsw,M2,as,bs,
     1  RHOd,CmVAL,RTP,SAS,PSEL,IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR,p,q,
     2                                                  NS,NL,NPHI,SAP
c-----------------------------------------------------------------------
      SAVE JFXRe,JFXDe,JFXVMIN, AREF,AREFp,Rep,betaINF,AA,BB, ULRe,
     1  dULRedRe
c=======================================================================
      IF(ABS(IDAT).LE.1) THEN
          JFXRe= IFXP(NPHI+1)
          JFXDe= IFXP(NPHI+2)
          JFXVMIN= IFXP(NPHI+3)
          ENDIF
      RDIST= RMSR
      IF(IDAT.GT.0) RDIST= RTP(IDAT)
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
                  dULRedRe= -RE3*(0.25d0*CmVAL(1)*(9.d0 + T0P)
     1               + RE3*(C6adj*(3.d0 + T0P) + CmVAL(4)*(3.d0 - T0P)
     2                                           + RE3*9.d0*C9adj))/Re
                  IF(NCMM.GT.4) THEN
                      dULRedRe= dULRedRe -RE8*4.d0*(CmVAL(5) 
     1                *(3.d0 + T0P) + CmVAL(6)*(3.d0 - T0P))/(3.d0*Re)
                      ENDIF
                ELSE
c** For normal inverse-power sum  ULR(r)  ....
                  ULRe= 0.d0
                  dULRedRe= 0.d0
                  DO  j= 1,NCMM
                      AA= CmVAL(j)/Re**MMLR(j)                  
                      ULRe= ULRe+ AA
                      dULRedRe= dULRedRe - MMLR(j)*AA/Re
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
              IF(Asw.LE.0.d0) THEN
                  yPOW= 1.d0 - yp
                ELSE
                  fsw= 1.d0/(dexp(Asw*(RDIST-Rsw)) + 1.d0)
                  yPOW= fsw
                ENDIF
              SUM= PV(1)*yPOW
              DSUM= 0.d0
              IF(NPOW.GE.2) THEN
                  DO  j= 2,NPOW
                      IF(RREF.LE.0.d0) DSUM= DSUM + PV(j)*(j-1)*yPOW
                      yPOW= yPOW*yq
                      SUM= SUM+ yPOW*PV(j)
                      ENDDO
                  ENDIF
              IF(Asw.LE.0.d0) THEN
                  XP= SUM + betaINF*yp
                ELSE
                  XP= SUM + betaINF*(1.d0 - fsw)
                ENDIF
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
c** For normal inverse-power sum  ULR(r)  ....
              ULR= 0.d0
              DO  J= 1,NCMM
                  ULR= ULR+ CmVAL(j)/RDIST**MMLR(j)
                  ENDDO
            ENDIF
          XPW= DEXP(-XP*ype) * ULR/ULRe
          YC= De*(1.d0 - XPW)**2 + VMIN
          DER= 2.d0*De*(1.d0- XPW)*XPW
          IF(Asw.LE.0.d0) THEN
              yPOW= DER*ype*(1.d0- yp)
            ELSE
              yPOW= DER*ype*fsw
            ENDIF
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
          IF(JFXDe.LE.0) THEN
              IF(Asw.LE.0.d0) THEN
                  PD(NPHI+2)= (1.d0- XPW)**2 + DER*ype*yp/De
                ELSE
                  PD(NPHI+2)= (1.d0- XPW)**2 + DER*ype*(1.d0- fsw)/De
                ENDIF
              ENDIF
          IF(JFXVMIN.LE.0) PD(NPHI+3)= 1.d0
c** If appropriate, also get partial derivative w.r.t. Re
          IF(JFXRe.LE.0) THEN
              dype= -0.5d0*(p/RE)*(1.d0 - yp**2)
              IF(Asw.LE.0.d0) THEN
c ... either for Huang exponent function ...
                  IF(RREF.LE.0.d0) THEN
                      DSUM= betaINF - SUM/(1.d0-yp) + DSUM 
                    ELSE
                      DSUM= 0.d0
                    ENDIF
                  PD(NPHI+1)= DER*(dype*(XP + ype*DSUM) 
     1                               + (1.d0 - ype*yp)*dULRedRe/ULRe )
                ELSE
c ... or for Hajigeorgiou exponent function ...
                  IF(RREF.GT.0.d0) DSUM= 0.d0
                  PD(NPHI+1)= DER*((ype*(1.d0- fsw)*MMLR(1)
     1                           - MMLR(1))/Re + dype*(XP + ype*DSUM))
                ENDIF
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
c** Evaluate ULR & its first 2 deriv. at  Re ... 
              IF(IDF.EQ.2) THEN
c... if using Scoles' damping function ...
                  DO  j= 1,NCMM
                      FCT= dexp(-3.97d0*RHOd*Re/MMLR(j)
     1                   - 0.39d0*(RHOd*Re)**2/DSQRT(DFLOAT(MMLR(j))))
                      VCN= CmVAL(j)*((1.d0- FCT)/Re)**MMLR(j)
                      ULRe= ULRe+ VCN
                      dULRe= dULRe- MMLR(j)*VCN/Re
                      DDER= FCT*(3.97d0*RHOd/MMLR(j)
     1                     + 0.78d0*Re*RHOd**2/DSQRT(DFLOAT(MMLR(j))))
                      dULRe= dULRe+  MMLR(j)*DDER*VCN/(1.d0- FCT)
                      d2ULRe= d2ULRe + MMLR(j)*(MMLR(j)+1)*VCN/Re**2
     1                    - 2.d0*MMLR(j)**2 *DDER*VCN/((1.d0- FCT)*Re)
     2                   + ((MMLR(j)-1)*DDER**2 - DDER*(1.d0- FCT)/FCT
     3        + (1.d0- FCT)*FCT*0.78d0*RHOd**2/DSQRT(DFLOAT(MMLR(j))))
     4                                     *MMLR(j)*VCN/(1.d0- FCT)**2
                      ENDDO
                  ENDIF
              IF(IDF.EQ.1) THEN
c... Using Tang-Toennies damping function ...
                  TT(0)= 1.d0
                  yPOW= 1.d0
                  RHOdR= RHOd*Re
                  DO  j= 1,MMLR(NCMM)
                      yPOW= yPOW*RHOdR/DFLOAT(J)
                      TT(J)= TT(J-1)+ yPOW
                      ENDDO
                  yPOW= DEXP(-RHOdR)
                  ULRe= 0.d0
                  dULRe= 0.d0
                  d2ULRe= 0.d0
                  DO  j= 1,NCMM
                      TTM= (1.d0- yPOW*TT(MMLR(j)))*CmVAL(j)/Re**MMLR(j)
                      ULRe= ULRe+ TTM
                      TTMM= yPOW*RHOd*(TT(MMLR(j)) - TT(MMLR(j)-1))
     1                                           *CmVAL(j)/Re**MMLR(j)
                      dULRe= dULRe+ TTMM - MMLR(j)*TTM/Re
                      d2ULRe= d2ULRe + MMLR(j)*(MMLR(j)+1)*TTM/Re**2
     1               + yPOW*RHOd**2*(-TT(MMLR(j)) + 2.d0*TT(MMLR(j)-1)
     2                         -TT(MMLR(j)-2))  - 2.d0*MMLR(j)*TTMM/Re
                      ENDDO
                  ENDIF
c-----------------------------------------------------------------------
              AA= De - ULRe - dULRe/PV(1)
              BB= 2.d0*(De - ULRe) - dULRe/PV(1)
              ENDIF
          RTPp = RDIST**p
          yp= (RTPp - AREFp)/(RTPp + AREFp)
          ULR= 0.d0
c ... evaluate ULR at the actual distance ...
          IF(IDF.EQ.2) THEN
c... if using Scoles' damping function ...
              DO  j= 1,NCMM
                  FCT= dexp(-3.97d0*RHOd*RDIST/MMLR(j)
     1             - 0.39d0*(RHOd*RDIST)**2/DSQRT(DFLOAT(MMLR(j))))
                  ULR= ULR+ CmVAL(j)*((1.d0- FCT)/RDIST)**MMLR(j)
                  ENDDO
              ENDIF
          IF(IDF.EQ.1) THEN
c... if using Tang-Toennies damping function ...
              RHOdR= RHOd*RDIST
              yPOW= DEXP(-RHOdR)
              TT(0)= yPOW
              DO  j= 1,MMLR(NCMM)
                  yPOW= yPOW*RHOdR/DFLOAT(J)
                  TT(J)= TT(J-1)+ yPOW
                  ENDDO
              DO  j=1,NCMM
                  ULR=ULR+CmVAL(j)*(1.d0-TT(MMLR(j)))/RDIST**MMLR(j)
                  ENDDO
              ENDIF
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
                  DSUM= -DSUM*0.5d0*(p/Re)*(1.d0-yp**2)
                ENDIF
              PD(NPHI+1)= DER*(DSUM - SUM) - XP*((XP- 2.d0)*dULRe 
     1                                       - (XP-1.d0)*D2ULRe/PV(1))
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
