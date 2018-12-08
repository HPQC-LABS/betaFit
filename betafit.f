
c***********************************************************************
c********  Program  betaFIT_2.1  version of  19 March 2013  **********
c***********************************************************************
c* Program to fit NTP read-in potential fx. values {RTP(i),VTP(i)} 
c  (with or without individual weights) to a chosen analytic form.
c***********************************************************************
      INTEGER MXDATA, MXPARM, MXMLR, NxPSE
      PARAMETER (MXDATA=1501, MXPARM=43, MXMLR= 15)
      INTEGER I,J,INFL,ITER,IROUND,JROUND,LPPOT,ROBUST,prFIT,prNLL,
     1  prDIFF,M,NPARM,NPR,NTP,NLIN,NS,NL,NDGF,IFXP(MXPARM)
      REAL*8 BETA(0:MXPARM),PV(MXPARM),PU(MXPARM),PS(MXPARM),
     1 CM(MXPARM,MXPARM),DYDP(MXDATA,MXPARM),VTP(MXDATA),
     2 uVTP(MXDATA),betay(MXDATA),Ubetay(MXDATA),YD(MXDATA),
     3 ypPSE(MXDATA),xPSE(MXDATA),rPSE(MXDATA),rKL(1:MXDATA,1:MXDATA),
     4 DM(MXMLR),DMP(MXMLR),DMPP(MXMLR),
     3 betaINF,UNC,yPOW,DSE,TSTPS,TSTPU,DSEB, TCM,UM,diff,
     4 Rep,AREF,AREFp,AREFq,RTPp,RTPq, AA,BB,ULR,dULR,FCT,RAT,UMAX,
     5 XX,YY,YH,yp,yq,ypRE,ReDE, ReIN,DeIN,VMINin,ULRe,dULRe,RE3,RE6,
     6 RE8,T0,T1,C6adj,C9adj,RTP3,RTP6,RTP8,RH,RR,RB,RBB,VV,VB,VBB,
     7 SCALC,yMIN,DRMSD ,BETA0,BETAN, RPR1,dRPR, SL,SLB,dSL,ycalc,dd
      REAL*8  DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),DEIGRe(1,1)
     1  ,DEIGDe(1,1)
      CHARACTER*4  NNAME,NAME(4)
      DATA NAME/' EMO',' MLR','DELR','GPEF'/
c-----------------------------------------------------------------------
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,sVSR2,IDSTT,KDER,NCMM,p,q,Nbeta,
     1                             APSE,MMLR(MXMLR),IFXCm(MXMLR)
      REAL*8 Re,De,VMIN,RREF,M2,as,bs,rhoAB, SUM,
     1  CmVAL(MXMLR),RTP(MXDATA),SAS(MXDATA,MXPARM)
      COMMON /DATABLK/Re,De,VMIN,RREF,M2,as,bs,rhoAB,CmVAL,RTP,SAS,PSEL,
     1 IFXRe,IFXDe,IFXVMIN,sVSR2,IDSTT,KDER,NCMM,MMLR,p,q,Nbeta,APSE,SUM
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
c  LPPOT specifies whether (LPPOT > 0) or not (LPPOT.le.0) the code should 
c             read instructions to generate an array of potential function 
c             values on some specified mesh and range
c       If LPPOT > 0, read a separate mesh and range for eash {N,p,q,Rref} case
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
      READ(5,*) PSEL, NTP, UNC, IROUND, LPPOT, prFIT, prDIFF
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
c** For an MLR (PSEL=2) or DELR (PSEL=3) potential, read the number of 
c   long-range terms NCMM used to define long-range potential tail:  
c    V(r)= De - \sum{CmVAL/r^MMLR}
c** If rhoAB .LE. 0.0  have NO damping functions: all  Dm(R)= 1.0
c   If rhoAB > 0.0  it is the molecule-dependent radial scaling factor
c                    of Douketis et al. [JCP 76, 3057 (1982)]
c     rhoAB =  2*rhoA*rhoB/(rhoA+rhoB)   where  rhoA  is the ionization
c           potential ratio  (I_p^A/I_p^H)^{2/3}  for atom A vs. atomic H
c  For rhoAB > 0.0,  sVSR2 specifies damping s.th.  Dm(r)/r^m --> r^{sVSR2/2}
c                    IDSTT > 0  use generalized Douketis et al. damping fx.
c                    IDSTT > 1  use DS damping for s=-1/2 AND impose very
c                       short-range 1/r behaviour with  C1= IDSTT = Z1*Z2
c                    IDSTT .LE. 0  use generalized Tang-Toennies damping fx
c
c** APSE is an integer: > 0  to use A.Pashov natural Spline in MLR Exponent
c                    : .le. 0  for normal constrained polynomial exponent
c*  For an MLR potential with  APSE > 0,  yMIN  is a negative{!} number
c      (-1 .le. yMIN < 0) which specifies the lower bound on the yp values 
c      selected to define the MLR exponent spline.  
c** For each long-range term read power  MMLR(i)  & coefficient CmVAL(i)
c** For special Aubert-Frecon case,  NCMM= 6,  MMLR= {3,0,6,6,8,8} and 
c  coefficients are: CmVAL(1)= C3(sig), CmVAL(2)= ASO, CmVAL(3)= C6(Sig),
c  CmVAL(4)= C6(pi), CmVAL(5)= C8(sig) & CmVAL(6)= C8(pi).
c=======================================================================
      IF((PSEL.EQ.2).OR.(PSEL.EQ.3)) THEN
          READ(5,*) NCMM, rhoAB, sVSR2, IDSTT, APSE, yMIN
          DO m= 1, NCMM
cc            READ(5,*) MMLR(m), CmVAL(m), IFXCm(m)
              READ(5,*) MMLR(m), CmVAL(m)
              ENDDO
c=======================================================================
          IF((NCMM.GE.4).AND.(MMLR(2).LE.0)) THEN
c** For Lyon 2x2 or 3x3 treatment of A- or b-state alkali dimers ...
              IF(MMLR(2).EQ.0)
     1   WRITE(6,602) 'A-State',CmVAL(2),CmVAL(1),(CmVAL(i),i=3,NCMM)
              IF(MMLR(2).EQ.-2)
     1   WRITE(6,602) 'b-State',CmVAL(2),CmVAL(1),(CmVAL(i),i=3, NCMM)
              IF(MMLR(2).EQ.-1)
     1              WRITE(6,6022) CmVAL(2),CmVAL(1),(CmVAL(m),m=3,NCMM)
            ELSE
c** For 'normal' inverse-power sum uL, with/without damping
              IF(rhoAB.LE.0) THEN
                  WRITE(6,698) (MMLR(m),CmVAL(m),m= 1,NCMM)
                ELSE
                 IF(IDSTT.GT.0) WRITE(6,696) 'DS',rhoAB,sVSR2,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
                 IF(IDSTT.LE.0) WRITE(6,696) 'TT',rhoAB,sVSR2,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
                ENDIF
            ENDIF
          ENDIF
c-----------------------------------------------------------------------
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
      IF((PSEL.EQ.2).AND.(NCMM.GE.4).AND.(MMLR(2).LE.0)) THEN
c** Shift of asymptote for special Li2 2x2 A-state case
          DO  I= 1,NTP
              VTP(I)= VTP(I) - 0.5d0*CmVAL(2)
              ENDDO
          ENDIF
c=======================================================================
c
  690 FORMAT(' uLR(r)  is a damped TT inverse-power sum with',6x,
     1'   rhoAB=',  F7.4,'   sVSR2=',i3/(49x,'C',i2,'=',1PD15.8:))
  696 FORMAT(' uLR(r) inverse-power terms incorporate ',a2,' damping wit
     1h  rhoAB=',f10.7/8x ,'defined to give very short-range damped uLR-
     2term behaviour  r^{',i2,'/2}'/(49x,'C',i2,'=',1PD15.8:)) 
  698 FORMAT(' uLR(r)  is a simple inverse-power sum with',6x,'C',I2,
     1  '=',1PD15.8:/(49x,'C',i2,'=',D15.8:))
  600 FORMAT(' Fit an ',A4,'p  potential function to the input points'/
     1  1x,26('==')/5x,'with initial   VMIN=',f13.4,'   Re=',f11.8:,
     2  '   De=',f11.4)
  601 FORMAT(7x,'in which   y_{p/q}(r)= [r^{p/q} -',f7.4,'^{p/q}]/[r^{p/
     1q} +',f7.4,'^{p/q}]' )
  602 FORMAT(5x,'Use Lyon 2x2  ',A7,'  uLR(r)  with   Aso=',F10.6/47x,
     1'C_3(1Sigma)='  ,1PD15.7:/47x,'C_6(Sigma) =',D15.7:/47x,'C_6(Pi)',
     2 '    =',D15.7:/  47x,'C_8(Sigma) =',1PD15.7:/47x,'C_8(Pi)    =',
     3  D15.7)
 6022 FORMAT(' Use Lyon 3x3  uLR(r)  with   Aso=',F10.6,'   C_3(Sigma)='
     1  ,1PD15.7:/47x,'C_6(Sigma)=',D15.7:/47x,'C_6(Pi)   =',D15.7:/
     2  47x,'C_8(Sigma)=',1PD15.7:/47x,'C_8(Pi)   =',D15.7)
  603 FORMAT(7x, 'in which   y_{p/q}(r)= [r^{p/q} - Re^{p/q}]/[r^{p/q} +
     1 Re^{p/q}]' )
  652 FORMAT('  & define beta(y(r)) as a natural spline through points a
     1t the',i4,'  yp values:'/(2x,7F11.7))
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
  611 FORMAT(3x,'Start with this long-range tail and   beta(0)=',f10.6)
  612 FORMAT(' *** Linearization Fails!!  for datum',I4,'  r =',f8.2/
     1  '  De too small or VMIN too high')
c=======================================================================
c** Now ... loop over different {p,q,NS,NL} combinations till end of data
c**  p and q are powers used to define radial variable in the exponent
c       beta(r)= yp*betaINF + [1 - yp]*sum{beta_i*yq^i}  where  
c       ya=(R^a - AREF^a)/(R^a + AREF^a)   for  a= p  or  q
c** NL= is the order of the\beta(y)  exponent expansion for EMO, DELR,
c    or MLR potentials with APSE.leq.0.
c   NS is a dummy variable for EMO, DELR & 'ordinary' MLR with APSE.le.0
c** For an MLR potential with  APSE > 0,  Nbeta=(NS+NL+1) is the number 
c   of yp values to be used to define the exponent spline used for beta(y):
c   NL specifies the number of yp values from  yp= 0  to  yp= 1.0
c   NS specifies the number of (equally spaced) yp values for  
c      yMIN .leq. yp = yMIN to 0.0  
c** For a GPEF potential, consider powers ranging from  NS(.ge.0) to NL
c* RREF   defines the reference distance in the expansion variable
c      - for  RREF.le.0 , define parameter  RREF = Re
c      - for  RREF.gt.0 , fix parameter  RREF   at its read-in value
c-----------------------------------------------------------------------
   20 READ(5,*,END= 999) p, q, NS, NL, RREF
c-----------------------------------------------------------------------
      IF(p.LE.0) GOTO 999
c=======================================================================
c** If desired, print the resulting potential for this case to channel-8
c     at for NPR  distances, starting at  RPR1 with a mesh of DPR
c*   IF(NPR.leq.0) omit printout for this case
c-----------------------------------------------------------------------
      IF(LPPOT.GT.0) READ(5,*) NPR, RPR1, dRPR
c-----------------------------------------------------------------------
      Nbeta= NL + 1
      IF(APSE.GE.0) Nbeta= NS+ NL+ 1
      NPARM= Nbeta
      Re= ReIN
      De= DeIN
      VMIN= VMINin
      IF(PSEL.LE.3) WRITE(6,600) NNAME, VMIN, Re, De
      IF(PSEL.EQ.4) THEN
          WRITE(6,600) NNAME, VMIN, Re
          WRITE(6,605) as,bs
          ENDIF
      IF((PSEL.EQ.1).OR.(PSEL.EQ.3)) THEN
          IF(RREF.GT.0.d0) THEN
              AREF= RREF
              WRITE(6,632) p,p,RREF,p,p,RREF,p
            ELSE
              AREF= RE
              IF(RREF.LE.0.d0) WRITE(6,633) p,p,p,p,p
            ENDIF
          ENDIF
      IF(PSEL.EQ.2) THEN
          IF(APSE.GT.0) WRITE(6,650) NS, NL
          IF(APSE.LE.0) WRITE(6,634) 
          IF(RREF.gt.0.d0) THEN
              AREF= RREF
              WRITE(6,601) AREF,AREF
            ELSE
              AREF= Re
              WRITE(6,603)
            ENDIF
          IF(APSE.GT.0) THEN
c** Set up yp abscissa values for APSpline Exponent case
              YH= -yMIN/NS
              ypPSE(1)= yMIN
              DO  I= 2,NS
                  ypPSE(I)= ypPSE(I-1) + YH
                  ENDDO
              NPARM= NS+ 1
              ypPSE(NS+1)= 0.d0
              YH= 1.d0/NL
              DO  I= 1,NL-1
                  ypPSE(NPARM+I)= ypPSE(NPARM+I-1) + YH
                  ENDDO
              NPARM= NPARM+ NL
              ypPSE(NPARM)= 1.d0
              WRITE(6,652) NPARM,(ypPSE(i), i= 1,NPARM)
              ENDIF
          ENDIF
      AREFp= AREF**p
      AREFq= AREF**q
      IF((PSEL.EQ.2).AND.(NCMM.GE.4).AND.(MMLR(2).LE.0)) THEN
c** Prepare uLR for Lyon 2x2 or 3x3 treatment of 2S + 2p Li2 dimers ...
          IF(MMLR(2).EQ.0) WRITE(6,602) 'A-state',CmVAL(2),
     1                                    CmVAL(1),(CmVAL(m),m=3,NCMM)
c... For Lyon treatment of b-state alkali dimers ...
          IF(MMLR(2).EQ.-2) WRITE(6,602) 'b-state',CmVAL(2),
     1                                    CmVAL(1),(CmVAL(m),m=3,NCMM)
          IF(MMLR(2).EQ.-1)
     1              WRITE(6,6022) CmVAL(2),CmVAL(1),(CmVAL(m),m=3,NCMM)
          ENDIF
c... Prepare uLR as 'conventional' (damped or non-damped) inverse-power sum
      IF((PSEL.EQ.2).OR.(PSEL.EQ.3)) THEN
          IF(rhoAB.LE.0) THEN
              WRITE(6,698) (MMLR(m),CmVAL(m),m= 1,NCMM)
            ELSE
              IF(IDSTT.GT.0) WRITE(6,696) 'DS',rhoAB,sVSR2,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
              IF(IDSTT.LE.0) WRITE(6,696) 'TT',rhoAB,sVSR2,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
            ENDIF
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
c     IF((PSEL.EQ.1).OR.(PSEL.EQ.3).AND.(beta(0).le.0)) THEN
      IF(PSEL.EQ.1) THEN
          q=0
c ... first define ordinate array
          DO  i= 1,NTP
              RTPp= RTP(i)**p
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              IF(RTP(i).GT.Re) THEN
                  T0= (1.d0 - DSQRT((VTP(i)-VMIN)/De))
                      IF(T0.LT.0.d0) THEN
                          WRITE(6,612) i,RTP(i)
                          STOP
                          ENDIF
                  betay(i)= - DLOG(T0)
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
              DO  j= 1, NPARM
                  DYDP(i,j)= yPOW
                  yPOW= yPOW*yp
                  ENDDO
c%%  elective printout for testing
c%%           if(i.eq.1) write(8,700) 
c%%           write(8,702) rtp(i),yp,ypRE,vtp(i),betay(i)
c%%  1                                           (dydp(i,j),j=1,Nbeta)
c 702 format(f6.3,f8.4,f9.2,1P2d13.5:/(14x,5d13.5))
c%%
              ENDDO
          CALL LLSQF(NTP,NPARM,MXDATA,MXPARM,betay,Ubetay,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
          IF(prFIT.GT.0) WRITE(6,620) NNAME,p,NL,DSE,
     1                        ('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPARM)
          ENDIF

c=======================================================================
      IF(PSEL.EQ.2) THEN
c-----------------------------------------------------------------------
c*** Preliminary linearized fit for an  MLR  potential ...
          IF(p.LE.(MMLR(NCMM)-MMLR(1))) WRITE(6,616) p,NCMM,MMLR(NCMM),
     1                                                        MMLR(1)
c** First define array of exponent values with uncertainties defined by 
c  the assumption that all potential values have equal uncertainties UNC
c
c... Begin by determining uLR(Re) and  betaINF
          IF((NCMM.GE.4).AND.(MMLR(2).LE.0)) THEN
c** For special Aubert-Frecon 2x2 or 3x3 Li2 {3,0,6,6,8,8} cases ...
              RE3= 1.d0/Re**3
              RE6= RE3*RE3
              RE8= RE6/Re**2
              C6adj= CmVAL(3) + 0.25D0*CmVAL(1)**2/De
              C9adj= 0.5d0*CmVAL(1)*C6adj/De
              IF((MMLR(2).EQ.0).OR.(MMLR(2).EQ.-2)) THEN
c ... for Aubert-Frecon 2x2 case ...
                  T1= (0.5d0*CmVAL(1)+ (C6adj-CmVAL(4))*RE3)*RE3/3.d0
                  T1= T1+ (CmVAL(5)- CmVAL(6))*RE8/3.d0
                  T0= DSQRT((T1- CmVAL(2))**2 + 8.d0*T1**2)
                  ULRe= 0.5d0*( -CmVAL(2)+ RE3*(1.5d0*CmVAL(1) 
     1           + RE3*(C6adj + CmVAL(4)))) + 0.5d0*T0 + C9adj*RE6*RE3
c   ...  now add C8 terms ...
                  ULRe= ULRe+ 0.5d0*(CmVAL(5)+ CmVAL(6))*RE8
c ** Now option for the b ^3\Pi_u state
                  IF(MMLR(2).EQ.-2) ULRe= ULRe - T0
                  ENDIF
              IF(MMLR(2).EQ.-1) THEN
c ... for Aubert-Frecon 3x3 Li2(c) {3,0,6,6,8,8} case ...
                  CALL AF3x3potret(Re,CmVAL(2),CmVAL(1),C6adj,CmVAL(5),
     1                      De,ULRe,DEIGM1,DEIGM3,DEIGM5,DEIGR,DEIGDe)
ccc             c9adj= 0.0
                  ULRe= ULRe + C9adj*RE6*RE3
                  ENDIF
            ELSE
c** For normal inverse-power sum MLR case, with or without damping
              KDER=0
              IF(rhoAB.GT.0.d0) CALL dampF(Re,rhoAB,NCMM,MMLR,sVSR2,
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
c.... completed determination of uLR(Re) and  betaINF
          Rep= RE**p
c.... Now loop to prepare linearized MLR arguments & derivatives
          DO  i= 1, NTP
              RTPp= RTP(i)**p
              RTPq= RTP(i)**q
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              yq= (RTPq - AREFq)/(RTPq + AREFq)
              ypRE= (RTPp - Rep)/(RTPp + Rep)
              xPSE(i)= yp
              IF((NCMM.GE.4).AND.(MMLR(2).LE.0)) THEN
c ... extension for Aubert-Frecon 2x2 & 3x3 {3,0,6,6,8,8} Li2 cases ...
                  RTP3= 1.d0/RTP(i)**3
                  RTP6= RTP3*RTP3
                  RTP8= RTP6/RTP(i)**2
                  IF((MMLR(2).EQ.0).OR.(MMLR(2).EQ.-2)) THEN
c ... for Aubert-Frecon Li2(A) 2x2 {3,0,6,6,8,8} case ...
                      T1= (0.5d0*CmVAL(1)+ (C6adj- CmVAL(4))*RTP3)
     1                                                      *RTP3/3.d0
c... now add C8 terms ...
                      T1= T1+ (CmVAL(5)- CmVAL(6))*RTP8/3.d0
                      T0= DSQRT((T1- CmVAL(2))**2 + 8.d0*T1**2)
                      ULR= 0.5d0*( -CmVAL(2) + RTP3*(1.5d0*CmVAL(1)
     1                       + RTP3*(CmVAL(3) + CmVAL(4)))) + 0.5d0*T0
     2                                               + C9adj*RTP3*RTP6
c... now add C8 terms ...
                      ULR= ULR+ 0.5d0*(CmVAL(5)+CmVAL(6))*RTP8
c** Now - for the b ^3\Pi_u case ....
                      IF(MMLR(2).EQ.-2) ULR= ULR - T0
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
     1                                     sVSR2,IDSTT,KDER,DM,DMP,DMPP)
                  ULR= 0.d0
                  DO  m= 1,NCMM
                      T0= CmVAL(m)/RTP(i)**MMLR(m)
                      IF(rhoAB.GT.0.d0) T0= T0*DM(m)
                      ULR= ULR + T0
                      ENDDO
                ENDIF
              IF(RTP(i).GT.Re) THEN
                  T0= (1.d0 - DSQRT((VTP(i)-VMIN)/De))
                  IF(T0.LT.0.d0) THEN
                      WRITE(6,612) i,RTP(i)
                      STOP
                      ENDIF
                  betay(i)= - DLOG(T0*ULRe/ULR)
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
              IF(APSE.LE.0) THEN
c** Subtract the \beta_\infty term to yield polynomial for fitting
c... For Huang MLR polynomial exponent function
                  betay(i)= betay(i)- betaINF*yp*ypRE
                  yPOW= ypRE*(1.d0- yp)
c... then create partial derivative array for linearized fit ...
                  DO  j= 1, NPARM
                      DYDP(i,j)= yPOW
                      yPOW= yPOW*yq
                      ENDDO
                  ENDIF
c%%%% printout for testing  polynomial exponent
c%%            if(i.eq.1) write(8,700) 
c%700 format('  RTP        yp       yp^{eq}       VTP(i)          uLR',
c%   1 10x,' beta(i)' )
c%%            write(8,702) rtp(i),yp,ypRE,vtp(i),ULR,betay(i),
c%%  1                                           (dydp(i,j),j=1,Nbeta)
c%             write(8,706) rtp(i),yp,ypRE,VTP(i),Ulr,betay(i)
c%706 Format( f7.4,2f12.8,4(1Pd15.7))
c%%%%
              ENDDO
          IF(APSE.LE.0) THEN
              CALL LLSQF(NTP,Nbeta,MXDATA,MXPARM,betay,Ubetay,DYDP,YD,
     1                                                PV,PU,PS,CM,DSE)
              IF(prFIT.GT.0) THEN
                  IF(RREF.GE.0.d0) WRITE(6,621) NNAME,p,q,RREF,NL,
     1                    DSE,('beta',j-1,PV(j),PU(j),PS(j),j= 1,Nbeta)
                  IF(RREF.LT.0.d0) WRITE(6,623) NNAME,p,q,NL,
     1                    DSE,('beta',j-1,PV(j),PU(j),PS(j),j= 1,Nbeta)
                  ENDIF
              ENDIF
          IF(APSE.GT.0) THEN
c** For Pashov spline-exponent, use spline through linearized input 
c  exponent values to define initial trial spline values of \beta(i) 
c  at the specified  xPSE(i) 0:qvalues !
              NxPSE= NTP
              DO  I= 1,NTP
                  IF(DABS(betay(I)).LT.1.d-04) THEN
c... If  xPSE\approx 0.0, collapse the xPES array to omit this point ... 
                      NxPSE= NxPSE-1
                      NTP= NTP-1
                      DO  J= I,NxPSE
                          xPSE(J)= xPSE(J+1)
                          betay(J)= betay(J+1)
                          ENDDO
                      ENDIF
                  ENDDO
              NxPSE= NxPSE+ 1
              xPSE(NxPSE)= 1.d0
              betay(NxPSE)= betaINF
              rtp(NxPSE)= 9.d99
              CALL Lkoef(NxPSE,xPSE,rKL,MXDATA)
c**  Elective printout for testing
cc            write(6,699) (rtp(i),xPSE(i),betay(i), betay(i)/xPSE(i),
cc   1                                                      I= 1,NxPSE)
cc699 FORMAT('  R=', f10.6,'   xPSE=',F12.9,'   beyay=', f12.9, f15.9)
              DO  I= 1, NPARM
                  XX= ypPSE(I)
                  IF(I.LT.NPARM) rPSE(I)= AREF*((1.d0 + XX)/
     1                                  (1.d0 - XX))**(1.d0/DFLOAT(p))
c... Now, use a spline through the exponent values defined by the input 
c    points to generate values of that exponent at the desired 
c    spline-definition points
                  PV(I)= 0.d0
                  DO  m= 1, NxPSE
                      PV(I)= PV(I) + 
     1               Scalc(XX,m,NxPSE,xPSE,rKL,MXDATA)*betay(m)/xPSE(m)
                      ENDDO
                  ENDDO
              NLIN= (NPARM+1)/2
              rPSE(NPARM)= 9.d99
              IF(prFIT.GT.0) WRITE(6,653) NNAME,p,
     1              ((ypPSE(I),rPSE(I),PV(I),I= J,NPARM,NLIN),J=1,NLIN)
  653 FORMAT(/' Linearized ',A4,'{p=',i1,'}-APSE treatment yields:'/
     1  2(4x,'ypSE',8x,'rPSE',6x,'betay  ')/(2(F12.6,f10.6,f10.6) ))
c
c... Finally, create the fixed array of S(m,x) to define the exponent
c   and its partial derivatives in subsequent fits, and compare predictions
c  of our preliminary spline function with the input data ...
              CALL Lkoef(NPARM,ypPSE,rKL,MXDATA)
              dd= 0.d0
              DO  I= 1,NxPSE
                  ycalc= 0.d0
                  DO  m= 1,NPARM
                      SAS(I,m)= Scalc(xPSE(I),m,NPARM,ypPSE,rKL,MXDATA)
                      ycalc= ycalc+ SAS(I,m)*PV(m)
                      ENDDO
                  dd= dd+ (ycalc - betay(I)/xPSE(I))**2
                  diff= ycalc - betay(I)/xPSE(I)
                  WRITE(6,655) xPSE(I),ycalc,betay(I)/xPSE(I),diff
                  ENDDO
              dd= DSQRT(dd/NxPSE)
              WRITE(6,656) NS, NL, RREF, dd
              ENDIF
          ENDIF
  655 FORMAT('   At  X=',F8.4,'   Y_{calc}=',F11.7,'   Y_{input}=',
     1  F11.7,F13.8)
  656 FORMAT('  SE-MLR Linearization:  NS=',I2,',  NL=',I3,
     1   ',  R_{ref}=', F7.3,'  yields   dd=',F9.5/)

c=======================================================================
c** Preliminary linearized fit for a  DELR  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.3) THEN
          BETA(0)= 1.d0
          WRITE(6,611) BETA(0)
c ... NOTE need iteration to determine self-consistent  beta(0)  value 
c    First generate  A & B  from input Re, De and trial beta(0)
          ITER= 0
   40     ULRe= 0.d0
          dULRe= 0.d0
          KDER= 1
          CALL dampF(Re,rhoAB,NCMM,MMLR,sVSR2,IDSTT,KDER,DM,DMP,DMPP)
          DO  m= 1,NCMM
              FCT= dexp(-3.95d0*rhoAB*Re/MMLR(m)
     1               - 0.39*(rhoAB*Re)**2/DSQRT(DFLOAT(MMLR(m))))
              AA= CmVAL(m)/Re**MMLR(m)
              ULRe= ULRe+ AA*DM(m)
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
              CALL dampF(RTP(i),rhoAB,NCMM,MMLR,sVSR2,IDSTT,KDER,DM,DMP,
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
              DO  j= 1, NPARM
                  DYDP(i,j)= yPOW
                  yPOW= yPOW*yp
                  ENDDO
              ENDDO
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
          IF(prFIT.GT.0) WRITE(6,620) NNAME,p,NL,DSE,
     1                        ('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPARM)
          ENDIF

      IF(PSEL.LE.3)  THEN
c======================================================================
c** Now ... do direct non-linear fits to potential values ... first with
c      Re and/or VMIN and De fixed, and then freeing them up too ...
c=======================================================================
c* FIRST optimize  BETA(j)'s (and VMIN) with  Re and De held fixed!
          DO  j= 1,NPARM
              BETA(j-1)= PV(j)
              IFXP(j)= 0
              ENDDO
c**?? Is this how I implement fixed limiting  C1  for MLR?
          IF(IDSTT.GT.1) IFXP(Nbeta)= IDSTT
c** Fix parameter value for  yp= 1  in SE-MLR
          IF(APSE.GT.0) IFXP(NPARM)= 1
c** For initial run, fix the  VMIN, De and Re values
          IFXP(NPARM+1)= 1
          IFXP(NPARM+2)= 1
          IFXP(NPARM+3)= 1
          NPARM= NPARM+ 3
          NDGF= NTP- NPARM
          IF(IFXVMIN.LE.0) THEN
              IFXP(NPARM)= 0
              NDGF= NDGF-1
              ENDIF
          PV(Nbeta+1)= Re
          PV(Nbeta+2)= De
          PV(Nbeta+3)= VMIN
c ..... On first run, free VMIN and the  \beta_i
          CALL NLLSSRR(NTP,NPARM,MXPARM,JROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
          IF(prFIT.GT.0) THEN
              DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
              IF(APSE.LE.0) THEN
                  IF(RREF.GT.0.d0) THEN
                      WRITE(6,622) NNAME,p,q,AREF,NL,DRMSD,0,PV(1),
     1                PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,Nbeta)
                      ENDIF
                  IF(RREF.LE.0.d0) THEN
                      WRITE(6,624) NNAME,p,q,NL,DRMSD,0,PV(1),PU(1),
     1                      PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,Nbeta)
                      ENDIF
                ELSE
                  IF(RREF.GT.0.d0) WRITE(6,626) NNAME,p,AREF,NS,
     1              NL,DRMSD,(j,ypPSE(j),j,PV(j),PU(j),PS(j),j=1,Nbeta)
                  IF(RREF.LE.0.d0) WRITE(6,628) NNAME,p,NS,NL,
     1                 DRMSD,(j,ypPSE(j),j,PV(j),PU(j),PS(j),j=1,Nbeta)
                ENDIF
              IF(IFXVMIN.LE.0) THEN
                  WRITE(6,660) PV(Nbeta+3),PU(Nbeta+3),PS(Nbeta+3)
                  VMIN= PV(Nbeta+3)
                  ENDIF
              ENDIF
c ... the, if appropriate, set  Re  free too ...
          IF(IFXRe.LE.0) THEN
              NDGF= NDGF- 1 
              IFXP(Nbeta+1)= 0
              CALL NLLSSRR(NTP,NPARM,MXPARM,JROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              Re= PV(Nbeta+1)
              DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
              IF(prFIT.GT.0) THEN
                  IF(APSE.LE.0) THEN
                      IF(RREF.GT.0.d0) THEN
                          WRITE(6,622) NNAME,p,q,AREF,NL,DRMSD,0,
     1          PV(1),PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,Nbeta)
                          ENDIF
                      IF(RREF.LE.0.d0) THEN
                          WRITE(6,624) NNAME,p,q,NL,DRMSD,0,PV(1),
     1                PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,Nbeta)
                          ENDIF
                    ELSE
                      IF(RREF.GT.0.d0) WRITE(6,626) NNAME,p,AREF,NS,NL,
     1                 DRMSD,(j,ypPSE(j),j,PV(j),PU(j),PS(j),j=1,Nbeta)
                      IF(RREF.LE.0.d0) WRITE(6,628) NNAME,p,NS,NL,DRMSD,
     1                       (j,ypPSE(j),j,PV(j),PU(j),PS(j),j=1,Nbeta)
                    ENDIF
                  WRITE(6,662) PV(Nbeta+1),PU(Nbeta+1),PS(Nbeta+1)
                  IF(IFXVMIN.LE.0) WRITE(6,660) PV(Nbeta+3),PU(Nbeta+3),
     1                                                     PS(Nbeta+3)
                  ENDIF
              ENDIF
c ... then with Re fixed again, free De & VMIN (as well as the beta's)
          IF(IFXDe.LE.0) THEN
              DSEB= DSE
              IFXP(Nbeta+1)= 1
              NDGF= NTP - Nbeta
              IF(IFXDe.LE.0) THEN
                  IFXP(Nbeta+2)= 0
                  NDGF= NDGF-1
                  ENDIF
              IF(IFXVMIN.LE.0) THEN
                  IFXP(Nbeta+3)= 0
                  NDGF= NDGF-1
                  ENDIF
              CALL NLLSSRR(NTP,NPARM,MXPARM,JROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              IF(IFXDe.LE.0) De= PV(Nbeta+2)
              IF(IFXVMIN.LE.0) VMIN= PV(Nbeta+3)
              IF((prFIT.GT.0).OR.(DSE.GT.DSEB*1.01)) THEN
                  IF(DSE.GT.DSEB*1.01) WRITE(6,654) DSEB,DSE

                  DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
                  IF(APSE.LE.0) THEN
                      IF(RREF.GT.0.d0) THEN
                          WRITE(6,622) NNAME,p,q,AREF,NL,DRMSD,0,PV(1),
     1               PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,Nbeta)
                          ENDIF
                      IF(RREF.LE.0.d0) THEN
                          WRITE(6,624) NNAME,p,q,NL,DRMSD,0,PV(1),
     1               PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=1,Nbeta)
                          ENDIF
                    ELSE
                      IF(RREF.GT.0.d0) WRITE(6,626) NNAME,p,AREF,NS,NL,
     1                 DRMSD,(j,ypPSE(j),j,PV(j),PU(j),PS(j),j=1,Nbeta)
                      IF(RREF.LE.0.d0) WRITE(6,628) NNAME,p,NS,NL,DRMSD,
     1                       (j,ypPSE(j),j,PV(j),PU(j),PS(j),j=1,Nbeta)
                    ENDIF
                  IF(IFXRe.LE.0)
     1                WRITE(6,662) PV(Nbeta+1),PU(Nbeta+1),PS(Nbeta+1)
                  IF(IFXDe.LE.0)
     1                WRITE(6,630) PV(Nbeta+2),PU(Nbeta+2),PS(Nbeta+2)
                  IF(IFXVMIN.LE.0)
     1                WRITE(6,660) PV(Nbeta+3),PU(Nbeta+3),PS(Nbeta+3)
                  ENDIF
              ENDIF
c ... and finally ... fit to all three of  VMIN, De and Re
          IFXP(Nbeta+1)= IFXRe
          IFXP(Nbeta+2)= IFXDe
          IFXP(Nbeta+3)= IFXVMIN
          PV(Nbeta+1)= Re
          PV(Nbeta+2)= De
          PV(Nbeta+3)= VMIN
          NDGF= NTP- Nbeta
          IF(IFXVMIN.LE.0) NDGF= NDGF-1
          IF(IFXDE.LE.0) NDGF= NDGF-1
          IF(IFXRE.LE.0) NDGF= NDGF-1
          CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
          DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
          IF(APSE.LE.0) THEN
              IF(RREF.GT.0.d0) THEN
                  IF((PSEL.EQ.1).OR.(PSEL.EQ.3)) THEN
                      WRITE(6,618) NNAME,p,AREF,NL,DRMSD,0,PV(1),PU(1),
     1                      PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,Nbeta)
                      WRITE(7,726) NNAME,p,AREF,NL,DRMSD,De,IFXDe,Re,
     1                        IFXRe,APSE,NL,p,q,AREF,(PV(j),j=1,Nbeta)
                      ENDIF
                  IF(PSEL.EQ.2) THEN
                      WRITE(6,622) NNAME,p,q,AREF,NL,DRMSD,0,PV(1),
     1               PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,Nbeta)
                      WRITE(7,722) NNAME,p,q,AREF,NL,DRMSD,De,IFXDe,Re,
     1                        IFXRe,APSE,NL,p,q,AREF,(PV(j),j=1,Nbeta)
                      ENDIF
                  ENDIF
              IF(RREF.LE.0.d0) THEN
                   IF((PSEL.EQ.1).OR.(PSEL.EQ.3)) THEN
                       WRITE(6,617) NNAME,p,NL,DRMSD,0,PV(1),PU(1),PS(1)
     1                          ,DSE,(j-1,PV(j),PU(j),PS(j),j=2,Nbeta)
                       WRITE(7,728) NNAME,p,NL,DRMSD,De,IFXDe,Re,IFXRe,
     1                              APSE,NL,p,q,AREF,(PV(j),j=1,Nbeta)
                       ENDIF
                   IF(PSEL.EQ.2) THEN
                       WRITE(6,624) NNAME,p,q,NL,DRMSD,0,PV(1),PU(1),
     1                     PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,Nbeta)
                       WRITE(7,724) NNAME,p,q,NL,DRMSD,De,IFXDe,Re,
     1                        IFXRe,APSE,NL,p,q,AREF,(PV(j),j=1,Nbeta)
                       ENDIF
                  ENDIF
            ELSE
              IF(RREF.GT.0.d0) WRITE(6,626) NNAME,p,AREF,NS,NL,DRMSD,
     1                       (j,ypPSE(j),j,PV(j),PU(j),PS(j),j=1,Nbeta)
              IF(RREF.LE.0.d0) WRITE(6,628) NNAME,p,NS,NL,DRMSD,
     1                       (j,ypPSE(j),j,PV(j),PU(j),PS(j),j=1,Nbeta)
            ENDIF
          IF(PSEL.EQ.3) BETA(0)= PV(1) 
          WRITE(6,662) PV(Nbeta+1),PU(Nbeta+1),PS(Nbeta+1) 
          WRITE(6,630) PV(Nbeta+2),PU(Nbeta+2),PS(Nbeta+2) 
          WRITE(6,660) PV(Nbeta+3),PU(Nbeta+3),PS(Nbeta+3)
ccc Print [calc.-obs.]
          IF(prDIFF.gt.0) WRITE(6,730) (RTP(I),YD(I),YD(I)/uVTP(I),
     1                                                       I= 1,NTP)
  730 FORMAT(1x,39('==')/3x,3(3x,'RTP',4x,'[c-o] [c-o]/unc')/
     1  1x,39('--')/(1x,3(f10.5,f8.4,f7.2)))
ccc
          IF(IFXRe.LE.0) Re= PV(Nbeta+1) 
          IF(IFXDe.LE.0) De= PV(Nbeta+2)
          IF(IFXVMIN.LE.0) VMIN= PV(Nbeta+3) 
          IF(PSEL.EQ.3) beta(0)= PV(1)
          ENDIF

c=======================================================================
c*** For case of a  GPEF  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.4) THEN
          DO  i= 1,NTP
ccc           uVTP(i)= uVTP(i)
              betay(i)= VTP(i)
              ENDDO
          DO  Nbeta= NS+1,NL+1
c*** Loop over expansion orders from NS to NL
              IFXP(Nbeta+1)= IFXVMIN
              IFXP(Nbeta+2)= IFXRe
              IFXP(Nbeta+3)= 1
c.... first, do fully linearized fit to get trial expansion coefficients
              Rep= Re**p
              NPARM= Nbeta
              IF(IFXVMIN.LE.0) NPARM= NPARM+ 1
cc            IF(IFXRe.LE.0) NPARM= NPARM+ 1
              DO  i= 1, NTP
                  RTPp= RTP(i)**p
                  yp= (RTPp - Rep)/(as*RTPp + bs*Rep)
                  yPOW= yp
                  DO  j=1,Nbeta
                      yPOW= yPOW*yp
                      DYDP(i,j)= yPOW
                      ENDDO
                  IF(IFXVMIN.LE.0) DYDP(i,NPARM)= 1.d0
                  ENDDO
              CALL LLSQF(NTP,NPARM,MXDATA,MXPARM,betay,uVTP,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
              IF(IFXVMIN.LE.0) VMIN= PV(Nbeta+1)
              IF(prFIT.GT.0) THEN
                  WRITE(6,620) NNAME,p,NL,DSE,
     1                         ('  c',j-1,PV(j),PU(j),PS(j),j= 1,Nbeta)
                  IF(IFXVMIN.LE.0) WRITE(6,660) PV(Nbeta+1),PU(Nbeta+1),
     1                                                      PS(Nbeta+1)
                  ENDIF
c.... then, proceed with fit to non-linear form
              NPARM= Nbeta+2 
              PV(Nbeta+2)= Re
              IFXDe= 1
              IF(IFXRe.LE.0) THEN
c... If Re is to be free, first optimize it in fit to initial form
                  CALL NLLSSRR(NTP,NPARM,MXPARM,JROUND,ROBUST,prNLL,
     1                   IFXP,VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
                  IF(prFIT.GT.0) THEN
                      WRITE(6,636)
                      WRITE(6,658) p,Nbeta,Rref,DSE,(i-1,PV(i),PU(i),
     1                                                PS(i),i=1,Nbeta)
                      IF(IFXVMIN.LE.0) WRITE(6,660) PV(Nbeta+1),
     1                                           PU(Nbeta+1),PS(Nbeta+1)
                      WRITE(6,662) PV(Nbeta+2),PU(Nbeta+2),PS(Nbeta+2)
                      ENDIF

                  ENDIF
              IF(Nbeta.GE.2) THEN
                  DO  j=2, Nbeta
                      PV(j)= PV(j+1)/PV(1)
                      ENDDO
                  ENDIF
c ... IFXDe is a flag indicating fit to final  c0*y**2(1 + c1*y + ... )
              IFXDe= 0 
              CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              WRITE(6,658) p,nbeta,DSE,(i-1,PV(i),PU(i),PS(i),
     1                                                      i=1,Nbeta)
              IF(IFXVMIN.LE.0) WRITE(6,660) PV(Nbeta+1),PU(Nbeta+1),
     1                                                     PS(Nbeta+1)
              IF(IFXRE.LE.0) WRITE(6,662) PV(Nbeta+2),PU(Nbeta+2),
     1                                                     PS(Nbeta+2)
              ENDDO
          NPARM= NPARM+ 1
          ENDIF
c======================================================================
c*** Step inward and check whether the repulsive inner potential wall 
c    has an inflection point or turnover.
      IF(APSE.LE.0) THEN
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
      IF((LPPOT.GT.0).AND.(NPR.GT.0)) THEN
c** If desired, print fitted potential at NPR points starting at r=RPR1
          WRITE(8,800)  NNAME, NL, p, q, AREF
          NPARM= Nbeta+3
          VB= 0.d0
          J= NTP+ 1
          DO  I= 1, NPR
              RTP(J)= RPR1 + (I-1)*dRPR 
              CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VV,PV,PU,PS,RR)
              SL= (VV-VB)/dRPR
              dSL= (SL-SLB)/dRPR
              SLB= SL
              VB= VV
              IF(I.EQ.1) WRITE(8,802) RTP(J),VV, SUM
              IF(I.EQ.2) WRITE(8,802) RTP(J),VV, SUM,SL
              IF(I.GT.2) WRITE(8,802) RTP(J),VV, SUM,SL,dSL
              ENDDO
          ENDIF
  800 FORMAT(/2x,A4,' potential for    NL=',i3,'  p=',I2,'  q=',i2,
     1  '   Rref=',f8.5/'    R       V(r)',10x,'beta(r)      dV(r)',
     2 8x,'d2V(r)'/2x,29('--') )
  802 FORMAT(  f9.4,1pd14.6,3d12.4)
      GOTO 20
c-----------------------------------------------------------------------
  616 FORMAT(//' !!! WARNING !!! Should set   p=',i2,' .GT. [MMLR(',i1,
     1 ')=',I2,' - MMLR(1)=',i2,']  CAUTION !!!!'/)
  619 FORMAT(' Linearized fit uses    beta(INF)=',f12.8) 
  620 FORMAT(/' Linearized ',A4,'{p=',i1,'; NL=',i2,'} fit yields   DSE=
     1',1Pd9.2/(3x,a4,'_{',i2,'}=',d19.11,' (+/-',d8.1,')   PS=', 
     2 d8.1))
  621 FORMAT(/' Linearized ',A4,'{p=',i1,', q=',i2,'; Rref=',f5.2,
     1 '; NL=',i2,'} fit yields   DSE=',1Pd9.2/(3x,a4,'_{',i2,'}=',
     2  d19.11,' (+/-',d8.1,')   PS=',d8.1))
  623 FORMAT(/' Linearized ',A4,'{p=',i1,', q=',i2,';  Rref=Re;  NL=',
     1 i2,'} fit yields   DSE=',1Pd9.2/(3x,a4,'_{',i2,'}=',d19.11, 
     2 ' (+/-',d8.1,')   PS=',d8.1))
  617 FORMAT(/' Direct fit to ',A4,'{p=',i1,'; Rref= Re  ; NL=',I2,
     1  '} potential:',10x,'dd=',1Pd12.5/ 
     2 '  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,3x,'DSE=',
     3  D9.2/('  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  618 FORMAT(/' Direct fit to ',A4,'{p=',i1,'; Rref=',f5.2,'; NL=',I2,
     1 '} potential:',10x,'dd=',1Pd12.5/ 
     2 '  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,3x,'DSE=',
     3   D12.5/('  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  622 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',i2,'; Rref=',f5.2,
     1 '; NL=',I2,'} potential:    dd=',1Pd12.5/ 
     2 '  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,3x,'DSE=', 
     3  D12.5/('  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  722 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',i2,'; Rref=',f5.2,
     1 ';  NL=',I2,'} potential:   dd=',1Pd12.5/d20.12,I3,9x,  
     2 '% De IFXDe'/d20.12,I3,9x,'% Re IFXRe'//2I3,2I4,D11.2,7x, 
     3 '% APSE Nbeta nPB nQB RREF'/(d20.12,'  0'))
  726 FORMAT(/' Direct fit to ',A4,'{p=',i1,',  Rref=',f5.2,';  NL=',
     1  I2,'} potential:',8x,'dd=',1Pd12.5/d20.12,I3,9x,'% De IFXDe'/
     1  d20.12,I3,9x,'% Re IFXRe'//2I3,2I4,D11.2,7x,'% APSE Nbeta nPB nQ
     3B RREF'/(d20.12,'  0'))
  624 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',I2,'; Rref= Re  ; NL=
     1',I2,'} potential:    dd=',1Pd12.5/ 
     2 '  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,3x'DSE=', 
     3  D9.2/('  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  724 FORMAT(/' Direct fit to ',A4,'{p=',i1,', q=',I2,'; Rref= Re  ;  NL
     1=',I2,'} potential:   dd=',1Pd12.5/d20.12,I3,9x,'% De IFXDe'/ 
     2  d20.12,I3,9x,'% Re IFXRe'//2I3,2I4,D11.2,7x,'% APSE Nbeta nPB nQ
     3B RREF'/(d20.12,'  0'))
  728 FORMAT(/' Direct fit to ',A4,'{p=',i1,',  Rref= Re  ;  NL=',I2,
     1  '} potential:',8x,'dd=',1Pd12.5/d20.12,I3,9x,'% De IFXDe'/ 
     2  d20.12,I3,9x,'% Re IFXRe'//2I3,2I4,D11.2,7x,'% APSE Nbeta nPB nQ
     3B RREF'/(d20.12,'  0'))
  626 FORMAT(/' Direct fit to ',A4,'{p=',i1,'; Rref=',f5.2,' ; NS=',i2,
     1  ', NL=',I2,'}  potential:    dd=',1Pd9.2/(' ypPSE{',i2,'}=',
     20PF11.7,'   beta_{',i2,'}=',1Pd18.10,'(+/-',d8.1,')   PS=',d8.1))
  628 FORMAT(/' Direct fit to ',A4,'{p=',i1,'; Rref= Re ; NS=',i2,
     1  ', NL=',I2,'}  potential:   DSE=',1Pd9.2/(' ypPSE{',i2,'}=',
     20PF11.7,'   beta_{',i2,'}=',1Pd18.10,'(+/-',d8.1,')   PS=',d8.1))
  630 FORMAT(8x,'De =',f13.6,' (+/-',f12.6,')      PS=',1pd8.1) 
  632 FORMAT('  Use exponent expansion variable:   y_',I1,'(r)= [r^',I1,
     1 ' -',f7.4,'^',I1,']/[r^',I1,' +',f7.4,'^',I1,']' )
  633 FORMAT(' Use exponent expansion variable:   y_',I1,'(r)= [r^',I1,
     1 ' - Re^',I1,']/[r^',I1,' + Re^',I1,']' )
  634 FORMAT(' Polynomial exponent fx is:  beta(R)= betaINF*y_p + (1-y_p
     1)*Sum{beta_i*[y_q]^i}')
  650 FORMAT(' Use Pashov natural spline exponent based on', i4,'  yp va
     1lues for  y < 0'/41x,'and',i4,'  yp values for  y > 0')
  636 FORMAT(/' First perform full non-linear GPEF fit without taking ou
     1t common factor of c_0')
  644 FORMAT(' Update  beta_0  from',f11.6,'   to',f11.6,'   by',
     1  1Pd9.1,' :   DSE=',1PD8.1)
  646 FORMAT(' !!! CAUTION !!! Iteration to optimize  beta(0)  not conve
     1rged after',i3,' tries')
  648 FORMAT('  Converge on   beta_0=',f11.6,'   Next change=',1Pd9.1)
  654 FORMAT(/' *** PROBLEM *** freeing De makes DSE increase from',
     1  1PD9.2,' to',D9.2)
  658 FORMAT(/' Fit to a GPEF{p=',i1,';  N=',I2,'} potential yields:',
     1  23x,'DSE=',1Pd9.2/ (5x,'c_{',i2,'} =',d18.10,' (+/-',d8.1,')',
     2  4x,'PS=',d8.1))
  660 FORMAT(6x,'VMIN =',f13.5,' (+/-',f12.6,')      PS=',1pd8.1)
  662 FORMAT(8x,'Re =',f13.9,' (+/-',f12.9,')      PS=',1pd8.1) 
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
c internal potential variables, while if  IDAT > 1  use SAVED values
c... [Must ensure that calculations based on the current UPDATED PV(j)]
c------------------------------------------------------------------------
      INTEGER MXDATA, MXPARM, MXMLR 
      PARAMETER (MXDATA=1501, MXPARM=43, MXMLR= 15) 
      INTEGER  i,j,m,IDAT,NPARM,NDATA, IFXP(MXPARM),JFXRe,JFXDe,JFXVMIN 
      REAL*8  YC,PV(NPARM),PD(NPARM),PS(NPARM),DM(MXMLR),DMP(MXMLR),
     1 DMPP(MXMLR),RMSR,RTPp,RTPq,Rep,AREF,AREFp,AREFq,ype,dype,Scalc,
     2 betaINF,yp,yq,yPOW,XP,XPW,DER,TCM,UM,TTMM,DERP,SUM,DSUM,AA,BB,
     3 FCT,FCT2,ULR,ULRe,dULRe,d2ULRe,DDER,T0,T0P,T1,RE3,RE6,RE8,RTP3,
     4 RTP6,RTP8,RDIST,C6adj,C9adj,dAAdRe,dBBdRe ,BETAN,beta0,C1tst
      REAL*8  DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),DEIGRe(1,1)
     1  ,DEIGDe(1,1)
c-----------------------------------------------------------------------
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,sVSR2,IDSTT,KDER,NCMM,p,q,Nbeta,
     1                             APSE,MMLR(MXMLR),IFXCm(MXMLR)
      REAL*8 Re,De,VMIN,RREF,M2,as,bs,rhoAB, SUM,
     1  CmVAL(MXMLR),RTP(MXDATA),SAS(MXDATA,MXPARM)
      COMMON /DATABLK/Re,De,VMIN,RREF,M2,as,bs,rhoAB,CmVAL,RTP,SAS,PSEL,
     1 IFXRe,IFXDe,IFXVMIN,sVSR2,IDSTT,KDER,NCMM,MMLR,p,q,Nbeta,APSE,SUM
c-----------------------------------------------------------------------
      SAVE JFXRe,JFXDe,JFXVMIN, AREF,AREFp,AREFq,Rep,C6adj,C9adj,
     1  betaINF,BETA0, AA,BB,ULRe,dULRe
c=======================================================================
      IF(ABS(IDAT).LE.1) THEN
          JFXRe= IFXP(Nbeta+1) 
          JFXDe= IFXP(Nbeta+2) 
          JFXVMIN= IFXP(Nbeta+3)
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
              IF(JFXRe.LE.0) Re= PV(Nbeta+1) 
              IF(JFXDe.LE.0) De= PV(Nbeta+2)
              IF(JFXVMIN.LE.0) VMIN= PV(Nbeta+3) 
              AREF= RREF 
              IF(RREF.LE.0) AREF= Re
              AREFp= AREF**p 
              ENDIF
          RTPp= RDIST**p
          yp= (RTPp - AREFp)/(RTPp + AREFp) 
          yPOW= 1.d0 
          SUM= PV(1) 
          DSUM= 0.d0 
          IF(Nbeta.GE.2) THEN
              DO  j= 2,Nbeta
                  IF(RREF.LE.0.d0) DSUM= DSUM+ (j-1)*PV(j)*yPOW 
                  yPOW= yPOW*yp
                  SUM= SUM+ yPOW*PV(j) 
                  ENDDO
              ENDIF
          XP= DEXP(-SUM*(RDIST- Re)) 
          YC= De*(1.d0 - XP)**2 + VMIN 
          DER= 2.d0*De*(1.d0- XP)*XP 
          DERP= DER*(RDIST- Re) 
          DO  j= 1, Nbeta
              PD(j)= DERP
              DERP= DERP*yp
              ENDDO
c** If appropriate, also get partial derivative w.r.t. Re
          IF(JFXRE.LE.0) THEN
              IF(RREF.LE.0.d0) SUM= SUM +(RDIST- Re)*DSUM*0.5d0*
     1                                             (p/Re)*(1.d0-yp**2)
              PD(Nbeta+1)= -DER*SUM
              ENDIF
c** If appropriate, also get partial derivative w.r.t. De & VMIN
          IF(JFXDe.LE.0) PD(Nbeta+2)= (1.d0- XP)**2
          IF(JFXVMIN.LE.0) PD(Nbeta+3)= 1.d0
          ENDIF

c=======================================================================
c  For the case of an  MLR_{p}  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.2) THEN
          IF(ABS(IDAT).LE.1) THEN
              IF(JFXRe.LE.0) Re= PV(Nbeta+1) 
              IF(JFXDe.LE.0) De= PV(Nbeta+2) 
              IF(JFXVMIN.LE.0) VMIN= PV(Nbeta+3) 
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
                  IF(rhoAB.GT.0.d0) CALL dampF(Re,rhoAB,NCMM,MMLR,sVSR2,
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
                  BETAN= -0.5D0*BETAN*(-1)**Nbeta 
                  T0= 1.d0 
                  DO  j= Nbeta-1, 1, -1
                      BETAN= BETAN + T0*PV(j) 
                      T0= -T0 
                      ENDDO
                  PV(Nbeta)= BETAN
                  IFXP(Nbeta)= 1 
                  WRITE(6,666) C1tst/116140.97d0, Nbeta,PV(Nbeta)
  666 FORMAT(/'  Test  C1(sr)=',1Pd10.3,'    PV(',i3,')=',1PD16.8/)
                  ENDIF
              ENDIF
          RTPp= RDIST**p 
          RTPq= RDIST**q 
          yp= (RTPp - AREFp)/(RTPp + AREFp)
          yq= (RTPq - AREFq)/(RTPq + AREFq)
          ype= (RTPp - Rep)/(RTPp + Rep)
          IF(APSE.GT.0) THEN
c*** Case of Pashov natural spline exponent ....  
c... Now, use a spline through the exponent values defined by the input 
c    points to generate values of that exponent at the desired 
c    spline-definition points
              XP= 0.d0
              DO  J= 1, Nbeta
                  PD(J)= SAS(IDAT,J)
                  XP= XP + PV(J)*PD(J)
                  ENDDO
            ELSE
c... For conventional case of a constrained polynomial exponent function
              yPOW= 1.d0 - yp
              SUM= PV(1)*yPOW
              DSUM= 0.d0
              IF(Nbeta.GE.2) THEN
                  DO  j= 2,Nbeta
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
              IF(rhoAB.GT.0.d0) CALL dampF(RDIST,rhoAB,NCMM,MMLR,sVSR2,
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
          IF(APSE.GT.0) THEN
c... finalize derivative w.r.t. exponent beta-function spline points ...
              DO  J= 1,Nbeta
                  PD(J)= PD(J)*DER*ype
                  ENDDO
            ELSE
c... finalize derivative w.r.t. exponent polynomial coefficient ....
              DO  j= 1,Nbeta
                  PD(j)= yPOW
                  yPOW= yPOW*yq
                  ENDDO
            ENDIF
c** If appropriate, also get partial derivative w.r.t. De & VMIN
          IF(JFXDe.LE.0) PD(Nbeta+2)= (1.d0- XPW)**2 + DER*ype*yp/De
          IF(JFXVMIN.LE.0) PD(Nbeta+3)= 1.d0
c** If appropriate, also get partial derivative w.r.t. Re
          IF(JFXRe.LE.0) THEN
              dype= -0.5d0*(p/RE)*(1.d0 - yp**2)
              IF(RREF.LE.0.d0) THEN
                  DSUM= betaINF - SUM/(1.d0-yp) + DSUM
                ELSE
                  DSUM= 0.d0
                ENDIF
              PD(Nbeta+1)= DER*(dype*(XP + ype*DSUM)
     1                               + (1.d0 - ype*yp)*dULRe/ULRe )
              ENDIF
          ENDIF

c=======================================================================
c** For the case of a  DELR_p  potential
c-----------------------------------------------------------------------
      IF(PSEL.EQ.3) THEN
          IF(ABS(IDAT).LE.1) THEN
              IF(JFXRe.LE.0) Re= PV(Nbeta+1)
              IF(JFXDe.LE.0) De= PV(Nbeta+2)
              IF(JFXVMIN.LE.0) VMIN= PV(Nbeta+3)
              AREF= RREF
              IF(RREF.LE.0) AREF= Re
              AREFp= AREF**p
              ULRe= 0.d0
              dULRe= 0.d0
              d2ULRe= 0.d0
c-----------------------------------------------------------------------
c** Evaluate uLR & its first 2 deriv. at  Re ... 
              KDER=2
              IF(rhoAB.GT.0.d0) CALL dampF(Re,rhoAB,NCMM,MMLR,sVSR2,
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
          CALL dampF(RDIST,rhoAB,NCMM,MMLR,sVSR2,IDSTT,KDER,DM,DMP,DMPP)
          ULR= 0.d0
          DO  m= 1,NCMM
              ULR= ULR + DM(m)*CmVAL(m)/RDIST**MMLR(m)
              ENDDO
          yPOW= 1.d0
          SUM= PV(1)
          DSUM= 0.d0
          IF(Nbeta.GE.2) THEN
              DO  j= 2,Nbeta
                  IF(RREF.LE.0.d0) DSUM= DSUM+ (j-1)*PV(j)*yPOW
                  yPOW= yPOW*yp
                  SUM= SUM+ yPOW*PV(j)
                  ENDDO
              ENDIF
          XP= DEXP(-SUM*(RDIST- Re))
          YC= (AA*XP - BB)*XP + De - ULR + VMIN
          DER= XP*(BB - 2.d0*AA*XP)
          DERP= DER*(RDIST- Re)
          DO  j= 1,Nbeta
              PD(j)= DERP
              DERP= DERP*yp
              ENDDO
          PD(1)= PD(1) + XP*(XP-1.d0)*dULRe/PV(1)**2
c** If appropriate, also get partial derivative w.r.t. De & VMIN
          IF(JFXDe.LE.0) PD(Nbeta+2)= (XP - 2.d0)*XP + 1.d0
          IF(JFXVMIN.LE.0) PD(Nbeta+3)= 1.d0
          IF(JFXRe.LE.0) THEN
c** If appropriate, also get partial derivative w.r.t. Re
              IF(RREF.GT.0.d0) THEN
                  DSUM= 0.d0
                ELSE
cc                DSUM= - DSUM*0.5d0*(p/Re)*(1.d0-yp**2)*(RDIST - Re)
                  DSUM= - DSUM*0.5d0*(p/Re)*(1.d0-yp**2)
                ENDIF
              PD(Nbeta+1)= DER*(DSUM - SUM) + XP*(dAAdRe*XP- dBBdRe)
              ENDIF
          ENDIF
c=======================================================================

c=======================================================================
c** For the case of a  GPEF(p,as,bs)  potential
c-----------------------------------------------------------------------
      IF(PSEL.EQ.4) THEN
          IF(ABS(IDAT).LE.1) THEN
              JFXVMIN= IFXP(Nbeta+1)
              IF(JFXVMIN.LE.0) VMIN= PV(Nbeta+1)
              JFXRe= IFXP(Nbeta+2)
              IF(JFXRe.LE.0) Re= PV(Nbeta+2)
              Rep= Re**p
              ENDIF
          RTPp= RDIST**p
          yp= (RTPp - Rep)/(as*RTPp + bs*Rep)
          yPOW= yp**2
          SUM= PV(1)*yPOW
          DSUM= 2.d0*PV(1)*yp
          IF(IFXDe.LE.0) yPOW= PV(1)*yPOW
          PD(1)= yp**2
          IF(Nbeta.GT.1) THEN
              DO  j= 2,Nbeta
                  DSUM= DSUM + (j+1)*PV(j)*yPOW
                  yPOW= yPOW* yp
                  PD(j)= yPOW
                  SUM= SUM+ PV(j)*yPOW
                  ENDDO
              ENDIF
          YC= SUM + VMIN
          IF(JFXVMIN.LE.0) PD(Nbeta+1)= 1.d0
          IF(JFXRe.LE.0) THEN
              PD(Nbeta+2)= -DSUM* (p/Re)*Rep*RTPp*(as+bs)/
     1                                            (as*RTPp +bs*Rep)**2
              ENDIF 
          ENDIF
c=======================================================================
c%%%%%%%  Optional printout for testing
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
      SUBROUTINE dampF(r,rhoAB,NCMM,MMLR,IDF,IDSTT,KDER,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. R of the chosen form of the damping
c    function, for  m= 1 to MMAX.
c---------------------- RJL Version of 19 January 2013 -----------------
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
     2  bpm(20,-4:0), cpm(20,-4:0),ZK
c------------------------------------------------------------------------
c  The following values for the numerical factors used in both TT and DS
c  were  normalized to the Hydrogen data presented
c  by Kreek and Meath in J.Chem.Phys. 50, 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c------------------------------------------------------------------------
       DATA bTT/2.10d0,2.44d0,2.78d0,3.13d0,3.47d0/
       DATA bDS/2.50d0,2.90d0,3.30d0,3.69d0,3.95d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0/
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
          Lsr= IDF/2
          IF((IDF.LT.-4).OR.(IDF.GT.4).OR.((2*LSR).NE.IDF)) THEN
                WRITE(6,600) IDSTT,IDF
                STOP
                ENDIF
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
                  DO  IDFF= -4,0
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

c=======================================================================
      SUBROUTINE AF3X3potRet(RDIST,DELTAE,C3val,C6val,C8val,De,ULR,
     1                              DEIGM1,DEIGM3,DEIGM5,DEIGR,DEIGDe)
c=======================================================================
      REAL*8  H(3,3),DM1(3,3),DM3(3,3),DM5(3,3),DR(3,3),
     1              DDe(3,3),Q(3,3)
      REAL*8  DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),
     1        DEIGDe(1,1), EIGVEC(3,1), RESID(3,1), W(3) 
      REAL*8  RDIST,RDIST2,RDIST3,DELTAE,C3val,C6val,C8val,De,ULR,
     1   RET,RETSig,RETPi,Modulus,M1,M3,M5,Z
      INTEGER          I,J,L,K
      M1= C3val
      M3= C6val
      M5= C8val
      RET= 9.36423830d-4*RDIST
      RETSig= DCOS(RET) + (RET)*DSIN(RET)
      RETPi= RETSig - RET**2 *DCOS(RET)
      RDIST2= RDIST**2
      RDIST3= RDIST*RDIST2
*      WRITE(25,*) 'Variables = "r", "U(r)","U(r)-U(r)^2/(4De)" ' 
*      WRITE(25,*) 'zone T = "U(r)"'
c  Initialize interaction matrix to 0.d0
      DO  I= 1,3
          H(I,I)=0.0D0
          ENDDO
ccccc Prepare interation matrix  H 
      H(1,1)= -(M1*RETSig+ M3/(RDIST3)+M5/(RDIST3*RDIST2))/(3.d0*RDIST3)
      H(1,2)= -(DSQRT(2.D0))*H(1,1)
      H(2,1)= H(1,2)
      H(1,3)= M1*RETPi/(DSQRT(6.D0)*RDIST3)
      H(3,1)= H(1,3)
      H(2,2)= 2*H(1,1) + DELTAE
      H(2,3)= H(1,3)/DSQRT(2.d0)
      H(3,2)= H(2,3)
      H(3,3)= DELTAE
cccccc Prepare radial derivative of interaction matrix (? is it needed ?)
      DR(1,1)= (3.d0*M1*RETSig + 6.d0*M3/RDIST3 
     1                  + 8.D0*M5/(RDIST3*RDIST2))/(3.d0*RDIST3*RDIST)
      DR(1,2)= -DSQRT(2.d0)*DR(1,1)
      DR(2,1)= DR(1,2)
      DR(2,2)= 2.d0*DR(1,1)
      DR(1,3)= -3.d0*H(1,3)/RDIST
      DR(3,1)= DR(1,3)
      DR(2,3)= -3.d0*H(2,3)/RDIST
      DR(3,2)= DR(2,3)
      DR(3,3)= 0.d0 
cccccc Partial derivative of interaction matric  H  w.r.t.  C3
      DM1(1,1)= -(RETSig + M1/(2.d0*De*RDIST3))/(3.d0*RDIST3)
      DM1(1,2)= -DSQRT(2.d0)*DM1(1,1)
      DM1(2,1)= DM1(1,2)
      DM1(2,2)= 2.d0*DM1(1,1)
      DM1(1,3)= RETPi/(DSQRT(6.d0)*RDIST3)
      DM1(3,1)= DM1(1,3)
      DM1(2,3)= DM1(1,3)/DSQRT(2.d0)
      DM1(3,2)= DM1(2,3)
      DM1(3,3)= 0.d0
cccccc Partial derivative of interaction matric  H  w.r.t.  C6
      DM3(1,1)= -1.d0/(3.d0*RDIST3**2)
      DM3(1,2)= -SQRT(2.d0)*DM3(1,1)
      DM3(1,3)= 0.D0
      DM3(2,1)= DM3(1,2)
      DM3(2,2)= 2.d0*DM3(1,1)
      DM3(2,3)= 0.D0
      DM3(3,1)= DM3(1,3)
      DM3(3,2)= DM3(2,3)
      DM3(3,3)= 0.D0
cccccc Partial derivative of interaction matric  H  w.r.t.  C8
      DM5(1,1)= DM3(1,1)/(RDIST2)
      DM5(1,2)= DM3(1,2)/(RDIST2)
      DM5(1,3)= 0.D0
      DM5(2,1)= DM3(1,2)
      DM5(2,2)= DM3(2,2)/(RDIST2)
      DM5(2,3)= 0.D0
      DM5(3,1)= DM5(1,3)
      DM5(3,2)= DM5(2,3)
      DM5(3,3)= 0.D0
cccccc Partial derivative of interaction matric  H  w.r.t.  De
      DDe(1,1)= M1**2/(12.D0*(RDIST3*De)**2)
      DDe(1,2)= -SQRT(2.D0)*DDe(1,1)
      DDe(1,3)= 0.D0
      DDe(2,1)= DDe(1,2)
      DDe(2,2)= 2.D0*DDe(1,1)
      DDe(2,3)= 0.d0
      DDe(3,1)= DDe(1,3)
      DDe(3,2)= DDe(2,3)
      DDe(3,3)= 0.D0
cccccc Call subroutine to prepare and invert interaction matrix  H
      CALL ZHEEVJ3(H,Q,W)
      L=1
ccc Nor - identify the lowest eigenvalue of  H  and label it  L
      DO J=2,3
          IF (W(J) .LT. W(L)) THEN
              L=J
              ENDIF
          ENDDO  
      ULR= -W(L)
      DO I=1,3      
          EIGVEC(I,1) = Q(I,L)
          ENDDO  
   30 DEIGM1= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM1,EIGVEC))
   40 DEIGM3= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM3,EIGVEC))
   50 DEIGM5= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM5,EIGVEC))           
   60 DEIGR = -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DR,EIGVEC))
   70 DEIGDe= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DDe,EIGVEC))
c     WRITE(25,600) RDIST ,ULR 
c 600 FORMAT(2D16.7)
c     WRITE(26,601) RDIST , DEIGM1, DEIGR ,DEIGDe
c 601 FORMAT(4D16.7)  
      Modulus = (Z)**2 
      RETURN
      end
**    CONTAINS
*=======================================================================
      SUBROUTINE ZHEEVJ3(H,Q,W)
*=======================================================================
c** Subroutine to setup and invert the matrix  H  and return 
c   eigenvalues W and eigenvector matric  Q
      INTEGER   N, I, X, Y, R
      PARAMETER (N=3)
      REAL*8    H(3,3),Q(3,3), W(3)
      REAL*8    SD,SO,S,T,C,G,B,Z,THRESH
cc    DOUBLE PRECISION FUNCTION SQRABS
* Initialize Q to the identitity matrix
* --- This loop can be omitted if only the eigenvalues are desired ---
      DO  X = 1, N
          Q(X,X) = 1.0D0
          DO  Y = 1, X-1
              Q(X, Y) = 0.0D0
              Q(Y, X) = 0.0D0
              ENDDO
          ENDDO
* Initialize W to diag(A)
      DO  X = 1, N
**        W(X) = DREAL(H(X, X))
          W(X) = H(X, X)
          ENDDO
* Calculate SQR(tr(A))
      SD= 0.0D0
      DO  X = 1, N
          SD= SD + ABS(W(X))
          ENDDO
      SD = SD**2
* Main iteration loop
      DO  I = 1, 50
* Test for convergence
          SO = 0.0D0
          DO  X = 1, N
              DO  Y = X+1, N
**                SO = SO + ABS(DREAL(H(X, Y)))
                  SO = SO + ABS(H(X, Y))
                  ENDDO
              ENDDO
          IF(SO.EQ.0.0D0) RETURN
          IF (I .LT. 4) THEN
              THRESH = 0.2D0 * SO / N**2
            ELSE
              THRESH = 0.0D0
            END IF
* Do sweep
          DO  X= 1, N
              DO  Y= X+1, N
**                G= 100.0D0*(ABS(DREAL(H(X, Y))) )
                  G= 100.0D0*(ABS(H(X, Y))) 
                  IF((I.GT.4).AND.((ABS(W(X))+G).EQ.ABS(W(X)))
     $                         .AND.((ABS(W(Y))+G).EQ.ABS(W(Y)))) THEN
                      H(X, Y)= 0.0D0
**                  ELSEIF(ABS(DREAL(H(X, Y))).GT.THRESH) THEN
                    ELSEIF(ABS(H(X, Y)).GT.THRESH) THEN
* Calculate Jacobi transformation
                      B= W(Y) - W(X)
                      IF((ABS(B)+G).EQ.ABS(B)) THEN
                          T= H(X, Y) / B
                        ELSE
                          IF(B .LE. 0.0D0) THEN
                              T= -2.0D0 * H(X, Y)
     $                       /(SQRT(B**2 + 4.0D0*(H(X, Y))**2) - B)
                            ELSE IF (B .EQ. 0.0D0) THEN
                              T= H(X, Y) * (1.0D0 / ABS(H(X, Y)))
                            ELSE
                              T= 2.0D0 * H(X, Y)
     $                       /(SQRT(B**2 + 4.0D0*(H(X, Y))**2) + B)
                            ENDIF
                        ENDIF
                      C= 1.0D0 / SQRT( 1.0D0 + (T)**2 )
                      S= T * C
                      Z= T * (H(X, Y))
cc                    Z= DREAL(T * (H(X, Y)))
* Apply Jacobi transformation
                      H(X, Y) = 0.0D0
                      W(X)    = W(X) - Z
                      W(Y)    = W(Y) + Z
                      DO  R = 1, X-1
                          T       = H(R, X)
                          H(R, X) = C * T - (S) * H(R, Y)
                          H(R, Y) = S * T + C * H(R, Y)
                          ENDDO
                      DO  R = X+1, Y-1
                          T       = H(X, R)
                          H(X, R) = C * T - S * (H(R, Y))
                          H(R, Y) = S * (T) + C * H(R, Y)
                          ENDDO
                      DO  R = Y+1, N
                          T       = H(X, R)
                          H(X, R) = C * T - S * H(Y, R)
                          H(Y, R) = (S) * T + C * H(Y, R)
                          ENDDO
* eigenvectors
* This loop can be omitted if only the eigenvalues are desired ---
                      DO  R = 1, N
                          T       = Q(R, X)
                          Q(R, X) = C * T - (S) * Q(R, Y)
                          Q(R, Y) = S * T + C * Q(R, Y)
                          ENDDO
                      ENDIF
                  ENDDO
              ENDDO
          ENDDO
      PRINT *, "ZHEEVJ3: No convergence."
**    END SUBROUTINE ZHEEVJ3
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

*=======================================================================
      DOUBLE PRECISION FUNCTION SQRABS(Z)
*=======================================================================
* Calculates the squared absolute value of a complex number Z
* ----------------------------------------------------------------------
*  Parameters ..
      REAL*8 Z
**    SQRABS = DREAL(Z)**2
      SQRABS = (Z)**2
      RETURN
      END
**      END FUNCTION SQRABS
*
**     END SUBROUTINE AF3X3potRet
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      double precision function Scalc(x,m,n,y,rKL,LMAX)
c** At the position 'x', Scalc is returned as the value of the m'th 
c  of the 'n' Sm(x) function defining a natural cubic spline through the
c  mesh points located at  x= y(x_i), for i=1,n.  LMAX specifies the 
c  maximum number of mesh points x= y(x_i) allowed by the calling program
c---------------------------------------------------------------------
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
c... Asen's original coding ...
cc       Scalc=ndirac(k,m)*A(x,y1,y2)+ndirac(k-1,m)*B(x,y1,y2)+
cc   +   C(x,y1,y2)*rKL(m,k)+D(x,y1,y2)*rKL(m,k-1)
cc       else
cc       Scalc=ndirac(k,m)*A(x,y1,y2)+ndirac(k-1,m)*B(x,y1,y2)
cc     A=(x1-z)/(x1-x2)
cc     B=(z-x2)/(x1-x2)
cc     C=((x1-z)*(((x1-z)/(x1-x2))**2-1)/6)*(x1-x2)
cc     D=((z-x2)*(((z-x2)/(x1-x2))**2-1)/6)*(x1-x2)
c... Asen's original coding ...
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      double precision function Sprime(x,m,n,y,rKL,LMAX)
c** At the position 'x', evaluate the derivative w.r.t. x of the m'th 
c  Sm(x) function contributing the definition of the the natural cubic
c  spline defined by function values at the  n  points  y(i) [i=1,n]
      INTEGER i,k,kk,m,n,LMAX
      REAL*8 x,del,y1,y2,y(1:LMAX),rKL(1:LMAX,1:LMAX)
      k=0
      kk=0
      do i=2,n
          if((x.gt.y(i-1)).and.(x.le.y(i)))  k=i
          enddo
      if(x.lt.y(1)) then
          k=2
          kk=1
          end if
      if (x.gt.y(n)) then
          k=n
          kk=1
          end if
      if (x.eq.y(1)) k=2
      y1=y(k-1)
      y2=y(k)
      del=y1-y2
      Sprime= 0.d0
      if(kk.eq.0) Sprime= (del-3.d0*(y1-x)**2/del)*rKL(m,k)/6.d0 +
     1                        (3.d0*(x-y2)**2/del-del)*rKL(m,k-1)/6.d0
      IF(k-1.eq.m) Sprime= Sprime + 1.d0/del 
      IF(k.eq.m) Sprime= Sprime - 1.d0/del 
ccc     if(kk.eq.0) then
ccc         Sprim=ndirac(k-1,m)/del-ndirac(k,m)/del+
ccc  +                    (del-3*(y1-x)**2/del)*rKL(m,k)/6+
ccc  +                    (3*(x-y2)**2/del-del)*rKL(m,k-1)/6
ccc       else
ccc         Sprim=ndirac(k-1,m)/del-ndirac(k,m)/del
ccc       end if
      end

c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      subroutine Lkoef(n,x,A,LMAX)   
c** Call this subroutine with list of the 'n' spline x_i values in array 
c   'x' with maximum dimension 'LMAX' and it will return the LMAX x LMAX
c   array of 'rKL' coefficients used for generating the 'n' S_n(x) 
c   spline coefficient functions
c-----------------------------------------------------------------------
c   
c*** Based on nespl subroutine          
      INTEGER LMAX
      INTEGER I,J,N,INDX(1:LMAX)
      REAL*8 X(1:LMAX),A(1:LMAX,1:LMAX),B(1:LMAX,1:LMAX), d
c
      DO  i= 1,LMAX
          DO  j= 1,LMAX
              A(i,j)= 0.d0
              B(i,j)= 0.d0
              ENDDO
          ENDDO
      A(1,1)= (x(3)-x(1))/3.d0
      A(1,2)= (x(3)-x(2))/6.d0
      do i= 2,n-3
          A(i,i-1)= (x(i+1)-x(i))/6.d0
          A(i,i)= (x(i+2)-x(i))/3.d0
          A(i,i+1)= (x(i+2)-x(i+1))/6.d0
          end do
      A(n-2,n-3)= (x(n-1)-x(n-2))/6.d0
      A(n-2,n-2)= (x(n)-x(n-2))/3.d0  
      do i= 1,n-2
          B(i,i)= 1.d0/(x(i+1)-x(i))
          B(i,i+1)= -1.d0/(x(i+2)-x(i+1))-1.d0/(x(i+1)-x(i))
          B(i,i+2)= 1.d0/(x(i+2)-x(i+1))
          end do  
      call ludcmp(A,n-2,LMAX,indx,d)
      do i= 1,n 
          call lubksb(A,n-2,LMAX,indx,B(1,i))
          end do 
      do i= 1,n-2
          do j= 1,n
              A(j,i+1)= B(i,j)
              end do
          end do 
      do i= 1,n
          A(i,1)= 0.0d0
          A(i,n)= 0.0d0
          end do
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      double precision d,a(np,np),TINY
      PARAMETER (NMAX= 500,TINY= 1.0e-20)
      INTEGER i,imax,j,k
      double precision aamax,dum,sum,vv(NMAX)
      d= 1.d0
      do  i= 1,n
          aamax= 0.d0
          do  j= 1,n
              if (abs(a(i,j)).gt.aamax) aamax= abs(a(i,j))
              enddo
          if (aamax.eq.0.) pause 'singular matrix in ludcmp'
          vv(i)= 1.d0/aamax
          enddo
      do  j= 1,n
          do  i= 1,j-1
              sum= a(i,j)
              do  k= 1,i-1
                  sum= sum-a(i,k)*a(k,j)
                  enddo
              a(i,j)= sum
              enddo
          aamax= 0.d0
          do  i= j,n
              sum= a(i,j)
              do  k= 1,j-1
                  sum= sum-a(i,k)*a(k,j)
                  enddo
              a(i,j)= sum
              dum= vv(i)*abs(sum)
              if (dum.ge.aamax) then
                  imax= i
                  aamax= dum
                  endif
              enddo
          if(j.ne.imax)then
              do  k= 1,n
                  dum= a(imax,k)
                  a(imax,k)= a(j,k)
                  a(j,k)= dum
                  enddo
              d= -d
              vv(imax)= vv(j)
              endif
          indx(j)= imax
          if(a(j,j).eq.0.)a(j,j)= TINY
              if(j.ne.n)then
                  dum= 1.d0/a(j,j)
                  do  i= j+1,n
                      a(i,j)= a(i,j)*dum
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
      ii= 0
      do  i= 1,n
          ll= indx(i)
          sum= b(ll)
          b(ll)= b(i)
          if (ii.ne.0)then
              do  j= ii,i-1
                  sum= sum-a(i,j)*b(j)
                  enddo
            else if (sum.ne.0.) then
              ii= i
            endif
          b(i)= sum
          enddo
      do  i= n,1,-1
          sum= b(i)
          do  j= i+1,n
              sum= sum-a(i,j)*b(j)
              enddo
          b(i)= sum/a(i,i)
          enddo
      return
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
c             COPYRIGHT 1998-2013  by  Robert J. Le Roy                +
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
c      which are held fixed [JFXP(j)=1], those PC(j) values are ignored.
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
c                COPYRIGHT 1998  by  Robert J. Le Roy                  +
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
c
c***********************************************************************
c     SUBROUTINE DYIDPJ(I,NDATA,NPTOT,IFXP,UU,PV,PD,PS,RMSR)
c** Illustrative dummy version of DYIDPJ for the case of a fit to a
c  power series of order (NPTOT-1) in X(i). ***  For datum number-i, 
c  calculate and return  PD(j)=[partial derivatives of datum-i] w.r.t. 
c  each of the free polynomial coefficients varied in the fit 
c  (for j=1 to NPTOT).  **  Elements of the integer array IFXP indicate
c  whether parameter j is being held fixed [IFXP(j) > 0] or varied in
c  the fit [IFXP(j).le.0].  If the former, the partial derivative 
c  for parameter j should be  PD(j)= 0.0. 
c* NOTE that  NDATA, PS and RMSR are useful for cases in which
c  derivatives-by-differences are generated (as for BCONT).
c=====================================================================
c** Use COMMON block(s) to bring in values of the independent variable 
c  [here XX(i)] and any other parameters or variables needeed to
c  calculate YC and the partial derivatives. 
c=====================================================================
c     INTEGER  I,J,NDATA,NPTOT,MXDATA,IFXP(NPTOT)
c     PARAMETER  (MXDATA= 501)
c     REAL*8  RMSR,YC,PV(NPTOT),PD(NPTOT),PS(NPTOT),POWER,XX(MXDATA)
c     COMMON /DATABLK/XX
c=====================================================================
c** NOTE BENE(!!) for non-linear fits, need to be sure that the
c  calculations of YC and PD(j) are based on the current UPDATED PV(j)
c  values.  If other (than PV) parameter labels are used internally
c  in the calculations, UPDATE them whenever (say)  I = 1 .
c=====================================================================
c     POWER= 1.D0
c     YC= PV(1)
c     PD(1)= POWER
c     DO 10 J= 2,NPTOT
c         POWER= POWER*XX(I)
c         YC= YC+ PV(J)*POWER
c         PD(J)= POWER
c  10     CONTINUE
c     RETURN
c     END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

