c========================================================================
c|..optional...Register...optional...Register...optional...Register...|c
c|--------------------------------------------------------------------|c
c| You have choosen to download the attached source code for our      |c
c| Fortran program betaFIT. I would appreciate it if you would please |c
c| go to the www address                                              |c
c|            http://scienide2.uwaterloo.ca/~rleroy/betaFIT16/   and  |c
c| fill in the registration form there if you wish to be accessible   |c
c| so that I can send you possible future updates and/or corrections  |c
c| for this code.  This address list will be held securely by me and  |c
c| used for no other purpose................. Robert J. Le Roy .......|c
c|..Register...optional....Register...optional...Register...optional..|c
c========================================================================


c***********************************************************************
      Program  betaFIT16 
c***********************************************************************
c+++++ COPYRIGHT 2016 by  R.J. Le Roy and Asen Pashov ++++++++++++++++++
c++ as per our paper in the 'Remote Sensing' issue of JQSRT (2016) +++++
c+++++++++++++++ version of  25 March 2016 +++++++++++++++++++++++++++++
c* Program to fit NTP read-in potential fx. values {RTP(i),VTP(i)} 
c  (with or without individual weights) to a chosen analytic form.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** This 'Block' Data Utility routine with the name 'arrsizes.h', that
c   governs array dimensioning in program betaFIT  MUST be installed
c   under this name in the same (sub)directory containing the FORTRAN
c   file(s) for this Program when it is being compiled, or incorporated
c   into the program wherever dimensioning is required (as it is below)
c-----------------------------------------------------------------------
      INTEGER MXDATA, LMAX, MXPARM, NCMMAX
      PARAMETER (MXDATA= 1001, LMAX= MXDATA, MXPARM=50, NCMMAX= 25)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=======================================================================
      INTEGER i,j,k,INFL,ITER,IROUND,JROUND,LPPOT,ROBUST,prFIT,prNLL,
     1 prDIFF,M,NPARM,NPR,NTP,NLIN,NS,NL,NDGF,NxPSE,CYCMAX,IFXP(MXPARM)
      REAL*8 PV(MXPARM),PU(MXPARM),PS(MXPARM),
     1 CM(MXPARM,MXPARM),DYDP(MXDATA,MXPARM),VTP(MXDATA),ypeq(MXDATA),
     2 uVTP(MXDATA),betay(MXDATA),Ubetay(MXDATA),YD(MXDATA),
     3 yqPSE(MXDATA),xPSE(MXDATA),rPSE(MXDATA),rKL(1:MXDATA,1:MXDATA),
     4 DM(NCMMAX),DMP(NCMMAX),DMPP(NCMMAX),ypRegrid(MXDATA),
     3 betaINF,UNC,yPOW,DSE,TSTPS,TSTPU,DSEB, TCM,UM,diff,
     4 Rep,Req,AREF,AREFp,AREFq,RTPp,RTPq, ULR,FCT,RAT,UMAX,
     5 XX,YY,YH,yp,yq,ypRE,ReDE, ReIN,DeIN,VMINin,ULRe,dULRedRe,
     6 dULRdR,RE3,RE6,RE8,T0,T1,C6adj,C9adj,RTP3,RTP6,RTP8,RH,RR,RB,RBB,
     7 VV,VB,VBB,SCALC,yMIN,DRMSD ,BETA0,BETAN, RPR1,dRPR, SL,SLB,dSL,
     8 ycalc,dd,VLIM,yqREFre
      REAL*8  DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),DEIGRe(1,1)
     1  ,DEIGDe(1,1)
      CHARACTER*4  NNAME,NAME(4)
      DATA NAME/' EMO',' MLR','DELR','GPEF'/,CYCMAX/25/
c-----------------------------------------------------------------------
      REAL*8 Re,De,VMIN,RREF,AA,BB,as,bs,rhoAB,BETA, CmVAL(NCMMAX),
     1  CmEFF(NCMMAX),dULRdCm(NCMMAX),RTP(MXDATA),SAS(MXDATA,MXPARM)
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,sVSR2,IDSTT,NCMM,MCMM,q,p,
     1  NPbeta,APSE,MMLR(NCMMAX),IFXCm(NCMMAX),updateCmADJ
      COMMON /DATABLK/Re,De,VMIN,RREF,AA,BB,as,bs,rhoAB,CmVAL,CmEFF,
     1 dULRdCm,RTP,SAS,BETA,PSEL,IFXRe,IFXDe,IFXVMIN,sVSR2,IDSTT,NCMM,
     2  MCMM,q,p,NPbeta,APSE,MMLR,IFXCm,updateCmADJ
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
c       If LPPOT > 0, read a separate mesh and range for eash {N,q,p,Rref} case
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
      READ(5,*)  PSEL, NTP, UNC, IROUND, LPPOT, prFIT, prDIFF
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
c** APSE is an integer: > 0  to use A.Pashov natural Spline in MLR Exponent
c                    : .le. 0  for normal constrained polynomial exponent
c*  For an MLR potential with  APSE > 0,  yMIN  is a negative{!} number
c      (-1 .le. yMIN < 0) which specifies the lower bound on the yp values 
c      selected to define the MLR exponent spline.  
c** For each long-range term read power  MMLR(i)  & coefficient CmVAL(i)
c** For special Aubert-Frecon 2x2 cases,  NCMM= 7,  MMLR= {x,3,3,6,6,8,8},
c   with x= 0 for the A state, x= -1 for the b state, and  CmVAL= {Aso, 
c   C3Sig, C3Pi, C6Sig, C6Pi, C8Sig, C8Pi}, 
c* while for the 3x3 diagonalization cases, NCMM=10, MMLR= {x,3,3,3,6,6,6,
c   8,8,8}  with x= -2 for the (lowest eigenvalue) c(1\,^3\Sigma_g^+ state,  
c   x= -3 for the (middle root) B^1\Pi_u state, and x=-4  for the 
c   highest-root state, while CmVal= {Aso, C3Sig, C3Pi1, C3Pi3, C6Sig, C6Pi1,
c    C6Pi3, C8Sig, C8Pi1, C8Pi3}
c=======================================================================
      IF((PSEL.EQ.2).OR.(PSEL.EQ.3)) THEN
          READ(5,*) NCMM, rhoAB, sVSR2, IDSTT, APSE, yMIN
          DO m= 1, NCMM
              READ(5,*) MMLR(m),CmVAL(m)
              CmEFF(m)= CmVAL(m)
              ENDDO
          ENDIF
          MCMM= NCMM
      IF(PSEL.NE.2) APSE= 0
c=======================================================================
c** For a GPEF potential, read coefficients to define expansion vble:
c       y = (r^q - Re^q)/(as*r^q + bs*Re^q)  where  q, as & bs all fixed
c=======================================================================
      IF(PSEL.EQ.4) THEN
          READ(5,*) as, bs 
          WRITE(6,604) as,bs
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
      IF((NCMM.GE.4).AND.((MMLR(1).EQ.0).OR.(MMLR(1).EQ.-1))) THEN
c** Shift of asymptote for special alkali dimer 2x2 cases
          DO  I= 1,NTP                         !! ???? huh??
              VTP(I)= VTP(I) - 0.5d0*CmVAL(1)  !! ???? huh??
              ENDDO                            !! ???? huh??
          ENDIF
c======================================================================= c
  690 FORMAT(' uLR(r)  is a damped TT inverse-power sum with',6x,
     1'   rhoAB=',  F7.4,'   sVSR2=',i3/(49x,'C',i2,'=',1PD15.8:))
  695 FORMAT(' uLR(r) inverse-power terms incorporate ',a2,' damping wit
     1h  rhoAB=',f10.7/8x ,'defined to give very short-range damped uLR-
     2term behaviour  r^{',i2,'/2}')
  696 FORMAT(' uLR(r) inverse-power terms incorporate ',a2,' damping wit
     1h  rhoAB=',f10.7/8x ,'defined to give very short-range damped uLR-
     2term behaviour  r^{',i2,'/2}'/(49x,'C',i2,'=',1PD15.8:))
  698 FORMAT(' uLR(r)  is a simple inverse-power sum with',6x,'C',I2,
     1  '=',1PD15.8:/(49x,'C',i2,'=',D15.8:))
  600 FORMAT(' Fit an ',A4,' potential function to the input points'/
     1  1x,57('=')/5x,'with initial   VMIN=',f13.4,'   Re=',f11.8:
     2 ,'   De=',f11.4)
  614 FORMAT(' Fit an ',A4,'(q=',I2,' p=',I2,')  potential function to t
     1he input points'/1x,57('=')/5x,'with initial   VMIN=',f13.4,
     2 '   Re=',f11.8:,'   De=',f11.4)
  615 FORMAT(' Fit an ',A4,'(q=',I2,') potential function to the input p
     1oints'/1x,57('=')/5x,'with initial   VMIN=',f13.4,'   Re=',f11.8:
     2 ,'   De=',f11.4)
  601 FORMAT(6x,'in which   y_{p/q}(r)= [r^{p/q} -',f7.4,'^{p/q}]/[r^{p
     1/q} +',f7.4,'^{p/q}]' )
  602 FORMAT(5x,'Use Lyon 2x2  ',A7,'  uLR(r)  with   Aso=',F11.6/
     1  47x,'C_3(^1Sig)=',1PD15.7:/47x,'C_3(^3Pi) =',D15.7:/
     1  47x,'C_6(^1Sig)=',1PD15.7:/47x,'C_6(^3Pi) =',D15.7:/
     1  47x,'C_8(^1Sig)=',1PD15.7:/47x,'C_8(^3Pi)  =',D15.7)
 6022 FORMAT(' Use Lyon 3x3 ',A7,'  uLR(r)  with   Aso=',F11.6 /
     1  47x,'C_3(^3Sig)=',D15.7:/47x,'C_3(^1Pi) =',D15.7:/
     2  47x,'C_3(^3Pi) =',D15.7:/
     3  47x,'C_6(^3Sig)=',D15.7:/47x,'C_6(^1Pi) =',D15.7:/
     4  47x,'C_6(^3Pi) =',D15.7:/
     5  47x,'C_8(^3Sig)=',D15.7:/47x,'C_8(^1Pi) =',D15.7:/
     6  47x,'C_8(^3Pi) =',D15.7)
  603 FORMAT(6x, 'in which   y_{p/q}(r)= [r^{p/q} - Re^{p/q}]/[r^{p/q} +
     1 Re^{p/q}]' )
  652 FORMAT( '  & define beta(y_q^{ref}(r)) as a natural spline through
     1 points at the',i4,'  yq^{ref} values:'/(2x,7F11.7))
  604 FORMAT(' GPEF expansion variable is:   y= (R^q - Re^q)//(',F9.5,
     1  '*R^q',1x,SP,F9.5,'*Re^q)' )
  605 FORMAT(' GPEF expansion variable is:   y= (R^',I1,' - Re^',i1,
     1  ')/(',F9.5,'*R^',I1,1x,SP,F9.5,'*Re^',SS,i1,')' )
  606 FORMAT(/' Fit to',I5,' input turning points assuming common uncert
     1ainty  u(VTP)=',1PD9.2/1x,39('--')/
     2  4('    RTP',7x,'VTP ')/1x,39('--')/(4(0pF9.5,f11.3)))
  607 FORMAT('      which yields initial values of   AA=',1PD14.7,
     1   '   BB=',D14.7)
  608 FORMAT(1x,39('--')/)
  609 FORMAT(/' Fit to',I5,' input turning points with initial energy mi
     1nimum   VMIN=',f11.4/1x,32('--')/2(5x,'RTP',8x,'VTP',6x,'unc',4x)/
     2  1x,32('--')/(2(0PF10.5,F12.4,1PD10.2)))
  610 FORMAT(' Shift read-in VMIN value below lowest input potential val
     1ue',f12.4)
  611 FORMAT(3x,'Start with this long-range tail and   beta(0)=',f10.6)
  612 FORMAT(' *** Linearization Fails **  for datum',I4,'  r =',f8.2/
     1  '  De too small or VMIN too high')
c=======================================================================
c** Now ... loop over different {q,p,NS,NL} combinations till end of data
c**  p and q are powers used to define radial variable in the exponent 
c       beta(r)= yp*betaINF + [1 - yp]*sum{beta_i*yq^i}  where
c       ya=(R^a - AREF^a)/(R^a + AREF^a)   for  a= p  or  q 
c** NL= is the order of the\beta(y)  exponent expansion for EMO, DELR, 
c    or MLR potentials with APSE.leq.0.  No. such param is NPbeta= NL+ 1 
c   NS is a dummy variable for EMO, DELR & 'ordinary' MLR with APSE.le.0 
c** For an MLR potential with  APSE > 0,  NPbeta=(NS+NL+1) is the number 
c of yp values to be used to define the exponent spline used for beta(y):
c   NL specifies the number of yp values from  yp= 0  to  yp= 1.0 
c NS specifies the number of (equally spaced) yp values for 
c      yMIN .leq. yp = yMIN to 0.0 
c** For a GPEF potential, consider powers ranging from  NS(.ge.0) to NL 
c* RREF   defines the reference distance in the expansion variable 
c      - for  RREF.le.0 , define parameter  RREF = Re 
c      - for  RREF.gt.0 , fix parameter  RREF   at its read-in value
c-----------------------------------------------------------------------
   20 READ(5,*,END= 999) q, p, NS, NL, RREF
c-----------------------------------------------------------------------
      IF(q.LE.0) GO TO 999
      IF(PSEL.NE.2) p = q
      NPR= 0
c=======================================================================
c** If desired, print the resulting potential for this case to channel-8 
c     at for NPR  distances, starting at  RPR1 with a mesh of DPR 
c*   IF(NPR.leq.0) omit printout for this case
c-----------------------------------------------------------------------
      IF(LPPOT.GT.0) READ(5,*) NPR, RPR1, dRPR
c-----------------------------------------------------------------------
      Re= ReIN
      De= DeIN
      VMIN= VMINin
      IF(PSEL.LE.3) NPbeta= NL + 1
      IF(PSEL.EQ.2) WRITE(6,614) NNAME, q, p, VMIN, Re, De
      IF((PSEL.EQ.1).OR.(PSEL.EQ.3)) WRITE(6,615) NNAME, q, VMIN, Re, De
      IF(PSEL.EQ.4) THEN
          WRITE(6,600) NNAME, VMIN, Re
          WRITE(6,605) q,q,as,q,bs,q
          ENDIF
      IF(PSEL.LE.3) THEN
          IF(RREF.GT.0.d0) THEN
              AREF= RREF
              WRITE(6,632) q,q,RREF,q,q,RREF,q
            ELSE
              AREF= RE
              IF(RREF.LE.0.d0) WRITE(6,633) q,q,q,q,q
            ENDIF
          ENDIF
      AREFp= AREF**p
      AREFq= AREF**q
      IF(PSEL.EQ.2) THEN
          IF(APSE.GT.0) WRITE(6,650) NS, NL
          IF(APSE.LE.0) WRITE(6,634) p, p, q
          IF(RREF.gt.0.d0) THEN
              AREF= RREF
              WRITE(6,601) AREF,AREF
            ELSE
              AREF= Re
              WRITE(6,603)
            ENDIF
          IF(APSE.GT.0) THEN
c** Set up yp abscissa values for APSpline Exponent case
c** first ... FIND    y_q^{ref}(r_e)
              yqREFre= (RE**q - AREFq)/(RE**q + AREFq)
              YH= (yqREFre -yMIN)/NS
              yqPSE(1)= yMIN
c** now .. select NS-1 points between  yqREFre and YMIN
              DO  I= 2,NS
                  yqPSE(I)= yqPSE(I-1) + YH
                  ENDDO
              NPbeta= NS+ 1
              yqPSE(NS+1)= yqREFre - 0.05d0
c** now .. select NL-1 points from  yqREFre to YMAX= 1.0
              YH= (1.d0 - yqREFre)/NL
              DO  I= 1,NL-1
                  yqPSE(NPbeta+I)= yqPSE(NPbeta+I-1) + YH
                  ENDDO
              NPbeta= NPbeta+ NL
              yqPSE(NPbeta)= 1.d0
              WRITE(6,652) NPbeta,(yqPSE(i), i= 1,NPbeta)
              DO i= 1,NPbeta
c*** Round yp values for 'portability' & to avoid precisely Re
                  J= 100* yqPSE(i)
                  yqPSE(i)= FLOAT(J)*1.d-2
                  ENDDO
              WRITE(6,652) NPbeta,(yqPSE(i), i= 1,NPbeta)
              ENDIF
          ENDIF
      IF((NCMM.GE.4).AND.(MMLR(1).LE.0)) THEN
c** uLR printout for Lyon 2x2 or 3x3 treatment of 2S + 2p alkali dimers ...
          IF(MMLR(1).EQ.0) WRITE(6,602) 'A-state',CmVAL(1),
     1                                    CmVAL(2),(CmVAL(m),m=3,NCMM)
c... For Lyon treatment of b-state alkali dimers ...
          IF(MMLR(1).EQ.-1) WRITE(6,602) 'b-state',CmVAL(1),
     1                                    CmVAL(2),(CmVAL(m),m=3,NCMM)
          IF(MMLR(1).EQ.-2) WRITE(6,6022) 'c-state',CmVAL(1) ,CmVAL(2),
     1             (CmVAL(m),m=3,NCMM)
          IF(MMLR(1).EQ.-3) WRITE(6,6022) 'B-state',CmVAL(1) ,CmVAL(2),
     1             (CmVAL(m),m=3,NCMM)
          IF(MMLR(1).EQ.-4) WRITE(6,6022) '1 ^3Pi',CmVAL(1) ,CmVAL(2),
     1             (CmVAL(m),m=3,NCMM)
          IF(IDSTT.GT.0) WRITE(6,695) 'DS',rhoAB,sVSR2
          IF(IDSTT.LE.0) WRITE(6,695) 'TT',rhoAB,sVSR2
          ENDIF
c... uLR printout for 'conventional' (damped or non-damped) inverse-power sum
      IF(((PSEL.EQ.3).OR.((PSEL.EQ.2))).AND.(MMLR(1).GT.0)) THEN
          IF(rhoAB.LE.0) THEN
              WRITE(6,698) (MMLR(m),CmVAL(m),m= 1,NCMM)
            ELSE
              IF(IDSTT.GT.0) WRITE(6,696) 'DS',rhoAB,sVSR2,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
              IF(IDSTT.LE.0) WRITE(6,696) 'TT',rhoAB,sVSR2,
     1                                    (MMLR(m),CmVAL(m),m= 1,NCMM)
            ENDIF
          ENDIF
c=======================================================================
       IF(PSEL.EQ.2) CALL quadCORR(NCMM,MCMM,NCMMAX,MMLR,De,CmVAL,CmEFF)
c=======================================================================
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
      IF(PSEL.EQ.1) THEN
c ... first define ordinate array
          DO  i= 1,NTP
              RTPp= RTP(i)**p
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              IF(RTP(i).GT.Re) THEN
                  T0= (1.d0 - DSQRT((VTP(i)-VMIN)/De))
                      IF(T0.LT.0.d0) THEN        ! This would mean we'd be about to do the ln() of a negative number.
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
              DO  j= 1, NPbeta
                  DYDP(i,j)= yPOW
                  yPOW= yPOW*yp
                  ENDDO
c%%  elective printout for testing
ccc           if(i.eq.1) write(8,700) 
ccc           write(8,702) rtp(i),yp,ypRE,vtp(i),betay(i),
ccc  1                                           (dydp(i,j),j=1,NPbeta)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              ENDDO
          CALL LLSQF(NTP,NPbeta,MXDATA,MXPARM,betay,Ubetay,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
          IF(prFIT.GT.0) WRITE(6,620) NNAME,q,NL,DSE,
     1                        ('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPbeta)
          ENDIF
c=======================================================================
c*** Preliminary linearized fit for an  MLR  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.2) THEN
          IF((MMLR(1).GT.0).AND.(p.LE.(MMLR(NCMM)-MMLR(1))))
     1                        WRITE(6,616) p,NCMM,MMLR(NCMM),1,MMLR(1) 
          IF((MMLR(1).LE.0).AND.(p.LE.(MMLR(NCMM)-MMLR(2))))
     1                        WRITE(6,616) p,NCMM,MMLR(NCMM),2,MMLR(2) 
c** First define array of exponent values with uncertainties defined by 
c  the assumption that all potential values have equal uncertainties UNC
c... Begin by determining uLR(Re) and  betaINF
          IF((NCMM.GE.4).AND.(MMLR(1).LE.0)) THEN
c** For special Aubert-Frecon 2x2 or 3x3 alkali dimer cases ...
              VLIM= VMIN+ De
              CALL AFdiag(Re,VLIM,NCMM,NCMMAX,MMLR,CmEFF,rhoAB,
     1                             sVSR2,IDSTT,ULRe,dULRdCm,dULRedRe)
            ELSE
c... for normal MLR/MLJ case ... with or without damping ---------------
              CALL dampF(Re,rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,DMP,
     1                                                           DMPP)
              ULRe= 0.d0 
              DO  m= 1,NCMM
                  ULRe= ULRe + DM(m)*CmEFF(m)/Re**MMLR(m)
                  ENDDO
            ENDIF
c*** Now, generate betaINF from the adjusted constants =====
          betaINF= DLOG(2.d0*De/ULRe) 
          WRITE(6,619) betaINF
c.... Have Completed determination of uLR(Re) and  betaINF  ============
c.... Now loop to prepare linearized MLR arguments & derivativesa------
          Rep= RE**p
          k= 0
          DO  i= 1, NTP
              RTPp= RTP(i)**p 
              RTPq= RTP(i)**q 
              yp= (RTPp - AREFp)/(RTPp + AREFp) 
              yq= (RTPq - AREFq)/(RTPq + AREFq) 
              ypRE= (RTPp - Rep)/(RTPp + Rep) 
              ypeq(i)= ypRE
cc            IF(DABS(yq).le.0.1d0) GO TO 66     !?? omit point?
              k=k+ 1
              xPSE(k)= yq      !!placed after the k counter
              IF((NCMM.GE.4).AND.(MMLR(1).LE.0)) THEN
c ... for Aubert-Frecon Li2(A) 2x2  or 3x3 cases
                  CALL AFdiag(RTP(i),VLIM,NCMM,NCMMax,MMLR,CmEFF,rhoAB,
     1                                 sVSR2,IDSTT,ULR,dULRdCm,dULRdR)
                ELSE
c... for normal inverse-power-sum MLR/MLJ case ... with or without damping 
                  CALL dampF(RTP(i),rhoAB,NCMM,NCMMAX,MMLR,
     1                                        sVSR2,IDSTT,DM,DMP,DMPP)
                  ULR= 0.d0 
                  DO  m= 1,NCMM
                      ULR= ULR + DM(m)*CmEFF(m)/RTP(i)**MMLR(m)
                      ENDDO
                ENDIF
              RAT= DSQRT((VTP(i)-VMIN)/De)
              IF(RTP(i).GT.Re) THEN
                  T0= (1.d0 - RAT)
                  IF(T0.LT.0.d0) THEN
                      WRITE(6,612) i,RTP(i)
                      STOP
                      ENDIF
                  betay(k)= - DLOG(T0*ULRe/ULR)       !! = beta*yp{eq}
                  Ubetay(k)= 0.5d0*uVTP(i)/(RAT*(1d0 - RAT))
                ELSE
                  betay(k)= - DLOG((1.d0 + RAT)*ULRe/ULR)
                  Ubetay(k)= 0.5d0*uVTP(i)/(RAT*(1.d0+ RAT))
                ENDIF
              IF(APSE.LE.0) THEN
c** Subtract the \beta_\infty term to yield polynomial for fitting 
c... For Huang MLR polynomial exponent function
                  betay(k)= betay(k)- betaINF*yp*ypRE 
                  yPOW= ypRE*(1.d0- yp)
c... then create partial derivative array for linearized fit ...
                  DO  j= 1, NPbeta
                      DYDP(i,j)= yPOW
                      yPOW= yPOW*yq
                      ENDDO
                  ENDIF
c%%%% printout for testing  polynomial exponent %%%%%%%%%%%%%%%%%%%%%%%%
c              if(i.eq.1) write(7,700) 
c 700 format('  RTP    yp       yp^{eq}       VTP(i)          uLR', 
c    1 10x,' beta(i)') 
c              write(7,702) rtp(i),yp,ypRE,vtp(i),ULR,betay(i),
c    1                                           (dydp(i,j),j=1,NPbeta)
c 702 format(f6.3,2f9.5,1P,3d13.5:/(14x,5d13.5))
c%%%% END of printout for testing  polynomial exponent %%%%%%%%%%%%%%%%%
   66         CONTINUE
              ENDDO
          IF(APSE.LE.0) THEN
              CALL LLSQF(NTP,NPbeta,MXDATA,MXPARM,betay,Ubetay,DYDP,YD,
     1                                                PV,PU,PS,CM,DSE)
              IF(prFIT.GT.0) THEN
                  IF(RREF.GE.0.d0) WRITE(6,621) NNAME,q,p,RREF,NL,
     1                    DSE,('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPbeta)
                  IF(RREF.LT.0.d0) WRITE(6,623) NNAME,q,p,NL,
     1                    DSE,('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPbeta)
                  ENDIF
              ENDIF
          IF(APSE.GT.0) THEN
c** For Pashov spline-exponent, use spline through linearized input
c  exponent values to define initial trial spline values of \beta(i)
c  at the specified  xPSE(i) 0:qvalues !   
c... First, add xPSE(NTP+1)= 1.0   where   betay= \betaINF   array point
              NxPSE= k + 1
              xPSE(NxPSE)= 1.d0
              betay(NxPSE)= betaINF
              ubetay(NxPSE)= 1.d0
              rtp(NxPSE)= 9.d99
              ypeq(NxPSE)= 1.d0
c**  Elective printout for linearized SE-MLR stuff *********************
cc            write(8,699) (rtp(i),xPSE(i),betay(i),betay(i)/ypeq(i),
cc   1                                  Ubetay(i),ypeq(i), I= 1,NxPSE)
cc699 FORMAT('  R=', 0p,f10.6,'   xPSE=',F12.9,'   betay=', f12.9,
cc   1   '   beta=',f15.9,'  ubetay=',1pd10.3,'  ypeq=',0p,f12.9)
c****** end of elective printout ***************************************
              CALL Lkoef(NxPSE,xPSE,rKL)
              DO  I= 1, NPbeta- 1 !! skip endpoint where rPSE= \infty
                  XX= yqPSE(I) 
                  rPSE(I)= AREF*((1.d0 + XX)/
     1                                  (1.d0 - XX))**(1.d0/DFLOAT(q))
c... Now, use a spline through the exponent values defined by the input 
c points to generate values of that exponent at the desired 
c spline-definition points
                  ypREgrid(I)= (rPSE(I)**p - Re**p)/(rPSE(I)**p+ Re**p)
                  PV(I)= 0.d0
                  DO  m= 1, NxPSE
                      PV(I)= PV(I) +
     1                      Scalc(XX,m,NxPSE,xPSE,rKL,MXDATA) * betay(m)
                      ENDDO
                  ENDDO
              rPSE(NPbeta)= 9.d99
              ypREgrid(NPbeta)= 1.d0
              PV(NPbeta)= betaINF     !!  now include known endpoint!!
              IF(prFIT.GT.0) WRITE(6,653) NNAME,q,p,
     1           (yqPSE(I),rPSE(I),PV(I),ypREgrid(i),PV(I)/ypREgrid(I),
     1                                                     I= 1,NPbeta)
c*********** Finally, convert  beta*ypRE  to beta
              DO I= 1, NPbeta
                  PV(I)= PV(I)/ypREgrid(i)
                  ENDDO
c c... Finally, create the fixed array of S(m,x) to define the exponent
c and its partial derivatives in subsequent fits, 
              CALL Lkoef(NPbeta,yqPSE,rKL)
c** Now Compare predictions of our preliminary spline function defined by 
c   these NPbeta points with the input data ...
              dd= 0.d0
              DO  I= 1,NxPSE
                  ycalc= 0.d0
                  DO  m= 1,NPbeta
                      SAS(I,m)= Scalc(xPSE(I),m,NPbeta,yqPSE,rKL,MXDATA)
                      ycalc= ycalc+ SAS(I,m)*PV(m)
                      ENDDO
c?? 12/01/15 hey ... maybe need to divide ycalc  by  xPSE
                  ycalc= ycalc*xPSE(I)      !! try this RJL ! 
c??  It appears to work ...  NOW Y_{calc} \approx Y_{input}  
c?? .. but is this  beta  or beta*yp  or ...??
c??  Now need to fix up  ubetay(i) !!!
                  diff= (ycalc - betay(I))/ubetay(I)
                  dd= dd+ diff**2
cc                WRITE(6,655) xPSE(I),ycalc,betay(I)/xPSE(I),diff 
                  ENDDO
              dd= DSQRT(dd/NxPSE) 
              WRITE(6,656) NS, NL, RREF, dd
              ENDIF
          ENDIF
cc655 FORMAT('   At  X=',F8.4,'   Y_{calc}=',F11.7,'   Y_{input}=',
cc   1  F11.7,1x,F13.8)
  656 FORMAT('  SE-MLR Linearization:  NS=',I2,',  NL=',I3,
     1   ',  R_{ref}=', F7.3,'  yields   dd=',F9.5/)
c=======================================================================
c** Preliminary linearized fit for a  DELR  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.3) THEN
          BETA0= 1.d0
          WRITE(6,611) BETA0
          IF((NCMM.GE.4).AND.(MMLR(1).LE.0)) THEN
c** For special Aubert-Frecon 2x2 or 3x3 alkali dimer cases ...
              VLIM= VMIN+ De
              CALL AFdiag(Re,VLIM,NCMM,NCMMAX,MMLR,CmEFF,rhoAB,
     1                             sVSR2,IDSTT,ULRe,dULRdCm,dULRedRe)
            ELSE
              CALL dampF(Re,rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,DMP,
     1                                                           DMPP)
              ULRe= 0.d0
              dULRedRe= 0.d0 
              DO  m= 1,NCMM
                  ULRe= ULRe+ AA*DM(m) 
                  dULRedRe= dULRedRe+ AA*(DMP(m) - MMLR(m)*DM(m)/Re) 
                  ENDDO
            ENDIF
c ... NOTE need iteration to determine self-consistent  beta(0)  value 
c First generate  A & B  from input Re, De and trial beta(0)
   40     CONTINUE
          AA= De - ULRe - dULRedRe/beta0 
          BB= 2.d0*(De - ULRe) - dULRedRe/beta0
          WRITE(6,607) AA,BB 
          RAT= 0.5d0*BB/AA 
          UMAX= DSQRT(RAT**2 + (uVTP(i) + ULRe - DE)/AA) 
          ReDE= Re- dlog(RAT)/beta0 
          DO  i= 1,NTP
              ULR= 0.d0 
              RTPp= RTP(i)**p 
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              IF((NCMM.GE.4).AND.(MMLR(1).LE.0)) THEN
c** For special Aubert-Frecon 2x2 or 3x3 alkali dimer cases ...
                  VLIM= VMIN+ De
                  CALL AFdiag(RTP(i),VLIM,NCMM,NCMMAX,MMLR,CmEFF,rhoAB,
     1                                  sVSR2,IDSTT,ULR,dULRdCm,dULRdR)
                ELSE
                  CALL dampF(RTP(i),rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,
     1                                                    DM,DMP,DMPP)
                  DO  m= 1,NCMM
                      ULR= ULR + DM(m)*CmEFF(m)/RTP(i)**MMLR(m) 
                      ENDDO
                ENDIF
              FCT= (VTP(i) - VMIN + ULR - De)/AA + RAT**2 
              IF(FCT.LT.0.d0) THEN
c** If estimate of Re/DE off a bit and  FCT < 0 , ignore & deweight point
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
              DO  j= 1, NPbeta
                  DYDP(i,j)= yPOW 
                  yPOW= yPOW*yp 
                  ENDDO
              ENDDO
          CALL LLSQF(NTP,NPbeta,MXDATA,MXPARM,betay,Ubetay,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
          IF(DABS(PV(1)-BETA0).GT.PS(1)) THEN
              WRITE(6,644) beta0,PV(1),PV(1)-beta0,DSE 
              beta0= PV(1)
              ITER= ITER+ 1 
              IF(ITER.LE.20) GOTO 40 
              WRITE(6,646) ITER
            ELSE
              WRITE(6,648) PV(1),PV(1)-beta0
            ENDIF
          IF(prFIT.GT.0) WRITE(6,620) NNAME,q,NL,DSE,
     1                        ('beta',j-1,PV(j),PU(j),PS(j),j= 1,NPbeta)
          ENDIF
      IF(PSEL.LE.3)  THEN
c======================================================================
c** Now ... do direct non-linear fits to potential values ... first
c     with Re and/or VMIN and De fixed, and then freeing them up too ...
c=======================================================================
c* FIRST optimize  BETA(j)'s (and VMIN) with  Re and De held fixed!
          DO  j= 1,NPbeta
              IFXP(j)= 0 
              ENDDO
c**?? Is this how I implement fixed limiting  C1  for MLR?
          IF(IDSTT.GT.1) IFXP(NPbeta)= IDSTT
c** Fix parameter value for  yp= 1  in SE-MLR
          IF(APSE.GT.0) IFXP(NPbeta)= 1
c** For initial run, fix the  VMIN, De and Re values
          NPARM= NPbeta+ 3
          IFXP(NPARM-2)= 1
          IFXP(NPARM-1)= 1
          IFXP(NPARM)= 1
          NDGF= NTP- NPARM
          IF(IFXVMIN.LE.0) THEN
              IFXP(NPARM)= 0
              NDGF= NDGF-1
              ENDIF
          PV(NPARM-2)= Re
          PV(NPARM-1)= De
          PV(NPARM)= VMIN
          updateCmADJ= 1
c ..... On first call to NLLSSRR, free only VMIN and the  \beta_i
          CALL NLLSSRR(NTP,NPARM,MXPARM,CYCMAX,JROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
          DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
          VMIN= PV(NPARM)
          IF(prFIT.GT.0) THEN
              IF(APSE.LE.0) THEN
                  IF(RREF.GT.0.d0) THEN
                      WRITE(6,622) NNAME,q,p,AREF,NL,DRMSD,0,PV(1),
     1                PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPbeta)
                      ENDIF
                  IF(RREF.LE.0.d0) THEN
                      WRITE(6,624) NNAME,q,p,NL,DRMSD,0,PV(1),PU(1),
     1                      PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPbeta)
                      ENDIF
                ELSE
                  IF(RREF.GT.0.d0) WRITE(6,626) NNAME,q,AREF,NS,
     1              NL,DRMSD,(j,yqPSE(j),j,PV(j),PU(j),PS(j),j=1,NPbeta)
                  IF(RREF.LE.0.d0) WRITE(6,628) NNAME,q,NS,NL,
     1                 DRMSD,(j,yqPSE(j),j,PV(j),PU(j),PS(j),j=1,NPbeta)
                ENDIF
              IF(IFXVMIN.LE.0) THEN 
                  WRITE(6,660) PV(NPARM),PU(NPARM),PS(NPARM)
                  ENDIF
              ENDIF
c ... then, if appropriate, repeat NLLSSRR call with  Re  free too ...
          IF(IFXRe.LE.0) THEN
              NDGF= NDGF- 1 
              IFXP(NPARM-2)= 0
              CALL NLLSSRR(NTP,NPARM,MXPARM,CYCMAX,JROUND,ROBUST,prNLL,
     1                   IFXP,VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              Re= PV(NPARM-2)
              DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
              IF(prFIT.GT.0) THEN
                  IF(APSE.LE.0) THEN
                      IF(RREF.GT.0.d0) THEN
                          WRITE(6,622) NNAME,q,p,AREF,NL,DRMSD,0,
     1          PV(1),PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPbeta)
                          ENDIF
                      IF(RREF.LE.0.d0) THEN
                          WRITE(6,624) NNAME,q,p,NL,DRMSD,0,PV(1),
     1                PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPbeta)
                          ENDIF
                    ELSE
                      IF(RREF.GT.0.d0) WRITE(6,626) NNAME,q,AREF,NS,NL,
     1                 DRMSD,(j,yqPSE(j),j,PV(j),PU(j),PS(j),j=1,NPbeta)
                      IF(RREF.LE.0.d0) WRITE(6,628) NNAME,q,NS,NL,DRMSD,
     1                       (j,yqPSE(j),j,PV(j),PU(j),PS(j),j=1,NPbeta)
                    ENDIF
                  WRITE(6,662) PV(NPbeta+1),PU(NPbeta+1),PS(NPbeta+1)
                  IF(IFXVMIN.LE.0) WRITE(6,660) PV(NPbeta+3),
     1                                        PU(NPbeta+3),PS(NPbeta+3)
                  ENDIF
              ENDIF
c ... then fix Re again, while freeing De & VMIN (as well as the beta's)
          IF(IFXDe.LE.0) THEN
              updateCmADJ= -30
              DSEB= DSE
              IFXP(NPARM-2)= 1
              NDGF= NTP - NPbeta
              IF(IFXDe.LE.0) THEN
                  IFXP(NPARM-1)= 0
                  NDGF= NDGF-1
                  ENDIF
              IF(IFXVMIN.LE.0) THEN
                  IFXP(NPARM)= 0
                  NDGF= NDGF-1
                  ENDIF
              CALL NLLSSRR(NTP,NPARM,MXPARM,CYCMAX,JROUND,ROBUST,prNLL,
     1                   IFXP,VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              IF(IFXDe.LE.0) De= PV(NPbeta+2)
              IF(IFXVMIN.LE.0) VMIN= PV(NPbeta+3)
              IF((prFIT.GT.0).OR.(DSE.GT.DSEB*1.01)) THEN
                  IF(DSE.GT.DSEB*1.01) WRITE(6,654) DSEB,DSE
                  DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
                  IF(APSE.LE.0) THEN
                      IF(RREF.GT.0.d0) THEN
                          IF((PSEL.EQ.1).OR.(PSEL.EQ.3)) THEN
                           WRITE(6,618) NNAME,q,AREF,NL,DRMSD,0,PV(1),
     1              PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPbeta)
                          IF(PSEL.EQ.3) WRITE(6,625) AA,BB
                              ENDIF
                          IF(PSEL.EQ.2) THEN
                         WRITE(6,622) NNAME,q,p,AREF,NL,DRMSD,0,PV(1),
     1              PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPbeta)
                              ENDIF
                          ENDIF
                      IF(RREF.LE.0.d0) THEN
                          IF((PSEL.EQ.1).OR.(PSEL.EQ.3)) THEN
                          WRITE(6,617) NNAME,q,NL,DRMSD,0,PV(1),PU(1),
     1                    PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPbeta)
                          IF(PSEL.EQ.3) WRITE(6,625) AA,BB
                              ENDIF
                          IF(PSEL.EQ.2) THEN
                              WRITE(6,624) NNAME,q,p,NL,DRMSD,0,PV(1),
     1              PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=1,NPbeta)
                              ENDIF
                          ENDIF
                    ELSE                     !! for Pashov Spline exponent
                      IF(RREF.GT.0.d0) WRITE(6,626) NNAME,q,AREF,NS,NL,
     1               DRMSD,(j,yqPSE(j),j,PV(j),PU(j),PS(j),j=1,NPbeta)
                      IF(RREF.LE.0.d0) WRITE(6,628) NNAME,q,NS,NL,DRMSD,
     1                     (j,yqPSE(j),j,PV(j),PU(j),PS(j),j=1,NPbeta)
                    ENDIF
                  IF(IFXRe.LE.0)
     1                WRITE(6,662) PV(NPARM-2),PU(NPARM-2),PS(NPARM-2)
                  IF(IFXDe.LE.0)
     1                WRITE(6,630) PV(NPARM-1),PU(NPARM-1),PS(NPARM-1)
                  IF(IFXVMIN.LE.0)
     1                      WRITE(6,660) PV(NPARM),PU(NPARM),PS(NPARM)
                  ENDIF
              ENDIF
c ... and finally ... fit to all three of  VMIN, De and Re
          IF(IFXDe.LE.0) updateCmADJ= -30
          IFXP(NPbeta+1)= IFXRe
          IFXP(NPbeta+2)= IFXDe
          IFXP(NPbeta+3)= IFXVMIN
          PV(NPbeta+1)= Re
          PV(NPbeta+2)= De
          PV(NPbeta+3)= VMIN
          NDGF= NTP- NPbeta
          IF(IFXVMIN.LE.0) NDGF= NDGF-1
          IF(IFXDE.LE.0) NDGF= NDGF-1
          IF(IFXRE.LE.0) NDGF= NDGF-1
          CALL NLLSSRR(NTP,NPARM,MXPARM,CYCMAX,IROUND,ROBUST,prNLL,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
          DRMSD= DSE*DSQRT(DFLOAT(NDGF)/DFLOAT(NTP))
          IF(APSE.LE.0) THEN
              IF(RREF.GT.0.d0) THEN
                  IF((PSEL.EQ.1).OR.(PSEL.EQ.3)) THEN
                      WRITE(6,618) NNAME,q,AREF,NL,DRMSD,0,PV(1),PU(1),
     1                      PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPbeta)
                      IF(PSEL.EQ.3) WRITE(6,625) AA,BB
                      WRITE(7,726) NNAME,q,AREF,NL,DRMSD,De,IFXDe,Re,
     1                        IFXRe,APSE,NL,q,p,AREF,(PV(j),j=1,NPbeta)
                      ENDIF
                  IF(PSEL.EQ.2) THEN
                      WRITE(6,622) NNAME,q,p,AREF,NL,DRMSD,0,PV(1),
     1               PU(1),PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPbeta)
                      WRITE(7,722) NNAME,q,p,AREF,NL,DRMSD,De,IFXDe,Re,
     1                        IFXRe,APSE,NL,q,p,AREF,(PV(j),j=1,NPbeta)
                      ENDIF
                  ENDIF
              IF(RREF.LE.0.d0) THEN
                  IF((PSEL.EQ.1).OR.(PSEL.EQ.3)) THEN
                      WRITE(6,617) NNAME,q,NL,DRMSD,0,PV(1),PU(1),PS(1)
     1                          ,DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPbeta)
                      IF(PSEL.EQ.3) WRITE(6,625) AA,BB
                      WRITE(7,728) NNAME,q,NL,DRMSD,De,IFXDe,Re,IFXRe,
     1                              APSE,NL,q,p,-AREF,(PV(j),j=1,NPbeta)
                      ENDIF
                  IF(PSEL.EQ.2) THEN
                      WRITE(6,624) NNAME,q,p,NL,DRMSD,0,PV(1),PU(1),
     1                     PS(1),DSE,(j-1,PV(j),PU(j),PS(j),j=2,NPbeta)
                      WRITE(7,724) NNAME,q,p,NL,DRMSD,De,IFXDe,Re,
     1                        IFXRe,APSE,NL,q,p,-AREF,(PV(j),j=1,NPbeta)
                      ENDIF
                  ENDIF
          ELSE        !! now print for Pashov Spline Exponent case
              IF(RREF.GT.0.d0) WRITE(6,626) NNAME,q,AREF,NS,NL,DRMSD,
     1                       (j,yqPSE(j),j,PV(j),PU(j),PS(j),j=1,NPbeta)
              IF(RREF.LE.0.d0) WRITE(6,628) NNAME,q,NS,NL,DRMSD,
     1                       (j,yqPSE(j),j,PV(j),PU(j),PS(j),j=1,NPbeta)
          ENDIF
          IF(PSEL.EQ.3) BETA0= PV(1) 
          WRITE(6,662) PV(NPbeta+1),PU(NPbeta+1),PS(NPbeta+1) 
          WRITE(6,630) PV(NPbeta+2),PU(NPbeta+2),PS(NPbeta+2) 
          WRITE(6,660) PV(NPbeta+3),PU(NPbeta+3),PS(NPbeta+3)
ccc Print [calc.-obs.]
          IF(prDIFF.gt.0) WRITE(6,730) (RTP(I),VTP(I),YD(I),
     1                                          YD(I)/uVTP(I),I= 1,NTP)
  730 FORMAT(1x,38('==')/3x,2(3x,'RTP',7x,'VTP',5x,'[c-o]  [c-o]/unc')/
     1  1x,38('--')/(1x,2(0P,f10.5,F10.2,f8.4,1P,D10.2)))   
ccc
          IF(IFXRe.LE.0) Re= PV(NPbeta+1) 
          IF(IFXDe.LE.0) De= PV(NPbeta+2)
          IF(IFXVMIN.LE.0) VMIN= PV(NPbeta+3) 
          ENDIF
c=======================================================================
c*** For case of a  GPEF  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.4) THEN
          DO  i= 1,NTP
ccc           uVTP(i)= uVTP(i)
              betay(i)= VTP(i)
              ENDDO
          DO  NPbeta= NS+1,NL+1
c*** Loop over expansion orders from NS to NL
              IFXP(NPbeta+1)= IFXVMIN
              IFXP(NPbeta+2)= IFXRe
              IFXP(NPbeta+3)= 1
c.... first, do fully linearized fit to get trial expansion coefficients
              Req= Re**q
              NPARM= NPbeta
              IF(IFXVMIN.LE.0) NPARM= NPARM+ 1
cc            IF(IFXRe.LE.0) NPARM= NPARM+ 1
              DO  i= 1, NTP
                  RTPq= RTP(i)**q
                  yq= (RTPq - Req)/(as*RTPq + bs*Req)
                  yPOW= yq
                  DO  j=1,NPbeta
                      yPOW= yPOW*yq
                      DYDP(i,j)= yPOW
                      ENDDO
                  IF(IFXVMIN.LE.0) DYDP(i,NPARM)= 1.d0
                  ENDDO
              CALL LLSQF(NTP,NPARM,MXDATA,MXPARM,betay,uVTP,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
              IF(IFXVMIN.LE.0) VMIN= PV(NPbeta+1)
              IF(prFIT.GT.0) THEN
                  WRITE(6,620) NNAME,q,NL,DSE,
     1                         ('  c',j-1,PV(j),PU(j),PS(j),j= 1,NPbeta)
                  IF(IFXVMIN.LE.0) WRITE(6,660) PV(NPbeta+1),
     1                                       PU(NPbeta+1),PS(NPbeta+1)
                  ENDIF
c.... then, proceed with fit to non-linear form
              NPARM= NPbeta+2 
              PV(NPbeta+2)= Re
              IFXDe= 1
              IF(IFXRe.LE.0) THEN
c... If Re is to be free, first optimize it in fit to initial form
                  CALL NLLSSRR(NTP,NPARM,MXPARM,CYCMAX,JROUND,ROBUST,
     1             prNLL,IFXP,VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
                  IF(prFIT.GT.0) THEN
                      WRITE(6,636)
                      WRITE(6,658) q,NPbeta,Rref,DSE,(i-1,PV(i),PU(i),
     1                                                PS(i),i=1,NPbeta)
                      IF(IFXVMIN.LE.0) WRITE(6,660) PV(NPbeta+1),
     1                                       PU(NPbeta+1),PS(NPbeta+1)
                      WRITE(6,662) PV(NPbeta+2),PU(NPbeta+2),
     1                                                    PS(NPbeta+2)
                      ENDIF
                  ENDIF
              IF(NPbeta.GE.2) THEN
                  DO  j=2, NPbeta
                      PV(j)= PV(j+1)/PV(1)
                      ENDDO
                  ENDIF
c ... IFXDe is a flag indicating fit to final  c0*y**2(1 + c1*y + ... )
              IFXDe= 0 
              CALL NLLSSRR(NTP,NPARM,MXPARM,CYCMAX,IROUND,ROBUST,prNLL,
     1                   IFXP,VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              WRITE(6,658) q,NPbeta-1,DSE,(i-1,PV(i),PU(i),PS(i),
     1                                                      i=1,NPbeta)
              IF(IFXVMIN.LE.0) WRITE(6,660) PV(NPbeta+1),PU(NPbeta+1),
     1                                                     PS(NPbeta+1)
              IF(IFXRE.LE.0) WRITE(6,662) PV(NPbeta+2),PU(NPbeta+2),
     1                                                     PS(NPbeta+2)
              WRITE(7,659) 2,PV(1),(i+1,PV(1)*pv(I),i=2,nPbeta)
              ENDDO           !! end of loop over polynomial order NPbeta
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
          CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VV,PV,PU,PS)
          RB= RR + RH
          RTP(J)= RB
          CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VB,PV,PU,PS)
          INFL= 0
          DO  I= 1,110
              RBB= RB
              VBB= VB
              RB= RR
              VB= VV
              RR= RR - RH
              RTP(J)= RR
              IF(RTP(J).LT.0.01d0) EXIT
              IF(APSE.GT.0) THEN
c*** Define SAS value required by dyidpj for this additional point
                  DO  m= 1,NPbeta
c ... gotta generate xPSE(J) for r= RTP(J) to alow SAS value to be generated
                      XPSE(J)= (RTP(J)**q - AREFq)/(RTP(J)**q + AREFq)
                      SAS(J,m)= Scalc(xPSE(J),m,NPbeta,yqPSE,rKL,MXDATA)
                      ENDDO
                  ENDIF
              CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VV,PV,PU,PS)
cc !!! temporary writes for testing.
cc    write(6,666) rr,vv, (vv-vb)/RH, (vv - 2.d0*vb + vbb)/RH**2
cc666 format(f8.4,3f16.4) 
cc  !!! end of temporary writes for testing.
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
          WRITE(8,800)  NNAME, NL, q, p, AREF
          NPARM= NPbeta+3
          VB= 0.d0 
          SLB= 0.d0
          J= NTP+ 1
          DO  I= 1, NPR
              RTP(J)= RPR1 + (I-1)*dRPR 
              IF(APSE.GT.0) THEN
c*** Define SAS values for this additional point
                  DO  m= 1,NPbeta
c ... gotta generate xPSE(J) for r= RTP(J) to alow SAS value to be generated
                      XPSE(J)= (RTP(J)**q - AREFq)/(RTP(J)**q + AREFq)
                      SAS(J,m)= Scalc(xPSE(J),m,NPbeta,yqPSE,rKL,MXDATA)
                      ENDDO
                  ENDIF
              CALL DYIDPJ(J,MXDATA,NPARM,IFXP,VV,PV,PU,PS)
              SL= (VV-VB)/dRPR
              dSL= (SL-SLB)/dRPR
              SLB= SL
              VB= VV
              IF(I.EQ.1) WRITE(8,802) RTP(J),VV
              IF(I.EQ.2) WRITE(8,802) RTP(J),VV, SL
              IF(I.GT.2) WRITE(8,802) RTP(J),VV, SL,dSL
              ENDDO
          ENDIF
  800 FORMAT(/2x,A4,' potential for    NL=',i3,'  q=',I2,'  p=',i2,
     1  '   Rref=',f8.5/'    R       V(r)',10x,'   dV(r)',
     2 8x,'d2V(r)'/2x,29('--') )
  802 FORMAT(  f9.4,1pd14.6,2d12.4)
      GOTO 20
c-----------------------------------------------------------------------
  616 FORMAT(//'!!! WARNING !!! Should set   p=',i2,' > [{MMLR(',i2,
     1 ')=',I2,'} - {MMLR(',i1,')=',i2,'}]  CAUTION !!!'/)
  619 FORMAT(' Linearized fit uses    beta(INF)=',f12.8) 
  620 FORMAT(/' Linearized ',A4,'{q=',i2,'; NL=',i2,'} fit yields   DSE=
     1',1Pd9.2/(3x,a4,'_{',i2,'}=',d19.11,' (+/-',d8.1,')   PS=', 
     2 d8.1))
  621 FORMAT(/' Linearized ',A4,'{q=',i2,', p=',i2,'; Rref=',f5.2,
     1 '; NL=',i2,'} fit yields   DSE=',1Pd9.2/(3x,a4,'_{',i2,'}=',
     2  d19.11,' (+/-',d8.1,')   PS=',d8.1))
  623 FORMAT(/' Linearized ',A4,'{q=',i2,', p=',i2,';  Rref=Re;  NL=',
     1 i2,'} fit yields   DSE=',1Pd9.2/(3x,a4,'_{',i2,'}=',d19.11, 
     2 ' (+/-',d8.1,')   PS=',d8.1))
  617 FORMAT(/' Direct fit to ',A4,'{q=',i2,'; Rref= Re   ; NL=',I2,
     1  '} potential:',8x,'dd=',1Pd12.5/ 
     2 '  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,3x,'DSE=',
     3  D9.2/('  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  618 FORMAT(/' Direct fit to ',A4,'{q=',i2,'; Rref=',f6.3,'; NL=',I2,
     1 '} potential:',8x,'dd=',1Pd12.5/ 
     2 '  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,3x,'DSE=',
     3   D12.5/('  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  622 FORMAT(/' Direct fit to ',A4,'{q=',i2,', p=',i2,'; Rref=',f5.2,
     1 '; NL=',I2,'} potential:   dd=',1Pd12.5/ 
     2 '  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,3x,'DSE=', 
     3  D12.5/('  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  625 FORMAT(8x,'AA=',1Pd20.12,'     BB=',d20.12)
  722 FORMAT(/' Direct fit to ',A4,'{q=',i2,', p=',i2,'; Rref=',f5.2,
     1 ';  NL=',I2,'} potential:  dd=',1Pd12.5/d20.12,I3,9x,  
     2 '% De IFXDe'/d20.12,I3,9x,'% Re IFXRe'//2I3,2I4,D11.2,7x, 
     3 '% APSE Nbeta nQB nPB RREF'/(d20.12,'  0'))
  726 FORMAT(/' Direct fit to ',A4,'{q=',i2,',  Rref=',f5.2,';  NL=',
     1  I2,'} potential:',7x,'dd=',1Pd12.5/d20.12,I3,9x,'% De IFXDe'/
     1  d20.12,I3,9x,'% Re IFXRe'//2I3,2I4,D11.2,7x,'% APSE Nbeta nQB nP
     3B RREF'/(d20.12,'  0'))
  624 FORMAT(/' Direct fit to ',A4,'{q=',i2,', p=',I2,'; Rref= Re  ; NL=
     1',I2,'} potential:   dd=',1Pd12.5/ 
     2 '  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1,3x'DSE=', 
     3  D9.2/('  beta_{',i2,'}=',d20.12,' (+/-',d8.1,')   PS=',d8.1))
  724 FORMAT(/' Direct fit to ',A4,'{q=',i2,', p=',I2,'; Rref= Re  ;  NL
     1=',I2,'} potential:  dd=',1Pd12.5/d20.12,I3,9x,'% De IFXDe'/ 
     2  d20.12,I3,9x,'% Re IFXRe'//2I3,2I4,D11.2,7x,'% APSE Nbeta nQB nP
     3B RREF'/(d20.12,'  0'))
  728 FORMAT(/' Direct fit to ',A4,'{q=',i2,',  Rref= Re  ;  NL=',I2,
     1  '} potential:',7x,'dd=',1Pd12.5/d20.12,I3,9x,'% De IFXDe'/ 
     2  d20.12,I3,9x,'% Re IFXRe'//2I3,2I4,D11.2,7x,'% APSE Nbeta nQB nP
     3B RREF'/(d20.12,'  0'))
  626 FORMAT(/' Direct fit to ',A4,'{q=',i2,'; Rref=',f5.2,' ; NS=',i2,
     1  ', NL=',I2,'}  potential:    dd=',1Pd9.2/(' yqPSE{',i2,'}=',
     20PF11.7,'   beta_{',i2,'}=',1Pd18.10,'(+/-',d8.1,')   PS=',d8.1))
  628 FORMAT(/' Direct fit to ',A4,'{q=',i2,'; Rref= Re ; NS=',i2,
     1  ', NL=',I2,'}  potential:   DSE=',1Pd9.2/(' yqPSE{',i2,'}=',
     20PF11.7,'   beta_{',i2,'}=',1Pd18.10,'(+/-',d8.1,')   PS=',d8.1))
  630 FORMAT(8x,'De =',f13.6,' (+/-',f12.6,')      PS=',1pd8.1) 
  632 FORMAT('  Use exponent expansion variable:   y_',I1,'(r)= [r^',I1,
     1 ' -',f7.4,'^',I1,']/[r^',I1,' +',f7.4,'^',I1,']' )
  633 FORMAT(' Use exponent expansion variable:   y_',I1,'(r)= [r^',I1,
     1 ' - Re^',I1,']/[r^',I1,' + Re^',I1,']' )
  634 FORMAT(' MLR polynomial exponent function is:'/29x,'beta(R)= betaI
     1NF*y_',I1,' + (1-y_',I1,')*Sum{beta_i*[y_',I1,'q]^i}')
  636 FORMAT(/' First perform full non-linear GPEF fit without taking ou
     1t common factor of c_0')
  644 FORMAT(' Update  beta_0  from',f11.6,'   to',f11.6,'   by',
     1  1Pd9.1,' :   DSE=',1PD8.1)
  646 FORMAT(' !!! CAUTION !!! Iteration to optimize  beta(0)  not conve
     1rged after',i3,' tries')
  648 FORMAT('  Converge on   beta_0=',f11.6,'   Next change=',1Pd9.1)
  650 FORMAT(' Use Pashov natural spline exponent based on', i4,' yq^{re
     1f} values for  r < r_e'/41x,'and',i4,' yq^{ref} values for  r > r_
     2e')
  653 FORMAT(/' Linearized ',A4,'{q=',i2,' p=',i2,'}-APSE treatment yiel
     1ds:'/4x,'qpSE',8x,'rPSE',6x,'ypRE'4x,'beta*ypRe   beta'/
     2  (F12.6,4f10.6))
  654 FORMAT(/' *** PROBLEM *** freeing De makes DSE increase from',
     1  1PD9.2,' to',D9.2)
  658 FORMAT(/' Fit to a GPEF{q=',i1,';  N=',I2,'} potential yields:',
     1  23x,'DSE=',1Pd9.2/ (5x,'c_{',i2,'} =',d18.10,' (+/-',d8.1,')',
     2  4x,'PS=',d8.1))
  659 FORMAT(/(5x,'a_{',i2,'} =',d18.10))
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

c=======================================================================
      SUBROUTINE quadCORR(NCMM,MCMM,NCMMAX,MMLR,De,CmVAL,CmEFF)
c=======================================================================
c** subroutine to generate and print MLR CmEFF values incorporating
c  quadratic 'Dattani' corrections to Cm values for both standard 'linear'
c  and A-F diagonalized uLR(r) functions for MLR potentials
c** Return MCMM= NCMM+1  for C9{adj} term for m_1= 3 potentials
c=======================================================================
      INTEGER NCMM,MCMM,NCMMAX,MMLR(NCMMAX)
      REAL*8 De,CmVAL(NCMMAX),CmEFF(NCMMAX)
c-----------------------------------------------------------------------
      IF(MMLR(1).GT.0) THEN
c** For 'normal' inverse-power sum MLR case, with or without damping,
c   set up Dattani's 'Quadratic-corrected' effective Cm values 
          IF((MMLR(1).EQ.6).AND.(NCMM.GE.4)) THEN
c... First, consider C6/C12adj(C14adj) for MMLR(m)={6,8,10,(11),12,14} case
              IF(MMLR(4).EQ.12) THEN             ! explicitly MMLR(4)=12
                  CmEFF(4)= CmVAL(4)+ 0.25D0*CmVAL(1)**2/De
                  WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
                  ENDIF
              IF(NCMM.GE.5) THEN
                 IF(MMLR(4).EQ.11) THEN         ! implicitly MMLR(5)=12
                     CmEFF(5)= CmVAL(5) + 0.25D0*CmVAL(1)**2/De
                     WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                     IF(NCMM.GE.6) THEN         ! implicitly MMLR(6)=14
                         CmEFF(6)= CmVAL(6)+ 0.5D0*CmVAL(1)*CmVAL(2)/De
                         WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                         ENDIF
                     ENDIF
                 IF(MMLR(4).EQ.12) THEN           ! assuming MMLR(5)=14
                     CmEFF(5)= CmVAL(5) + 0.5D0*CmVAL(1)*CmVAL(2)/De
                     WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                     ENDIF
                 ENDIF
              ENDIF
          IF((MMLR(1).EQ.5).AND.(NCMM.GE.4)) THEN
c... Then, consider C5/C10adj + C12adj for MMLR(m)={5,6,8,10,12,14} cases
              CmEFF(4)= CmVAL(4) + 0.25D0*CmVAL(1)**2/De
              WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
              IF(NCMM.GE.5) THEN                 ! introduce C12^{adj}
                  CmEFF(5)= CmVAL(5) + 0.25D0*CmVAL(2)**2/De
                  WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                  IF(NCMM.GE.6) THEN             ! introduce C14^{adj}
                      CmEFF(6)= CmVAL(6) + 0.5D0*CmVAL(2)*CmVAL(3)/De
                      WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                      ENDIF
                  ENDIF
              ENDIF
          IF((MMLR(1).EQ.4).AND.(NCMM.GE.3)) THEN
c... Then, consider C4/C8adj + C12adj for MMLR(m)={4,6,8,10,12,14} cases
              CmEFF(3)= CmVAL(3) + 0.25D0*CmVAL(1)**2/De
              WRITE(6,712) MMLR(3),MMLR(3),CmEFF(3)
              IF(NCMM.GE.4) THEN                 ! implicitly MMLR(4)=10
                  CmEFF(4)= CmVAL(4) + 0.5D0*CmVAL(1)*CmVAL(2)/De
                  WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
                  IF(NCMM.GE.5) THEN             ! implicitly MMLR(5)=12
                      CmEFF(5)= CmVAL(5) + 0.5D0*CmVAL(1)*CmVAL(3)/De
     1                                       + 0.25D0*CmVAL(2)**2/De
                      WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                      IF(NCMM.GE.6) THEN         ! implicitly MMLR(6)=14
                          CmEFF(6)= CmVAL(6)+ 0.5D0*CmVAL(2)*CmVAL(3)/De
     1                                      + 0.5D0*CmVAL(1)*CmVAL(4)/De
                          WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                          ENDIF
                      ENDIF
                  ENDIF
                ENDIF
          IF((MMLR(1).EQ.3).AND.(NCMM.GE.2)) THEN
c... Then, consider C3/C6adj + C9adj for MMLR(m)={3,6,8,(9),10,12,14} cases
              CmEFF(2)= CmVAL(2) + 0.25D0*CmVAL(1)**2/De 
              WRITE(6,712) MMLR(2),MMLR(2),CmEFF(2)
              IF(NCMM.GE.3) THEN              ! introduce C9adj & MMLR=9
                  MCMM= NCMM + 1
                  MMLR(MCMM)= 9 
                  CmEFF(MCMM)= 0.5d0*CmVAL(1)*CmEFF(2)/De
                  WRITE(6,714) MMLR(MCMM),CmEFF(MCMM)
                  IF(NCMM.GE.5) THEN             ! implicitly MMLR(5)=12
                      CmEFF(5)= CmVAL(5) + 0.5D0*CmVAL(1)*CmEFF(MCMM)/De
     1                                         + 0.25D0*CmEFF(2)**2/De
                      WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                      IF(NCMM.GE.6) THEN         ! implicitly MMLR(6)=14
                          CmEFF(6)= CmVAL(6)+ 0.5D0*CmEFF(2)*CmVAL(3)/De
                          WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF
c======================================================================= c
c** End of  CmEFF= Cm + CmADJ  setup for non-AF case ===================
  710 Format("  'Quadratic correction' for   C",I2,'(MLR)   yields',
     1  6x,'C',I2,'{adj}=',1PD15.8)
  712 Format("  'Quadratic correction' for   C",I1,'(MLR)    yields',
     1  7x,'C'I1,'{adj}=',1PD15.8)
  714 Format("  'Quadratic corrn' for  MLR(m_1=3)  introduces   C",
     1    I1,'(',A4,',adj) =',1PD15.8)
  716 Format("  'Quadratic correction' for  C",I1,'(Sigma)  yields   C',
     1    I1,'(Sigma,adj)=',1PD15.8)
  718 Format("  'Quadratic correction' for  C",I1,'(^3Pi)   yields   C',
     1    I1,'(^3Pi,adj) =',1PD15.8)
  720 Format("  'Quadratic correction' for  C",I1,'(^1Pi)   yields  C',
     1    I1,'(^1Pi,adj) =',1PD15.8)
c=========================================================================      
      IF(MMLR(1).LE.0) THEN
c** implement Quadratic 'Dattani' MLR corrections for AF cases         
          IF(MMLR(1).GE.-1) THEN         !! first for the 2x2 cases ...
              CmEFF(4)= CmVAL(4) + 0.25*CmVAL(2)**2/De
              CmEFF(5)= CmVAL(5) + 0.25*CmVAL(3)**2/De
              WRITE(6,716) MMLR(4),MMLR(4),CmEFF(4)
              WRITE(6,718) MMLR(5),MMLR(5),CmEFF(5)
c*  prepare C9{adj} coefficients for addition to chosen root
              MMLR(8)= 9               !! These terms added just
              MMLR(9)= 9               !! before exit from  AFdiag
              Cmeff(8)= 0.5*CmVAL(2)*CmEFF(4)/De   
              WRITE(6,714) MMLR(8),'Sigm',CmEFF(8)
              Cmeff(9)= 0.5*CmVAL(3)*CmEFF(5)/De
              WRITE(6,714) MMLR(9),'^3Pi',CmEFF(9)
              ENDIF
          IF(MMLR(1).LE.-2) THEN         !! now for the 3x3 cases ...
              CmEFF(5)= CmVAL(5) + 0.25*CmVAL(2)**2/De
              WRITE(6,716) MMLR(5),MMLR(5),CmEFF(5)
              CmEFF(6)= CmVAL(6) + 0.25*CmVAL(3)**2/De
              WRITE(6,720) MMLR(6),MMLR(6),CmEFF(6)
              CmEFF(7)= CmVAL(7) + 0.25*CmVAL(4)**2/De
              WRITE(6,718) MMLR(7),MMLR(7),CmEFF(7)
c*  prepare C9{adj} coefficients for addition to chosen root
              MMLR(11)= 9               !! These terms added just
              MMLR(12)= 9               !! before exit from  AFdiag
              MMLR(13)= 9
              Cmeff(11)= 0.5*CmVAL(2)*CmEFF(5)/De   
              IF(MMLR(1).EQ.-2) WRITE(6,714) MMLR(11),'Sigm',CmEFF(11)
              Cmeff(12)= 0.5*CmVAL(3)*CmEFF(6)/De
              IF(MMLR(1).EQ.-3) WRITE(6,714) MMLR(12),'^3Pi',CmEFF(12)
              Cmeff(13)= 0.5*CmVAL(4)*CmEFF(7)/De
              IF(MMLR(1).EQ.-4) WRITE(6,714) MMLR(13),'^1Pi',CmEFF(13)
              ENDIF
          ENDIF
      RETURN
      END
c23456789012345678901234567890123456789012345678901234567890123456789012

c***********************************************************************
      SUBROUTINE DYIDPJ(IDAT,NDATA,NPTOT,IFXP,YC,PV,PD,PS)
c** Subroutine to calculate potential function value YC at distance
c  RDIST= RTP(IDAT), and its partial derivatives w.r.t. the various
c  potential parameters.  If  IDAT.LE.1  generate a new set of 
c internal potential variables, while if  IDAT > 1  use SAVED values
c... [Must ensure that calculations based on the current UPDATED PV(j)]
c-----------------------------------------------------------------------
c********* Introduced AFdiag) 22 December 2015 **********
c********* Last updated 27 January 2015 **********
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** This 'Block' Data Utility routine with the name 'arrsizes.h', that
c   governs array dimensioning in program betaFIT  MUST be installed
c   under this name in the same (sub)directory containing the FORTRAN
c   file(s) for this Program when it is being compiled, or incorporated
c   into the program wherever dimensioning is required (as it is below)
c-----------------------------------------------------------------------
      INTEGER MXDATA, LMAX, MXPARM, NCMMAX
      PARAMETER (MXDATA= 1001, LMAX= MXDATA, MXPARM=50, NCMMAX= 25)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER  i,j,m,IDAT, IFXP(MXPARM),JFXRe,JFXDe,JFXVMIN,NPARM
     1      ,NDATA,NPTOT     !! 2 param unused here but reqd. by NNLSSRR
      REAL*8 YC,PV(MXPARM),PD(MXPARM),PS(MXPARM),DM(NCMMAX),DMP(NCMMAX),
     1 DMPP(NCMMAX),BETAsm,RTPp,RTPq,Rep,AREF,AREFp,AREFq,ype,dype,dyqe,
     2 Scalc,betaINF,yp,yq,yPOW,XP,XPW,DER,TCM,UM,TTMM,DERP,DERPe,DSUM,
     3 ypRe,dypRedRe,betaRe,yRePOW,FCT,FCT2,ULR,ULRe,dULRdr,
     4 dULRedRe,d2ULRe,DbetaRe,dAAdbi,dAAdRe,dBBdRe,DDER,T0,
     5 T0P,T1,RE3,RE6,RE8,RTP3,RTP6,RTP8,RDIST,C6adj,C9adj,BETAN,beta0,
     6 C1tst,VLIM
      REAL*8  DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),DEIGRe(1,1)
c-----------------------------------------------------------------------
      REAL*8 Re,De,VMIN,RREF,AA,BB,as,bs,rhoAB,BETA, CmVAL(NCMMAX),
     1  CmEFF(NCMMAX),dULRdCm(NCMMAX),RTP(MXDATA),SAS(MXDATA,MXPARM)
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,sVSR2,IDSTT,NCMM,MCMM,q,p,
     1  NPbeta,APSE,MMLR(NCMMAX),IFXCm(NCMMAX),updateCmADJ
      COMMON /DATABLK/Re,De,VMIN,RREF,AA,BB,as,bs,rhoAB,CmVAL,CmEFF,
     1 dULRdCm,RTP,SAS,BETA,PSEL,IFXRe,IFXDe,IFXVMIN,sVSR2,IDSTT,NCMM,
     2  MCMM,q,p,NPbeta,APSE,MMLR,IFXCm,updateCmADJ
c-----------------------------------------------------------------------
      SAVE JFXRe,JFXDe,JFXVMIN, AREF,AREFp,AREFq,Rep,C6adj,C9adj,
     1  betaINF,BETA0, ULRe,dULRedRe
c=======================================================================
      IF(ABS(IDAT).LE.1) THEN
          JFXRe= IFXP(NPbeta+1) 
          JFXDe= IFXP(NPbeta+2) 
          JFXVMIN= IFXP(NPbeta+3)
          ENDIF
      RDIST= RTP(IDAT) 
      DO  j=1,MXPARM
          PD(j)= 0.d0
          ENDDO
c=======================================================================
c** For case of an  EMO(q)  potential
c-----------------------------------------------------------------------
      IF(PSEL.EQ.1) THEN
          IF(ABS(IDAT).LE.1) THEN
              IF(JFXRe.LE.0) Re= PV(NPbeta+1) 
              IF(JFXDe.LE.0) De= PV(NPbeta+2)
              IF(JFXVMIN.LE.0) VMIN= PV(NPbeta+3) 
              AREF= RREF 
              IF(RREF.LE.0) AREF= Re
              AREFp= AREF**p 
              ENDIF
          RTPp= RDIST**p
          yp= (RTPp - AREFp)/(RTPp + AREFp) 
          yPOW= 1.d0 
          BETAsm= PV(1) 
          DSUM= 0.d0 
          IF(NPbeta.GE.2) THEN
              DO  j= 2,NPbeta
                  IF(RREF.LE.0.d0) DSUM= DSUM+ (j-1)*PV(j)*yPOW 
                  yPOW= yPOW*yp
                  BETAsm= BETAsm+ yPOW*PV(j) 
                  ENDDO
              ENDIF
          XP= DEXP(-BETAsm*(RDIST- Re)) 
          YC= De*(1.d0 - XP)**2 + VMIN 
          DER= 2.d0*De*(1.d0- XP)*XP 
          DERP= DER*(RDIST- Re) 
          DO  j= 1, NPbeta
              PD(j)= DERP
              DERP= DERP*yp
              ENDDO
c** If appropriate, also get partial derivative w.r.t. Re
          IF(JFXRE.LE.0) THEN
              IF(RREF.LE.0.d0) BETAsm= BETAsm +(RDIST- Re)*DSUM*0.5d0*
     1                                             (p/Re)*(1.d0-yp**2)
              PD(NPbeta+1)= -DER*BETAsm
              ENDIF
c** If appropriate, also get partial derivative w.r.t. De & VMIN
          IF(JFXDe.LE.0) PD(NPbeta+2)= (1.d0- XP)**2
          IF(JFXVMIN.LE.0) PD(NPbeta+3)= 1.d0
          ENDIF
c=======================================================================
c  For the case of an  MLR(q,p)  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.2) THEN
          IF(ABS(IDAT).LE.1) THEN
c*** Begin with calculations at  r=Re  to get uLR(Re) and betaINF
              IF(JFXRe.LE.0) Re= PV(NPbeta+1) 
              IF(JFXDe.LE.0) De= PV(NPbeta+2) 
              IF(JFXVMIN.LE.0) VMIN= PV(NPbeta+3) 
              AREF= RREF
              IF(RREF.LE.0.d0) AREF= Re
              AREFp= AREF**p 
              AREFq= AREF**q 
              Rep= Re**p 
c!!  update Cm{adj} values at beginning of each iteration
              IF(updateCmADJ.LE.0) CALL quadCORR(NCMM,MCMM,NCMMAX,MMLR,
     1                                                 De,CmVAL,CmEFF)
              updateCmADJ= updateCmADJ+ 1
C!!!
              IF(MMLR(1).LE.0) THEN
c** for Aubert-Frecon 2x2 or 3x3 diagonalization eigenvalues
                  VLIM= VMIN+ De                        !!!!   0.d0
                  CALL AFdiag(Re,VLIM,NCMM,NCMMax,MMLR,CmEFF,rhoAB,
     1                             sVSR2,IDSTT,ULRe,dULRdCm,dULRedRe)
                ELSE
c** For 'normal' inverse-power sum MLR/MLJ case, with or without damping
                  CALL dampF(Re,rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,
     1                                                        DMP,DMPP)
                  ULRe= 0.d0
                  dULRedRe= 0.d0 
                  DO  m= 1,MCMM
                      T0= CmEFF(m)/Re**MMLR(m) 
                      IF(rhoAB.GT.0.d0) T0= T0*DM(m)
                      ULRe= ULRe + T0 
                      dULRedRe= dULRedRe - MMLR(m)*T0/Re    !! corrected!!
                      IF(rhoAB.GT.0.d0)  dULRedRe= dULRedRe + 
     1                                                 T0*DMP(m)/DM(m)
                      ENDDO
                ENDIF
              betaINF= DLOG(2.d0*De/ULRe)
              IF((IDSTT.GT.1).AND.(sVSR2.EQ.-1)) THEN
c** For MLR with DS(s=-1/2) damping constrained to have C1= Z1*Z2= IDSTT
c   then fix the  geta(N) value ...
                  BETA0= 0.d0 
                      DO  m= 1,MCMM
                          BETA0= BETA0 + CmEFF(m)*DSQRT(3.69d0*rhoAB/
     1                                         MMLR(m))**(2*MMLR(m)-1)
                      ENDDO
                  C1tst= BETA0/uLre 
                  BETAN= DLOG(DSQRT(4.d0*De*IDSTT*116140.97d0)/BETA0) 
                  BETA0= DLOG(DSQRT(IDSTT*116140.97d0/De)*uLRe/BETA0) 
                  C1tst= De*(DEXP(BETA0)*C1tst)**2
c*** For constrained C1/r case, now define  beta(N)
                  BETAN= -0.5D0*BETAN*(-1)**NPbeta 
                  T0= 1.d0 
                  DO  j= NPbeta-1, 1, -1
                      BETAN= BETAN + T0*PV(j) 
                      T0= -T0 
                      ENDDO
                  PV(NPbeta)= BETAN
                  IFXP(NPbeta)= 1 
                  WRITE(6,666) C1tst/116140.97d0, NPbeta,PV(NPbeta)
  666 FORMAT(/'  Test  C1(sr)=',1Pd10.3,'    PV(',i3,')=',1PD16.8/)
                  ENDIF
              ENDIF
c*** END of calculations at  r=Re for  IDAT=1 to get uLR(Re) and betaINF
c--------- now proceed with calcn. for  IDAT.ge.1  points --------------
c... First, generate the MLR exponent variables & function .............
          RTPp= RDIST**p 
          RTPq= RDIST**q 
          yp= (RTPp - AREFp)/(RTPp + AREFp)
          yq= (RTPq - AREFq)/(RTPq + AREFq)
          ype= (RTPp - Rep)/(RTPp + Rep)
          IF(APSE.GT.0) THEN
c*** Case of Pashov natural spline exponent ............................  
c... use a spline through the exponent values defined by the input 
c    points to generate values of that exponent at the desired 
c    spline-definition points
              XP= 0.d0
              DO  J= 1, NPbeta
                  PD(J)= SAS(IDAT,J)
                  XP= XP + PV(J)*PD(J)
                  ENDDO
              XP= XP*ype     !!  hey ... I had dropped the  ype here  !!!
            ELSE
c... For conventional constrained polynomial exponent function
              yPOW= 1.d0 - yp
              BETAsm= PV(1)*yPOW
              DSUM= 0.d0
              IF(NPbeta.GE.2) THEN
                  DO  j= 2,NPbeta
                      IF(RREF.LE.0.d0) DSUM= DSUM + PV(j)*(j-1)*yPOW
                       yPOW= yPOW*yq
                      BETAsm= BETAsm+ yPOW*PV(j)
                      ENDDO
                  ENDIF
              XP= (BETAsm + betaINF*yp)*ype         !! total exponent
            ENDIF
c... Now, generate  uLR(r) for actual distance ... RDIST=RTP(i) ........
          IF(MMLR(1).LE.0) THEN
c** For Aubert-Frecon 2x2 or 3x3 diagonalization  uLR(r)
              CALL AFdiag(RDIST,VLIM,NCMM,NCMMax,MMLR,CmEFF,rhoAB,
     1                                 sVSR2,IDSTT,ULR,dULRdCm,dULRdR)
            ELSE
c** For normal inverse-power sum MLR/MLJ case, with or without damping
              CALL dampF(RDIST,rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,
     1                                                       DMP,DMPP)
              ULR= 0.d0
              DO  m= 1,MCMM
                  T0= CmEFF(m)/RDIST**MMLR(m)
                  IF(rhoAB.GT.0.d0) T0= T0*DM(m)
                  ULR= ULR + T0
                  ENDDO
            ENDIF
          XPW= DEXP(-XP) * ULR/ULRe 
          YC= De*(1.d0 - XPW)**2 + VMIN
          DER= 2.d0*De*(1.d0- XPW)*XPW
          yPOW= DER*ype*(1.d0- yp)
          IF(APSE.GT.0) THEN
c... finalize derivative w.r.t. exponent beta-function spline points ...
              DO  J= 1,NPbeta
                  PD(J)= PD(J)*DER * ype      !! I had  forgotten  ype
                  ENDDO
            ELSE
c... finalize derivative w.r.t. exponent polynomial coefficient ....
              DO  j= 1,NPbeta
                  PD(j)= yPOW
                  yPOW= yPOW*yq
                  ENDDO
            ENDIF
c** If appropriate, also get partial derivative w.r.t. De & VMIN
          IF(JFXDe.LE.0) PD(NPbeta+2)= (1.d0- XPW)**2 + DER*ype*yp/De
          IF(JFXVMIN.LE.0) PD(NPbeta+3)= 1.d0
c** If appropriate, also get partial derivative w.r.t. Re
          IF(JFXRe.LE.0) THEN
              dype= -0.5d0*(p/RE)*(1.d0 - yp**2)
              dyqe= -0.5d0*(q/RE)*(1.d0 - yq**2)
              PD(NPbeta+1)= DER*(dype*XP/ype 
     1                                + (1.d0 - ype*yp)*dULRedRe/ULRe)
              IF(RREF.LE.0.d0) PD(NPbeta+1)= PD(NPbeta+1) + 
     1                DER*(dype*(betaINF - BETAsm/(1.d0-yp))+ dyqe*DSUM)
              ENDIF
          ENDIF
c=======================================================================
c** For the case of a  DELR_p  potential
c-----------------------------------------------------------------------
      IF(PSEL.EQ.3) THEN
          IF(ABS(IDAT).LE.1) THEN
              IF(JFXRe.LE.0) Re= PV(NPbeta+1)
              IF(JFXDe.LE.0) De= PV(NPbeta+2)
              IF(JFXVMIN.LE.0) VMIN= PV(NPbeta+3)
              AREF= RREF
              IF(RREF.LE.0.d0) AREF= Re
              AREFp= AREF**p
              ULRe= 0.d0
              dULRedRe= 0.d0
              d2ULRe= 0.d0
c-----------------------------------------------------------------------
c** Evaluate uLR and its first 2 derivs. w,r,t, r,  at  r=Re ... 
              IF(MMLR(1).LE.0) THEN
c** For Aubert-Frecon 2x2 or 3x3 diagonalization  uLR(r)
                  CALL AFdiag(Re,VLIM,NCMM,NCMMax,MMLR,CmEFF,rhoAB,
     1                             sVSR2,IDSTT,ULRe,dULRdCm,dULRedRe)
                ELSE
                  CALL dampF(Re,rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,
     1                                                       DMP,DMPP)
                  ULRe= 0.d0
                  dULRedRe= 0.d0
                  d2ULRe= 0.d0
                  DO  m= 1,NCMM
                      T0= CmEFF(m)/Re**MMLR(m)
                      ULRe= ULRe + T0*DM(m)
                      dULRedRe= dULRedRe + T0*(DMP(m)- DM(m)*MMLR(m)/RE)
                      d2ULRe= d2ULRe+T0*(DMPP(m)- 2.d0*MMLR(m)*DMP(m)/RE
     1                           + MMLR(m)*(MMLR(m)+1.d0)*DM(m)/RE**2)
                      ENDDO
                  Rep= Re**p
                  ypRe= (Rep - AREFp)/(Rep + AREFP)
                  dypRedRe= 2.d0*p*Rep*AREFp/(Re*(Rep + AREFP)**2)
                  IF(RREF.LE.0.d0) dypRedRe= 0.d0
                ENDIF  
c** now get AA & BB their derivs w.r.t. Re  & parially w.r.t. beta_i  
              betaRe= PV(1)
              DbetaRe= 0.d0            !! this is d{beta}/d{y}  at r= Re
              yPOW= 1.d0
              IF(NPbeta.GE.1) THEN
                  DO j= 1,NPbeta
                      DbetaRe= DbetaRe + j*PV(j+1)*yPOW
                      yPOW= yPOW*ypRe
                      betaRe= betaRe + PV(J+1)*yPOW   !! corrn found by AEM
c                     betaRe= betaRe + yPOW   !! error found by AEM
                      ENDDO
                  ENDIF
              AA= De - ULRe - dULRedRe/betaRe
              BB= 2.d0*(De - ULRe) - dULRedRe/betaRe
              dAAdbi= dULRedRe/betaRe**2
              dAAdRe= -dULRedRe -d2ULRe/betaRe +dAAdbi*DbetaRe*dypRedRe
              dBBdRe= dAAdRe - dULRedRe
              ENDIF
c===== end of calcn. for properties at Re performed, for 1'st point ====
          RTPp = RDIST**p
          yp= (RTPp - AREFp)/(RTPp + AREFp)
          ULR= 0.d0
c ... evaluate uLR at the actual distance ...
          IF(MMLR(1).LE.0) THEN
c** For Aubert-Frecon 2x2 or 3x3 diagonalization  uLR(r)
              CALL AFdiag(RDIST,VLIM,NCMM,NCMMax,MMLR,CmEFF,rhoAB,
     1                                 sVSR2,IDSTT,ULR,dULRdCm,dULRdR)
            ELSE
              CALL dampF(RDIST,rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,
     1                                                       DMP,DMPP)
              ULR= 0.d0
              DO  m= 1,NCMM
                  ULR= ULR + DM(m)*CmEFF(m)/RDIST**MMLR(m)
                  ENDDO
            ENDIF
          yPOW= 1.d0
          yRePOW= 1.d0
          BETAsm= PV(1)
          DSUM= 0.d0
          IF(NPbeta.GE.2) THEN
              DO  j= 2,NPbeta
                  IF(RREF.LE.0.d0) DSUM= DSUM+ (j-1)*PV(j)*yPOW
                  yPOW= yPOW*yp
                  BETAsm= BETAsm+ yPOW*PV(j)
                  ENDDO
              ENDIF
          XP= DEXP(-BETAsm*(RDIST- Re))
          YC= (AA*XP - BB)*XP + De - ULR + VMIN
          DER= XP*(BB - 2.d0*AA*XP)
          DERP= DER*(RDIST- Re)
          DERPe= XP*(XP- 1.d0)*dAAdbi
          DO  j= 1,NPbeta
              PD(j)= DERP + DERPe
              DERP= DERP*yp
              DERPe= DERPe*yp             !! another AEM correction
              ENDDO
c** If appropriate, also get partial derivative w.r.t. De & VMIN
          IF(JFXDe.LE.0) PD(NPbeta+2)= (XP - 2.d0)*XP + 1.d0
          IF(JFXVMIN.LE.0) PD(NPbeta+3)= 1.d0
          IF(JFXRe.LE.0) THEN
c??? For the old sitn. in which  r_{ref} \equiv r_e ...............
c** If appropriate, also get partial derivative w.r.t. Re
              PD(NPbeta+1)=  BETAsm*(2.d0*XP*AA - BB)*XP 
     1                                        + XP*(dAAdRe*XP- dBBdRe)
              ENDIF
          ENDIF
c=======================================================================
c=======================================================================
c** For the case of a  GPEF(p,as,bs)  potential
c-----------------------------------------------------------------------
      IF(PSEL.EQ.4) THEN
          IF(ABS(IDAT).LE.1) THEN
              JFXVMIN= IFXP(NPbeta+1)
              IF(JFXVMIN.LE.0) VMIN= PV(NPbeta+1)
              JFXRe= IFXP(NPbeta+2)
              IF(JFXRe.LE.0) Re= PV(NPbeta+2)
              Rep= Re**p
              ENDIF
          RTPp= RDIST**p
          yp= (RTPp - Rep)/(as*RTPp + bs*Rep)
          yPOW= yp**2
          BETA= PV(1)*yPOW
          DSUM= 2.d0*PV(1)*yp
          IF(IFXDe.LE.0) yPOW= PV(1)*yPOW
          PD(1)= yp**2
          IF(NPbeta.GT.1) THEN
              DO  j= 2,NPbeta
                  DSUM= DSUM + (j+1)*PV(j)*yPOW
                  yPOW= yPOW* yp
                  PD(j)= yPOW
                  BETA= BETA+ PV(j)*yPOW
                  ENDDO
              ENDIF
          YC= BETA + VMIN
          IF(JFXVMIN.LE.0) PD(NPbeta+1)= 1.d0
          IF(JFXRe.LE.0) THEN
              PD(NPbeta+2)= -DSUM* (p/Re)*Rep*RTPp*(as+bs)/
     1                                            (as*RTPp +bs*Rep)**2
              ENDIF 
          ENDIF
c=======================================================================
c%%%%%%%  Optional printout for testing
cc        if(IDAT.eq.1) then
cc              NPARM= NPbeta
cc              IF(IFXDe.LE.0) NPARM= NPARM+ 1
cc              IF(IFXRe.LE.0) NPARM= NPARM+ 1
cc              IF(IFXVMIN.LE.0) NPARM= NPARM+ 1
cc              write(8,700) nparm,(PV(i),i=1,nparm)
cc              write(8,702)  (i,i=1,min0(5,nparm))
cc              ENDIF
cc    write(8,704) i,yc,(pd(i),i=1,nparm)
cc700 FORMAT(///' Partial derivatives for',i3,' input parameters:'/
cc   1  (1P5D16.8))
cc702 FORMAT('  I    YC   ',5('     PD(',i1,')   ':))
cc704 format(i3,f9.2,1P5D13.5:/(12x,5d13.5:))
c%%%%%%%
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c** Asen Pashov's subroutines for constructing spline functions and
c   their derivatives.
      double precision function Scalc(x,m,n,XGRID,rKL,LMAX)
c** At the position 'x', Scalc is returned as the value of the m'th 
c  of the 'n' Sm(x) function defining a natural cubic spline through the
c  mesh points located at  x= XGRID(x_i), for i=1,n.  LMAX specifies the 
c  maximum number of mesh points x= XGRID(x_i) allowed by the calling program
c---------------------------------------------------------------------
      INTEGER  LMAX,I,K,KK,M,N
      REAL*8  x,y1,y2,XGRID(LMAX),rKL(LMAX,LMAX)
      k= 0
      kk= 0
      do i=2,n
c... select interval
          if ((x.gt.XGRID(i-1)).and.(x.le.XGRID(i)))  k=i
          end do
      if (x.lt.XGRID(1)) then
          k=2
          kk=1
          end if
      if (x.gt.XGRID(n)) then
          k=n
          kk=1
          end if
      if(x.eq.XGRID(1)) k=2
      y1=XGRID(k-1)
      y2=XGRID(k)
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
      double precision function Sprime(x,m,n,XGRID,rKL,LMAX)
c** At the position 'x', evaluate the derivative w.r.t. x of the m'th 
c  Sm(x) function contributing the definition of the the natural cubic
c  spline defined by function values at the  n  points  XGRID(i) [i=1,n]
      INTEGER i,k,kk,m,n,LMAX
      REAL*8 x,del,y1,y2,XGRID(LMAX),rKL(LMAX,LMAX)
      k=0
      kk=0
      do i=2,n
          if((x.gt.XGRID(i-1)).and.(x.le.XGRID(i)))  k=i
          enddo
      if(x.lt.XGRID(1)) then
          k=2
          kk=1
          end if
      if (x.gt.XGRID(n)) then
          k=n
          kk=1
          end if
      if (x.eq.XGRID(1)) k=2
      y1=XGRID(k-1)
      y2=XGRID(k)
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
      subroutine Lkoef(NGRID,XGRID,rKL)   
c** Call this subroutine with list of the 'NGRID' spline x_i values in 
c   array 'XGRID' with maximum dimension 'LMAX', and it will return the 
c   LMAX x LMAX  array of 'rKL' coefficients used for generating the 
c   'NGRID' S_{NGRID}(x) spline coefficient functions
c----------------- Based on nespl subroutine ---------------------------
c** CAUTION .. must dimension internal arrays B, INDX & vv @ compilation
cc    INCLUDE 'arrsizes.h'                !! needed only to define  LMAX
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** This 'Block' Data Utility routine with the name 'arrsizes.h', that
c   governs array dimensioning in program betaFIT  MUST be installed
c   under this name in the same (sub)directory containing the FORTRAN
c   file(s) for this Program when it is being compiled, or incorporated
c   into the program wherever dimensioning is required (as it is below)
c-----------------------------------------------------------------------
      INTEGER MXDATA, LMAX, MXPARM, NCMMAX
      PARAMETER (MXDATA= 1001, LMAX= MXDATA, MXPARM=50, NCMMAX= 25)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c*** Subroutine by written A. Pashov (2000) ----------------------------
      INTEGER I,J,NGRID,INDX(1:LMAX)
      REAL*8 XGRID(LMAX),rKL(LMAX,LMAX),B(LMAX,LMAX),vv(LMAX), d
c ...  note vv dimensioned here, but only used in   ludcmp !!
      DO  i= 1,LMAX
          DO  j= 1,LMAX
              rKL(i,j)= 0.d0
              B(i,j)= 0.d0
              ENDDO
          ENDDO
      rKL(1,1)= (XGRID(3)-XGRID(1))/3.d0
      rKL(1,2)= (XGRID(3)-XGRID(2))/6.d0
      do i= 2,NGRID-3
          rKL(i,i-1)= (XGRID(i+1)-XGRID(i))/6.d0
          rKL(i,i)= (XGRID(i+2)-XGRID(i))/3.d0
          rKL(i,i+1)= (XGRID(i+2)-XGRID(i+1))/6.d0
          end do
      rKL(NGRID-2,NGRID-3)= (XGRID(NGRID-1)-XGRID(NGRID-2))/6.d0
      rKL(NGRID-2,NGRID-2)= (XGRID(NGRID)-XGRID(NGRID-2))/3.d0  
      do i= 1,NGRID-2
          B(i,i)= 1.d0/(XGRID(i+1)-XGRID(i))
          B(i,i+1)= -1.d0/(XGRID(i+2)-XGRID(i+1))-1.d0/
     1                                           (XGRID(i+1)-XGRID(i))
          B(i,i+2)= 1.d0/(XGRID(i+2)-XGRID(i+1))
          end do  
      call ludcmp(rKL,NGRID-2,LMAX,indx,vv,d)
      do i= 1,NGRID 
          call lubksb(rKL,NGRID-2,LMAX,indx,B(1,i))
          end do 
      do i= 1,NGRID-2
          do j= 1,NGRID
              rKL(j,i+1)= B(i,j)
              end do
          end do 
      do i= 1,NGRID
          rKL(i,1)= 0.0d0
          rKL(i,NGRID)= 0.0d0
          end do
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ludcmp(rKL,NGRID,LMAX,indx,vv,d)
c** Subroutine taken from "Numerical Recipies", by Press et al. (1986)
      INTEGER NGRID,LMAX,indx(LMAX),NMAX,i,imax,j,k
      double precision d,rKL(LMAX,LMAX),vv(LMAX),TINY,aamax,dum,sum
      PARAMETER (TINY= 1.0e-20)
      d= 1.d0
      do  i= 1,NGRID
          aamax= 0.d0
          do  j= 1,NGRID
              if (abs(rKL(i,j)).gt.aamax) aamax= abs(rKL(i,j))
              enddo
          if (aamax.eq.0.) WRITE(6,*) 'singular matrix in ludcmp'
          vv(i)= 1.d0/aamax
          enddo
      do  j= 1,NGRID
          do  i= 1,j-1
              sum= rKL(i,j)
              do  k= 1,i-1
                  sum= sum-rKL(i,k)*rKL(k,j)
                  enddo
              rKL(i,j)= sum
              enddo
          aamax= 0.d0
          do  i= j,NGRID
              sum= rKL(i,j)
              do  k= 1,j-1
                  sum= sum-rKL(i,k)*rKL(k,j)
                  enddo
              rKL(i,j)= sum
              dum= vv(i)*abs(sum)
              if (dum.ge.aamax) then
                  imax= i
                  aamax= dum
                  endif
              enddo
          if(j.ne.imax)then
              do  k= 1,NGRID
                  dum= rKL(imax,k)
                  rKL(imax,k)= rKL(j,k)
                  rKL(j,k)= dum
                  enddo
              d= -d
              vv(imax)= vv(j)
              endif
          indx(j)= imax
          if(rKL(j,j).eq.0.)rKL(j,j)= TINY
              if(j.ne.NGRID)then
                  dum= 1.d0/rKL(j,j)
                  do  i= j+1,NGRID
                      rKL(i,j)= rKL(i,j)*dum
                      enddo
                  endif
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE lubksb(rKL,NGRID,LMAX,indx,b)
c** Subroutine taken from "Numerical Recipies", by Press et al. (1986)
      INTEGER i,ii,j,ll, NGRID,LMAX,indx(LMAX)
      double precision rKL(LMAX,LMAX),b(LMAX), sum
      ii= 0
      do  i= 1,NGRID
          ll= indx(i)
          sum= b(ll)
          b(ll)= b(i)
          if (ii.ne.0)then
              do  j= ii,i-1
                  sum= sum-rKL(i,j)*b(j)
                  enddo
            else if (sum.ne.0.) then
              ii= i
            endif
          b(i)= sum
          enddo
      do  i= NGRID,1,-1
          sum= b(i)
          do  j= i+1,NGRID
              sum= sum-rKL(i,j)*b(j)
              enddo
          b(i)= sum/rKL(i,i)
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE NLLSSRR(NDATA,NPTOT,NPMAX,CYCMAX,IROUND,ROBUST,LPRINT,
     1                      IFXP,YO,YU,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c**  Program for performing linear or non-linear least-squares fits and
c  (if desired) automatically using sequential rounding and refitting 
c  to minimize the numbers of parameter digits which must be quoted [see
c  R.J. Le Roy, J.Mol.Spectrosc. 191, 223-231 (1998)].         25/03/16
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c             COPYRIGHT 1998-2016  by  Robert J. Le Roy                +
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
c    CYCMAX is the upper bound on the allowed number of iterative cycles 
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
      INTEGER MXPdim             !! internal limit on max # parameters
      PARAMETER (MXPdim=50)      !! must be .GE. external max # NPMAX
      INTEGER I,J,K,L,IDF,ITER,NITER,CYCMAX,IROUND,ISCAL,JROUND,LPRINT,
     1 NDATA,NPTOT,NPMAX,NPARM,NPFIT,JFIX,QUIT,ROBUST,
     2 IFXP(NPMAX),JFXP(MXPdim)
      REAL*8  YO(NDATA), YU(NDATA), YD(NDATA), PV(NPTOT), PU(NPTOT), 
     1 PS(NPTOT),PSS(MXPdim),PC(MXPdim),PCS(MXPdim),PX(MXPdim),
     2 PY(MXPdim),CM(NPMAX,NPMAX), F95(10),
     3 RMSR, RMSRB, DSE, TSTPS, TSTPSB, TSTPU, TFACT, S, YC, Zthrd
      DATA F95/12.7062D0,4.3027D0,3.1824D0,2.7764D0,2.5706D0,2.4469D0,
     1  2.3646D0,2.3060D0,2.2622D0,2.2281D0/
      IF((NPTOT.GT.NPMAX).OR.(NPTOT.GT.MXPdim).OR.(NPTOT.GT.NDATA)
     1                                     .OR.(NPMAX.GT.MXPdim)) THEN
c** If array dimensioning inadequate, print warning & then STOP
          WRITE(6,602) NPTOT,MXPdim,NPMAX,NDATA
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
      DO 50 ITER= 1, CYCMAX
          ISCAL= 0
          NITER= NITER+ 1
          DSE= 0.d0 
          TSTPSB= TSTPS
          RMSRB= RMSR
c** Zero out various arrays
   10     IF(NPARM.GT.0) THEN
              DO  I = 1,NPARM
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
c  values {PV} to generate predicted datum # I [y(calc;I)=YC] and its
c  partial derivatives w.r.t. each of the parameters, returning the 
c  latter in 1-D array PC.  See dummy sample version at end of listing.
c* NOTE 1: if more convenient, DYIDPJ could prepare the y(calc) values 
c     and derivatives for all data at the same time (when I=1), but only
c     returned the values here one datum at a time (for I > 1).]
c* NOTE 2: the partial derivative array PC returned by DYIDPJ must have
c     an entry for every parameter in the model, though for parameters 
c     which are held fixed [JFXP(j)=1], those PC(j) values are ignored.
              CALL DYIDPJ(I,NDATA,NPTOT,JFXP,YC,PV,PC,PSS)
              IF(NPARM.LT.NPTOT) THEN
c** For constrained parameter or sequential rounding, collapse partial 
c   derivative array here
                  DO  J= NPTOT,1,-1
                      IF(JFXP(J).GT.0) THEN
c!! First ... move derivative for constrained-parameter POTFIT case
cc                        IF(JFXP(J).GT.1) THEN
cc                            write(6,666) I,J,PC(J),JFXP(J),PC(JFXP(J))
cc                            PC(JFXP(J))= PC(JFXP(J))+ PC(J)
cc666 FORMAT(' For  IDAT=',I5,'  add PC(',I3,') =',1pD15.8,
cc   1  '  to PC(',0pI3,') =',1pD15.8)
cc                            ENDIF
c  ... now continue collapsing partial derivative array
                          IF(J.LT.NPTOT) THEN
                              DO  K= J,NPTOT-1
                                  PC(K)= PC(K+1)
                                  ENDDO
                              ENDIF
                          PC(NPTOT)= 0.d0
                          ENDIF
                      ENDDO
                  ENDIF
              YD(I)= YC - YO(I)
              S = 1.D0/YU(I)
cc *** For 'Robust' fitting, adjust uncertainties here
              IF(Zthrd.GT.0.d0) S= 1.d0/DSQRT(YU(I)**2 + Zthrd*YD(I)**2)
              YC= -YD(I)*S
              DSE= DSE+ YC*YC
              IF(NPARM.GT.0) THEN
                  DO  J = 1,NPARM
                      PC(J) = PC(J)*S
                      PS(J) = PS(J)+ PC(J)**2
                      ENDDO
                  CALL QROD(NPARM,NPMAX,NPMAX,CM,PC,PU,YC,PX,PY)
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
                  YC = 0.D0
                  DO  K = J,NPARM
                      YC = YC + CM(I,K) * CM(J,K)
                      ENDDO
                  CM(I,J) = YC
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
          YC= DSE*0.1d0/DFLOAT(NPARM)
          S= DSE*TFACT
          DO  J = 1,NPARM
              PU(J)= S* PU(J)
              PS(J)= YC*DSQRT(NDATA/PS(J))
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
cc               IF(JFXP(J).GT.1) THEN       !! a special PotFit option
c** If this parameter constrained to equal some earlier parameter ....
cc                   PV(J)= PV(JFXP(J))      
cc                   WRITE(6,668) J,JFXP(J),PV(J),ITER
cc                   ENDIF
cc668 FORMAT(' Constrain  PV('i3,') = PV(',I3,') =',1pd15.8,
cc   1 '  on cycle',i3)
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
          IF(ITER.GT.1) THEN
c** New Convergence test requires  RMSD to be constant to 1 part in 10^7
c     in adjacent cycles (unlikely to occur by accident!)
c** Replaces less severe requirement that  TSTPS < 1.0
              IF(ABS((RMSR/RMSRB)-1.d0).LT.1.d-07) THEN
                  IF(LPRINT.GE.3) WRITE(6,607) ITER,
     1                                      ABS(RMSR/RMSRB-1.d0),TSTPS
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
          JFIX= NPTOT
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
      YC= PV(JFIX)
      CALL ROUND(JROUND,NPMAX,NPTOT,NPTOT,JFIX,PV,PU,PS,CM)
      JFXP(JFIX)= 1
      IF(LPRINT.GE.2)
     1       WRITE(6,614) JFIX,YC,PU(JFIX),PS(JFIX),JFIX,PV(JFIX),RMSR
      NPARM= NPARM-1
      IF(NPARM.EQ.0) THEN
c** After rounding complete, make one more pass with all non-fixed 
c  parameters set free to get full correct final correlation matrix, 
c  uncertainties & sensitivities.  Don't update parameters on this pass!
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
  602 FORMAT(/' *** NLLSSRR problem:  [NPTOT=',i4,'] > min{MXPdim=',
     1  i4,' NPMAX=',i4,', NDATA=',i6,'}')
  604 FORMAT(' After Cycle #',i2,':  DRMSD=',1PD14.7,'    test(PS)=',
     1  1PD8.1,'   test(PU)=',D8.1)
  606 FORMAT(/' Effective',i3,'-cycle Cgce:  MAX{|change/unc.|}=',
     1  1PD8.1,' < 0.01   DRMSD=',D10.3)
  607 FORMAT(/' Full',i3,'-cycle convergence:  {ABS(RMSR/RMSRB)-1}=',
     1  1PD9.2,'  TSTPS=',D8.1)
  610 FORMAT(/ ' !! CAUTION !! fit of',i5,' parameters to',I6,' data not
     1 converged after',i3,' Cycles'/5x,'DRMS(deviations)=',1PD10.3,
     2 '    test(PS) =',D9.2,'    test(PU) =',D9.2/1x,31('**'))
  612 FORMAT((3x,'PV(',i4,') =',1PD22.14,' (+/-',D8.1,')    PS=',d8.1,
     1  '   PC=',d9.1))
  614 FORMAT(' =',39('==')/' Round Off  PV(',i4,')=',1PD21.13,' (+/-',
     1 D9.2,')    PS=',d9.2/4x,'fix PV(',I4,') as ',D19.11,
     2 '  & refit:  DRMS(deviations)=',D12.5)
  616 FORMAT(/i6,' data fit to',i5,' param. yields  DRMS(devn)=',
     1 1PD14.7:'  tst(PS)=',D8.1)
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
          IF(DABS(Z(1)).LE.0.D0) GOTO 10
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

c***********************************************************************
c     SUBROUTINE DYIDPJ(I,NDATA,NPTOT,IFXP,YC,PV,PD.PS)
c** Illustrative dummy version of DYIDPJ for the case of a fit to a
c  power series of order (NPTOT-1) in X(i). ***  For datum number-i, 
c  calculate and return  PD(j)=[partial derivatives of datum-i] w.r.t. 
c  each of the free polynomial coefficients varied in the fit 
c  (for j=1 to NPTOT).  **  Elements of the integer array IFXP indicate
c  whether parameter j is being held fixed [IFXP(j) > 0] or varied in
c  the fit [IFXP(j).le.0].  If the former, the partial derivative 
c  for parameter j should be  PD(j)= 0.0. 
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
      SUBROUTINE AFdiag(RDIST,VLIM,NCMM,NCMMax,MMLR,Cm,rhoAB,sVSR2,
     1                                       IDSTT,ULR,dULRdCm,dULRdR)
c***********************************************************************
c**   Aubert-Frecon Potential Model for u_{LR}(r)
c***********************************************************************
c** Subroutine to generate, at the onee distance RDIST, an eigenvalue 
c  of the 2x2 or 3x3 long-range interaction matrix described by Eqs.1
c and 10, resp., of J.Mol.Spec.188, 182 (1998) (Aubert-Frecon et al)
c** and its derivatives w.r.t. the C_m long-range parameters.
c***********************************************************************
c==> Input:  r= RDIST, VLIM, NCMM, m=MMLR & Cm's, rhoAB, sVSR2, IDSTT
c==> Output: ULR, partial derivatives dULRdCm & radial derivative dULRdR
c-----------------------------------------------------------------------
c** Original Version from Nike Dattani in June 2011 for 3x3 case
c** Generalized to incorporate 2x2 case, removed retardation terms and
c   incorporate damping  ...  by Kai Slaughter:                July 2014
c* rj:  C6{adj} & C9{adj} included in CmEFF & fixed dampF call  Jan 2016
c-----------------------------------------------------------------------
      INTEGER NCMMax
c-----------------------------------------------------------------------
      REAL*8 RDIST,VLIM,Cm(NCMMax),ULR,dULRdCm(NCMMax),dULRdR,R2,R3,R5,
     1       R6,R8,R9,T1,T0,T2,T0P,T0P23,Dm(NCMMax),Dmp(NCMMax),
     2       Dmpp(NCMMax),rhoAB,A(3,3),DR(3,3),Q(3,3),DMx(NCMMax,3,3),
     3       DMtemp(3,3),DEIGMx(NCMMax,1,1),DEIGMtemp(1,1),DEIGR(1,1),
     4       EIGVEC(3,1),RESID(3,1),W(3),RPOW(NCMMax),DELTAE,Modulus,Z
      INTEGER H,I,J,K,L,M,X,NCMM,MMLR(NCMMax),sVSR2,IDSTT,MMtemp
c-----------------------------------------------------------------------
      DELTAE=Cm(1)
      R2= 1.d0/RDIST**2
      R3= R2/RDIST
      R5= R2*R3
      R6= R3*R3
      R8= R6*R2
c-----------------------------------------------------------------------
c....... for rhoAB.le.0.0   returns Dm(m)=1 & Dmp(m)=Dmpp(m)=0  
      CALL dampF(RDIST,rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,Dm,Dmp,Dmpp)
c-----------------------------------------------------------------------
      IF(MMLR(1).GE.-1) THEN           !!  For the A (0)  or b (-1) state
c***********************************************************************
c************* Aubert Frecon 2x2 case   NCMM= 7  and  ...
c***              Cm(1) = DELTAE
c***              Cm(2) = C3Sig
c***              Cm(3) = C3Pi
c***              Cm(4) = C6Sig
c***              Cm(5) = C6Pi
c***              Cm(6) = C8Sig
c***              Cm(7) = C8Pi
c***********************************************************************
          T1= R3*(Dm(2)*(Cm(2)-Cm(3)) + R3*Dm(4)*(Cm(4)-Cm(5)) + 
     1        R5*Dm(6)*(Cm(6)-Cm(7)))/3.d0
          T0= DSQRT((T1 - Cm(1))**2 + 8.d0*T1**2)
          ULR= 0.5d0*(-Cm(1) + R3*(Dm(2)*(Cm(2)+Cm(3)) + 
     1         R3*Dm(4)*(Cm(4)+Cm(5)) + R5*Dm(6)*(Cm(6)+Cm(7))) + T0)
c-----------------------------------------------------------------------
          IF(MMLR(1).EQ.0) THEN     
              ULR= ULR + Cm(8)*R3*R6          !! add C9{adj correction
              ENDIF
c...  adjustment for the b-state
          IF(MMLR(1).EQ.-1) THEN
              ULR=ULR-T0
              ULR= ULR + Cm(9)*R3*R6          !! add C9{adj correction
              ENDIF
c...  now get derivatives
          T0P= 0.5d0*(9.d0*T1 - Cm(1))/T0
          T0P23= 0.5d0 + T0P/3.d0
c...  another adjustment for the b-state
          IF(MMLR(1).EQ.-1) T0P23=T0P23-2.d0*T0P/3.d0
          dULRdCm(1)= 0.d0
          dULRdCm(2)= R3*(T0P23)
          dULRdCm(3)= R3*(1.d0-T0P23)
          dULRdCm(4)= R6*(T0P23)
          dULRdCm(5)= R6*(1.d0 - T0P23)
          dULRdCm(6)= R8*T0P23
          dULRdCm(7)= R8*(1.d0-T0P23)
          T2        =-T0P*R3*((Dm(2)*(Cm(2)-Cm(3))+R3*(Dm(4)*2.d0*(Cm(4)
     1                -Cm(5))+R2*Dm(6)*8.d0/3.d0*(Cm(6)-Cm(7))))/RDIST
     2                +(Dmp(2)*(Cm(2)-Cm(3))+R3*Dmp(4)*(Cm(4)-Cm(5))+
     3                R2*R3*Dmp(6)*(Cm(6)-Cm(7)))/3.d0)
          dULRdR    = -R3*((1.5d0*Dm(2)*(Cm(2)+Cm(3)) + R3*(Dm(4)*3.d0*
     1                (Cm(4)+Cm(5))+4.d0*Dm(6)*R2*(Cm(6)+Cm(7))))/RDIST
     2                + 0.5d0*(Dmp(2)*(Cm(2)+Cm(3)) + Dmp(4)*R3*(Cm(4)+
     3                Cm(5)) + Dmp(6)*R3*R2*(Cm(6)+Cm(7)))) + T2
c... and a final adjustment for the b-state
          IF(MMLR(1).EQ.-1) dULRdR= dULRdR- 2.d0*T2
c-----------------------------------------------------------------------
      ELSE
c***********************************************************************
c********* Aubert Frecon 3x3 case   NCMM= 10  and ...
c*********        Cm(1) = DELTAE
c*********        Cm(2) = C3Sig
c*********        Cm(3) = C3Pi1
c*********        Cm(4) = C3Pi3
c*********        Cm(5) = C6Sig
c*********        Cm(6) = C6Pi1
c*********        Cm(7) = C6Pi3
c*********        Cm(8) = C8Sig
c*********        Cm(9) = C8Pi1
c*********        Cm(10)= C8Pi3
c***********************************************************************      
c...      Initialize interaction matrix to 0.d0
          DO  I= 1,3
              DO J= 1,3
                  A(I,J)=0.0D0
                  DR(I,J)=0.d0
                  DO  K= 1,NCMMax
                      DMx(K,I,J)=0.d0
                      ENDDO
                  ENDDO
              ENDDO
c...      Prepare interaction matrix  A
          DO  I= 2,NCMM,3
              RPOW(I)= RDIST**MMLR(I)
           A(1,1)= A(1,1) - Dm(I)*(Cm(I)+Cm(I+1)+Cm(I+2))/(3.d0*RPOW(I))
           A(1,2)= A(1,2) - Dm(I)*(Cm(I+2)+Cm(I+1)-2.d0*Cm(I))/(RPOW(I))
           A(1,3)= A(1,3) - Dm(I)*(Cm(I+2)-Cm(I+1))/(RPOW(I))
           A(2,2)= A(2,2) - Dm(I)*(Cm(I+2)+Cm(I+1)+4.d0*Cm(I))
     1                             /(6.d0*RPOW(I))
           A(3,3)= A(3,3) - Dm(I)*(Cm(I+2)+Cm(I+1))/(2.d0*RPOW(I))
           ENDDO
          A(1,1) = A(1,1) + VLIM
          A(1,2) = A(1,2)/(3.d0*DSQRT(2.d0))
          A(2,1) = A(1,2)
          A(2,2) = A(2,2) + VLIM + DELTAE
          A(2,3) = A(1,3)/(2.d0*DSQRT(3.d0))
          A(1,3) = A(1,3)/(DSQRT(6.d0))
          A(3,1) = A(1,3)
          A(3,2) = A(2,3)
          A(3,3) = A(3,3) + VLIM + DELTAE
c...      Prepare radial derivative of interaction matrix (? is it needed ?)
          DO  I= 2,NCMM,3
              DR(1,1)= DR(1,1) + Dm(I)*MMLR(I)*(Cm(I)+Cm(I+1)+Cm(I+2))
     1                             /(3.d0*RPOW(I)*RDIST)
     2                    -Dmp(I)*(Cm(I)+Cm(I+1)+Cm(I+2))/(3.d0*RPOW(I))
              DR(1,2)= DR(1,2) + Dm(I)*MMLR(I)*(Cm(I+2)+Cm(I+1)-2.d0*
     1                           Cm(I))/(RPOW(I)*RDIST)
     2                    -Dmp(I)*(Cm(I+2)+Cm(I+1)-2.d0*Cm(I))/(RPOW(I))
              DR(1,3)= DR(1,3) + Dm(I)*MMLR(I)*(Cm(I+2)-Cm(I+1))
     1                            /(RPOW(I)*RDIST)
     2                        -Dmp(I)*(Cm(I+2)-Cm(I+1))/(RPOW(I))
              DR(2,2)= DR(2,2) + Dm(I)*MMLR(I)*(Cm(I+2)+Cm(I+1)+
     1                           4.d0*Cm(I))/(6.d0*RPOW(I)*RDIST)
     2                        -Dmp(I)*(Cm(I+2)+Cm(I+1)+4.d0*Cm(I))
     3                            /(6.d0*RPOW(I))
              DR(3,3)= DR(3,3) + Dm(I)*MMLR(I)*(Cm(I+2)+Cm(I+1))
     1                            /(2.d0*RPOW(I)*RDIST)
     2                        -Dmp(I)*(Cm(I+2)+Cm(I+1))/(2.d0*RPOW(I)) 
              ENDDO
          DR(1,2) = DR(1,2)/(3.d0*DSQRT(2.d0))
          DR(2,1) = DR(1,2)
          DR(2,3) = DR(1,3)/(2.d0*DSQRT(3.d0))
          DR(1,3) = DR(1,3)/(DSQRT(6.d0))
          DR(3,1) = DR(1,3)
          DR(3,2) = DR(2,3)
c...      Partial derivatives of interaction matrix A  w.r.t.  Cm's
          DO  I= 2,NCMM,3 
              DMx(I,1,1)= -Dm(I)/(3.d0*RPOW(I))
              DMx(I+1,1,1)= DMx(I,1,1) 
              DMx(I+2,1,1)= DMx(I,1,1)
              DMx(I,1,2)= 2.d0*Dm(I)/(3.d0*DSQRT(2.d0)*RPOW(I))
              DMx(I+1,1,2)= -DMx(I,1,2)/2.d0
              DMx(I+2,1,2)= DMx(I+1,1,2)
              DMx(I,2,1)= DMx(I,1,2)
              DMx(I+1,2,1)= DMx(I+1,1,2)
              DMx(I+2,2,1)= DMx(I+2,1,2)
              DMx(I,1,3)= 0.d0
              DMx(I,3,1)= 0.d0
              DMx(I+1,1,3)= Dm(I)/(DSQRT(6.d0)*RPOW(I))
              DMx(I+1,3,1)= DMx(I+1,1,3)
              DMx(I+2,1,3)= -DMx(I+1,1,3)
              DMx(I+2,3,1)= DMx(I+2,1,3)
              DMx(I,2,2)= 2.d0*Dm(I)/(3.d0*RPOW(I))
              DMx(I+1,2,2)= DMx(I,2,2)/4.d0
              DMx(I+2,2,2)= DMx(I+1,2,2)
              DMx(I,2,3)= 0.d0
              DMx(I,3,2)= 0.d0
              DMx(I+1,2,3)= Dm(I)/(2.d0*DSQRT(3.d0)*RPOW(I))
              DMx(I+1,3,2)= DMx(I+1,2,3)
              DMx(I+2,2,3)= -DMx(I+1,2,3)
              DMx(I+2,3,2)= DMx(I+2,2,3)
              DMx(I,3,3)= 0.d0
              DMx(I+1,3,3)= Dm(I)/(2.d0*RPOW(I))
              DMx(I+2,3,3)= DMx(I+1,3,3)
              ENDDO
c...      Call subroutine to prepare and invert interaction matrix  A
          CALL DSYEVJ3(A,Q,W)
          L=1
c...      Now - identify the lowest eigenvalue of  A  and label it  L
          DO J=2,3
              IF (W(J) .LT. W(L)) THEN
                  L=J
                  ENDIF
              ENDDO
c...      Identifiy the highest eigenvalue of A and label it H
          H=1 
          DO J=2,3
              IF(W(J).GT.W(H)) THEN
                  H=J
                  ENDIF
              ENDDO
c...      Identify the middle eigenvalue of A and label it M
          M=1 
          DO J=2,3
              IF((J.NE.L).AND.(J.NE.H)) M= J
              ENDDO
c...      Select which eigenvalue to use based on user input
          IF(MMLR(1).EQ.-2) THEN 
              X = L
          ELSEIF(MMLR(1).EQ.-3) THEN
              X = M
          ELSE         
              X = H
              ENDIF
c...      determine ULR and eigenvectors
          ULR= -W(X)
          IF(MMLR(1).EQ.-2) ULR= ULR+ Cm(11)*R3*R6        !! C9adj term
          IF((MMLR(1).EQ.-3).OR.(MMLR(1).EQ.-4)) ULR = ULR + DELTAE
          IF(MMLR(1).EQ.-3) ULR= ULR+ Cm(12)*R3*R6        !! C9adj term
          IF(MMLR(1).EQ.-4) ULR= ULR+ Cm(13)*R3*R6        !! C9adj term
          DO I=1,3      
              EIGVEC(I,1) = Q(I,X)
              ENDDO 
cc  loop over values of m to determine partial derivatives w.r.t. each Cm
          DO I=2,NCMM
             DMtemp(1:3,1:3) = DMx(I,1:3,1:3) 
             DEIGMtemp= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DMtemp,EIGVEC))
             dULRdCm(I)= DEIGMtemp(1,1)
             ENDDO
          DEIGR = -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DR,EIGVEC))
          dULRdR= DEIGR(1,1)    !! radial derivative w.r.t. r (I think!)
c------------------------------------------------------------------------
          ENDIF
c------------------------------------------------------------------------
      RETURN
      CONTAINS
c=======================================================================
      SUBROUTINE DSYEVJ3(A, Q, W)
c ----------------------------------------------------------------------
c** Subroutine to setup and diagonalize the matrix  A  and return 
c   eigenvalues W and eigenvector matrix  Q
      INTEGER N, I, X, Y, R
      PARAMETER (N=3)
      REAL*8 A(3,3), Q(3,3), W(3)
      REAL*8 SD, SO, S, C, T, G, H, Z, THETA, THRESH
c     Initialize Q to the identitity matrix
c --- This loop can be omitted if only the eigenvalues are desired ---
      DO  X = 1, N
          Q(X,X) = 1.0D0
          DO  Y = 1, X-1
              Q(X, Y) = 0.0D0
              Q(Y, X) = 0.0D0
              ENDDO
          ENDDO
c Initialize W to diag(A)
      DO  X = 1, N
          W(X) = A(X, X)
          ENDDO
c Calculate SQR(tr(A))  
      SD = 0.0D0
      DO  X = 1, N
          SD = SD + ABS(W(X))
          ENDDO
      SD = SD**2
c Main iteration loop
      DO 40 I = 1, 50
c Test for convergence
          SO = 0.0D0
          DO  X = 1, N
              DO  Y = X+1, N
                  SO = SO + ABS(A(X, Y))
                  ENDDO
              ENDDO
          IF(SO .EQ. 0.0D0)  RETURN
          IF(I .LT. 4) THEN
              THRESH = 0.2D0 * SO / N**2
            ELSE
              THRESH = 0.0D0
            END IF
c Do sweep
          DO 60 X = 1, N
              DO 61 Y = X+1, N
                  G = 100.0D0 * ( ABS(A(X, Y)) )
                  IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X))
     1                          .AND. ABS(W(Y)) + G .EQ. ABS(W(Y))) THEN
                      A(X, Y) = 0.0D0
                    ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
c Calculate Jacobi transformation
                      H = W(Y) - W(X)
                      IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                          T = A(X, Y) / H
                        ELSE
                          THETA = 0.5D0 * H / A(X, Y)
                          IF (THETA .LT. 0.0D0) THEN
                              T= -1.0D0/(SQRT(1.0D0 + THETA**2)-THETA)
                            ELSE
                              T= 1.0D0/(SQRT(1.0D0 + THETA**2) + THETA)
                            END IF
                        END IF
                      C = 1.0D0 / SQRT( 1.0D0 + T**2 )
                      S = T * C
                      Z = T * A(X, Y)
c Apply Jacobi transformation
                      A(X, Y) = 0.0D0
                      W(X)    = W(X) - Z
                      W(Y)    = W(Y) + Z
                      DO  R = 1, X-1
                          T       = A(R, X)
                          A(R, X) = C * T - S * A(R, Y)
                          A(R, Y) = S * T + C * A(R, Y)
                          ENDDO
                      DO  R = X+1, Y-1
                          T       = A(X, R)
                          A(X, R) = C * T - S * A(R, Y)
                          A(R, Y) = S * T + C * A(R, Y)
                          ENDDO
                      DO  R = Y+1, N
                          T       = A(X, R)
                          A(X, R) = C * T - S * A(Y, R)
                          A(Y, R) = S * T + C * A(Y, R)
                          ENDDO
c Update eigenvectors
c --- This loop can be omitted if only the eigenvalues are desired ---
                      DO  R = 1, N
                          T       = Q(R, X)
                          Q(R, X) = C * T - S * Q(R, Y)
                          Q(R, Y) = S * T + C * Q(R, Y)
                          ENDDO
                    END IF
   61             CONTINUE
   60         CONTINUE
   40     CONTINUE
      WRITE(6,'("DSYEVJ3: No convergence.")')
      END SUBROUTINE DSYEVJ3
      END SUBROUTINE AFdiag
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE dampF(r,rhoAB,NCMM,NCMMAX,MMLR,sVRS2,IDSTT,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. r of the chosen form of the damping
c    function, for  m= 1 to MMAX.
c---------------------- RJL Version of 30 January 2016 ------------------
c-----------------------------------------------------------------------
c                 Upon Input
c* r - the radial distance in Angsroms (!) 
c* RHOab  'universal' scaling coefficient used for systems other than H_2
c       RHOab= 2*(RHOa*RHOb)/(RHOa+RHOb) where RHOa = (I_p^A/I_p^H)^0.66
c              where I_p^A is the ionization potential of atom A
c              and I_p^H is the ionization potential of atomic hydrogen
c* NCMM  the number of inverse-power terms to be considered
c* MMLR  are the powers of the NCMM inverse-power terms
c* sVRS2 defines damping s.th.  Dm(r)/r^m --> r^{sVRS2/2} as r --> 0
c* IDSTT specifies damping function type:  > 0  use Douketis et al. form 
c                               if  IDSTT .LE. 0  use Tang-Toennies form
c-----------------------------------------------------------------------
c                 Upon Output
c  DM(m) - The value of the damping function for the long range term 
c          C_MMLR(m)/r^MMLR(m)    {m= 1, NCMM}
c  DMP(m): first derivative of the damping function  DM(m) w.r.t. r
c  DMPP(m): second derivative of the damping function  DM(m) w.r.t. r
c  IF(rhoAB.LE.0.0) return w. DM(m)= 1.0 & DMP(m)=DMPP(m)=0.0 for all m
c-----------------------------------------------------------------------
      INTEGER NCMM,NCMMAX,MMLR(NCMMAX),sVRS2,IDSTT,sVRS2F,FIRST, Lsr,m,
     1  MM,MMAX,MMTEMP
      REAL*8 r,rhoAB,bTT(-2:2),cDS(-4:4),bDS(-4:4),aTT,br,XP,YP,
     1  TK, DM(NCMMAX),DMP(NCMMAX),DMPP(NCMMAX),SM(-3:25),
     2  bpm(20,-4:0), cpm(20,-4:0),ZK
c------------------------------------------------------------------------
c  The following values for the numerical factors used in both TT and DS
c  were  normalized to the Hydrogen data presented
c  by Kreek and Meath in J.Chem.Phys. 50, 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c------------------------------------------------------------------------
      DATA bTT/2.10d0,2.44d0,2.78d0,3.13d0,3.47d0/
      DATA bDS/2.50d0,2.90d0,3.30d0,3.69d0,3.95d0,0.d0,4.53d0,0.d0,
     1         4.99d0/
      DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0,0.d0,0.360d0,
     1           0.d0,0.340d0/
      DATA FIRST/ 1/
      SAVE FIRST, bpm, cpm
c-----------------------------------------------------------------------
      MMTEMP = MMLR(1)
      IF(MMLR(1).LE.0) MMLR(1) = 1
      IF(RHOab.LE.0) THEN
          DO  m=1,NCMMax
              DM(m)=1.d0
              DMP(m)= 0.d0
              DMPP(m)= 0.d0
              ENDDO
          RETURN
          ENDIF
      IF(IDSTT.LE.0) THEN
c===========================================
c** For Tang-Toennies type damping functions
c===========================================
          Lsr= sVRS2/2
          IF((sVRS2.LT.-4).OR.(sVRS2.GT.4).OR.((2*LSR).NE.sVRS2)) THEN
                WRITE(6,600) 'TT',sVRS2
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
                  DMP(m)= aTT*XP*(SM(MM) - SM(MM-1))
                  DMPP(m)= -aTT*aTT*XP*(SM(MM) 
     1                                     - 2.d0*SM(MM-1) + SM(MM-2))
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
                  DMP(m)= aTT*XP*SM(m-1)
                  DMPP(m)= aTT*aTT*XP*(SM(m-2)-SM(m-1))
                  ENDDO
            ENDIF
          ENDIF
c
      IF(IDSTT.GT.0) THEN
c=======================================================================
c** For Douketis-Scoles-Marchetti-Zen-Thakkar type damping function ...
c=======================================================================
          IF((sVRS2.LT.-4).OR.(sVRS2.GT.4).OR.(sVRS2.EQ.1).OR.
     1                                              (sVRS2.EQ.3)) THEN
              WRITE(6,600) 'DS',sVRS2
              STOP
              ENDIF
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  sVRS2F= -4,0
                      bpm(m,sVRS2F)= bDS(sVRS2F)/DFLOAT(m)
                      cpm(m,sVRS2F)= cDS(sVRS2F)/DSQRT(DFLOAT(m))
                      ENDDO
                  ENDDO
              FIRST= 0 
              ENDIF
          br= rhoAB*r
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,sVRS2) + cpm(MM,sVRS2)*br)*br)
              YP= 1.d0 - XP
              ZK= MM + 0.5d0*sVRS2
              DM(m)= YP**ZK
              TK= (bpm(MM,sVRS2) + 2.d0*cpm(MM,sVRS2)*br)*rhoAB
              DMP(m) = ZK*XP*TK*DM(m)/YP
c ... calculate second derivative [for DELR case] {check this!}
              DMPP(m)= (ZK-1.d0)*DMP(m)*(XP*TK)/YP
     1               - DMP(m)*TK + DMP(m)*2.d0*cpm(MM,sVRS2)*rhoAB**2/TK
              ENDDO   
          ENDIF  
      MMLR(1) = MMTEMP
      RETURN
  600 FORMAT(/,' *** ERROR ***  For  ',A2,'-damping functions not yet de
     1fined for   sVRS2=',i3)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

