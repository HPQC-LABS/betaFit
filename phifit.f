c***********************************************************************
c************  Program  PHIFIT_1.2  dated  25 April 2007  **************
c***********************************************************************
c* Program to fit NTP read-in potential fx. values {RTP(i),VTP(i)} to
c  a chosen analytic form, to determine realistic initial estimates of 
c  exponent expansion coefficients  phi_i  , for other purposes.
c** See  http://leroy.uwaterloo.ca/programs/   for further documentation
c***********************************************************************
      INTEGER MXDATA, MXPARM, MXMLR
      PARAMETER (MXDATA=1501, MXPARM=30, MXMLR= 8)
      INTEGER i,j,ITER,IROUND,ROBUST,LPRINT,IWR,NPARM,NTP,IFXP(MXPARM)
      REAL*8 PHI(0:MXPARM),PV(MXPARM),PU(MXPARM),PS(MXPARM),
     1  CM(MXPARM,MXPARM),DYDP(MXDATA,MXPARM),VTP(MXDATA),
     2  uVTP(MXDATA),phiy(MXDATA),Uphiy(MXDATA),YD(MXDATA),
     3  phiINF,UNC,yPOW,DSE,TSTPS,TSTPU,DSEB,TT(0:20),RHOdR,RHOp,TTM,
     4  Rep,AREF,AREFp,RTPp, AA,BB,VLR,dVLR,FCT,RAT,UMAX,
     5  yp,fsw,ypRE,ReDE, ReIN,DeIN,VMINin ,VLRe,RE3,T1,RTP3
      CHARACTER*4  NNAME,NAME(5)
      DATA NAME/' EMO',' MLJ',' MLR','DELR','GPEF'/
c-----------------------------------------------------------------------
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR(MXMLR),p,NS,NL,
     1                                                            NPHI
      REAL*8 Re,De,VMIN,RREF,Asw,Rsw,M2,ASO,R01,R12,as,bs,RHOd,
     1  CmVAL(MXMLR),RTP(MXDATA)
      COMMON /DATABLK/Re,De,VMIN,RREF,Asw,Rsw,M2,ASO,R01,R12,as,bs,
     1  RHOd,CmVAL,RTP,PSEL,IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR,p,NS,NL,
     2                                                            NPHI
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
c   of exponent coefft  phi(y)= [1-fsw] phi_inf + fsw \sum {phi_i y^i} 
c   read positive values of Asw and Rsw; otherwise, set them .LE. 0 and
c        phi(r)= yp phi_inf + [1 - yp] Sum{ phi_i yp^i }
c** For each long-range term read power  MMLR(i)  & coefficient CmVAL(i)
c** For special Aubert-Frecon case,  NCMM= 4,  MMLR= {3,0,6,6} and the 
c  coefficients are:  CmVAL(1)= M^2, CmVAL(2)= ASO, CmVAL(3)= R10, and
c  CmVAL(4)= R12
c=======================================================================
      IF(PSEL.EQ.2) THEN
          READ(5,*) NCMM, Asw, Rsw
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
c Also ...  PHI(0) is the initial trial value of  \phi_0  used to
c    generate initial trial values of the A & B potential parameters.
c  * If the read-in  PHI(0) \leq 0.0 , the program uses a preliminary 
c    EMO_{p} fit to generate an estimate of  PHI(0); this usually works
c    if  u_{LR} is attractive, so that the potential has no barrier.
c  * If the potential has a barrier, one must determine an initial trial
c    value another way.  A one way would be to do an EMO_{p} fit to the
c    potential points, treating the barrier maximum as dissociation
c=======================================================================
          READ(5,*) NCMM, RHOd, IDF, PHI(0)
          READ(5,*) (MMLR(i), CmVAL(i), i= 1, NCMM)
c=======================================================================
          PHI(0)= 1.d0
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
          IF((NCMM.EQ.4).AND.(MMLR(2).EQ.0)) THEN
c** For Lyon treatment of A-state alkali dimers ...
              WRITE(6,618)
              M2= CmVAL(1)
              ASO= CmVAL(2)
              R01= CmVAL(3)
              R12= CmVAL(4)
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
      IF((NCMM.EQ.4).AND.(MMLR(2).EQ.0)) THEN
          DO  I= 1,NTP
              VTP(I)= VTP(I) - 0.5d0*ASO
              ENDDO
          ENDIF
      IF((PSEL.EQ.2).AND.(NCMM.EQ.4).AND.(MMLR(2).EQ.0)) THEN
          ENDIF
c
  600 FORMAT(' Determine  EMOp  exponent expansion coefficients'/
     1  1x,24('==')/' Start with   Re=',f11.8,'   De=',f11.4)
  601 FORMAT(' Using exponent expansion variable  yp(r)= [r^p -',f6.2,
     `  '^p]/[r^p +',f6.2,'^p]' )
  602 FORMAT(' Determine ',A4,'p  exponent expansion coefficients'/1
     1  x,24('==')/' Start with   Re=',f11.8,'   De=',f11.4,'    C',
     2  i2,'=',1PD15.8:/(49x,'C',i2,'=',D15.8:))
  618 FORMAT(4x,'Use Lyon  uLR(r) with   C_0= ASO   C_6(1)= R01',
     1  '   C_6(2)= R12')
  603 FORMAT(' Use exponent expansion variable  yp(r)= [r^p - Re^p]/[r^p
     1 + Re^p]' )
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
c** Now ... loop over different {p,NS,NL} combinations till end of data
c*  p  is power in expansion variable  yp=(R^p - AREF^p)/(R^p + AREF^p)
c** Read powers NS used in  \phi(y)  expansion for  r > Re  and 
c               NL used in  \phi(y)  expansion for  r = Re  or  > Re 
c** For a GPEF potential, consider powers ranging from  NS(.ge.0) to NL
c* RREF   defines the reference distance in the expansion variable
c      - for  RREF.le.0 , define parameter  RREF = Re
c      - for  RREF.gt.0 , fix parameter  RREF   at its read-in value
c-----------------------------------------------------------------------
   10 READ(5,*, END= 999) p, NS, NL, RREF
c-----------------------------------------------------------------------
      Re= ReIN
      De= DeIN
      VMIN= VMINin
      IF((p.LE.0).OR.(NS.LE.0).OR.(NL.LE.0)) GOTO 999
      IF(PSEL.EQ.1) WRITE(6,600) Re, De
      IF(PSEL.EQ.2) THEN
          WRITE(6,602) NNAME, Re, De, (MMLR(i),CmVAL(i),i= 1,NCMM)
          IF(Asw.GT.0.d0) WRITE(6,632) Asw, Rsw
          IF(Asw.LE.0.d0) WRITE(6,634) 
          ENDIF
      IF(PSEL.EQ.3) THEN
          WRITE(6,604) Re, De, RHOp, NCMM,(MMLR(j),CmVAL(j),j=1, NCMM)
          IF(IDF.EQ.1) WRITE(6,614)
          IF(IDF.EQ.2) WRITE(6,612)
          ENDIF
      IF(PSEL.LE.3) THEN
          IF(RREF.gt.0.d0) THEN
              AREF= RREF
              WRITE(6,601) AREF,AREF
              ENDIF
          IF(RREF.LE.0.d0) THEN
              AREF= Re
              WRITE(6,603)
              ENDIF
          AREFp= AREF**p
          ENDIF
      NPHI= MAX(NS,NL)+ 1
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
      IF((PSEL.EQ.1).OR.(PSEL.EQ.3).AND.(phi(0).le.0)) THEN
c ... first define ordinate array
          NNAME= NAME(1)
          DO  i= 1,NTP
              RTPp= RTP(i)**p
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              IF(RTP(i).GT.Re) THEN
                  phiy(i)= - DLOG(1.d0 - DSQRT((VTP(i)-VMIN)/De))
                  IF(VTP(i).GT.UNC) THEN
                      Uphiy(i)= 0.5d0*UNC/(DSQRT((VTP(i)-VMIN)*De) - 
     1                                                  (VTP(i)-VMIN))
                    ELSE
                      Uphiy(i)= DSQRT(UNC/De)
                    ENDIF
                ELSE
                  phiy(i)= - DLOG(1.d0 + DSQRT((VTP(i)-VMIN)/De))
                  IF(VTP(i).GT.UNC) THEN
                      Uphiy(i)= 0.5d0*UNC/(DSQRT((VTP(i)-VMIN)*De) +
     1                                                  (VTP(i)-VMIN))
                    ELSE
                      Uphiy(i)= DSQRT(UNC/De)
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
cc        write(8,702) rtp(i),yp,vtp(i),phiy(i),Uphiy(i)
cc   1                                 ,(dydp(i,j),j=1,nphi)
cc700 format('  RTP     yp      VTP       phi*y      unc(phi*y) :',
cc   1  ' {dY/dp}')
cc702 format(f6.3,f8.4,f9.2,1P2d13.5:/(14x,5d13.5))
c%%
              ENDDO
          NPARM= NPHI
          CALL LLSQF(NTP,NPARM,MXDATA,MXPARM,phiy,Uphiy,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
          IF(IWR.GT.0) WRITE(6,620) NNAME,p,NS,NL,DSE,
     1                        ('phi',j-1,PV(j),PU(j),PS(j),j= 1,NPARM)
          ENDIF

c=======================================================================
      IF(PSEL.EQ.2) THEN
c-----------------------------------------------------------------------
c*** Preliminary linearized fit for an  MLRp  potential ...
c** First define array of exponent values with uncertainties defined by 
c  the assumption that all potential values have equal uncertainties UNC
          NNAME= NAME(3)
          IF((NCMM.EQ.4).AND.(MMLR(2).EQ.0)) THEN
c** Aubert-Frecon based VLR(r)
              RE3= Re**3
              T1= M2/(9.d0*Re3) + (5.d0*R01 + R12)/(45.d0*RE3**2)
              VLRe= 0.5d0*(M2/RE3 - ASO) + (5.d0*R01 + 8.2d0*R12)
     1    /(18.d0*RE3**2) + 0.5d0*DSQRT((T1- ASO)**2 + 8.d0*T1**2)
              WRITE(6,618) ASO,R01,R12
            ELSE
c** For normal inverse-power sum MLR/MLJ case
              IF(NCMM.EQ.1) NNAME= NAME(2)
              IF(p.LE.(MMLR(NCMM)-MMLR(1))) THEN
                  WRITE(6,616) p, NCMM,MMLR(NCMM)-MMLR(1)
                  GOTO 10
                  ENDIF
              VLRe= 0.d0
              DO  i= 1,NCMM
                  VLRe= VLRe + CmVAL(i)/Re**MMLR(i)
                  ENDDO
            ENDIF
          phiINF= DLOG(2.d0*De/VLRe)
          WRITE(6,619) phiINF
          Rep= RE**p
          DO  i= 1, NTP
              RTPp= RTP(i)**p
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              ypRE= (RTPp - Rep)/(RTPp + Rep)
              IF((NCMM.EQ.4).AND.(MMLR(2).EQ.0)) THEN
c... for Aubert-Frecon {3,0,6,6} case ...
                  RTP3= RTP(i)**3
                  T1= M2/(9.d0*RTP3)+ (5.d0*R01+R12)/(45.d0*RTP3**2)
                  VLR= 0.5d0*(M2/RTP3 - ASO) + (5.d0*R01+ 8.2d0*R12)
     1   /(18.d0*RTP3**2) + 0.5d0*DSQRT((T1- ASO)**2 + 8.d0*T1**2)
                ELSE
c... for normal MLR/MLJ case ...
                  VLR= 0.d0
                  DO  j= 1, NCMM
                      VLR= VLR+ CmVAL(j)/RTP(i)**MMLR(j)
                      ENDDO
                ENDIF
              IF(RTP(i).GT.Re) THEN
                  phiy(i)= - DLOG((1.d0 - DSQRT((VTP(i)-VMIN)/De))
     1                                                      *VLRe/VLR)
                  IF((VTP(i)-VMIN).GT.UNC) THEN
                      Uphiy(i)= 0.5d0*UNC
     1                      /(DSQRT((VTP(i)-VMIN)*De) - (VTP(i)-VMIN))
                    ELSE
                      Uphiy(i)= DSQRT(UNC/De)
                    ENDIF
                ELSE
                  phiy(i)= - DLOG((1.d0 + DSQRT((VTP(i)-VMIN)/De))
     1                                                      *VLRe/VLR)
                  IF((VTP(i)-VMIN).GT.UNC) THEN
                      Uphiy(i)= 0.5d0*UNC
     1                      /(DSQRT((VTP(i)-VMIN)*De) + (VTP(i)-VMIN))
                    ELSE
                      Uphiy(i)= DSQRT(UNC/De)
                    ENDIF
                ENDIF
c** Subtract the \phi_\infty term to yield polynomial for fitting
              IF(Asw.LE.0.d0) THEN
c... For Huang's MLR exponent function
                  phiy(i)= phiy(i)- phiINF*yp*ypRE
                  yPOW= ypRE*(1.d0- yp)
                ELSE
c... For Photos' origonal MLJ exponent switching function
                  fsw= 1.d0/(DEXP(Asw*(RTP(i)- Rsw)) + 1.d0)
                  phiy(i)= phiy(i)- phiINF*ypRE*(1.d0 - fsw)
                  yPOW= ypRE*fsw
                  ENDIF
              uVTP(i)= UNC
c... then create partial derivative array for linearized fit ...
              DO  j= 1, NPHI
                  DYDP(i,j)= yPOW
                  IF((RTP(i).GT.Re).AND.(j.GT.NL+1)) DYDP(i,j)= 0.d0
                  IF((RTP(i).LE.Re).AND.(j.GT.NS+1)) DYDP(i,j)= 0.d0
                  yPOW= yPOW*yp
                  ENDDO
c%%%
cc    if(i.eq.1) write(8,700) 
cc            write(8,702) rtp(i),yp,vtp(i),phiy(i),Uphiy(i)
cc   1                                 ,(dydp(i,j),j=1,nphi)
c             write(8,800) rtp(i),yp,ypRE,vlr,phiy(i)
c 800 Format( f7.4,2f12.8,4(1Pd15.7))
c%%%
              ENDDO
          CALL LLSQF(NTP,NPHI,MXDATA,MXPARM,phiy,Uphiy,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
          IF(IWR.GT.0) WRITE(6,620) NNAME,p,NS,NL,DSE,
     1                        ('phi',j-1,PV(j),PU(j),PS(j),j= 1,NPHI)
          ENDIF

c=======================================================================
c** Preliminary linearized fit for a  DELR  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.3) THEN
c ... NOTE need iteration to determine self-consistent  phi(0)  value
c  First generate  A & B  from input Re, De and trial phi(0)
          NNAME= NAME(4)
          ITER= 0
   40     VLR= 0.d0
          dVLR= 0.d0
          IF(IDF.EQ.2) THEN
c... Using Scoles-type damping function ...
              DO  j= 1,NCMM
                  FCT= dexp(-3.97d0*RHOd*Re/MMLR(j)
     1                    - 0.39d0*(RHOd*Re)**2/DSQRT(DFLOAT(MMLR(j))))
                  AA= CmVAL(j)*((1.d0- FCT)/Re)**MMLR(j)
                  VLR= VLR+ AA
                  dVLR= dVLR- MMLR(j)*AA/Re
                  BB= FCT*(3.97d0*RHOd/MMLR(j) 
     1                      + 0.78d0*Re*RHOd**2/DSQRT(DFLOAT(MMLR(j))))
                  dVLR= dVLR+  MMLR(j)*BB*AA/(1.d0- FCT)
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
              VLR= 0.d0
              dVLR= 0.d0
              DO  j= 1,NCMM
                  TTM= (1.d0- yPOW*TT(MMLR(j)))*CmVAL(j)/Re**MMLR(j)
                  VLR= VLR+ TTM
		      dVLR= dVLR+ yPOW*RHOd*(TT(MMLR(j)) - TT(MMLR(j)-1))
     1                          *CmVAL(j)/Re**MMLR(j) - MMLR(j)*TTM/Re
                  ENDDO
              ENDIF
          AA= De - VLR - dVLR/phi(0)
          BB= 2.d0*(De - VLR) - dVLR/phi(0)
          WRITE(6,607) AA,BB
          RAT= 0.5d0*BB/AA
          UMAX= DSQRT(RAT**2 + (UNC + VLR - DE)/AA)
          ReDE= Re- dlog(RAT)/phi(0)
          DO  i= 1,NTP
              VLR= 0.d0
              RTPp= RTP(i)**p
              yp= (RTPp - AREFp)/(RTPp + AREFp)
              IF(IDF.EQ.2) THEN
c... if using Scoles damping function ...
                  DO  j= 1,NCMM
                      FCT= dexp(-3.97d0*RHOd*RTP(i)/MMLR(j)
     1               - 0.39d0*(RHOd*RTP(i))**2/DSQRT(DFLOAT(MMLR(j))))
                      VLR= VLR+ CmVAL(j)*((1.d0-FCT)/RTP(i))**MMLR(j)
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
                      VLR= VLR+CmVAL(j)*(1.d0-TT(MMLR(j)))/
     1                                                 RTP(i)**MMLR(j)
                      ENDDO
                  ENDIF
              FCT= (VTP(i) - VMIN + VLR - De)/AA + RAT**2
              IF(FCT.LT.0.d0) THEN
c** If estimate of ReDE off a bit and  FCT < 0 , ignore & deweight point
                    phiy(i)= 0.d0
                    Uphiy(i)= 9.d99
                    GOTO 44
                    ENDIF
              FCT= DSQRT(FCT)
              IF(RTP(i).GT.ReDE) THEN
                  IF(RAT.GT.FCT) THEN
                      phiy(i)= - DLOG(RAT - FCT)
                      IF((VTP(i)-VMIN).GT.UNC) THEN
                          Uphiy(i)= 0.5d0*UNC/(AA*(RAT- FCT)*FCT)
                        ELSE
                          Uphiy(i)= UNC/(AA*UMAX*(RAT- UMAX))
                        ENDIF
                    ELSE
c ... deweight away points for which \ln argument would be negative
                      phiy(i)= 0.d0
                      Uphiy(i)= 9.d9
                    ENDIF
                ELSE
                  phiy(i)= - DLOG(RAT + FCT)
                  IF((VTP(i)-VMIN).GT.UNC) THEN
                      Uphiy(i)= 0.5d0*UNC/(AA*(RAT+ FCT)*FCT)
                    ELSE
                      Uphiy(i)= UNC/(AA*UMAX*(RAT+ UMAX))
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
          CALL LLSQF(NTP,NPARM,MXDATA,MXPARM,phiy,Uphiy,DYDP,YD,PV,
     1           PU,PS,CM,DSE)
          IF(DABS(PV(1)-PHI(0)).GT.PS(1)) THEN
              WRITE(6,644) phi(0),PV(1),PV(1)-phi(0),DSE
              phi(0)= PV(1)
              ITER= ITER+ 1
              IF(ITER.LE.10) GOTO 40
              WRITE(6,646) ITER
            ELSE
              WRITE(6,648) PV(1),PV(1)-phi(0)
            ENDIF
          IF(IWR.GT.0) WRITE(6,620) NNAME,p,NS,NL,DSE,
     1                        ('phi',j-1,PV(j),PU(j),PS(j),j= 1,NPARM)
          ENDIF 

      IF(PSEL.LE.3)  THEN
c======================================================================
c** Now ... do direct non-linear fit to potential values ... first with
c      Re and/or VMIN and De fixed, and then freeing them up too ...
c=======================================================================
c* FIRST optimize  PHI(j)'s (and VMIN) with  Re and De held fixed!
          DO  j= 1,NPHI
              PHI(j-1)= PV(j)
              IFXP(j)= 0
              ENDDO
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
              IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,AREF,NS,NL,DSE,
     1                                (j-1,PV(j),PU(j),PS(j),j=1,NPHI)
              IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,NS,NL,DSE,
     1                                (j-1,PV(j),PU(j),PS(j),j=1,NPHI)
              IF(IFXVMIN.LE.0)
     1                   WRITE(6,660) PV(NPHI+3),PU(NPHI+3),PS(NPHI+3)
              ENDIF
c ... the, if appropriate, set  Re  free too ...
          IF(IFXRe.LE.0) THEN
              IFXP(NPHI+1)= 0
              CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,LPRINT,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              Re= PV(NPHI+1)
              IF(IWR.GE.1) THEN
                  IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,AREF,NS,NL,DSE,
     1                               (j-1,PV(j),PU(j),PS(j),j=1,NPHI)
                  IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,NS,NL,DSE,
     1                               (j-1,PV(j),PU(j),PS(j),j=1,NPHI)
                  WRITE(6,662) PV(NPHI+1),PU(NPHI+1),PS(NPHI+1)
                  IF(IFXVMIN.LE.0)
     1                   WRITE(6,660) PV(NPHI+3),PU(NPHI+3),PS(NPHI+3)
                  ENDIF
              ENDIF
c ... then with Re fixed again, free De & VMIN (as well as the phi's)
          IF(IFXDe.LE.0) THEN
              DSEB= DSE
              IFXP(NPHI+1)= 1
              IF(IFXDe.LE.0) IFXP(NPHI+2)= 0
              IF(IFXVMIN.LE.0) IFXP(NPHI+3)= 1
              CALL NLLSSRR(NTP,NPARM,MXPARM,IROUND,ROBUST,LPRINT,IFXP,
     1                        VTP,uVTP,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
              IF(IFXDe.LE.0) De= PV(NPHI+2)
              IF(IFXVMIN.LE.0) VMIN= PV(NPHI+3)
              IF((IWR.GE.1).OR.(DSE.GT.DSEB*1.01)) THEN
                  IF(DSE.GT.DSEB*1.01) WRITE(6,654) DSEB,DSE
                  IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,AREF,NS,NL,DSE,
     1                               (j-1,PV(j),PU(j),PS(j),j=1,NPHI)
                  IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,NS,NL,DSE,
     1                               (j-1,PV(j),PU(j),PS(j),j=1,NPHI)
                  IF(IFXRe.LE.0) 
     1                WRITE(6,662) PV(NPHI+1),PU(NPHI+1),PS(NPHI+1)
                  IF(IFXDe.LE.0)
     1                WRITE(6,628) PV(NPHI+2),PU(NPHI+2),PS(NPHI+2)
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
          IF(RREF.GT.0.d0) WRITE(6,622) NNAME,p,AREF,NS,NL,DSE,
     1                               (j-1,PV(j),PU(j),PS(j),j=1,NPHI)
          IF(RREF.LE.0.d0) WRITE(6,624) NNAME,p,NS,NL,DSE,
     1                               (j-1,PV(j),PU(j),PS(j),j=1,NPHI)
          IF(PSEL.EQ.3) PHI(0)= PV(1)
          WRITE(6,662) PV(NPHI+1),PU(NPHI+1),PS(NPHI+1)
          WRITE(6,628) PV(NPHI+2),PU(NPHI+2),PS(NPHI+2)
          WRITE(6,660) PV(NPHI+3),PU(NPHI+3),PS(NPHI+3)
          IF(IFXRe.LE.0) Re= PV(NPHI+1)
          IF(IFXDe.LE.0) De= PV(NPHI+2)
          IF(IFXVMIN.LE.0) VMIN= PV(NPHI+3)
          IF(PSEL.EQ.3) phi(0)= PV(1)
          ENDIF

c=======================================================================
c*** For case of a  GPEF  potential ...
c-----------------------------------------------------------------------
      IF(PSEL.EQ.4) THEN
          NNAME= NAME(5)
          DO  i= 1,NTP
              uVTP(i)= UNC
              phiy(i)= VTP(i)
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
              CALL LLSQF(NTP,NPARM,MXDATA,MXPARM,phiy,uVTP,DYDP,YD,PV,
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
      WRITE(6,608)
      GOTO 10
c-----------------------------------------------------------------------
  616 FORMAT(//' *** Since  p=',i2,' .LE. [MMLR(',i1,')-MMLR(1)]=',i2,
     1  '  STOP !!!!!! ******')
  619 FORMAT(' Linearized fit uses    phi(INF)=',f10.6)
  620 FORMAT(/' Linearized ',A4,'{p=',i1,'} fit with   NS=',i2,'   NL=',
     1 i2,'   yields   DSE=',1Pd9.2/(4x,a3,'_{',i2,'} =',d17.9,
     2 ' (+/-',d8.1,')   PS=',d8.1))
  622 FORMAT(/' Direct fit to ',A4,'{p=',i1,'; Rref=',f5.2,' ; NS=',i2,
     1  ', NL=',I2,'}  potential:   DSE=',1Pd9.2/
     2  ('    phi_{',i2,'} =',d17.9,' (+/-',d8.1,')   PS=',d8.1))
  624 FORMAT(/' Direct fit to ',A4,'{p=',i1,'; Rref= Re ; NS=',i2,
     1  ', NL=',I2,'}  potential:   DSE=',1Pd9.2/
     2  ('    phi_{',i2,'} =',d17.9,' (+/-',d8.1,')   PS=',d8.1))
  628 FORMAT(10x,'De =',f13.6,' (+/-',f12.6,')   PS=',1pd8.1)
  632 FORMAT(' Use exponent switching function with   Asw=',F9.6,
     1  '   Rsw=',F9.6)
  634 FORMAT(' Use Huang exponent function:  phi(R)= phiINF*y_p + (1-y_p
     1)* Sum{phi_i*[y_p]^i}')
  636 FORMAT(/' First perform full non-linear GPEF fit without taking ou
     1t common factor of c_0')
  644 FORMAT('   Update  phi_0  from',f9.6,'   to',f9.6,'   by',
     1  1Pd9.1,' :   DSE=',1PD8.1)
  646 FORMAT(' !!! CAUTION !!! Iteration to optimize  phi(0)  not conver
     1ged after',i3,' tries')
  648 FORMAT('   Converge on   phi_0=',f9.6,'   Next change=',1Pd9.1)
  654 FORMAT(/' *** PROBLEM *** freeing De makes DSE increase from',
     1  1PD9.2,' to',D9.2)
  658 FORMAT(/' Fit to GPEF{p=',i1,'} potential with   As=',F5.2,
     1  '   Bs=',F5.2,'   yields   DSE=',1Pd9.2/
     2  (6x,'c_{',i2,'} =',d17.9,' (+/-',d8.1,')   PS=',d8.1))
  660 FORMAT(8x,'VMIN =',f13.7,' (+/-',f12.6,')   PS=',1pd8.1)
  662 FORMAT(10x,'Re =',f13.9,' (+/-',f12.9,')   PS=',1pd8.1)
  999 STOP
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DYIDPJ(IDAT,NDATA,NPARM,IFXP,YC,PV,PD,PS,RMSR)
      INTEGER MXDATA, MXPARM, MXMLR
      PARAMETER (MXDATA=1501, MXPARM=30, MXMLR= 8)
      INTEGER  j,IDAT, NPOW,NPARM,NDATA, IFXP(MXPARM),JFXRe,JFXDe,
     1  JFXVMIN
      REAL*8  YC,PV(NPARM),PD(NPARM),PS(NPARM),TT(0:20),RHOdR,RMSR,RTPp,
     1  Rep,AREF,AREFp,ype,dype,phiINF,yp,fsw,yPOW,XP,XPW,DER,TTM,TTMM,
     2  DERP,SUM,DSUM,AA,BB,FCT,VLR,VLRe,dVLRe,d2VLRe,VCN,DDER, T0,T1,
     3  RE3,RTP3,dVLRedRe
c-----------------------------------------------------------------------
      INTEGER PSEL,IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR(MXMLR),p,NS,NL,
     1                                                            NPHI
      REAL*8 Re,De,VMIN,RREF,Asw,Rsw,M2,ASO,R01,R12,as,bs,RHOd,
     1  CmVAL(MXMLR),RTP(MXDATA)
      COMMON /DATABLK/Re,De,VMIN,RREF,Asw,Rsw,M2,ASO,R01,R12,as,bs,
     1  RHOd,CmVAL,RTP,PSEL,IFXRe,IFXDe,IFXVMIN,IDF,NCMM,MMLR,p,NS,NL,
     2                                                            NPHI
c-----------------------------------------------------------------------
      SAVE JFXRe,JFXDe,JFXVMIN, AREF,AREFp,Rep,phiINF,AA,BB, VLRe,
     1  dVLRedRe
c=======================================================================
c** NOTE BENE(!!) for non-linear fits, need to be sure that the
c  calculations of YC and PD(j) are based on the current UPDATED PV(j)
c  values.  If other (than PV) parameter labels are used internally
c  in the calculations, UPDATE them whenever (say)  IDAT = 1 .
      IF(IDAT.EQ.1) THEN
          JFXRe= IFXP(NPHI+1)
          JFXDe= IFXP(NPHI+2)
          JFXVMIN= IFXP(NPHI+3)
          ENDIF
      DO  j=1,NPARM
          PD(j)= 0.d0
          ENDDO
c=======================================================================
c** For case of an  EMO_p  potential
c-----------------------------------------------------------------------
      IF(PSEL.EQ.1) THEN
          IF(IDAT.EQ.1) THEN
              IF(JFXRe.LE.0) Re= PV(NPHI+1)
              IF(JFXDe.LE.0) De= PV(NPHI+2)
              IF(JFXVMIN.LE.0) VMIN= PV(NPHI+3)
              AREF= RREF
              IF(RREF.LE.0) AREF= Re
              AREFp= AREF**p
              ENDIF
          RTPp= RTP(IDAT)**p
          yp= (RTPp - AREFp)/(RTPp + AREFp)
          NPOW= NS+1
          IF(RTP(IDAT).GE.Re) NPOW= NL+1
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
          XP= DEXP(-SUM*(RTP(IDAT)- Re))
          YC= De*(1.d0 - XP)**2 + VMIN
          DER= 2.d0*De*(1.d0- XP)*XP
          DERP= DER*(RTP(IDAT)- Re)
          DO  j= 1,NPOW
              PD(j)= DERP
              DERP= DERP*yp
              ENDDO
c** If appropriate, also get partial derivative w.r.t. Re
          IF(JFXRE.LE.0) THEN
              IF(RREF.LE.0.d0) SUM= SUM +(RTP(IDAT)- Re)*DSUM*0.5d0*
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
          IF(IDAT.EQ.1) THEN
              IF(JFXRe.LE.0) Re= PV(NPHI+1)
              IF(JFXDe.LE.0) De= PV(NPHI+2)
              IF(JFXVMIN.LE.0) VMIN= PV(NPHI+3)
              AREF= RREF
              IF(RREF.LE.0.d0) AREF= Re
              AREFp= AREF**p
              Rep= Re**p
              IF((NCMM.EQ.4).AND.(MMLR(2).EQ.0)) THEN
c** For Aubert-Frecon based  uLR(r)
                  RE3= Re**3
                  T1= M2/(9.d0*Re3)+ (5.d0*R01+ R12)/(45.d0*RE3**2)
                  T0= DSQRT((T1- ASO)**2 + 8.d0*T1**2)
                  VLRe= 0.5d0*(M2/RE3 - ASO) + 0.5d0*T0  
     1                         + (5.d0*R01 + 8.2d0*R12)/(18.d0*RE3**2)
                  dVLRedRe= (-1.5d0*M2 - (5.d0*R01 + 8.2d0*R12)
     1                        /(3.d0*RE3) - 0.5d0*((9.d0*T1- ASO)/T0)*
     2        (M2/3.d0 + (10.d0*R01 + 2.d0*R12)/(15.d0*Re3)))/(Re3*Re)
                ELSE
c** For normal inverse-power sum  uLR(r)  ....
                  VLRe= 0.d0
                  dVLRedRe= 0.d0
                  DO  j= 1,NCMM
                      AA= CmVAL(j)/Re**MMLR(j)                  
                      VLRe= VLRe+ AA
                      dVLRedRe= dVLRedRe - MMLR(j)*AA/Re
                      ENDDO
                ENDIF
              phiINF= DLOG(2.d0*De/VLRe)
              ENDIF
          RTPp= RTP(IDAT)**p
          yp= (RTPp - AREFp)/(RTPp + AREFp)
          ype= (RTPp - Rep)/(RTPp + Rep)
          NPOW= NS+1
          IF(RTP(IDAT).GE.Re) NPOW= NL+1
          IF(Asw.LE.0.d0) THEN
              yPOW= 1.d0 - yp
            ELSE
              fsw= 1.d0/(dexp(Asw*(RTP(IDAT)-Rsw)) + 1.d0)
              yPOW= fsw
            ENDIF
          SUM= PV(1)*yPOW
          DSUM= 0.d0
          IF(NPOW.GE.2) THEN
              DO  j= 2,NPOW
                  IF(RREF.LE.0.d0) DSUM= DSUM + PV(j)*(j-1)*yPOW
                  yPOW= yPOW*yp
                  SUM= SUM+ yPOW*PV(j)
                  ENDDO
              ENDIF
          IF(Asw.LE.0.d0) THEN
              XP= SUM + phiINF*yp
            ELSE
              XP= SUM + phiINF*(1.d0 - fsw)
            ENDIF
          IF((NCMM.EQ.4).AND.(MMLR(2).EQ.0)) THEN
c** For Aubert-Frecon based  uLR(r)
              RTP3= RTP(IDAT)**3
              T1= M2/(9.d0*RTP3) + (5.d0*R01+R12)/(45.d0*RTP3**2)
              VLR= 0.5d0*M2/RTP3 + (5.d0*R01 + 8.2d0*R12)
     1       /(18.d0*RTP3**2) + 0.5d0*DSQRT((T1- ASO)**2 + 8.d0*T1**2)
            ELSE
c** For normal inverse-power sum  uLR(r)  ....
              VLR= 0.d0
              DO  J= 1,NCMM
                  VLR= VLR+ CmVAL(j)/RTP(IDAT)**MMLR(j)
                  ENDDO
            ENDIF
          XPW= DEXP(-XP*ype) * VLR/VLRe
          YC= De*(1.d0 - XPW)**2 + VMIN
          DER= 2.d0*De*(1.d0- XPW)*XPW
          IF(Asw.LE.0.d0) THEN
              yPOW= DER*ype*(1.d0- yp)
            ELSE
              ypow= DER*ype*fsw
            ENDIF
          DO  j= 1,NPOW
              PD(j)= yPOW
              yPOW= yPOW*yp
              ENDDO
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
                      DSUM= phiINF - SUM/(1.d0-yp) + DSUM 
                    ELSE
                      DSUM= 0.d0
                    ENDIF
                  PD(NPHI+1)= DER*(dype*(XP + ype*DSUM) 
     1                               + (1.d0 - ype*yp)*dVLRedRe/VLRe )
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
          IF(IDAT.EQ.1) THEN
              IF(JFXRe.LE.0) Re= PV(NPHI+1)
              IF(JFXDe.LE.0) De= PV(NPHI+2)
              IF(JFXVMIN.LE.0) VMIN= PV(NPHI+3)
              AREF= RREF
              IF(RREF.LE.0) AREF= Re
              AREFp= AREF**p
              VLRe= 0.d0
              dVLRe= 0.d0
              d2VLRe= 0.d0
c-----------------------------------------------------------------------
c** Evaluate VLR & its first 2 deriv. at  Re ... 
              IF(IDF.EQ.2) THEN
c... if using Scoles' damping function ...
                  DO  j= 1,NCMM
                      FCT= dexp(-3.97d0*RHOd*Re/MMLR(j)
     1                   - 0.39d0*(RHOd*Re)**2/DSQRT(DFLOAT(MMLR(j))))
                      VCN= CmVAL(j)*((1.d0- FCT)/Re)**MMLR(j)
                      VLRe= VLRe+ VCN
                      dVLRe= dVLRe- MMLR(j)*VCN/Re
                      DDER= FCT*(3.97d0*RHOd/MMLR(j)
     1                     + 0.78d0*Re*RHOd**2/DSQRT(DFLOAT(MMLR(j))))
                      dVLRe= dVLRe+  MMLR(j)*DDER*VCN/(1.d0- FCT)
                      d2VLRe= d2VLRe + MMLR(j)*(MMLR(j)+1)*VCN/Re**2
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
                  VLRe= 0.d0
                  dVLRe= 0.d0
                  d2VLRe= 0.d0
                  DO  j= 1,NCMM
                      TTM= (1.d0- yPOW*TT(MMLR(j)))*CmVAL(j)/Re**MMLR(j)
                      VLRe= VLRe+ TTM
                      TTMM= yPOW*RHOd*(TT(MMLR(j)) - TT(MMLR(j)-1))
     1                                           *CmVAL(j)/Re**MMLR(j)
                      dVLRe= dVLRe+ TTMM - MMLR(j)*TTM/Re
                      d2VLRe= d2VLRe + MMLR(j)*(MMLR(j)+1)*TTM/Re**2
     1               + yPOW*RHOd**2*(-TT(MMLR(j)) + 2.d0*TT(MMLR(j)-1)
     2                         -TT(MMLR(j)-2))  - 2.d0*MMLR(j)*TTMM/Re
                      ENDDO
                  ENDIF
c-----------------------------------------------------------------------
              AA= De - VLRe - dVLRe/PV(1)
              BB= 2.d0*(De - VLRe) - dVLRe/PV(1)
              ENDIF
          RTPp = RTP(IDAT)**p
          yp= (RTPp - AREFp)/(RTPp + AREFp)
          VLR= 0.d0
c ... evaluate VLR at the actual distance ...
          IF(IDF.EQ.2) THEN
c... if using Scoles' damping function ...
              DO  j= 1,NCMM
                  FCT= dexp(-3.97d0*RHOd*RTP(IDAT)/MMLR(j)
     1             - 0.39d0*(RHOd*RTP(IDAT))**2/DSQRT(DFLOAT(MMLR(j))))
                  VLR= VLR+ CmVAL(j)*((1.d0- FCT)/RTP(IDAT))**MMLR(j)
                  ENDDO
              ENDIF
          IF(IDF.EQ.1) THEN
c... if using Tang-Toennies damping function ...
              RHOdR= RHOd*RTP(IDAT)
              yPOW= DEXP(-RHOdR)
              TT(0)= yPOW
              DO  j= 1,MMLR(NCMM)
                  yPOW= yPOW*RHOdR/DFLOAT(J)
                  TT(J)= TT(J-1)+ yPOW
                  ENDDO
              DO  j=1,NCMM
                  VLR=VLR+CmVAL(j)*(1.d0-TT(MMLR(j)))/RTP(IDAT)**MMLR(j)
                  ENDDO
              ENDIF
          NPOW= NS+1
          IF(RTP(IDAT).GE.Re) NPOW= NL+1
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
          XP= DEXP(-SUM*(RTP(IDAT)- Re))
          YC= (AA*XP - BB)*XP + De - VLR + VMIN
          DER= XP*(BB - 2.d0*AA*XP)
          DERP= DER*(RTP(IDAT)- Re)
          DO  j= 1,NPOW
              PD(j)= DERP
              DERP= DERP*yp
              ENDDO
          PD(1)= PD(1) + XP*(XP-1.d0)*dVLRe/PV(1)**2
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
              PD(NPHI+1)= DER*(DSUM - SUM) - XP*((XP- 2.d0)*dVLRe 
     1                                       - (XP-1.d0)*D2VLRe/PV(1))
              ENDIF
          ENDIF
c=======================================================================

c=======================================================================
c** For the case of a  GPEF(p,as,bs)  potential
c-----------------------------------------------------------------------
       IF(PSEL.EQ.4) THEN
          IF(IDAT.EQ.1) THEN
              JFXVMIN= IFXP(NPHI+1)
              IF(JFXVMIN.LE.0) VMIN= PV(NPHI+1)
              JFXRe= IFXP(NPHI+2)
              IF(JFXRe.LE.0) Re= PV(NPHI+2)
              Rep= Re**p
              ENDIF
          RTPp= RTP(IDAT)**p
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
cc    if(IDAT.eq.1) then
cc          write(7,700) nparm,(PV(i),i=1,nparm)
cc          write(7,702)  (i,i=1,min0(5,nparm))
cc          ENDIF
cc    write(7,704) i,yc,(pd(i),i=1,nparm)
cc700 FORMAT(///' Partial derivatives for',i3,' input parameters:'/
cc   1  (1P5D16.8))
cc702 FORMAT('  I    YC   ',5('     PD(',i1,')   ':))
cc704 format(i3,f9.2,1P5D13.5:/(12x,5d13.5:))
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
c             COPYRIGHT 1998-2004  by  Robert J. Le Roy                +
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
c     which are held fixed [JFXP(j)=1], those PC(j) values are ignored.
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
          IF(Z(1) .EQ. 0.D0) GOTO 10
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

