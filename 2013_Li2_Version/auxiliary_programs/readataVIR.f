c***********************************************************************
      SUBROUTINE READATA(NSTATES,PASok,UCUTOFF,JTRUNC,EFSEL,VMIN,VMAX,
     1                                            NDAT,NOWIDTHS,PRINP)
c***********************************************************************
c** Subroutine to read, do book-keeping for, and print summary of
c  experimental data used in fits to spectroscopic data for one or more
c  electronic states and one or more isotopomers. 
c             ********* Version of 6 June 2006 *********
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++  COPYRIGHT 1997-2006 by  Robert J. Le Roy & Dominique R.T. Appadoo +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The present program version can treat seven types of experimental
c   experimental data, for up to NISTPMX isotopomers of a given species.
c   The data are read in grouped as "bands", as (fluorescence) series, 
c   as binding energies (from photoassociation spectroscopy), as a set
c   of Bv values for a given electronic state, and [in a potential-fit
c   aanalysis] as tunneling predissociation level widths.  The types are
c   identified by the values of the 'electronic state label' parameters
c   IEP & IEPP.  They are:
c (i)  microwave transitions within a given electronic state;
c (ii)  infrared bands among the vibrational levels a given state;
c (iii) fluorescence series from some initial excited state level into 
c    vibration-rotation levels of a given electronic state
c (iv)  visible (electronic) absorption or emission bands between vib.
c    levels of two electronic state.
c (v)  binding energies - as from photoassociation spectroscopy
c (vi) "experimental" B_v values for vibrational levels of one of the
c    electronic states.
c (vii) Widths of tunneling predissociation quasibound levels (this 
c    option only meaningful for program DSPotFit).  
c-----------------------------------------------------------------------
c** On Entry:
c  NSTATES is the number of electronic states involved in the data set
c    considered (don't count states giving rise to fluorescence series).
c  PASok indicates how photoassociation data to be treated in analysis:
c    If(PASok(ISTATE).GE.1) treat it as proper PA binding energy data.
c    If(PASok(ISTATE).LE.0) treat PAS data as fluorescence series.
c    Set PASok= 0 if potential model has no explicit Dissoc. Energy
c  Data cutoffs:  for levels of electronic state  s , neglect data with:
c     J(s) > JTRUNC(s),  or vibrational levels lying outside the range
c     VMIN(s)  to  VMAX(s),  AND  NEGLECT any data for which the read-
c     in uncertainty is  > UCUTOFF (cm-1).  EFSEL(s) > 0 causes f-parity
c     levels to be neglected, EFSEL(s) < 0 omits e-parity levels
c     while  EFSEL(s) = 0  allows both types of parity to be included.
c  NOWIDTHS > 0  causes the program to ignore any tunneling widths in
c            the data set.
c  PRINP > 0  turns on the printing of a summary description of the data.
c** On Return:
c  UCUTOFF (cm-1)  is the smallest uncertainty in the (accepted) data
c  NDAT(v,i,s)  is the number of transitions associated with 
c    vibrational level-v of isotopomer-i of state-s [for NDEGB < 0 case]
c** This subroutine reads in the experimental data on channel-4
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
c
      INTEGER I,IBB,NTRANS,COUNT,IBAND,JMAX(NPARMX),JMIN(NPARMX),
     1  VMX(NSTATEMX),ISOT,NBND,ESP,ESPP,ISTATE,ISTATEE,MN1,MN2,PRINP,
     2  FSOMIT,VMAXesp,VMINesp,VMAXespp,VMINespp,JTRUNCesp,JTRUNCespp
      INTEGER NSTATES,NOWIDTHS,JTRUNC(NSTATEMX),EFSEL(NSTATEMX),
     1  VMIN(NSTATEMX),VMAX(NSTATEMX),NDAT(0:NVIBMX,NISTPMX,NSTATEMX),
     2  PASok(NSTATES)
      REAL*8 UCUTOFF,UMIN,TOTUFREQ
      CHARACTER*3 NEF(-1:1)
      CHARACTER*2 LABLP,LABLPP

      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKTYPE.h'
c
      DATA NEF/'  f','   ','  e'/
c-----------------------------------------------------------------------
      WRITE(6,603) UCUTOFF 
      DO  ISTATE= 1,NSTATES
          IF(JTRUNC(ISTATE).GE.0) THEN
              WRITE(6,607) SLABL(ISTATE),JTRUNC(ISTATE),VMIN(ISTATE),
     1                                                    VMAX(ISTATE)
            ELSE
              WRITE(6,605) SLABL(ISTATE),-JTRUNC(ISTATE),VMIN(ISTATE),
     1                                                    VMAX(ISTATE)
            ENDIF 
          IF(EFSEL(ISTATE).GT.0) WRITE(6,601) NEF(-1)
          IF(EFSEL(ISTATE).LT.0) WRITE(6,601) NEF(1)
          ENDDO
      UMIN= UCUTOFF
c** Initialize counters for book-keeping on input data
      COUNT= 0
      DO  ISOT= 1,NISTP
          DO  ISTATE= 1,NSTATES
              NTRANSFS(ISOT,ISTATE)= 0
              NTRANSIR(ISOT,ISTATE)= 0
              NTRANSMW(ISOT,ISTATE)= 0
              NBANDFS(ISOT,ISTATE)= 0
              NBANDVIS(ISOT,ISTATE)= 0
              NBANDIR(ISOT,ISTATE)= 0
              NBANDMW(ISOT,ISTATE)= 0
              NBVPP(ISOT,ISTATE)= 0
              NWIDTH(ISOT,ISTATE)= 0
              NEBPAS(ISOT,ISTATE)= 0
              NVIRIAL(ISOT,ISTATE)= 0
              DO  I= 1,NSTATES
                  NTRANSVIS(ISOT,ISTATE,I)= 0
                  NBANDEL(ISOT,ISTATE,I)= 0
                  ENDDO
              ENDDO
          NBANDS(ISOT)= 0
          ENDDO
      DO  ISTATE= 1,NSTATES
          VMX(ISTATE)= 0
          ENDDO
      NFSTOT= 0
      FSOMIT= 0
c========================================================================
c** Begin loop to read in data, band(or series)-by-band(or series).
c  STOP when run out of bands or when encounter a negative vibrational
c  quantum number.
c** Read all data for each isotopomer at one time.
      IBAND= 0
   10 CONTINUE
      IBAND= IBAND+1
      IF(IBAND.GT.NPARMX) THEN
            IF(PRINP.GT.0) WRITE(6,609) IBAND,NPARMX
            IBAND= IBAND-1
            GOTO 20
            ENDIF
c
c For each "band", read in:  (i) upper/lower vibrational quantum numbers
c   VP & VPP,  (ii) a two-character electronic-state alphameric label 
c   {enclosed in single quotes; e.g., 'X0' or 'A1'} for the upper
c   (LABLP) and lower (LABLP) state, and  (iii) integers NM1 & NM2 are
c   the mass numbers [corresponding to input atomic numbers AN(1) & 
c   AN(2)] identifying the particular isotopomer.  Note that LABLP also
c   identifies the type of data in the 'band' or data-group (see below).
c
c** LABLP = LABLPP  and  VP = VPP  for a microwave band
c   LABLP = LABLPP  and  VP.ne.VPP  for an infrared band 
c   LABLP = 'FS'  identifies this data group/band as a fluorescence 
c           series from a single emitting level into vibrational levels
c           of electronic state LABLPP.  In this case: VP is the quantum
c           number v' for the emitting level, while VPP is actually the 
c           rotational quantum number J' for the emitting level and JP
c           [see below] the lower state vibrational quantum number v".
c   LABLP = 'PA'  identifies this data group/band as a set of binding
c           energies [D-E(v,J,p)] for a given state.  Labels as for 'FS'
c   LABLP = 'BV'  identifies this data group/band as a set of Bv values
c           for electronic state LABLPP.  In this case, parameters  VP
c           & VPP are dummy variables, as are EFP, JPP and EFPP [see
c           below],  JP is actually the vibrational quantum number v",
c           FREQ the Bv value & UFREQ its uncertainty
c   LABLP = 'WI'  identifies this data group/band as a set of tunneling 
c           predissociation widths for electronic state LABLPP.  In this
c           case, parameters VP, VPP and EFP are dummy variables, while
c           the predissociating level is identified as: v"=JP, J"=JPP,
c           and parity p"=EFPP.
c   LABLP = 'VR' identifies this data group/band as a set of virial 
c           coefficients for electronic state LABLPP.  In this case, 
c           parameters VP, VPP are dummy variables.
c   NOTE: !!!!!!!!!!! This option is ignored by DSParFit !!!!!!!!!!!!!!!
c** STOP reading when run out of bands OR when read-in VPP is negative   
c-----------------------------------------------------------------------
      READ(4,*,END=20) VP(IBAND), VPP(IBAND), LABLP, LABLPP, MN1,MN2
c-----------------------------------------------------------------------
      IF(VP(IBAND).LT.0) GO TO 20
      IEP(IBAND)= -99
      IEPP(IBAND)= -99
      DO  I= -4,NSTATES
          IF(LABLP.EQ.SLABL(I)) IEP(IBAND)= I
          IF(LABLPP.EQ.SLABL(I)) IEPP(IBAND)= I
          ENDDO
c** Check that this isotopomer is one of those chosen to be fitted ...
      ISOT= 0
      DO  I= 1,NISTP
          IF((MN1.EQ.MN(1,I)).AND.(MN2.EQ.MN(2,I))) ISOT= I
          ENDDO
      IF(IEP(IBAND).EQ.-4) THEN
c** Now ... consider the case of virial coefficients ...
          COUNT= COUNT+ 1
          IFIRST(IBAND)= COUNT
   14     READ(4,*) TEMP(COUNT),FREQ(COUNT),UFREQ(COUNT)
          IF(TEMP(COUNT).GT.0.d0) THEN
              COUNT= COUNT+1
              GOTO 14
              ENDIF
          ILAST(IBAND)= COUNT
          NTRANS= ILAST(IBAND) - IFIRST(IBAND) + 1
          IF(ISOT.LE.0) THEN
              COUNT= COUNT- NTRANS
              IBAND= IBAND- 1
              ENDIF
          GOTO 10
          ENDIF
c... now ... for the case of spectroscopic data ...
      ISTP(IBAND)= ISOT
      TOTUFREQ= 0.D0
      MAXUFREQ(IBAND)= 0
      JMAX(IBAND)= 0
      JMIN(IBAND)= 9999
      COUNT= COUNT+1
      IF(COUNT.GT.NDATAMX) THEN
          WRITE(6,640) COUNT,NDATAMX
          STOP
          ENDIF
      NTRANS= 0
      IFIRST(IBAND)= COUNT
      ESP= IEP(IBAND)
      ESPP= IEPP(IBAND)
      VMAXespp= VMAX(ESPP)
      VMINespp= VMIN(ESPP)
      JTRUNCespp= JTRUNC(ESPP)
      IF(ISOT.GT.1) THEN
          VMAXespp= INT((VMAX(ESPP)+0.5d0)/RSQMU(ISOT)-0.5d0)
          VMINespp= INT((VMIN(ESPP)+0.5d0)/RSQMU(ISOT)-0.5d0)
          JTRUNCespp= INT(JTRUNC(ESPP)/RSQMU(ISOT))
          ENDIF
      VMAXesp= VMAX(ESPP)
      IF(ESP.GT.0) THEN
          VMAXesp= VMAX(ESP)
          VMINesp= VMIN(ESP)
          JTRUNCesp= JTRUNC(ESP)
          IF(ISOT.GT.1) THEN
              VMAXesp= INT((VMAX(ESP)+ 0.5d0)/RSQMU(ISOT) - 0.5d0)
              VMINesp= INT((VMIN(ESP)+ 0.5d0)/RSQMU(ISOT) - 0.5d0)
              JTRUNCesp= INT(JTRUNC(ESP)/RSQMU(ISOT))
              ENDIF
          ENDIF
c** For each of the lines in a given band/series, read upper level
c  rotational quantum number (JP) and e/f parity [EFP= +1 for e, -1 for
c  f, and  0 if e/f splitting unresolved and to  be ignored], and lower
c  level rotational quantum number (JPP) and parity [EFPP, as above],
c  the transition frequency  FREQ, and its uncertainty UFREQ.
c** For PAS or Tunneling Width data,  JP(COUNT)=v", JPP(COUNT)=J", 
c  EFPP(COUNT)=p", FREQ is the observable (a positive No.), while 
c  EFP(COUNT), VP(IBAND) & VPP(IBAND) are dummy variables.
c** For Bv values, JP(COUNT)=v" while JPP(COUNT), EFP(COUNT) and
c   EFPP(COUNT) as well as VP(IBAND) & VPP(IBAND) are dummy variables.
c-----------------------------------------------------------------------
   15 READ(4,*) JP(COUNT), EFP(COUNT), JPP(COUNT), EFPP(COUNT), 
     1                                       FREQ(COUNT), UFREQ(COUNT)
c-----------------------------------------------------------------------
c=======================================================================
c   Sample IR band data of HF for the '.4' file:
c   --------------------------------------------                          
c   1 0  'X0' 'X0'  1 19             % VP VPP LABLP LABLPP MN1 MN2
c   8 1   9 1  266.0131002  0.005    % JP EFP JPP EFPP FREQ UFREQ
c   9 1  10 1  265.8885896  0.003
c  10 1  11 1  265.7716591  0.002
c   .    .      .            .
c   .    .      .            .
c   [end of a band indicated by -ve JP and/or JPP value(s)]
c  -1 1  -1 1  -1.1         -1.1
c=======================================================================
      IF(EFP(COUNT).GT.1) EFP(COUNT)= 1
      IF(EFP(COUNT).LT.-1) EFP(COUNT)= -1
      IF(EFPP(COUNT).GT.1) EFPP(COUNT)= 1
      IF(EFPP(COUNT).LT.-1) EFPP(COUNT)= -1
c** At end of a band, exit from implicit loop
      IF((JPP(COUNT).LT.0).OR.(JP(COUNT).LT.0)) GOTO 18
c** If this band is not for one of the isotopomers chosen to be fitted,
c  omit its data from the fit
      IF(ISOT.EQ.0) GO TO 15
c** If this band involves electronic states other than those chosen to 
c   be treated, omit its data from the fit
      IF((ESP.EQ.-99).OR.(ESPP.EQ.-99)) GO TO 15
c** If a datum uncertainty of zero is accidentally read in, STOP
      IF(DABS(UFREQ(COUNT)).LE.0.d0) THEN
          WRITE(6,600) COUNT,FREQ(COUNT),IBAND
          STOP
          ENDIF
c** Omit data  with uncertainties outside specified limit UCUTOFF
      IF(UFREQ(COUNT).GT.UCUTOFF) GOTO 15
c** Require that datum lies within specified J & v ranges
      IF(ESP.GE.-2) THEN
          IF(((JTRUNCespp.GE.0).AND.(JPP(COUNT).GT.JTRUNCespp)).OR.
     1       ((JTRUNCespp.LT.0).AND.(JPP(COUNT).LT.-JTRUNCespp)))
     2                                                         GOTO 15
          IF((EFPP(COUNT)*EFSEL(ESPP)).LT.0) GOTO 15
          ENDIF
      IF(ESP.GT.0) THEN
          IF(VPP(IBAND).GT.VMAXespp) GOTO 15
          IF(VPP(IBAND).LT.VMINespp) GOTO 15
          IF(VP(IBAND).GT.VMAXesp) GOTO 15
          IF(VP(IBAND).LT.VMINesp) GOTO 15
          IF((JTRUNCesp.GE.0).AND.(JP(COUNT).GT.JTRUNCesp)) GOTO 15
          IF((JTRUNCesp.LT.0).AND.(JP(COUNT).LT.-JTRUNCesp)) GOTO 15
          IF((EFP(COUNT)*EFSEL(ESP)).LT.0) GOTO 15
        ELSE
          IF(JP(COUNT).GT.VMAXespp) GOTO 15
          IF(JP(COUNT).LT.VMINespp) GOTO 15
        ENDIF
c** If NOWIDTHS > 0  omit any tunneling width data from the fit.
      IF((ESP.EQ.-2).AND.(NOWIDTHS.GT.0)) GOTO 15
c
c** End of tests for datum inclusion.  Now count/sort data
c=======================================================================
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%% Convert  MHz  to  cm-1
c     freq(count)=freq(count)/2.99792458d+4
c     ufreq(count)=ufreq(count)/2.99792458d+4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      TVUP(COUNT)= 0
      TVLW(COUNT)= 0
      IF(ESP.GE.-1) UMIN= MIN(UMIN,UFREQ(COUNT))
c** Determine actual v & J range of data & count data for each v
c  JMIN & JMAX needed for printout summary & data-count for testing
c  no. parameters allowed in Band Constant fit.
c??? This segment imperfect & needs re-examination ?????????????
      IF(ESP.GT.0) THEN
          IF(JPP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JPP(COUNT) 
          IF(JPP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JPP(COUNT)
          IF(JP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JP(COUNT) 
          IF(JP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JP(COUNT)
          VMX(ESP)= MAX(VMX(ESP),VP(IBAND))
          VMX(ESPP)= MAX(VMX(ESPP),VPP(IBAND))
c
c** Accumulate count of data associated with each vibrational level ...
          NDAT(VPP(IBAND),ISTP(IBAND),ESPP)=
     1                            NDAT(VPP(IBAND),ISTP(IBAND),ESPP)+ 1
          NDAT(VP(IBAND),ISTP(IBAND),ESP)=
     1                              NDAT(VP(IBAND),ISTP(IBAND),ESP)+ 1
        ELSEIF((ESP.LE.0).OR.(ESP.GE.-2)) THEN
          IF(JP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JP(COUNT) 
          IF(JP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JP(COUNT)
          VMX(ESPP)= MAX(VMX(ESPP),JP(COUNT))
          NDAT(JP(COUNT),ISTP(IBAND),ESPP)=
     1                      NDAT(JP(COUNT),ISTP(IBAND),ESPP)+ 1
        ELSEIF(ESP.EQ.-3) THEN
c... and for Bv data ...
          IF(JPP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JPP(COUNT) 
          IF(JPP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JPP(COUNT)
          NDAT(JPP(COUNT),ISTP(IBAND),ESPP)=
     1                     NDAT(JPP(COUNT),ISTP(IBAND),ESPP)+ 1
        ENDIF
      DFREQ(COUNT)= 0.d0
      IB(COUNT)= IBAND 
      TOTUFREQ= TOTUFREQ+UFREQ(COUNT) 
      IF(UFREQ(COUNT).GT.MAXUFREQ(IBAND)) MAXUFREQ(IBAND)= UFREQ(COUNT)
      COUNT= COUNT+1 
      IF(COUNT.GT.NDATAMX) THEN
          WRITE(6,640) COUNT,NDATAMX
          STOP
          ENDIF
      GOTO 15 
c** End of loop reading data for a given band/series 
c
c** Tidy up at end of reading for a given band
   18 COUNT= COUNT-1
      ILAST(IBAND)= COUNT 
      NTRANS= ILAST(IBAND)-IFIRST(IBAND)+1
      IF(NTRANS.GT.0) THEN
c** Treat PAS data as Fluorescence series unless  PASok > 0
          IF((IEP(IBAND).EQ.-1).AND.(PASok(IEPP(IBAND)).LE.0)) 
     1                                                    IEP(IBAND)=0
          IF((NTRANS.EQ.1).AND.(LABLP.EQ.'FS')) THEN
c** Ignore any fluorescence series consisting of only one datum
              COUNT= COUNT-1
              IBAND= IBAND-1
              FSOMIT= FSOMIT+1
              GOTO 10
              ENDIF
          AVEUFREQ(IBAND)= TOTUFREQ/NTRANS
          NBANDS(ISTP(IBAND))= NBANDS(ISTP(IBAND))+1
        ELSE
          IBAND= IBAND-1
          GOTO 10
        ENDIF
c=======================================================================
c** Accumulate counters for bands/series of different types
      IF(ESP.EQ.0) THEN
c** For Fluorescence Series ... first enumerate the No. of bands & lines
          NFSTOT= NFSTOT+1
          FSBAND(NFSTOT)= IBAND
c** Define counter to label which f.s. is associated with band IBAND 
          NFS(IBAND)= NFSTOT
          NBANDFS(ISOT,ESPP)= NBANDFS(ISOT,ESPP)+1
          NBND= NBANDFS(ISOT,ESPP)
          NTRANSFS(ISOT,ESPP)= NTRANSFS(ISOT,ESPP)+NTRANS
c ... and then set up labels/ranges/properties for each band
          YPR(ISOT,ESPP,1,1,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,1,2,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,1,3,NBND)= NTRANS
          YPR(ISOT,ESPP,1,4,NBND)= IBAND
          YPR(ISOT,ESPP,1,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,1,6,NBND)= JMAX(IBAND)
          ENDIF
c
      IF((ESP.GT.0).AND.(ESP.NE.ESPP)) THEN
c** For vibrational band of a normal 2-state electronic transition
c ... count bands and transitions in visible (electronic) spectrum
          NBANDEL(ISOT,ESP,ESPP)= NBANDEL(ISOT,ESP,ESPP)+ 1
          NBANDVIS(ISOT,ESPP)= NBANDVIS(ISOT,ESPP)+ 1
          NBND= NBANDVIS(ISOT,ESPP)
          NTRANSVIS(ISOT,ESP,ESPP)= NTRANSVIS(ISOT,ESP,ESPP)+NTRANS
c ... and then set up labels/ranges/properties for each of them
          YPR(ISOT,ESPP,2,1,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,2,2,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,2,3,NBND)= NTRANS
          YPR(ISOT,ESPP,2,4,NBND)= IBAND
          YPR(ISOT,ESPP,2,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,2,6,NBND)= JMAX(IBAND)
          ENDIF 
c
      IF((ESP.EQ.ESPP).AND.(VP(IBAND).NE.VPP(IBAND))) THEN
c** For an Infrared band of electronic state  s=ESPP=ESP
c** First cumulatively count the number of IR bands & transitions
          NBANDIR(ISOT,ESPP)= NBANDIR(ISOT,ESPP)+1
          NBND= NBANDIR(ISOT,ESPP)
          NTRANSIR(ISOT,ESPP)= NTRANSIR(ISOT,ESPP)+NTRANS 
c ... and then set up labels/ranges/properties for each of them
          YPR(ISOT,ESPP,3,1,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,3,2,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,3,3,NBND)= NTRANS
          YPR(ISOT,ESPP,3,4,NBND)= IBAND
          YPR(ISOT,ESPP,3,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,3,6,NBND)= JMAX(IBAND)
          ENDIF
c
      IF((ESP.EQ.ESPP).AND.(VP(IBAND).EQ.VPP(IBAND))) THEN
c** For Microwave transitions in electronic state  s=ESPP=ESP
c** First cumulatively count the number of MW bands & transitions
          NBANDMW(ISOT,ESPP)= NBANDMW(ISOT,ESPP)+1
          NBND= NBANDMW(ISOT,ESPP)
          NTRANSMW(ISOT,ESPP)= NTRANSMW(ISOT,ESPP)+NTRANS
c ... and then set up labels/ranges/properties for each of them
          YPR(ISOT,ESPP,4,1,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,4,2,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,4,3,NBND)= NTRANS
          YPR(ISOT,ESPP,4,4,NBND)= IBAND
          YPR(ISOT,ESPP,4,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,4,6,NBND)= JMAX(IBAND)
          ENDIF
c
c** NOTE ... in YPR array a last index counts bands of this type for 
c  this isotopomer of this electronic state ... and put all Bv's, 
c  Tunneling Widths or PAS binding energies in one group.
      IF(ESP.EQ.-3) THEN
c** Data are not transition energies, but rather the values of Bv in
c  electronic state s=IEPP  [As in the published IBr(A-X) analysis].
ccc       IF((NBVPP(ISOT,ESPP).GT.0).AND.(NTRANS.GT.0)) THEN
              WRITE(6,612) ESPP,ISOT
ccc           STOP
ccc           ENDIF
          NBVPP(ISOT,ESPP)= NTRANS
          YPR(ISOT,ESPP,5,3,1)= NTRANS
          YPR(ISOT,ESPP,5,4,1)= IBAND
          YPR(ISOT,ESPP,5,5,1)= JMIN(IBAND)
          YPR(ISOT,ESPP,5,6,1)= JMAX(IBAND)
          ENDIF
c
      IF(ESP.EQ.-2) THEN
c** Data are tunneling predissociation linewidths (in cm-1) for levels
c  of electronic state IEPP=ESPP
ccc       IF((NWIDTH(ISOT,ESPP).GT.0).AND.(NTRANS.GT.0)) THEN
              WRITE(6,626) ESPP,ISOT
ccc           STOP
ccc           ENDIF
          NWIDTH(ISOT,ESPP)= NTRANS
          YPR(ISOT,ESPP,6,3,1)= NTRANS
          YPR(ISOT,ESPP,6,4,1)= IBAND
          YPR(ISOT,ESPP,6,5,1)= JMIN(IBAND)
          YPR(ISOT,ESPP,6,6,1)= JMAX(IBAND)
          ENDIF
c
      IF(ESP.EQ.-1) THEN
c** Data are PhotoAssociation Binding Energies (in cm-1) for levels
c  of electronic state IEPP=ESPP
          WRITE(6,636) LABLPP,ISOT
          NEBPAS(ISOT,ESPP)= NTRANS
          YPR(ISOT,ESPP,7,3,1)= NTRANS
          YPR(ISOT,ESPP,7,4,1)= IBAND
          YPR(ISOT,ESPP,7,5,1)= JMIN(IBAND)
          YPR(ISOT,ESPP,7,6,1)= JMAX(IBAND)
          ENDIF
c** Now return to read the next band
      GOTO 10
c========================================================================
c** Now, write a summary of the input data to the output file
   20 COUNTOT= COUNT
      NBANDTOT= 0
      DO  I= 1,NISTP
          NBANDTOT= NBANDTOT+ NBANDS(I)
          ENDDO
      ISOT= 1
      UCUTOFF= UMIN
      IF(FSOMIT.GT.0) WRITE(6,650) FSOMIT
      IF(PRINP.LE.0) RETURN
c** Print a summary of the data, one isotopomer at a time.
   26 WRITE(6,602) NBANDS(ISOT), (NAME(I),MN(I,ISOT),I=1,2)
c
      DO 50 ISTATE= 1,NSTATES
c ... For internal use, may wish to update VMAX(ISTATE) to the actual 
c  highest v in the data set for this state. ** Reactivate as needed.
c      VMAX(ISTATE)= VMX(ISTATE)
c ... and separately list data for each (lower) electronic state in turn
      IF(NTRANSMW(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for Micowave data
          WRITE(6,604) NTRANSMW(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDMW(ISOT,ISTATE)
          DO  I= 1,NBANDMW(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,4,4,I)
              WRITE(6,606) YPR(ISOT,ISTATE,4,2,I),
     1                     YPR(ISOT,ISTATE,4,1,I),
     2                  YPR(ISOT,ISTATE,4,3,I),YPR(ISOT,ISTATE,4,5,I),
     3                  YPR(ISOT,ISTATE,4,6,I), 
     3                  AVEUFREQ(YPR(ISOT,ISTATE,4,4,I)),
     4                  MAXUFREQ(YPR(ISOT,ISTATE,4,4,I))
              ENDDO
	    ENDIF
c
      IF(NTRANSIR(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Infrared data
          WRITE(6,608) NTRANSIR(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDIR(ISOT,ISTATE)
          DO  I= 1,NBANDIR(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,3,4,I)
              WRITE(6,606) YPR(ISOT,ISTATE,3,2,I),
     1                     YPR(ISOT,ISTATE,3,1,I),
     2                  YPR(ISOT,ISTATE,3,3,I),YPR(ISOT,ISTATE,3,5,I),
     3                  YPR(ISOT,ISTATE,3,6,I), 
     4                  AVEUFREQ(YPR(ISOT,ISTATE,3,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,3,4,I))
              ENDDO
          ENDIF
c
c** Book-keeping for electronic vibrational band data
      DO  ISTATEE= 1,NSTATES
          IF((ISTATEE.NE.ISTATE).AND.
     1                 (NTRANSVIS(ISOT,ISTATEE,ISTATE).GT.0)) THEN
c ... for ISTATEE{upper}-ISTATE{lower} electronic vibrational bands
              WRITE(6,610) NTRANSVIS(ISOT,ISTATEE,ISTATE),
     1         (NAME(I),MN(I,ISOT),I=1,2),SLABL(ISTATEE),SLABL(ISTATE),
     2                                    NBANDEL(ISOT,ISTATEE,ISTATE)
              DO  I= 1,NBANDVIS(ISOT,ISTATE)
                  IBB= YPR(ISOT,ISTATE,2,4,I)
                  IF(IEP(IBB).EQ.ISTATEE) THEN
                      WRITE(6,606) YPR(ISOT,ISTATE,2,2,I),
     1                            YPR(ISOT,ISTATE,2,1,I),
     2                  YPR(ISOT,ISTATE,2,3,I),YPR(ISOT,ISTATE,2,5,I),
     3                  YPR(ISOT,ISTATE,2,6,I), 
     4                  AVEUFREQ(YPR(ISOT,ISTATE,2,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,2,4,I))
                      ENDIF
                  ENDDO
              ENDIF
          ENDDO
      IF(NTRANSFS(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Fluorescence data
          WRITE(6,614) NTRANSFS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDFS(ISOT,ISTATE)
          DO  I= 1,NBANDFS(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,1,4,I)
              WRITE(6,616) YPR(ISOT,ISTATE,1,1,I),
     1                   YPR(ISOT,ISTATE,1,2,I),NEF(EFP(IFIRST(IBB))),
     2                  YPR(ISOT,ISTATE,1,3,I),YPR(ISOT,ISTATE,1,5,I),
     3                  YPR(ISOT,ISTATE,1,6,I), 
     4                  AVEUFREQ(YPR(ISOT,ISTATE,1,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,1,4,I))
              ENDDO
          ENDIF
      IF(NBVPP(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for  Bv  data
          WRITE(6,618) NBVPP(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBB= YPR(ISOT,ISTATE,5,4,1)
          WRITE(6,620) YPR(ISOT,ISTATE,5,3,1),YPR(ISOT,ISTATE,5,5,1),
     1       YPR(ISOT,ISTATE,5,6,1),AVEUFREQ(YPR(ISOT,ISTATE,5,4,1)),
     2                               MAXUFREQ(YPR(ISOT,ISTATE,5,4,1))
          ENDIF
      IF(NWIDTH(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Tunneling Width  data
          WRITE(6,628) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBB= YPR(ISOT,ISTATE,6,4,1)
          WRITE(6,630) YPR(ISOT,ISTATE,6,3,1),
     1                  YPR(ISOT,ISTATE,6,5,1),YPR(ISOT,ISTATE,6,6,1),
     2                               AVEUFREQ(YPR(ISOT,ISTATE,6,4,1)),
     3                                MAXUFREQ(YPR(ISOT,ISTATE,6,4,1))
          ENDIF
      IF(NEBPAS(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  PAS Binding Energy  data
          WRITE(6,632) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBB= YPR(ISOT,ISTATE,6,4,1)
          WRITE(6,630) YPR(ISOT,ISTATE,7,3,1),
     1                  YPR(ISOT,ISTATE,7,5,1),YPR(ISOT,ISTATE,7,6,1),
     2                               AVEUFREQ(YPR(ISOT,ISTATE,7,4,1)),
     3                                MAXUFREQ(YPR(ISOT,ISTATE,7,4,1))
          ENDIF
   50 CONTINUE
      IF(ISOT.LT.NISTP) THEN
c** If NISTP > 1, return to print data summaries for other isotopomers
          ISOT= ISOT+1
          GO TO 26
          ENDIF 
      WRITE(6,622)
      RETURN
  600 FORMAT(/' *** INPUT ERROR ***  Datum   FREQ(',i5,')=',f12.4,
     1 '  in   IBAND=',i4,'   has zero uncertainty!!!')
  601 FORMAT(23x,'or with',A3,'-parity.')
  603 FORMAT(/' Neglect data with:  Uncertainties > UCUTOFF=',G12.3,
     1  ' (cm-1)')
  605 FORMAT(7x,'and State ',A2,' data with  J < JTRUNC=',I4,
     1  '  or  v  outside range',i3,'  to',i4)
  607 FORMAT(7x,'and State ',A2,' data with  J > JTRUNC=',I4,
     2  '  or  v  outside range',i3,'  to',i4)
  602 FORMAT(/1x,20('===')/'  *** Input data for',i5,' bands/series of '
     1  ,A2,'(',I3,')-',A2,'(',I3,') ***'/1x,20('==='))
  604 FORMAT(1x,28('--')/I5,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') MW transitions in',i4,' sets'/1x,28('--')/"   v'  ",
     1 'v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/1x,25('--'))
  606 FORMAT(I4,I4,3I7,1x,1P2D10.1)
  608 FORMAT(1x,32('--')/I6,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') InfraRed transitions in',I4,' bands'/1x,32('--')/
     2 "   v'  ",'v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/
     3 1x,25('--'))
  609 FORMAT(/' *** ERROR *** Dimension allocated for number of bands ex
     1ceeded:'/' (IBAND=',i4,') > (NBANDMX=',i4,')   so truncate input a
     2nd TRY to continue ...')
  610 FORMAT(/1x,35('==')/I6,1x,A2,'(',I3,')-',A2,'(',i3,')  {State ',
     1  A2,'}--{State ',A2,'} Transitions in',i4,' Bands'/1x,35('--')/
     2 "   v'",'  v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/
     3 1x,25('--'))
  612 FORMAT(/" NOTE that all read-in Bv's for   ISTATE=",i2,'   ISOT=',
     1  i2/32x,' must be input as a single "band" or data group')
cc612 FORMAT(/" *** STOP INPUT *** and put all read-in Bv's for   ISTATE
cc   1=",i2,'   ISOT=',i2/ 10x,'into one "band" or data group.')
  614 FORMAT(1x,38('==')/I5,' Fluorescence transitions into State ',
     1 A2,2x,A2,'(',I3,')-',A2,'(',I3,')  in',i5,' series'/
     2 1x,38('==')/"   v'  j' p' ",'#data  v"min  v"max  Avge.Unc.  Max.
     3Unc.'/1x,51('-'))
  616 FORMAT(2I4,A3,I6,2I7,1x,1P2D10.1)
  618 FORMAT(1x,65('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') Bv values treated as independent data'/1x,24('--')/
     2 '  #values   v(min)  v(max)  Avge.Unc.   Max.Unc.'/
     3 1x,24('--'))
  620 FORMAT(I7,I9,I8,3x,1P2D11.1)
  622 FORMAT(1x,25('===')/1x,25('==='))
  626 FORMAT(/" NOTE that all read-in Tunneling Widths for   ISTATE=",
     1 i2,'   ISOT=',i2/10x,' must be in a single "band" or data group')
cc626 FORMAT(/" *** STOP INPUT *** and put all read-in Tunneling Widths'
cc   1  '  for   ISTATE=",i2,'   ISOT=',i2/ 
cc   2  10x,'into one "band" or data group.')
  628 FORMAT(1x,61('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') Tunneling Widths included as data'/
     2 1x,61('-')/'  #values   v(min)  v(max)   Avge.Unc.   Max.Unc.'/
     3 1x,24('--'))
  630 FORMAT(I7,I9,I8,2x,1P2D11.1)
  632 FORMAT(1x,70('=')/I4,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') PAS Binding Energies included in data set'/
     2 1x,70('-')/'  #values   v(min)  v(max)   Avge.Unc.   Max.Unc.'/
     3 1x,24('--'))
  636 FORMAT(/' NOTE that all read-in PAS Binding Energies for   ISTATE=
     1 ',a2,'  ISOT=',i2/10x,' must be in a single "band" or data group'
     2 )
  640 FORMAT(/' *** Input Data Count reaches',i6,' which EXCEEDS ARRAY L
     1IMIT of',i6)
  650 FORMAT(/' Data input IGNORES',i4,' fluorescence series consisting'
     1 ,' of only  onee  line!')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

