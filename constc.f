       SUBROUTINE CONSTC
 
c compute commonly used constants to be saved

      INCLUDE 'zeq.common'
 
      BETA = 2.28E-11
      SQRT2=SQRT(2.e0)
      PI =ACOS(-1.e0)
 
c the various dimensional scaling factors.....
c if hequiv=0, inputs are assumed in nondimensional form

      IF(HEQUIV.GT.0.0) GO TO 5
      CWAVE=1.0
      ELEQ=1.00
      TEQ=1.00
      TEQM=1.00
      WSCALE=1.00
      GO TO 6

c eleq is in degrees latitude, teq is in days
    5 CWAVE=SQRT(9.81E-2*HEQUIV)
      ELEQ= SQRT(CWAVE/BETA)/111.E3
      TEQ=1.E0/ SQRT(CWAVE*BETA)/86400.
      TEQM=TEQ*12.E0/365.25

c scale the wind stress forcing so that H will come out in meters
c it is assumed the wind stress is specified in dy/cm/cm
      WSCALE =1.E-4*ELEQ*111.E3/CWAVE/CWAVE
    6 WRITE(6,7) HEQUIV,CWAVE,ELEQ,TEQ,TEQM,WSCALE
    7 FORMAT(//' HEQUIV =',F4.0,' CM.  CWAVE =',F5.2,' M/S.  LEQ =',
     A         F5.2,' DEGREES.  TEQ =',F5.2,' DAYS =',F7.4,' MONTHS.'/
     B         ' WSCALE =',F5.2,' M/(CM/SEC)**2'/)
 
 
      XE=(XED-XWD)/ELEQ
      YS=YSD/ELEQ
      YN=YND/ELEQ
            NTIMES=(TENDD-TZERO)/DTD+0.1
            DX = XE/FLOAT(NX)
            DY = (YN-YS)/FLOAT(NY)
            NXP = NX+1
            NYP = NY+1
      DT=DTD/TEQM
      RFRIC=TEQM/TDECAY
      DXD=DX*ELEQ
      DYD=DY*ELEQ
 
 
      WRITE(6,8) XWD,XED,YSD,YND,DXD,DYD,DTD,TDECAY,XE,YS,YN,DX,DY,DT
     A           ,RFRIC
    8 FORMAT(25X,'XW        XE        YS        YN        DX       DY'
     A   ,'       DT       TDECAY'
     A ,/' DIMENSIONAL      ',4F10.2,3F10.5,F10.2
     A ,/' NONDIMENSIONAL   ','       0.0',3F10.2,4F10.5/)
      WRITE(6,9) NX,NY,TZERO,TENDD,NTIMES
    9 FORMAT(/'  NX =',I5,'  NY =',I5,'  TZERO = ',F7.2,
     A '  TENDD = ',F10.2,'  NTIMES = ',I7/)
C**********************************************************************
 
      DYH=0.5*DY
      DYINV = 2.0/DY
 
       DO 10 J=1,NY
      YY(J)=YS+(J-.5)*DY
10    C(J)=YY(J)*YY(J)
 
c compute phik, the meridional structure function of kelvin mode
c anorm is the analytic normalization factor for kelvin amplitudes
c phik is the kelvin mode; ie., anorm*exp(-.5*y*y)
 
       ENORM = SQRT(PI)
       ANORM = ( ENORM*( ERF(YN)+ERF(-YS) ) ) ** (-0.5)
      JNOT = 1.55-YY(1)/DY
      IF(JNOT.LT.2) JNOT=2
      JNOT1=JNOT+1
      PHIK(JNOT)=ANORM*EXP(-0.5*YY(JNOT)*YY(JNOT))
      DO 11 J=JNOT1,NY
11    PHIK(J)=PHIK(J-1)*(DYINV-YY(J-1))/(DYINV+YY(J))
      DO 12 J=2,JNOT
      K=JNOT1-J
12    PHIK(K)=PHIK(K+1)*(DYINV+YY(K+1))/(DYINV-YY(K))
 
c For normalization of phik, pone is projection of 1 on kelvin mode
c and pnorm is sum of phik*phik over all y.
      PONE(1)=0.0
      PNORM(1)=0.0
      EPSLON=1.E-10*PHIK(JNOT)
      DO 21 J=1,NY
      PH=PHIK(J)
      PONE(J+1)=PONE(J)+PH
      IF(PH.LT.EPSLON) PH=0.0
      PNORM(J+1)=PNORM(J)+PH*PH
21    CONTINUE
 
c set up segmented geometry. Each segment represents a box that describes
c the north and south boundaries within a specified longitude range.  The
c full basin is made up by placing the longitudinal segments side by side.
c Thus, it is not possible to have more than one north and south boundary
c at any longitude - excludes most general geometries.

c initialize the default values (rectangular ocean basin)
      XWEST(1) = XWD
      IF (NSEG .EQ. 1) YSOUTH(1) = YSD
      IF (NSEG.EQ.1) YNORTH(1)=YND
 
      DO 30 K=1,NSEG
          IWEST(K)=(XWEST(K)-XWD)/DXD+1.5
          JSOUTH(K)=(YSOUTH(K)-YSD)/DYD+1.5
          JNORTH(K)=(YNORTH(K)-YSD)/DYD+1.5
          SONEK(K)=PONE(JNORTH(K))-PONE(JSOUTH(K))
          SNORM(K)=2.0*(PNORM(JNORTH(K))-PNORM(JSOUTH(K)))
30    CONTINUE
      IWEST(NSEG+1) = NXP
 
c find the boundary type and calculate appropriate transmission coefficiet
 
      WRITE(6,60)
60    FORMAT(//30X,'BASIN CONFIGURATION'//' SEG NO    TYPE    WEST     '
     A,'    SOUTH',8X,'NORTH',10X,'ONEK',6X,'NORM',8X,'TK',5X,'NSTAR'/)
 
      KTYPE(NSEG+1)=1
      KTYPE(1)=6
      STK(1)=0.0
      DO 40 K=1,NSEG
      IF (K.EQ.1) GO TO 62
          IA=JNORTH(K)-JNORTH(K-1)
          IB=JSOUTH(K-1)-JSOUTH(K)
 
C    IA=0 & IB>0,OR, IA>0 & IB=0 IMPLY WESTERN BOUNDARY TYPE
C    IA=0 & IB<0,OR, IA<0 & IB=0 IMPLY EASTERN BOUNDARY TYPE
 
          IF (IA) 41,42,43
41        IF(IB) 44,45,46
44    GO TO 68
45      KTYPE(K)=2
          NSTAR(K)=JNORTH(K)-1
          OFD(K)=1.0-DYH*YY(NSTAR(K))
        GO TO 48
46      GO TO 68
42        IF(IB)49,68,47
49      KTYPE(K)=3
          NSTAR(K)=JSOUTH(K)
          OFD(K)=1.0+DYH*YY(NSTAR(K))
48        DFD(K)=SNORM(K)+OFD(K)*PHIK(NSTAR(K))*(SONEK(K-1)-SONEK(K))
        GO TO 39
47      KTYPE(K)=5
          NSTAR(K)=JSOUTH(K-1)
          OFD(K)=1.0+DYH*YY(NSTAR(K))
        GO TO 38
43        IF(IB) 50,51,52
51      KTYPE(K)=4
          NSTAR(K)=JNORTH(K-1)-1
          OFD(K)=1.0-DYH*YY(NSTAR(K))
38        DFD(K)=SNORM(K-1)+OFD(K)*PHIK(NSTAR(K))*(SONEK(K)-SONEK(K-1))
39        STK(K)=SNORM(K-1)/DFD(K)
        GO TO 62
50      GO TO 68
52    CONTINUE

c Error conditions
68    WRITE(6,69)K,XWEST(K),YSOUTH(K),YSOUTH(K+1),YNORTH(K),YNORTH(K+1)
69    FORMAT(//,' MUST EXIT...ILLEGAL SEGMENT CONFIGURATION'/
     A2X,'AT WEST END OF SEGMENT ',I2,'  X=',F6.2,' YSOUTH OF K
     A,K+1..',2F6.2,'YNORTH OF K,K+1..',2F6.2)
      STOP
668   WRITE(6,669)
669   FORMAT(//'ERROR IN RPEQ INPUT .... SEGMENT OUT OF BOUNDS ...'/)
      STOP
 
62    WRITE(6,61)K,KTYPE(K),IWEST(K),XWEST(K),JSOUTH(K),YSOUTH(K),JNORTH
     A(K),YNORTH(K),SONEK(K),SNORM(K),STK(K),NSTAR(K)
61    FORMAT(I6,I8,3(I6,F7.2),3X,3F10.4,5X,I3)

c check the bounds
       IF ( YSOUTH(K).LT.YSD ) GO TO 668
       IF ( YNORTH(K).GT.YND ) GO TO 668
       if ( K.GT.1 ) THEN
       IF ( XWEST(K).LE.XWEST(K-1) ) GO TO 668
       ENDIF
 
40    CONTINUE
 
c precalculation of frequently used quantities...
 
       ALPHA= DT/DX
       IALPHA = ALPHA
       RALPHA = ALPHA - IALPHA
       IDTX = NX - IALPHA
       ASQP = ( 1+ALPHA*ALPHA )/DT
       ASQM = ( 1-ALPHA*ALPHA )/DT
       ALPHDT = 2.0*ALPHA/DT
      ADY2 = 2.0*ALPHA/DY
      EFRIC = EXP(-RFRIC*DT)
       ASQMI=1.0/ASQM
 
       WRITE(6,220) ADY2,ALPHA,ALPHDT,DYINV,ASQM,ASQMI,ASQP,
     A                EFRIC,PI,SQRT2
220    FORMAT(' ADY2= ',1PG10.4,'  ALPHA= ',1PG10.4,'  ALPHDT=',1PG10.4,
     A  '  DYINV= ',1PG10.4,'  ASQM= ',1PG10.4,'  ASQMI= ',G11.4,
     A  '  ASQP= ',1PG10.4,/,'  EFRIC= ',1PG10.4,
     A  '  PI= ',1PG10.7,'  SQRT2= ',1PG10.7)
 
c calculate the factors a,b,r1 
      DYISQ4 = DYINV*DYINV
      DYISQ8=2.0*DYISQ4
       DO 80 J = 2,NY
          A(J) = C(J-1)-DYISQ4
          R1(J) = C(J-1)+C(J)+DYISQ8
          B(J) = R1(J)  +  4. * ALPHA
80        CONTINUE

c zero out nyp value of a,r1,b, which are not used but show up in the code
      A(NYP)=0.
      R1(NYP)=0.
      B(NYP)=0.
 
210   DO 215 I=1, NXP
      BK(I)=0.E0
      DO 215 J=1, NYP
      F(J,I)=0.E0
      G(J,I)=0.E0
CQ    Q(J,I)=0.E0
215   CONTINUE


      RETURN
       END
