       SUBROUTINE CFORCE
 
c this routine supplies wind stress forcing to the ocean dynamics model
c For coupled ocean/atmosphere model the wind stress is derived from
c the atmosphere model, and is saved in the array HTAU.

      INCLUDE 'zeq.common'
 
      DIMENSION  HT(34,30),XOFF(2),YOFF(2)

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34)

c following are the wind field parameters
      DATA   WXW/+101.25/,WXE/+286.875/,WYS/-29./,WYN/+29./
      DATA    IW/34/,JW/30/

c  offsets for the staggering of F and G
      DATA XOFF/0.5E0,1.E0/,YOFF/0.5E0,1.E0/
      DATA LP/1/

      SAVE

      L1=LP
      IF(ISTART.NE.0) GO TO 100

c......................................................................
c some initialization computations. Gridding factors...

      WDX=(WXE-WXW)/(IW-1)
      AWX=1.E0+(XWD-WXW)/WDX
      WDXS=ELEQ*DX/WDX
      WDY=(WYN-WYS)/(JW-1)
      AWY=1.E0+(YSD-WYS)/WDY
      WDYS=ELEQ*DY/WDY
 
c on the initial time step, only G is needed, unless restarting.
      L1 = 2
      IF(NSTART.GT.0 .AND. NATMR .GT. 0) L1=1
 
      ISM=1
      IEM=NX
      JSM=1
      JEM=NY
      LP=1
      L2=2

   50 WRITE(6,51)
      WRITE(6,54) LP,L2,ISM,IEM,JSM,JEM
   51 FORMAT(/,30X,'********** WIND FIELD PARAMETERS ********'/)
   54 FORMAT(/
     A 'L1 =',I1,' L2 =',I1,' WEST,EAST,SOUTH,NORTH LIMITS ...
     A',4I8//)
 
c....................................................................
c this is the main loop section
100   CONTINUE

c zero out bk to start
      DO 105 I=1,NXP
      BK(I)=0.E0
105   CONTINUE

      IF(L1.GT.L2) RETURN
 
c get the wind stress by accessing the atmosphere model array HTAU...
c the stress data are interpolated onto ocean grids...
c the code assumes the wind data domain includes the entire ocean domain

      DO 130 L=L1,L2
!      ASSIGN 119 TO LFG
!      IF(L.EQ.2) ASSIGN 118 TO LFG
      LFG=2 !liang
      IF(L.EQ.2) LFG=1 !liang
      AAWX=AWX-XOFF(L)*WDXS
      AAWY=AWY-YOFF(L)*WDYS
 
      DO 110 J=1,JW
      J13=JW+1-J
      DO 111 I=1, IW
      HT(I,J)=WSCALE*HTAU(I,J,L)
111   CONTINUE
110   CONTINUE
 
      K = 1
      IEG = IEM + (L-1)

      DO 125 I=ISM,IEG
      RX=AAWX+WDXS*I
      IX=RX
      RX=RX-IX

c ocean basin segment boundary check
112   IF(I.LT.IWEST(K)) GO TO 116

C     CROSSED A SEGMENT BOUNDARY
      IF(L.EQ.1) GO TO 114

C     A G POINT. USE THE LONGER PIECE
      IF(I.GT.IWEST(K)) GO TO 114
      IF(KTYPE(K).LE.3) GO TO 116
114   JSS = MAX0(JSM,JSOUTH(K))
      JEG = MIN0(JEM,JNORTH(K)-1) +L-1
      K = K+1
      GO TO 112

116   DO 120 J=JSS,JEG
      RY=AAWY+WDYS*J
      JY=RY
      RY=RY-JY
      FA=HT(IX  ,JY)+RY*(HT(IX  ,JY+1)-HT(IX  ,JY))
      FB=HT(IX+1,JY)+RY*(HT(IX+1,JY+1)-HT(IX+1,JY))
      FA = FA + RX * ( FB - FA )
!      GO TO LFG,(118,119)
      GO TO (118,119) LFG
118   G(J,I)=FA
      GO TO 120
119   F(J,I)=FA

c compute bk, the projection of F and Q on kelvin mode
      BK(I)=BK(I)+FA*PHIK(J)
120   CONTINUE

125   CONTINUE
130   CONTINUE
 
c  normalize bk
 
      K=1
      DO 250 I=ISM,IEM
212   IF(I.LT.IWEST(K)) GO TO 216
      K=K+1
      GO TO 212
216   BK(I)=BK(I)/SNORM(K-1)
250   CONTINUE
      RETURN
      END
