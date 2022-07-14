      SUBROUTINE SSTA(T)

c routine for updating ssta; grid is same as atmosphere model, and different
c from ocean dynamics grid.... computations are done dimensionally with
c SSTA expressed in degrees C.  The code is specific to a timestep of 10days
c and other time steps cannot be used correctly without changes to various
c parameters....

      DIMENSION HBAR(34),TD(30,34),DTZM(34),DTZP(30,34)
      DIMENSION DT(30,34,7)

      COMMON/TIME/IT,MP,TY

      COMMON/SSTPAR/GAM1,GAM2,TDA1,TDB1,TDA2,TDB2,TLOSS

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34)
 
      COMMON/ZDAT2/H1(30,34),U1(30,34),V1(30,34)

      DATA DTZM/15*.2,.175,.15,.125,.2,.25,.3,.35,.4,.45,.5,9*.55/

      DATA HBAR/15*1.65,1.7,3*1.75,1.5,1.25,1.,.932,.864,.796,.728,
     + .66,.592,.524,5*.5/
 
      SAVE
 
1234  DO 6663 I=1, 30
      DO 6663 J=1, 34
      DO 6663 N=1, 7
      DT1(I,J,N)=0.
6663  CONTINUE
 
      UFACTR=5.4
      R=.23
 
c calculate the total temperature and mean currents:
      DO 61 I=5, 26
      DO 62 J=6, 32
      WM1(I,J)=WEM(I,J,IT)+TY*(WEM(I,J,MP)-WEM(I,J,IT))
      UV1(I,J)=UV(I,J,1,IT)+TY*(UV(I,J,1,MP)-UV(I,J,1,IT))
      UV2(I,J)=UV(I,J,2,IT)+TY*(UV(I,J,2,MP)-UV(I,J,2,IT))
      TT(I,J)=TO(I,J)+SSTM(I,J,IT)+TY*(SSTM(I,J,MP)-SSTM(I,J,IT))
62    CONTINUE
61    CONTINUE

c calculate the vertical temperature gradient from h1:
      DO 20 I=6, 25
      DO 25 J=6, 32
      IF(H1(I,J) .LT. 0.) GO TO 27
29    TD(I,J)=TDA1*(TANH(TDB1*(HBAR(J)+.01*H1(I,J)))-
     + TANH(TDB1*HBAR(J)))
      GO TO 28
27    TD(I,J)=TDA2*(TANH(TDB2*(HBAR(J)-.01*H1(I,J)))-
     + TANH(TDB2*HBAR(J)))
28    DTZP(I,J)=.1*(TO(I,J)-TD(I,J))
25    CONTINUE
20    CONTINUE

c  compute the surface current anomalies.....constants are such as to give
c  results in cm/sec.  ufactr assumes MLD and upper layer depths of 50m and
c  150m, and the factor of 2. multiplying u1,v1 assumes oceanic equivalent
c  depth is 86 cm.

22    DO 30 I=5, 26
      I1=31-I
      Y=-FLOAT(I-1)*.2+2.9
      A1=UFACTR/(Y*Y+R*R)
      DO 35 J=6, 32
      US(I,J)=U1(I,J)*2.+A1*(R*HTAU(J,I1,1)+Y*HTAU(J,I1,2))
      VS(I,J)=V1(I,J)*2.+A1*(R*HTAU(J,I1,2)-Y*HTAU(J,I1,1))
35    CONTINUE
30    CONTINUE
 
c next, calculate the anomaly in upwelling velocity
c computation assumes currents are calculated on the 2 X 5.635 degree grid.
 
      DO 42 I=1, 30
      DO 42 J=1, 34
      WP(I,J)=0.
42    CONTINUE
      DX=.5625*2.
      DY=.2*2.
      DO 43 I=6, 25
      DO 44 J=7, 31
      WP(I,J)=.045*(US(I,J+1)-US(I,J-1))/DX+.045*(VS(I-1,J)-
     A VS(I+1,J))/DY
44    CONTINUE
      WP(I,32)=.045*(-US(I,31))/DX+.045*(VS(I-1,32)-VS(I+1,32))/DY
      WP(I,6)=WP(I,7)
43    CONTINUE
 
c calculate the change in sst, using variable sub-timestep to insure
c numerical stability
 
      CALL GETN(NSST)
      WRITE(6,3554) NSST
3554  FORMAT(1X, 'NUMBER OF SST LOOPS IS: ' I2)

      DO 50 N=1, NSST
      DO 51 I=6, 25
      DO 52 J=6, 32
      IF(J.EQ.32.OR.J.EQ.6) GO TO 45

c zonal advection
      DT(I,J,1)=.008*ADV(1,US,TT,I,J)
      DT(I,J,1)=DT(I,J,1)+.08*ADV(1,UV1,TO,I,J)
      DT(I,J,2)=0.
      GO TO 46
45    DT(I,J,1)=0.
      DT(I,J,2)=0.

c meridional advection
46    DT(I,J,3)=.008*ADV(2,VS,TT,I,J)
      DT(I,J,3)=DT(I,J,3)+.08*ADV(2,UV2,TO,I,J)
      DT(I,J,4)=0.
      A1=WM1(I,J)
      A2=WP(I,J)

c upwelling
      DT(I,J,5)=-.864*HF(A1)*DTZP(I,J)*GAM1
      DT(I,J,6)=-.864*GF(A1,A2)*(DTZM(J)+DTZP(I,J))*GAM2

c surface heat flux
!Newton's law of cooling
      DT(I,J,7)=(TLOSS-1.)*TO(I,J)*4.
52    CONTINUE
51    CONTINUE

c add up the different tendencies to update ssta
!this is Eq A11 in ZC1987,  liang
      DO 53 J=6, 32
      DO 54 I=6, 25
      TO(I,J)=TO(I,J)+(DT(I,J,1)+DT(I,J,3)+DT(I,J,5)+DT(I,J,6)+
     + DT(I,J,7))/FLOAT(NSST)
54    CONTINUE
53    CONTINUE

c update total sst and do not allow anomalies that cause total sst to
c exceed 30 degrees C.
!SSTM is clmt SST, read from rcsstmn.data in initdat.f
      DO 70 J=6, 32
      DO 71 I=5, 26
      TT(I,J)=TO(I,J)+SSTM(I,J,IT)+TY*(SSTM(I,J,MP)-SSTM(I,J,IT))
      IF(TT(I,J) .LE. 30.) GO TO 71
      TO(I,J)=TO(I,J)+30.-TT(I,J)
      TT(I,J)=30.
71    CONTINUE

c update the vertical temp. gradient
!this is (T-Te)/H1 in A11, and Te in A12. ZC1987
      DO 72 I=6, 25
      DTZP(I,J)=.1*(TO(I,J)-TD(I,J))

c save the temperature tendency terms
      DO 73 M=1, 7
      DT1(I,J,M)=DT1(I,J,M)+DT(I,J,M)/FLOAT(NSST)
73    CONTINUE
72    CONTINUE
70    CONTINUE
50    CONTINUE
 
 
c  save total currents in uv1,uv2,wm1
      DO 81 J=6, 32
      DO 82 I=6, 25
      WM1(I,J)=WM1(I,J)+WP(I,J)
82    CONTINUE
      DO 83 I=5, 26
      UV1(I,J)=10.*UV1(I,J)+US(I,J)
      UV2(I,J)=10.*UV2(I,J)+VS(I,J)
83    CONTINUE
81    CONTINUE

C     NOW GET THE U1,V1,H1 FIELDS FOR THE NEXT SST ITERATION
 
215   CALL ZAVG(H1,U1,V1)

      RETURN
      END


      SUBROUTINE GETN(NSST)

c routine to determine number of subintervals to divide 10 days into in
c order to assure numerical stability by CFL type criterion

      INCLUDE 'zeq.common'

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34)

      COMMON/SGRID/XWRS,DXRS,YNRS,DYRS,NYPS
      SAVE

      DX1=DXRS*1.11E7
      DY1=DYRS*1.11E7
      DTT=DTD*2.628E6
      NSST=2

      DO 10 I=6, 25
      DO 11 J=6, 32
      AAA=ABS(US(I,J))*DTT/DX1
      BBB=ABS(VS(I,J))*DTT/DY1
      IAA=1+AAA
      IBB=1+BBB
      NSST=MAX0(NSST,IAA)
      NSST=MAX0(NSST,IBB)
11    CONTINUE
10    CONTINUE

      RETURN
      END


      SUBROUTINE ZAVG(H1,U1,V1)

c routine to average fields from ocean dynamics grid onto the coarser sst grid
      INCLUDE 'zeq.common'
      DIMENSION H1(30,34),U1(30,34),V1(30,34)
      DO 10 I=6, 25
      I1=120-4*I
      IS=I1-2
      IE=I1+2
      DO 20 J=6, 32
      X=FLOAT(J-1)*.5625
      J1=(10.*X+.25-20.)/2.
      JS=J1-1
      JE=J1+1
      AA=0.
      AA1=0.
      AA2=0.
      DO 21 II=IS, IE
      DO 22 JJ=JS, JE
      AA=AA+HB(II,JJ)+AKB(JJ)*PHIK(II)
      AA1=AA1+UB(II,JJ)+AKB(JJ)*PHIK(II)
      AA2=AA2+V(II,JJ)
22    CONTINUE
21    CONTINUE
      H1(I,J)=AA/15.
      U1(I,J)=AA1/15.
      V1(I,J)=AA2/15.
20    CONTINUE
      H1(I,33)=H1(I,32)
      H1(I,34)=H1(I,32)
      U1(I,33)=0.
      U1(I,34)=0.
      V1(I,33)=0.
      V1(I,34)=0.
      DO 15 J=1, 5
      H1(I,J)=H1(I,6)
15    CONTINUE
10    CONTINUE
      DO 12 J=1, 34
      DO 13 I=1, 5
      H1(I,J)=0.
      H1(31-I,J)=0.
      U1(I,J)=0.
      U1(31-I,J)=0.
      V1(I,J)=0.
      V1(31-I,J)=0.
13    CONTINUE
12    CONTINUE
 
      RETURN
      END


 
      REAL FUNCTION HF(X)
      IF(X .LE. 0.) HF=0.
      IF(X .GT. 0.) HF=X
      RETURN
      END
 
 
      REAL FUNCTION GF(X,Y)
      Z=X+Y
      IF(X .LE. 0. .AND. Z .LE. 0.) GF=0.
      IF(X .LE. 0. .AND. Z .GT. 0.) GF=Z
      IF(X .GT. 0. .AND. Z .LE. 0.) GF=-X
      IF(X .GT. 0. .AND. Z .GT. 0.) GF=Y
      RETURN
      END
 
 
      REAL FUNCTION ADV(ID1,U,T,I,J)
c routine for evaluating advection by modified upwind differencing scheme
      DIMENSION U(30,34),T(30,34)
      DX=2.*.5625
      DY=2.*.2
      IF(ID1.EQ.2) GO TO 11
      IF(U(I,J).GE.0.) ADV=(U(I,J-1)+U(I,J))*(T(I,J-1)-T(I,J))/DX
      IF(U(I,J).LT.0.) ADV=(U(I,J)+U(I,J+1))*(T(I,J)-T(I,J+1))/DX
      RETURN
11    IF(U(I,J).GE.0.) ADV=(U(I+1,J)+U(I,J))*(T(I+1,J)-T(I,J))/DY
      IF(U(I,J).LT.0.) ADV=(U(I-1,J)+U(I,J))*(T(I,J)-T(I-1,J))/DY
      RETURN
      END

