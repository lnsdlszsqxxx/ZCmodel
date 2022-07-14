         SUBROUTINE BNDARY(K,I)
 
c This routine calculates dynamical variables at boundaries....

      INCLUDE 'zeq.common'
 
      KTYP=KTYPE(K)
      GO TO (10,20,30,40,50,60),KTYP
 
c calculate U and H along far eastern boundary
 
10    JTOP=JNORTH(K-1)-1
      JBOT=JSOUTH(K-1)
      JBOT1=JBOT+1
      GINT(JBOT)=0.0
      GG=0.0
      DO 12 J=JBOT1,JTOP
        GINT(J)=GINT(J-1)+G(J,NXP)*DY
12    GG=GG+GINT(J)*PHIK(J)
 
      HKG=(AK(NXP)*SNORM(K-1)-GG)/SONEK(K-1)
      DO 15 J=JBOT,JTOP
      U(J,NXP)=-AK(NXP)*PHIK(J)
15    H(J,NXP)=HKG+U(J,NXP)+GINT(J)
      RETURN
 
c handle eastern boundary types
c                                        **
c    (ktype=2)                                    ***
c                                        ********
20    LBOT=JSOUTH(K)
      LTOP=NSTAR(K)
      JBOT=LTOP+1
      JTOP=JNORTH(K-1)-1
      HRNUM=H(LTOP,I)-DYH*YY(LTOP)*U(LTOP,I)
      GO TO 25
C
C                                                **********
C    (ktype=3)                                   ***
C                                                **
30    LBOT=NSTAR(K)
      LTOP=JNORTH(K)-1
      JBOT=JSOUTH(K-1)
      JTOP=LBOT-1
      HRNUM=H(LBOT,I)+DYH*YY(LBOT)*U(LBOT,I)
25    GG=0.0
      GINT(JBOT-1)=0.0
      JBOT1=JBOT+1
      DO 26 J=JBOT,JTOP
      GINT(J)=GINT(J-1)+G(J,I)*DY
26    GG=GG+GINT(J)*PHIK(J)
      IF (KTYP.EQ.3) GG=-GG
 
      TR(K)=-(HRNUM*(SONEK(K-1)-SONEK(K))+GG)/DFD(K)
 
c    calculate U and H for these cases
 
      D1=HRNUM+(AK(I)*STK(K)+TR(K))*OFD(K)*PHIK(NSTAR(K))
      DO 27 J=JBOT,JTOP
      H(J,I)=D1-AK(I)*PHIK(J)+GINT(J)
27    U(J,I)=-AK(I)*PHIK(J)
 
      DO 28 J=LBOT,LTOP
      H(J,I)=H(J,I)+(TR(K)+AK(I)*(STK(K)-1.0))*PHIK(J)
28    U(J,I)=U(J,I)+(TR(K)+AK(I)*(STK(K)-1.0))*PHIK(J)
 
      RETURN
 
 
c     handle the western boundary types
c                                        **
c    (ktype=4)                          ***
c                                **********
40    LTOP=NSTAR(K)
      LBOT=JSOUTH(K-1)
      JTOP=JNORTH(K)-1
      JBOT=LTOP+1
      JVBOT=LTOP
      JVTOP=JTOP
      GO TO 45
c                              **********
c    (ktype=5)                        ***
c                                      **
50    LTOP=JNORTH(K-1)-1
      LBOT=NSTAR(K)
      JTOP=LBOT-1
      JBOT=JSOUTH(K)
      JVBOT=JBOT
      JVTOP=LBOT
45    SHU=0.0
      DO 51 J=LBOT,LTOP
51    SHU=SHU+(H(J,I)+U(J,I))*PHIK(J)
      SU=0.0
      GINT(JBOT)=U(JBOT,I)
      JBOT1=JBOT+1
      DO 52 J=JBOT1,JTOP
52    GINT(J)=GINT(J-1)+U(J,I)
 
      TR(K)=-(SHU+OFD(K)*PHIK(NSTAR(K))*GINT(JTOP))/DFD(K)
 
      ROSS(K)=AK(I)*STK(K)+TR(K)
 
c  Save values without transmission coefficient for the next time step
 
      J1=JSOUTH(K)
      J2=JNORTH(K)-1
      DO 205 J=J1, J2
      HBNDY(J,K)= H(J,I)*EFRIC
      UBNDY(J,K)= U(J,I)*EFRIC
205   CONTINUE
 
      DYODX=DY/DX
      UINT=U(JVBOT,I)+ROSS(K)*PHIK(JVBOT)
      JVBOT1=JVBOT+1
      JVTOP1=JVTOP-1
      DO 250 J=JVBOT1, JVTOP1
      V(J,I)=V(J,I)-DYODX*UINT
250   UINT=UINT+U(J,I)+ROSS(K)*PHIK(J)
      V(JVTOP,I)=V(JVTOP,I)-DYODX*UINT
      IP1=I+1
      DO 260 J=JBOT, JTOP
      U(J,I)=-AK(I)*PHIK(J)
      H(J,I)=HB(J,IP1)/EFRIC-DX*
     1 ((V(J,I)+V(J+1,I))*.5*YY(J)+F(J,I))+(AK(IP1)-AK(I))*PHIK(J)
260   CONTINUE

      DO 203 J=LBOT,LTOP
      U(J,I)=U(J,I)+(ROSS(K)-AK(I))*PHIK(J)
      H(J,I)=H(J,I)+(ROSS(K)-AK(I))*PHIK(J)
203   CONTINUE
 
      ROSS(K)=ROSS(K)*(SONEK(K)-SONEK(K-1))+GINT(JTOP)
      U(NSTAR(K),I) = U(NSTAR(K),I)+ROSS(K)
      IF (KTYP.EQ.4) ROSS(K)=-ROSS(K)
      H(NSTAR(K),I) = H(NSTAR(K),I)+ROSS(K)*DYH*YY(NSTAR(K))
 
      RETURN
 
c calculate the kelvin amplitude at the far western boundary due to 
c incoming rossby waves

60    J1=JSOUTH(1)
      J2=JNORTH(1)-1
      SUM=0.0
      DO 90 J=J1,J2
      SUM=SUM+U(J,1)
      HBNDY(J,K)=H(J,I)*EFRIC
      UBNDY(J,K)=U(J,I)*EFRIC
90    CONTINUE
      AK(1)=-SUM/SONEK(1)
 
      JBOT=J1
      JBOT1=JBOT+1
      JTOP=J2
      DYODX=DY/DX

      UINT=U(JBOT,I)+AK(I)*PHIK(JBOT)
      DO 350 J=JBOT1, JTOP
      V(J,I)=V(J,I)-DYODX*UINT
350   UINT=UINT+U(J,I)+AK(I)*PHIK(J)
      IP1=I+1
      DO 360 J=JBOT, JTOP-1
      U(J,I)=-AK(I)*PHIK(J)
      H(J,I)=HB(J,IP1)/EFRIC-DX*
     1 ((V(J,I)+V(J+1,I))*.5*YY(J)+F(J,I))+(AK(IP1)-AK(I))*PHIK(J)
360   CONTINUE
      J=JTOP
      U(J,I)=-AK(I)*PHIK(J)
      H(J,I)=HB(J,IP1)/EFRIC-DX*
     1 (V(J,I)*.5*YY(J)+F(J,I))+(AK(IP1)-AK(I))*PHIK(J)

      RETURN
      END
