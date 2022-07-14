       SUBROUTINE MLOOP
 
c this routine cycles through the processes of the ocean dynamics calculation
c for each time step

      INCLUDE 'zeq.common'
 
       IF (ISTART.EQ.0) THEN
c  model initialization;  initialize U,H,V in uhinit
 
          CALL UHINIT
          CALL AKCALC( NXP,0,NSEG+1 )
          CALL BNDARY( NSEG+1,NXP )
          CALL BNDARY( 1,1 )
          ISTART=1
          GO TO 110
        ENDIF
 
C********************* TOP OF THE MAIN LOOP ****************************
 
        K=NSEG+1
       I=NXP
 
50     IF (I.NE.1) CALL AKCALC( I,0,K )
       CALL BNDARY ( K,I )
       K=K-1
       IF ( K.LT.1 ) GO TO 30
 
40    I=I-1
      CALL UHCALC(I,K)
      IF (I-IWEST(K)) 50,50,40
 
 30    CONTINUE
 
c calculation of ak - kelvin wave amplitude
 
c first, get values near western boundary
      IF (IALPHA .EQ. 0) GO TO 102
      AST = (AKB(1) - AK(1))/ALPHA
      DO 15 I = 1,IALPHA
15       AK(I+1) = AK(I) + AST + BK(I)*DX
102   CONTINUE

c now the remaining values
      ILEF = IALPHA + 2
      K=2
      DO 101 IX = ILEF,NX
101     CALL AKCALC(IX,1,K)
 
 
c Save ak,u,h for next time step computations
 
  110 I2 = 0
 
      DO 130 K = 1,NSEG
      I1 = I2+1
      I2 = IWEST(K+1)
      IF (KTYPE(K+1).GE.4) I2=I2-1
      J1=JSOUTH(K)
      J2=JNORTH(K)-1
         DO 130 I=I1,I2
         AKB(I) = AK(I)*EFRIC
             DO 120 J=J1,J2
             UB(J,I) = U(J,I)*EFRIC
             HB(J,I) = H(J,I)*EFRIC
120          CONTINUE
130    CONTINUE

 
c write diagnostic output line to standard output....

      IF(MOD(NT,NPRINT).EQ.0) WRITE(6,135) NT,TD,AKB(1),AKB(NXP)
  135 FORMAT(' TIMESTEP',I6,' TIME =',F8.3,' MONTHS. AK(1),AKXE...',
     X 1PD12.4,D12.4)
 
       RETURN
       END
