       SUBROUTINE UHCALC(I,K)
 
c this routine produces u, h, and v for all but the initialization timestep

      INCLUDE 'zeq.common'
 
c *** GET R1 AND R3
 
      J1=JSOUTH(K)
      J2=JNORTH(K)-1
 
      IF (I.NE.IWEST(K)) GO TO 8
      IF(KTYPE(K) .LT. 4) GO TO 8
 
      DO 11 J = J1,J2
 
         T1 = F(J,I)+F(J,I)
CQ       T1 = T1+ 2.0*ALPHA*Q(J,I)
          T2 = ASQM*UB(J,I+1) - ASQP*( U(J,I+1)-UBNDY(J,K) )
          T3 = ALPHDT*( H(J,I+1)-HBNDY(J,K) )
          R1(J) = T1+T2-T3
 
         T1 = (ALPHA+ALPHA)*F(J,I)
CQ       T1 = T1+ Q(J,I)+Q(J,I)
          T2 = ASQM*HB(J,I+1) - ASQP*( H(J,I+1)-HBNDY(J,K) )
          T3 = ALPHDT*( U(J,I+1)-UBNDY(J,K) )
          R3(J) = T1+T2-T3
11        CONTINUE
       GO TO 12
 
  8      DO 10 J = J1,J2
 
         T1 = F(J,I)+F(J,I)
CQ       T1 = T1+ 2.0*ALPHA*Q(J,I)
          T2 = ASQM*UB(J,I+1) - ASQP*( U(J,I+1)-UB(J,I) )
          T3 = ALPHDT*( H(J,I+1)-HB(J,I) )
          R1(J) = T1+T2-T3
 
         T1 = (ALPHA+ALPHA)*F(J,I)
CQ       T1 = T1+ Q(J,I)+Q(J,I)
          T2 = ASQM*HB(J,I+1) - ASQP*( H(J,I+1)-HB(J,I) )
          T3 = ALPHDT*( U(J,I+1)-UB(J,I) )
          R3(J) = T1+T2-T3
10        CONTINUE
 
 12      CONTINUE
 
C *** GET V FOR COLUMN I
 
      V(J1,I)=0.0
      V(J2+1,I)=0.0
      J3=J1+1
       DO 20 J=J3,J2
          T1 = 2.0 * ASQM * G(J,I)
          T2 = YY(J-1)*R1(J-1)+YY(J)*R1(J)
          T3 = DYINV*( R3(J)-R3(J-1) )
          D(J) = T1-T2-T3
20     CONTINUE
 
c   solve the tridiagnonal system in tridag
      CALL TRIDAG(A,B,C,D,V(1,I),J1,J2)
 
c  now calculate u and h
 
      DO 30 J=J1,J2
         U(J,I) = ASQMI*( R1(J)+V(J,I)*( YY(J)+ADY2 )
     A             +V(J+1,I)*( YY(J)-ADY2 ) )
         H(J,I) = ASQMI*( R3(J)+V(J,I)*( DYINV+ALPHA*YY(J) )
     A             +V(J+1,I)*(ALPHA*YY(J)-DYINV) )
30       CONTINUE

 
c  remove the spurious kelvin wave component
 
60    AKELV=0.0
      DO 40 J = J1,J2
40    AKELV = AKELV+(U(J,I)+H(J,I))*PHIK(J)
      AKELV=AKELV/SNORM(K)
 
      DO 50 J = J1,J2
      U(J,I) = U(J,I)-AKELV*PHIK(J)
50    H(J,I) = H(J,I)-AKELV*PHIK(J)
 
      RETURN
      END
