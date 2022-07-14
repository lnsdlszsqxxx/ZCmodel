       SUBROUTINE UHINIT
 
c this routine produces u, h, and v for the initialization timestep

      INCLUDE 'zeq.common'
 
      J1=JSOUTH(1)
      J2=JNORTH(1)-1
      K=2
      DO 1 I=1,NX
 
c *** A,R1,AND D ARE INTERMEDIATE QUANTITIES USED TO SOLVE FOR V
 
      IF (I.LE.IWEST(K)) GO TO 6
          J2=JNORTH(K)-1
          J1=JSOUTH(K)
      K=K+1
6     DO 20 J=J1,J2
20        D(J) = G(J,I)+G(J,I)
 
C *** TRIDAG DOES ACTUAL SOLUTION FOR V
 
      V(J1,I)=0.0
      V(J2+1,I)=0.0
      CALL TRIDAG(A,R1,C,D,V(1,I),J1,J2)   
 
C *** COMPUTATION OF U AND H FOR COLUMN I
 
       DO 30 J=J1,J2
          U(J,I) = YY(J)*( V(J,I)+V(J+1,I) )
          H(J,I) = -DYINV*( V(J+1,I)-V(J,I) )
30     CONTINUE
1     CONTINUE
       RETURN
       END
