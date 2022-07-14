       SUBROUTINE AKCALC(IX,KBDL,K)

c This routine computes the kelvin wave amplitude by integrating along
c the characteristics.  Points near the westernmost boundary are special;
c that is, points within a distance dx = c dt where c and dt are the kelvin
c wave speed and timestep, respectively.  These points are treated differently.
 
      INCLUDE 'zeq.common'
 
      IMTX=IX-IALPHA
 
       SUM=RALPHA * BK(IMTX - 1)
       IF(IMTX.EQ.IX)GO TO 5
       IXM = IX - 1
       DO 10 I=IMTX,IXM
10          SUM = SUM+BK(I)
 
5      T1 = AKB(IMTX)+RALPHA*( AKB(IMTX-1)-AKB(IMTX))
      IF (KBDL.EQ.0) GO TO 6
       IF (IX.LE.IWEST(K)) GO TO 6
c      stop    
      IPLUS=IWEST(K)+IALPHA
      IF (IX-IPLUS)1,1,2
2     K=K+1
       GO TO 6
1     T1=TR(K)+STK(K)*T1
6      AK(IX) = T1+DX * SUM
       RETURN
       END
