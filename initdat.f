       SUBROUTINE INITDAT
 
c this routine reads in and sets initial values of data arrays for
c the coupled model.

      INCLUDE 'zeq.common'
 
      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34)

c read in the mean wind, sst, divergence, and currents data
 
      DO 302 L=1, 12
      READ(84) ((SSTM(I,J,L),I=1, 30),J=1, 34)
      READ(83) ((DIVM(I,J,L),I=1, 30),J=1, 34)
      READ(82) (((WM(I,J,K,L),I=1, 30),J=1, 34),K=1,2)
      READ(87) ((WEM(I,J,L),I=1, 30), J=1, 34)
      READ(89) (((UV(I,J,K,L),I=1, 30),J=1, 34),K=1,2)
302   CONTINUE
      REWIND 82
      REWIND 83
      REWIND 84
      REWIND 87
      REWIND 89
 
c if not restarting, zero out arrays q0o,uo,vo,do,to

      IF (NSTART .GT. 0 .AND. NATMR .GT. 0) GO TO 311
 
      DO 310 I=1, 30
      DO 310 J=1, 34
      Q0O(I,J)=0.E0
      UO(I,J)=0.E0
      VO(I,J)=0.E0
      DO(I,J)=0.E0
      TO(I,J)=0.E0
310   CONTINUE
311   CONTINUE

      RETURN
      END

