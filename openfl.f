      SUBROUTINE OPENFL
 
c routine to open appropriate files for coupled model.

      CHARACTER*45 FN61,FN62,FN51,FN52
 
      READ(5,111) FN61
111     FORMAT(A45)
      OPEN(UNIT=61, FILE=FN61, FORM='UNFORMATTED')
      REWIND 61
      READ(5,112) FN62
112     FORMAT(A45)
      OPEN(UNIT=62, FILE=FN62, FORM='UNFORMATTED')
      REWIND 62
      OPEN (UNIT=82,FILE='rcwindmn.data', FORM='UNFORMATTED')
      REWIND 82
      OPEN(UNIT=83,FILE='rcdivmn.data', FORM='UNFORMATTED')
      REWIND 83
      OPEN(UNIT=84,FILE='rcsstmn.data', FORM='UNFORMATTED')
      REWIND 84
      OPEN ( UNIT=87, FILE='wem.zeb', FORM='UNFORMATTED' )
      REWIND 87
      OPEN ( UNIT=89, FILE='uv.zeb', FORM='UNFORMATTED' )
      REWIND 89
 
      RETURN
      END
