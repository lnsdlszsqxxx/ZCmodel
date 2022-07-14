c****************************zeq1main.f*********************************
c  driver routine for coupled model
c note time convention in months; 0.5 corresponds to mid January 1960.
c***********************************************************************

      INCLUDE 'zeq.common'

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34)

      COMMON/TIME/IT,MP,TY
 
      COMMON/INIT/IC2

      DATA TFIND/1.E20/
c wo jia de
!no need to open file here, use them as follow:
!./ocn3.out < fctest.data >& test.out
!liang 2022.6.29
      OPEN(UNIT=5,FILE='fctest.data',FORM='FORMATTED',status='old')
      REWIND 5 
      OPEN(UNIT=6,FILE='test.out',FORM='FORMATTED',status='replace')
      REWIND 6 
c read case title information and write to standard output
      READ(5,450) HEADER 
450   FORMAT(80A1)                                                      
      WRITE(6,475) HEADER
475   FORMAT('  CASE TITLE: ',80A1/)
 
c first, open necessary files...........................................
      CALL OPENFL
 
c read control paramters from input data file on unit 5.................
      NATM=2
      READ(5,40,END=999) NSTART, TFIND
      IF(NSTART.LT.0) GO TO 999
      WRITE(6,50) NSTART, TFIND
40    FORMAT ( 8X, I15 / 8X, F15.0 )
50    FORMAT ( / I15, ' NSTART ' / E15.7, ' TFIND ' / )

c if history input is being used for startup, call rdhist...............
      IF(NSTART.GT.0) CALL RDHIST(0,TFIND)
      NATMR=NATM

c now pick up (possibly) new ocean model parameters unless nstart=1.....
      IF(NSTART.NE.1) CALL SETUP
 
c if nstart=1 or nstart=2, time is gotten from the history file.........
c otherwise (nstart=3) time is set to tzero from input data file........
      IF(NSTART .EQ. 3 .OR. NSTART .EQ. 0) THEN
      TD=TZERO
      ELSE
      TZERO=TD
      ENDIF

      NT=0
      IT1=TD+120.501
      IT=MOD(IT1, 12)
      MP=IT+1
      IF(IT.EQ.0) IT=12
      TY=TD+120.501-IT1

c read in ssta model and atmosphere model constants, etc...
       CALL SETUP2 

c read in mean field data...............................................
       
       CALL INITDAT

c do initialization of ocean dynamics model constants, etc..............
       CALL CONSTC

c initial stress field must be recalculated from winds saved in history.
      ISTART=0
      CALL STRESS(TD,ISTART)

c if one wanted to overwrite any fields with data, it would be done here
c based on the paramter IC2... routine INSERT would be created to do it.
c     IF(IC2 .GT. 0) CALL INSERT

c initial call to ocean force routine and dynamics update routine mloop
c are special if starting up from rest - ie, if not restarting..........
      CALL CFORCE
 
      IF(NSTART .GT. 0) THEN
      ISTART=1
      GO TO 400
      ELSE
      GO TO 350
      ENDIF

c.......................................................................
c this is the top of main loop for all but the first cycle

c update the ssta field.................................................
300   CALL SSTA(TD)

c update the winds......................................................
      CALL ZATMC(TD)

c compute the wind stress anomalies.....................................
      CALL STRESS(TD,ISTART)

c insert the stress forcing in ocean model..............................
      CALL CFORCE
 
c update the ocean dynamics.............................................
350   CALL MLOOP

400   CONTINUE

c///////////////////////////////////////////////////////////////////////
c  here one can add any desired output write statements.................
c///////////////////////////////////////////////////////////////////////

!TZERO the initial month for the run, 0.5 means middle of the month
!TENDD the final month for the run. 
!NTAPE=3 means write hist file every 3 steps, i.e. monthly output, not the
!monthly mean but the last step in each month.
!NTAPE=0 means only one output at the end
!NT is the current number of step, counter in main 
!NTIMES is the total steps, calculated from TZERO and TENDD
!TPLMIN: earliest time that history output will be written. see fctest.data

c write history output if conditions are met............................
      IF(NTAPE.NE.0.AND.MOD(NT,NTAPE).EQ.0.AND.NT.NE.0.AND.
     * TD .GE. TPLMIN) then
      
!status=NEW, can avoid replace the old one accidentaly
!remove the old file menually
!save hist in little endian format
      if(NT.eq.NTAPE) then !first time saving
        iwrite=1
        open(unit=1000,file='ZCCPD_SSTA_34x30_Deg.bin',status='NEW',
     &convert='little_endian',access='direct',recl=34*30*4)
!        open(unit=1001,file='ZCCPD_TAUA_34x30_dyne.bin',status='NEW',
!     &convert='little_endian',access='direct',recl=34*30*2*4)
      end if

!Note that the order is different from the model, J (x direction) goes first
!latitude and running from north to south, reverse it
!HTAU in a diff order
!it is 'direct', not 'sequential', easier for plotting
       write(1000,rec=iwrite) ((TO(I,J),J=1, 34),I=30,1,-1)
!       write(1001,rec=iwrite) HTAU !no need to change order
!       write(1001,rec=iwrite) ((UO(I,J),J=1, 34),I=30,1,-1)
!       write(1002,rec=iwrite) ((VO(I,J),J=1, 34),I=30,1,-1)
!       write(1002,rec=iwrite) ((H1(I,J),J=1, 34),I=30,1,-1)
       iwrite=iwrite+1

      END IF
 
c check to see if it's time to stop
      IF(NT.GE.NTIMES) GO TO 900

c if not, loop again
      NT=NT+1
      TD=TZERO+FLOAT(NT)*DTD
      IT1=TD+120.501
      IT=MOD(IT1, 12)
      MP=IT+1
      IF(IT.EQ.0) IT=12
      TY=TD+120.501-IT1
      GO TO 300

c time to finish up
900   CALL WRHIST
      close(1000)
!      close(1001)
 
       WRITE(6,998)
998    FORMAT(' **************** NORMAL END-OF-JOB ZEQ *************')
 
999   CONTINUE 
      CLOSE(5)
      CLOSE(6)
       END
 
