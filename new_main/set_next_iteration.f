	  IF(MYPE .EQ. 0)WRITE(6,*)'Start of end iteration section',MYPE,MAXCH,FIXED_T
!
!*****************************************************************************
!*****************************************************************************
!                END OF MAIN ITERATION LOOP
!*****************************************************************************
!*****************************************************************************
!
! Close units 2 and 16 to force writing of information. We check if OUTGEN
! has been opened --- if not we are writing to the terminal and nothing
! needs to be done.
!
	  INQUIRE(FILE='OUTGEN',OPENED=FILE_OPEN)
	  IF(FILE_OPEN)THEN
	    CLOSE(UNIT=LUER)
	    CALL GEN_ASCI_OPEN(LUER,'OUTGEN','OLD','APPEND',' ',IZERO,IOS)
	    CALL SET_LINE_BUFFERING(LUER)
	  END IF
	  CLOSE(UNIT=LU_SE)
	  IF(MYPE .EQ. 0)CALL GEN_ASCI_OPEN(LU_SE,'STEQ_VALS','OLD','APPEND',' ',IZERO,IOS)
!
! Adjust X-ray filling factors upwards if within a factor of 100 of convergence.
! We adjust MAXCH to ensure that NUM_ITS_TO_DO is not set to 1.
!
	  IF(XRAYS .AND. ADD_XRAYS_SLOWLY .AND. RD_LAMBDA .AND. MAXCH .LT. 100)THEN
	     IF(FILL_FAC_XRAYS_1 .NE. FILL_X1_SAV .OR.  FILL_FAC_XRAYS_2 .NE. FILL_X2_SAV)THEN
	       FILL_FAC_XRAYS_1=MIN(FILL_FAC_XRAYS_1*SLOW_XRAY_SCL_FAC,FILL_X1_SAV)
	       FILL_FAC_XRAYS_2=MIN(FILL_FAC_XRAYS_2*SLOW_XRAY_SCL_FAC,FILL_X2_SAV)
	       IF(MYPE .EQ. 0)THEN
	         WRITE(LUER,*)'Have adjusted X-ray values closer to the desired values'
	         WRITE(LUER,*)'Current filling factor is (1st component)',FILL_FAC_XRAYS_1
	         WRITE(LUER,*)'Current filling factor is (2nd component)',FILL_FAC_XRAYS_2
	         WRITE(LUER,*)'Current on desired (1st component)=',FILL_FAC_XRAYS_1/FILL_X1_SAV
	         WRITE(LUER,*)'Current on desired (2nd component)=',FILL_FAC_XRAYS_2/FILL_X2_SAV
	         CALL UPDATE_KEYWORD(FILL_FAC_XRAYS_1,'[XFI1_BEG]','VADAT',L_TRUE,L_FALSE,LUIN)
	         CALL UPDATE_KEYWORD(FILL_FAC_XRAYS_2,'[XFI2_BEG]','VADAT',L_FALSE,L_TRUE,LUIN)
	         MAXCH=100			!To force run to continue
	       END IF
	     ELSE
	       IF(MYPE .EQ. 0)CALL UPDATE_KEYWORD(L_FALSE,'[XSLOW]','VADAT',L_TRUE,L_TRUE,LUIN)
	       IF(DO_LAMBDA_AUTO)THEN
	         RD_LAMBDA=.FALSE.
	         LAMBDA_ITERATION=.FALSE.
	         IF(MYPE .EQ. 0)CALL UPDATE_KEYWORD(L_FALSE,'[DO_LAM_IT]','IN_ITS',L_TRUE,L_TRUE,LUIN)
	       END IF
	    END IF
	  ELSE IF(XRAYS .AND. MAXCH .LT. 100.0_LDP .AND. .NOT. ADD_XRAYS_SLOWLY)THEN
!
! We do not do the scaling if there is intrinsic X-ray emssion from the star.
! DESIRED_XRAY_LUM should be in units of LSTAR.
!
	    IF(OBS_XRAY_LUM_0P1 .LT. SUM(XRAY_LUM_0P1) .AND. SCALE_XRAY_LUM)THEN
	      IF( (LUM*DESIRED_XRAY_LUM/OBS_XRAY_LUM_0P1-1.0_LDP) .GT. ALLOWED_XRAY_FLUX_ERROR)THEN
	        T1=SQRT(LUM*DESIRED_XRAY_LUM/OBS_XRAY_LUM_0P1)
	        IF(T1 .GT. 10.0_LDP)T1=10.0_LDP
	        IF(T1 .LT. 0.1_LDP)T1=0.1_LDP
	        FILL_FAC_XRAYS_1=T1*FILL_FAC_XRAYS_1
	        FILL_FAC_XRAYS_2=T1*FILL_FAC_XRAYS_2
	        IF(MYPE .EQ. 0)THEN
	          CALL UPDATE_KEYWORD(FILL_FAC_XRAYS_1,'[FIL_FAC_1]','VADAT',L_TRUE,L_FALSE,LUIN)
	          CALL UPDATE_KEYWORD(FILL_FAC_XRAYS_2,'[FIL_FAC_2]','VADAT',L_FALSE,L_TRUE,LUIN)
	          WRITE(6,*)'Adjusted filling factors to match desired X-ray luminosity'
	          WRITE(6,*)'Adjustment factor is',T1
	        END IF
	      END IF
	    END IF
	  END IF
	  IF(INCL_ADVECTION .AND. ADVEC_RELAX_PARAM .LT. 1.0_LDP .AND. MAXCH .LT. 100)THEN
	    COMPUTE_BA=.TRUE.
	    ADVEC_RELAX_PARAM=MIN(1.0_LDP,ADVEC_RELAX_PARAM*2.0_LDP)
	    IF(MYPE .EQ. 0)THEN
	      WRITE(LUER,*)'Have adjusted advection relaxation parameter to:',ADVEC_RELAX_PARAM
	      CALL UPDATE_KEYWORD(ADVEC_RELAX_PARAM,'[ADV_RELAX]','VADAT',L_TRUE,L_TRUE,LUIN)
	    END IF
	    MAXCH=100
	  END IF
!
! Adjust non-thermal decay energy scale factor. This is option is useful when adding non-thermal ioizations
! to a thermal model. During the convergence process we typically increase DEC_NRG_SCL_FAC by a factor of 10.
! If only 2 iterations were done, befores changing, we increase it by a factor of 100.
!
	  IF(TREAT_NON_THERMAL_ELECTRONS .AND. ADD_DEC_NRG_SLOWLY .AND. RD_LAMBDA .AND. MAXCH .LT. 100)THEN
	    IF(DEC_NRG_SCL_FAC .NE. 1.0_LDP)THEN
	      DEC_NRG_SCL_FAC=MIN(DEC_NRG_SCL_FAC*10.0_LDP,1.0_LDP)
	      IF(GAMMA_ADD_SLOWLY_COUNTER .EQ. 2)THEN
	        DEC_NRG_SCL_FAC=MIN(DEC_NRG_SCL_FAC*10.0_LDP,1.0_LDP)
	      END IF
	      IF(MYPE .EQ. 0)THEN
	        WRITE(LUER,*)'Have adjusted radioactivity decay energy scale factor'
	        WRITE(LUER,*)'New scale factor is',DEC_NRG_SCL_FAC
	        CALL UPDATE_KEYWORD(DEC_NRG_SCL_FAC,'[DECNRG_SCLFAC_BEG]','VADAT',L_TRUE,L_TRUE,LUIN)
	      END IF
	      MAXCH=100
	    ELSE
	      IF(MYPE .EQ. 0)THEN
	        CALL UPDATE_KEYWORD(L_FALSE,'[GAMMA_SLOW]','VADAT',L_TRUE,L_TRUE,LUIN)
	      END IF
	    END IF
	    GAMMA_ADD_SLOWLY_COUNTER=1
	  ELSE IF(DEC_NRG_SCL_FAC .NE. 1.0_LDP)THEN
	    GAMMA_ADD_SLOWLY_COUNTER=GAMMA_ADD_SLOWLY_COUNTER+1
	  END IF
!
	  IF (INC_SHOCK_POWER .AND. ADD_SHOCK_POWER_SLOWLY .AND.  RD_LAMBDA .AND. MAXCH .LT. 100)THEN
	    IF (SHOCK_POWER_FAC .NE. 1.0_LDP)THEN
	      SHOCK_POWER_FAC = MIN(SHOCK_POWER_FAC*10._LDP,1.0_LDP)
	      IF(MYPE .EQ. 0)THEN
	        WRITE(LUER,*)'Have adjusted SHOCK POWERscale factor'
	        WRITE(LUER,*)'New scale factor is',SHOCK_POWER_FAC
	        CALL UPDATE_KEYWORD(SHOCK_POWER_FAC,'[SHOCK_POWER_FAC_BEG]','VADAT',L_TRUE,L_TRUE,LUIN)
	      END IF
	      MAXCH=100
	    ELSE
	      IF(MYPE .EQ. 0)THEN
	        CALL UPDATE_KEYWORD(L_FALSE,'[ADD_SHOCK_POWER_SLOWLY]','VADAT',L_TRUE,L_TRUE,LUIN)
	      END IF
	    END IF
	  END IF
!
	  IF(SN_MODEL .AND. MYPE .EQ. 0)THEN
	    CALL WRITE_SEQ_TIME_FILE_MPI_V1(SN_AGE_DAYS,ND,LUSCR)
	  END IF
!
! If we have reached desired convergence, we do one final loop
! so as to write out all relevant model data.
!
	  IF(DO_T_AUTO .AND. RD_FIX_T .AND. MAXCH .LT. 50.0_LDP .AND.
	1                                 (LAST_LAMBDA .NE. MAIN_COUNTER) )THEN
	       RD_FIX_T=.FALSE.
	       FIXED_T=RD_FIX_T
               COMPUTE_BA=.TRUE.
	       IF(MYPE .EQ. 0)CALL UPDATE_KEYWORD(L_FALSE,'[FIX_T]','VADAT',L_TRUE,L_TRUE,LUIN)
	  ELSE IF( ( (RD_LAMBDA .AND. .NOT. DO_LAMBDA_AUTO) .OR. (LAST_LAMBDA .NE. MAIN_COUNTER)) .AND.
	1      MAXCH .LT. EPS .AND. NUM_ITS_TO_DO .NE. 0 .AND. .NOT. DONE_HYDRO_REVISION)THEN
	      IF(MAIN_COUNTER .NE. LAST_NG)NUM_ITS_TO_DO=1
	  ELSE
!
! If we are USING a fixed J, autmatically switched to variable J when convergence achieved.
! We switch when the convergence is 20%.
!
	    IF(USE_FIXED_J .AND. DO_LAMBDA_AUTO .AND. RD_LAMBDA .AND. MAXCH .LT. 50.0_LDP)THEN
	       USE_FIXED_J=.FALSE.
	       IF(MYPE .EQ. 0)CALL UPDATE_KEYWORD(L_FALSE,'[USE_FIXED_J]','VADAT',L_TRUE,L_TRUE,LUIN)
	       COMPUTE_EDDFAC=.TRUE.
	       IF(DO_GREY_T_AUTO)THEN
	         CALL GREY_T_ITERATE_MPI_V1(POPS,Z_POP,NU,NU_EVAL_CONT,FQW,
	1               LUER,LUIN,NC,ND,NP,NT,NCF,N_LINE_FREQ,MAX_SIM)
                 MAIN_COUNTER=MAIN_COUNTER+1
	         WRITE(6,*)'GREY_T_ITEARTE needs to be fixed for SET_NEXT_ITERATION'; STOP
	         IF(MYPE .EQ. 0THEN
	           CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
	1               LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	           CALL MPI_BCAST(IREC,IONE,MPI_INTEGER,IZERO,MPI_COMM_WORLD,IERR)
	         END IF
	       END IF
	       I=WORD_SIZE*(NDEXT+1)/UNIT_SIZE
!
! NB: If not ACCURATE, NDEXT was set to ND. The +1 arises since we write
! NU on the same line as RJ. J is used to get the REC_LENGTH, while string
! will contain the date.
!
	       CALL OPEN_RW_EDDFACTOR(R,V,LANG_COORD,ND,
	1             REXT,VEXT,LANG_COORDEXT,NDEXT,
	1             ACCESS_F,L_TRUE,COMPUTE_EDDFAC,L_FALSE,'EDDFACTOR',LU_EDD)
	       CALL INIT_GET_J_FOR_TWO_PHOT
	       IF(MYPE .EQ. 0)THEN
	         T1=0.0D0; WRITE(LU_EDD,REC=FINISH_REC)T1
	       END IF
	       COHERENT_ES=.TRUE.
!
! If DONE_HYDRO_REVISION is TRUE, we do a LAMBDA iteration immediately afterwoods.
! LAMBDA_ITERATION has already been set.
!
	    ELSE IF(DONE_HYDRO_REVISION .OR. R_GRID_REVISED)THEN
	       LAMBDA_ITERATION=.TRUE.
!
! If RD_LABDA is TRUE., switch to full iteration when convergence has been achieved.
! We switch when the convergence is 20%.
!
	    ELSE IF(DO_LAMBDA_AUTO .AND. RD_LAMBDA .AND. MAXCH .LT. 50.0_LDP)THEN
	       IF(MYPE .EQ. 0)CALL UPDATE_KEYWORD(L_FALSE,'[DO_LAM_IT]','IN_ITS',L_TRUE,L_TRUE,LUIN)
	       RD_LAMBDA=.FALSE.
	       LAMBDA_ITERATION=.FALSE.
	       FIXED_T=RD_FIX_T
	    END IF
!
	    CALL SPECIFY_IT_CYCLE_V3(MAIN_COUNTER,COMPUTE_BA,LAMBDA_ITERATION,FIXED_T,NEXT_NG,NEXT_AV)
	    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	    CALL MPI_BCAST(NEWMOD,IONE,MPI_LOGICAL,IZERO,MPI_COMM_WORLD,IERR)
!	    CALL MPI_BCAST(RD_LAMBDA,IONE,MPI_LOGICAL,IZERO,MPI_COMM_WORLD,IERR)
!	    CALL MPI_BCAST(LAMBDA_ITERATION,IONE,MPI_LOGICAL,IZERO,MPI_COMM_WORLD,IERR)
!	    CALL MPI_BCAST(FIXED_T,IONE,MPI_LOGICAL,IZERO,MPI_COMM_WORLD,IERR)
!
! Check to see if the user has changed IN_ITS to modify the number of iterations being
! undertaken. If the file has not been modified, no action will be taken. The use may
! also change whether LAMBDA iterations are being done. At least one final iteration
! will be undertaken.
!
	      WRITE(STRING,'(I3.3)')MYPE; STRING='MODEL_SCR_'//STRING
	      CALL GEN_ASCI_OPEN(LUSCR,STRING,'UNKNOWN',' ',' ',IZERO,IOS)
	      IF(IOS .EQ. 0)CALL GEN_ASCI_OPEN(LUIN,'IN_ITS','OLD',' ','READ',IZERO,IOS)
	      IF(IOS .NE. 0)THEN
	         WRITE(LUER,*)'Error opening IN_ITS or '//TRIM(STRING)//'in CMFGEN, IOS=',IOS
	         WRITE(LUER,*)'Error occurs at the end of CMFGEN_SUB.'
	         WRITE(LUER,*)'Error will be ignored.'
	         GOTO 20000
	      END IF
	      CALL RD_OPTIONS_INTO_STORE(LUIN,LUSCR)
	      OLD_RD_LAMBDA=RD_LAMBDA
	      I=NUM_ITS_RD
	      CALL RD_STORE_INT(NUM_ITS_RD,'NUM_ITS',L_TRUE,'Number of iterations to perform')
	      CALL RD_STORE_LOG(RD_LAMBDA,'DO_LAM_IT',L_TRUE,'Do LAMBDA iterations ?')
	      CALL RD_STORE_LOG(DO_LAMBDA_AUTO,'DO_LAM_AUTO',L_FALSE,
	1                  'Start non-lambda iterations automatically?')
	      CALL RD_STORE_LOG(DO_GREY_T_AUTO,'DO_GT_AUTO',L_FALSE,
	1                  'Do a grey temperature iteration after revising USE_FIXED_J?')
	      CALL RD_STORE_LOG(DO_T_AUTO,'DO_T_AUTO',L_FALSE,
	1                  'Allow temperature to vary when sufficent convergence has been obtained?')
	      CALL RD_STORE_LOG(SET_POPS_D2_EQ_D1,'D2_EQ_D1',L_FALSE,
	1                  'Replace pops at depth 2 with those at depth 1 for non-LAMBDA it?')
	      CALL CLEAN_RD_STORE()
	    CLOSE(UNIT=LUIN)
	    NUM_ITS_TO_DO=NUM_ITS_TO_DO+(NUM_ITS_RD-I)
	    IF(NUM_ITS_TO_DO .LE. 0)NUM_ITS_TO_DO=1
	    IF(RD_LAMBDA)THEN
	      LAMBDA_ITERATION=.TRUE.
	      FIX_IMPURITY=.FALSE.
	      FIXED_T=.TRUE.
!
! Don't wish to change ITERATION cycle values unless we have to.
!
	    ELSE IF(OLD_RD_LAMBDA)THEN
	      FIX_IMPURITY=RD_FIX_IMP
	      FIXED_T=RD_FIX_T
	    END IF
	    CLOSE(UNIT=LUSCR,STATUS='DELETE')
	  END IF
	  CALL TUNE(ITWO,'GIT')
	  CALL TUNE(ITHREE,' ')
	  IF(MYPE .EQ. 0)WRITE(6,*)'End of iteration section',MYPE,MAXCH,FIXED_T
	  FLUSH(UNIT=6)
!
	  GOTO 20000				!Begin another iteration
