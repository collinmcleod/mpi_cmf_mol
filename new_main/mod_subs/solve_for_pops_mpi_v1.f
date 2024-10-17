!
! Routine designed to solve the linearized statistical and radiative equilibrium
! equations for the new population estimates. Various tricks (e.g., NG Acceleration)
! are done to help ensure convergence. Routine was cut from CMFGEN_SUB in order
! to facilitate latter changes.
!
	SUBROUTINE SOLVE_FOR_POPS_MPI_V1(POPS,NT,NION,ND,NC,NP,NUM_BNDS,DIAG_INDX,
	1            MAXCH,MAIN_COUNTER,IREC,LU_SE,LUSCR,LST_ITERATION)
	USE SET_KIND_MODULE
	USE MPI
	USE MOD_CMFGEN
        USE ANG_QW_MOD
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
! Altered:   26-Apr-2019 : Added "CALL SMOOTH_POPS_AS_WE_ITERATE" (done earlier on OSIRIS).
! Altered:   18-Nov-2016 : Before doing AV, check NG acceleration was not done recently.
! Altered:   01-Apr-2015 : Changed ESTAU call to ESTAU_V2 to take clumping to account.
! Altered:   21-Jan-2014 : Changed to CALL SOLVEBA_V12 - now pass MAIN_COUNTER.
! Aleterd:   31-Dec-2013 : No longer set COMPUTE_BA=T when T_MIN_BA_EXTRAP is true. We may need to
!                            change this.
! Altered:   26-Mar-2012 : Changed to CALL SOLVEBA_V11. Has POP_ATOM.
! Altered:   12-Mar-2012 : Changed to CALL SOLVEBA_V10. Allows populations at depth 2 to be
!                            equated to those at depth 1.
! Altered:   18-May-2010 : Change to allow BA to be held fixed after a LAMBDA iteration.
! Altered:   23-Feb-2007 : Call to SOLVEBA_V8 changed to SOLVEBA_V9; LAM_SCALE_OPT inserted
!                            into SOLVEBA_V9 call.
! Finalized: 15-Feb-2006
!
	INTEGER NT
	INTEGER NION
	INTEGER ND
	INTEGER NC
	INTEGER NP
	INTEGER NUM_BNDS
	INTEGER DIAG_INDX
!
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) MAXCH
!
	INTEGER MAIN_COUNTER
	INTEGER IREC
	INTEGER LU_SE
	INTEGER LUSCR
	LOGICAL LST_ITERATION
!
! Local data
!
	REAL(KIND=LDP), ALLOCATABLE :: SOL(:,:)
	REAL(KIND=LDP) LOC_SOL(NT,DST:DEND)
	REAL(KIND=LDP) R_OLD(ND)
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) ESEC(ND)
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP), SAVE :: MAXCH_SUM=0.0_LDP
	INTEGER, SAVE :: COUNT_FULL_LAM_BA=1
!
	INTEGER I,J
	LOGICAL FILE_OPEN
	LOGICAL BA_SUCCESSFULLY_OUTPUT
	LOGICAL SUCCESS
	CHARACTER*20 TEMP_CHAR
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!	include 'mpif.h'
!
	LUER=ERROR_LU()	
!
! We can adjust the BA/STEQ equations so that certain species are
! held fixed. For all species, this now done using FIXPOP_IN_BA_V2
! called by GENERATE_FULL_MATRIX.
!
	MOD_FIXED_NE=FIXED_NE
	MOD_FIX_IMPURITY=FIX_IMPURITY
	MOD_FIXED_T=.FALSE.
!
! Determine those depths where the temperature is to be held fixed.
! If either FIXED_T of RD_FIX_T is true, the temperature is fixed independent
! of VARFIXT. We zero all elements of the R.E. Eq. except the local variation
! with respect to T. For a LAMBDA iteration this element is zero, and hence
! must be set.
!
	MOD_TAU_SCL_T=TAU_SCL_T
	MOD_FIX_T_D_ST=0
	MOD_FIX_T_D_END=0
	IF(FIXED_T .OR. RD_FIX_T)THEN
	  MOD_FIXED_T=.TRUE.
	  MOD_FIX_T_D_ST=1
	  MOD_FIX_T_D_END=ND
	  CALL ESTAU_V2(TA,R,ED,CLUMP_FAC,TB,ND)
	ELSE IF(VARFIXT .AND. CON_SCL_T .NE. 0.0_LDP)THEN
	  CALL ESTAU_V2(TA,R,ED,CLUMP_FAC,TB,ND)
	  MOD_FIXED_T=.TRUE.
	  MOD_FIX_T_D_ST=1
	  DO I=1,ND
	    IF(TA(I) .LT. TAU_SCL_T)THEN
	      MOD_FIX_T_D_END=I
	    END IF
	  END DO
	  IF(MOD_FIX_T_D_END .LT. 5)MOD_FIX_T_D_END=5  	!avoid inconsistency at outer boundary.
	  WRITE(LUER,*)' '
	  WRITE(LUER,'(A,ES14.4)')' Electron scattering optical depth at outer boundary'//
	1                        ' in SOLVE_FOR_POPS is',TA(1)
	  WRITE(LUER,'(A,I3,A,I3)')' T will be held fixed over the following depths',
	1                           MOD_FIX_T_D_ST,' to ',MOD_FIX_T_D_END
	  WRITE(LUER,*)' '
	END IF
! 
!
	IF(LAMBDA_ITERATION .AND. MYPE .EQ. 0)THEN
	  WRITE(LUER,*)'LAMBDA iteration used'
	END IF
!
	IF(COMPUTE_BA .AND. MYPE .EQ. 0)THEN
	  WRITE(LUER,*)'BA matrix computed'
	ELSE IF(MYPE .EQ. 0)THEN
	  WRITE(LUER,*)'BA matrix **NOT** computed'
	END IF
	FLUSH(LUER)
!
	NLBEGIN=0	!Initialize for next iteration
	CONTINUE
	CALL WR_ASCI_STEQ_MPI_V1(NION,DST,DEND,ND,'STEQ ARRAY',LU_SE)
!
! If we are currently doing a LAMBDA iteration we allow for bigger
! changes - up to a factor of 100. Otherwise changes are limited to
! a factor of 10.0D0. (Installed 28-Dec-1989). Specified by parameter
! T1 in SOLVEBA call. TEMP_CHAR is used to utilize the BLK_DIAGONAL
! nature of the BA matrix when performing a LAMBDA iteration.
!
	IF(LAMBDA_ITERATION)THEN
	  T1=MAX_LAM_COR
	  TEMP_CHAR='DIAG'
	ELSE
	  T1=MAX_LIN_COR
	  TEMP_CHAR=METH_SOL
	END IF
	CALL SOLVEBA_MPI_V1(LOC_SOL,POPS,POP_ATOM,DIAG_INDX,
	1       NT,NION,NUM_BNDS,DST,DEND,ND,
	1       MAXCH,TEMP_CHAR,SUCCESS,SCALE_OPT,
	1       LAM_SCALE_OPT,T1,MAX_dT_COR,T_MIN,
	1       COMPUTE_BA,WR_BA_INV,WR_PART_OF_INV,LAMBDA_ITERATION,
	1       MAIN_COUNTER,SET_POPS_D2_EQ_D1)
!
	I=NT*ND
	IF(MYPE .EQ. 0)ALLOCATE(SOL(NT,ND))
	CALL GATHER_POPS_MPI_V1(LOC_SOL,SOL,NT,DST,DEND,ND)
	LOC_SOL(:,DST:DEND)=POPS(:,DST:DEND)
	CALL GATHER_POPS_MPI_V1(LOC_SOL,POPS,NT,DST,DEND,ND)
!
! Complicated algorithim to decide when to switch off BA computation.
! We only switch off BA computation in WRBAMAT_RDIN is TRUE.
!
! 1. We can switch off the BA computation, once we are no longer doing
!       LAMBDA iteartions. This is only done for N_ITS_TO_FIX_BA iterations.
!
! 2. We switch of the BA computation semi-permanently once the corrections
!       are less than VAL_FIX_BA.
!
! 3. We recompute BA if the accumulated sum of the maximum corrections
!       exceeds 3*VAL_FIX_BA.
!
	IF(N_ITS_TO_FIX_BA .GT. 0 .AND. .NOT. LAMBDA_ITERATION)THEN
	  IF(COMPUTE_BA)THEN
	    CNT_FIX_BA=0
	    MAXCH_SUM=0.0_LDP
	  END IF
	  MAXCH_SUM=MAXCH_SUM+MAXCH
	  CNT_FIX_BA=CNT_FIX_BA+1
	  IF(CNT_FIX_BA .GT. N_ITS_TO_FIX_BA .AND. MAXCH .GT. VAL_FIX_BA)THEN
	    COMPUTE_BA=L_TRUE
	    WRBAMAT=WRBAMAT_RDIN
	  ELSE
	    IF(WRBAMAT_RDIN)COMPUTE_BA=L_FALSE
	  END IF
	  IF(CNT_FIX_BA .GT. N_ITS_TO_FIX_BA .AND. MAXCH_SUM .GT. 3.0_LDP*VAL_FIX_BA)THEN
	    MAXCH_SUM=0.0_LDP
	    COMPUTE_BA=L_TRUE
	    WRBAMAT=WRBAMAT_RDIN
	  END IF
!
! This will force BA to be recomuted again, after N_ITS_TO_FIX_BA, because T was variable
! at some depths.
!
	  IF(CON_SCL_T .NE. 0.0_LDP)THEN
	    MAXCH_SUM=1000.0_LDP*VAL_FIX_BA
	  END IF
	ELSE
!
! Switch off BA computation of MAXCH if less than VAL_FIX_BA.
! We only do this provided the last iteration was not a LAMBDA
! iteration, and if T was not partially held fixed at some depths.
!
! We now only right out the BA matrix if we are nearing the correct solution,
! and assuming that it is a full solution matrix (i.e. not from a
! LAMBDA iteration and CON_SCL_T was not set).
!
	  COMPUTE_BA=L_TRUE
	  IF( MAXCH .LT. VAL_FIX_BA .AND. WRBAMAT
	1          .AND. .NOT. LAMBDA_ITERATION
	1          .AND. CON_SCL_T .EQ. 0.0_LDP)THEN
	    COMPUTE_BA=COMPUTE_BARDIN
	  END IF
	  IF(WRBAMAT_RDIN .AND. MAXCH .LT. 2.0_LDP*VAL_FIX_BA
	1         .AND. .NOT. LAMBDA_ITERATION
	1          .AND. CON_SCL_T .EQ. 0.0_LDP)THEN
	    WRBAMAT=.TRUE.
	  END IF
	END IF
!
	IF(MYPE .EQ. 0)THEN
          CALL WR2D_V2(SOL,NT,ND,'STEQ SOLUTION ARRAY','#',L_TRUE,LU_SE)
	END IF
!
	IF(AUTO_SMOOTH_POPS)THEN
	  CALL SMOOTH_POPS_AS_WE_ITERATE(POPS,SOL,ND,NT)
	END IF
!
! NB: The call to SUM_STEQ_SOL corrupts SOL. I is used for output --
!      the unit is closed on exit.
!
	IF(MYPE .EQ. 0)THEN
	  I=7; CALL SUM_STEQ_SOL(SOL,NT,ND,I)
	END IF
!
! Determine whether convergence is sufficient to consider using
! NG acceleration. The first NG acceleration is done 4 iterations after
! the iteration on which MAXCH < VAL_DO_NG. IST_PER_NG must be greater
! than, or equal to, 4. LAST_NG now always indicates the last iteration
! on which an NG acceleration was performed. NEXT_NG is used to indicate
! the next iteration on which an NG acceleration is to occur.
!
! NEXT_NG=1000 indicates a NEW model
! NEXT_NG=1500 indicates that MAXCH is again above VAL_DO_NG
! NEXT_NG=2000 indicates a CONTINUING model.
!
! A NG acceleration can be forced after 1 iteration if LAST_NG is set to
! some vale .LE. MAIN_COUNTER-ITE_PER_NG in the POINT1 file, and provide the
! change on that iteration is less than VAL_DO_NG.
!
	IF(MAIN_COUNTER .LT. IT_TO_BEG_NG-3)THEN
	  NEXT_NG=1000
	ELSE IF(MAIN_COUNTER .GE. IT_TO_BEG_NG-3 .AND. NEXT_NG .EQ. 1000)THEN
	  IF( MAXCH .LE. VAL_DO_NG )NEXT_NG=MAX(LAST_NG+ITS_PER_NG,MAIN_COUNTER+4)
	ELSE IF(MAXCH .LT. VAL_DO_NG .AND. NEXT_NG .EQ. 1500)THEN
	  NEXT_NG=MAX(LAST_NG+ITS_PER_NG,MAIN_COUNTER+4)
	ELSE IF(MAXCH .LT. VAL_DO_NG .AND. NEXT_NG .EQ. 2000)THEN
	  NEXT_NG=MAX(LAST_NG+ITS_PER_NG,MAIN_COUNTER+1)
	  IF(LAST_NG .EQ. -1000)NEXT_NG=MAIN_COUNTER+4
	ELSE IF( MAXCH .GE. VAL_DO_NG )THEN
	  NEXT_NG=1500
	END IF
!
! We switch between LINEARIZATION and LAMBDA iterations until the
! maximum percentage change is less  than VAL_DO_LAM. RD_CNT_LAM iterations
! are performed per full linearization. RD_LAMBDA overrides this section if
! TRUE. We only do a LAMBDA iteration with T fixed. NB --- FIX_IMPURITY
! is a soft option --- is removed when convergence nearly obtained.
!
	IF(LAMBDA_ITERATION)LAST_LAMBDA=MAIN_COUNTER
	IF(.NOT. RD_LAMBDA)THEN
	  IF( MAXCH .LT. VAL_DO_LAM )THEN
	     FIX_IMPURITY=RD_FIX_IMP .AND. LAMBDA_ITERATION
	     LAMBDA_ITERATION=.FALSE.
	     FIXED_T=RD_FIX_T
	     CNT_LAM=0
	  ELSE IF(LAMBDA_ITERATION)THEN
	     CNT_LAM=CNT_LAM+1
	     IF(CNT_LAM .GE. RD_CNT_LAM)THEN
	       LAMBDA_ITERATION=.FALSE.
	       FIX_IMPURITY=RD_FIX_IMP
	       FIXED_T=RD_FIX_T
	     END IF
	     IF(MAXCH .GT. 1.0E+05_LDP)THEN
	       LAMBDA_ITERATION=.TRUE.
	       COMPUTE_BA=.TRUE.
	       FIX_IMPURITY=.FALSE.
	       FIXED_T=.TRUE.
	     END IF
	  ELSE
	     LAMBDA_ITERATION=.TRUE.
	     COMPUTE_BA=.TRUE.
	     CNT_LAM=0
	     FIXED_T=.TRUE.
	     FIX_IMPURITY=.FALSE.
	  END IF
	END IF
!
! Automatically adjust R grid, so that grid is uniformally spaced on the
! FLUX optical depth scale. Used for SN models with very sharp ioinization
! fronts. By doing it before the output to SCRTEMP, we ensure that
! a continuuing model starts with the revised R grid.
!
! Write pointer file and output data necessary to begin a new
! iteration.
!
	IF(MYPE .EQ. 0)THEN
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
	1                LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	END IF
!
! Program to create a NEW_R_GRID which is equally spaced in LOG(Tau) where
! TAU is based on the FLUX mean opacity.
!
! In VADAT REVISE_R_GRID should be set to TRUE.
!
	IF(REVISE_R_GRID)THEN
	  R_OLD(1:ND)=R(1:ND)
	  CALL ESOPAC(ESEC,ED,ND)               !Electron scattering emission factor.
          CALL ADJUST_R_GRID_V4(POPS,ESEC,MAIN_COUNTER,R_GRID_REVISED,ND,NT)
	  IF(R_GRID_REVISED)THEN
	     MAIN_COUNTER=MAIN_COUNTER+1
	     IF(MYPE .EQ. 0)CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
	1                          LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  END IF
	ELSE
	  R_GRID_REVISED=.FALSE.
	END IF
! 
!
! Perform the acceleration. T1 and T2 return the percentage changes
! in the populations. I is used as an INTEGER error flag. We never do
! an NG acceleration on the last iteration.
!
! We perform a LAMBDA iteration after the NG acceleration in the
! case the where the NG acceleration has caused population changes
! > 20%. In this case, the BA matrix will also be recaluated on
! the next full iteration.
!
!	WRITE(6,*)MAIN_COUNTER,NEXT_AV,LAST_LAMBDA+NUM_OSC_AV+2
!	WRITE(6,*)NEXT_NG-4,LAST_AV+ITS_PER_AV
!	WRITE(6,*)AVERAGE_DO,LST_ITERATION
!
	AVERAGE_DONE=.FALSE.
	IF(MAIN_COUNTER .NE. LAST_LAMBDA .AND. UNDO_LAST_IT .AND. MYPE .EQ. 0)THEN
	  I=1; J=5
	  CALL  UNDO_IT(SOL,POPS,NT,ND,MAIN_COUNTER,I,J,AVERAGE_DONE)
	  IF(AVERAGE_DONE)THEN
	    WRITE(LUER,'(A,I4,A,I4)')'Last corrections undone for depths',I,' to',J
	  ELSE
	    WRITE(LUER,*)'Error undoing last corrections.'
	  END IF
	END IF
!
	AVERAGE_DONE=.FALSE.
	IF(AVERAGE_DO .AND. MAIN_COUNTER .GE. NEXT_AV
	1             .AND. MAIN_COUNTER .GT. LAST_NG+NUM_OSC_AV
	1             .AND. MAIN_COUNTER .GT. LAST_LAMBDA+NUM_OSC_AV+2
	1             .AND. (MAIN_COUNTER .LT. NEXT_NG-4 .OR. .NOT. NG_DO)
	1             .AND. MAIN_COUNTER .GE. LAST_AV+ITS_PER_AV
	1             .AND. MYPE .EQ. 0
	1             .AND. .NOT. LST_ITERATION)THEN
	  CALL AVE_FLIPS(SOL,POPS,NT,ND,MAIN_COUNTER,NUM_OSC_AV,L_TRUE,AVERAGE_DONE)
	  IF(AVERAGE_DONE)THEN
	    WRITE(LUER,*)'Averaging of oscilating populations',
	1                   ' performed on iteration',MAIN_COUNTER
	    LAST_AV=MAIN_COUNTER
	    NEXT_AV=MAIN_COUNTER+ITS_PER_AV
	  ELSE
	    NEXT_AV=MAIN_COUNTER+8
	    WRITE(LUER,*)'Error performing AV acceleration.'
	    WRITE(LUER,*)'Will try again in 8 iterations.'
	  END IF
	END IF
!
! Save the corrected poplatons to SCREMP.
!
	IF(AVERAGE_DONE .AND. MYPE .EQ. 0)THEN
	  MAIN_COUNTER=MAIN_COUNTER+1
	  POPS=SOL
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,
	1             RITE_N_TIMES,LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	END IF
!
	IF(NG_DO .AND. (.NOT. LST_ITERATION) .AND.
	1       MYPE .EQ. 0 .AND. 
	1       (NEXT_NG .EQ. MAIN_COUNTER) )THEN
	  CALL DO_NG_BAND_ACCEL_V4(POPS,R,V,SIGMA,R_OLD,
	1         NT,ND,NG_BAND_WIDTH,NG_DONE,T1,T2,
	1         DO_NG_VALIDITY_CHECK,LUSCR,LUER)
	  IF(NG_DONE)THEN
	    MAIN_COUNTER=MAIN_COUNTER+1
	    LAST_NG=MAIN_COUNTER
	    NEXT_NG=MAIN_COUNTER+ITS_PER_NG
	    CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,
	1             RITE_N_TIMES,LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	    IF(T1 .GT. 50.0_LDP .AND. T2 .GT. 50.0_LDP)THEN
	      LAMBDA_ITERATION=.TRUE.
	      CNT_LAM=0
	      FIXED_T=.TRUE.
	      FIX_IMPURITY=.FALSE.
	      COMPUTE_BA=.TRUE.
	    END IF
	  ELSE
	    NEXT_NG=MAIN_COUNTER+8
	    WRITE(LUER,*)'Error performing NG acceleration.'
	    WRITE(LUER,*)'Will try again in 8 iterations.'
	    WRITE(LUER,*)'Error flag=',I
	  END IF
	END IF
	FLUSH(LUER)
!
! Here we check whether we keep BA fixed. Previously we always computed the BA
! matrix after a LAMBDA iteration. This will allow us to keep it fixed for
! at least 1 iteration. Since the computation of the BA matrix takes considerable time,
! this may reduce computational effort.
!
	IF(R_GRID_REVISED)THEN
	   LAMBDA_ITERATION=.TRUE.
	   COMPUTE_BA=.TRUE.
	   COUNT_FULL_LAM_BA=1000
	ELSE IF(LAST_LAMBDA .EQ. MAIN_COUNTER .AND. .NOT. LAMBDA_ITERATION)THEN
	   IF(COUNT_FULL_LAM_BA .NE. 0)THEN
	      COUNT_FULL_LAM_BA=0
	   ELSE
	     OPEN(UNIT=7, FILE='BAMATPNT',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	     IF(IOS .EQ. 0)READ(7,*,IOSTAT=IOS)BA_SUCCESSFULLY_OUTPUT
	     INQUIRE(UNIT=7,OPENED=FILE_OPEN)
	     IF(FILE_OPEN)CLOSE(UNIT=7)
	     IF(IOS .EQ. 0 .AND. BA_SUCCESSFULLY_OUTPUT)THEN
	       COUNT_FULL_LAM_BA=COUNT_FULL_LAM_BA+1
	     END IF
	  END IF
	END IF
!
! Ensure consistency in parameters across all processes.
!
	I=ND*NT
	CALL MPI_BCAST(POPS,I,MPI_DOUBLE_PRECISION,IZERO,MPI_COMM_WORLD,IERR)
!
	CALL MPI_BCAST(LAST_NG,IONE,MPI_INTEGER,IZERO,MPI_COMM_WORLD,IERR)
	CALL MPI_BCAST(NEXT_NG,IONE,MPI_INTEGER,IZERO,MPI_COMM_WORLD,IERR)
	CALL MPI_BCAST(MAIN_COUNTER,IONE,MPI_INTEGER,IZERO,MPI_COMM_WORLD,IERR)
	CALL MPI_BCAST(CNT_LAM,IONE,MPI_INTEGER,IZERO,MPI_COMM_WORLD,IERR)
	CALL MPI_BCAST(COUNT_FULL_LAM_BA,IONE,MPI_INTEGER,IZERO,MPI_COMM_WORLD,IERR)
!
	CALL MPI_BCAST(FIXED_T,IONE,MPI_LOGICAL,IZERO,MPI_COMM_WORLD,IERR)
	CALL MPI_BCAST(FIX_IMPURITY,IONE,MPI_LOGICAL,IZERO,MPI_COMM_WORLD,IERR)
	CALL MPI_BCAST(LAMBDA_ITERATION,IONE,MPI_LOGICAL,IZERO,MPI_COMM_WORLD,IERR)
	CALL MPI_BCAST(COMPUTE_BA,IONE,MPI_LOGICAL,IZERO,MPI_COMM_WORLD,IERR)
!
	IF(MYPE .EQ. 0)DEALLOCATE(SOL)
	RETURN
	END
