	SUBROUTINE SOLVEBA_MPI_V1(SOL_MAT,POPS,POP_ATOM,DIAG_INDX,
	1              NT,NION,NUM_BNDS,DST,DEND,ND,
	1              MAXCH,METH_SOL,SUCCESS,SCALE_OPT,
	1              LAM_SCALE_OPT,CHANGE_LIM,MAX_dT_COR,T_MIN,
	1              BA_COMPUTED,WR_BA_INV,WR_PRT_INV,LAMBDA_IT,
	1              MAIN_COUNTER,SET_POPS_D2_EQ_D1)
	USE SET_KIND_MODULE
	USE MPI
	USE CONTROL_VARIABLE_MOD, ONLY: LTE_MODEL, TRI_SOL_OPTIONS
	IMPLICIT NONE
!
! Altered 08-Aug-2025:  POPS at depth 2 set to those at depth 1 even if LAMBDA operation
!                           (assuming SET_POPS_D2_EQ_D1 is set).
!
	INTEGER NT
	INTEGER NION
	INTEGER NUM_BNDS
	INTEGER DST,DEND,ND
	INTEGER DIAG_INDX
	INTEGER MAIN_COUNTER
!
	REAL(KIND=LDP) SOL_MAT(NT,DST:DEND)
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) POP_ATOM(ND)
	REAL(KIND=LDP) MAXCH
	REAL(KIND=LDP) CHANGE_LIM
	REAL(KIND=LDP) MAX_dT_COR
	REAL(KIND=LDP) T_MIN
	CHARACTER*(*) METH_SOL
	CHARACTER*(*) SCALE_OPT
	CHARACTER*(*) LAM_SCALE_OPT
	LOGICAL SUCCESS
	LOGICAL BA_COMPUTED
	LOGICAL WR_BA_INV
	LOGICAL WR_PRT_INV
	LOGICAL LAMBDA_IT
	LOGICAL SET_POPS_D2_EQ_D1
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
!
!
! Local variables.
!
	INTEGER, PARAMETER :: NV=10
	REAL(KIND=LDP) MAX_INC_VEC(NV)
	REAL(KIND=LDP) MAX_DEC_VEC(NV)
	REAL(KIND=LDP) ACTUAL_MAX_dT_COR
!
	REAL(KIND=LDP) SCALE,MINSCALE,T1,T2,T3
	REAL(KIND=LDP) INCREASE,DECREASE
	REAL(KIND=LDP) INCREASE_SAVE,DECREASE_SAVE
	REAL(KIND=LDP) BIG_LIM,LIT_LIM
	INTEGER I,J,K,IINC,IDEC,IOS
	INTEGER COUNT(7)
	INTEGER IERR
	LOGICAL LOC_WR_BA_INV
	LOGICAL DO_LEVEL_CHK
!
!	include 'mpif.h'
!
	LUER=ERROR_LU()
	CALL SET_CASE_UP(SCALE_OPT,IONE,IZERO)
	IF(SCALE_OPT(1:5) .NE. 'LOCAL' .AND. SCALE_OPT(1:4) .NE. 'NONE'
	1   .AND. SCALE_OPT(1:6) .NE. 'GLOBAL'
	1   .AND. SCALE_OPT(1:5) .NE. 'MAJOR')THEN
	  SCALE_OPT='MAJOR'
	  WRITE(LUER,*)'Warning - Invalid scale option in SOLVEBA',
	1          ' MAJOR scaling assumed'
	END IF
!
! Solve for the perturbations. ND: The scaling of BA is now done within CMF_XXX_BAND routines.
!
	IF(METH_SOL(1:4) .EQ. 'DIAG')THEN
!
	    CALL TUNE(IONE,'DIAG_BAND')
	    LOC_WR_BA_INV=WR_BA_INV
	    IF(LAMBDA_IT)LOC_WR_BA_INV=.FALSE.
	    CALL CMF_DIAG_BAND_MPI_V1(SOL_MAT,POPS,METH_SOL,SUCCESS,
	1              DIAG_INDX,NT,NION,NUM_BNDS,DST,DEND,ND,
	1              BA_COMPUTED,LOC_WR_BA_INV,WR_PRT_INV)
	    IF(.NOT. SUCCESS)THEN
	      WRITE(LUER,*)'Error in CMF_DIAG_BAND_MPI_V1 - shutting code down'
	      STOP
	    END IF
	    CALL TUNE(ITWO,'DIAG_BAND')
!
	ELSE IF( METH_SOL(1:3) .EQ. 'TRI')THEN
!
	    CALL TUNE(IONE,'TRI_BAND')
	    LOC_WR_BA_INV=WR_BA_INV
	    IF(INDEX(TRI_SOL_OPTIONS,'THOMAS') .NE. 0)THEN
	      CALL CMF_TRIBAND_THOMAS_MPI_V1(SOL_MAT,POPS,METH_SOL,SUCCESS,
	1                DIAG_INDX,NT,NION,NUM_BNDS,DST,DEND,ND,
	1                BA_COMPUTED,LOC_WR_BA_INV,WR_PRT_INV,TRI_SOL_OPTIONS)
	    ELSE
!
! We use an iterative technique.
!
	      CALL CMF_TRI_BAND_MPI_V2(SOL_MAT,POPS,METH_SOL,SUCCESS,
	1               DIAG_INDX,NT,NION,NUM_BNDS,DST,DEND,ND,
	1               BA_COMPUTED,LOC_WR_BA_INV,WR_PRT_INV,TRI_SOL_OPTIONS)
	      IF(.NOT. SUCCESS)THEN
	        WRITE(LUER,*)'Error in CMF_TRI_BAND_MPI_V1 - shutting code down'
	        STOP
	      END IF
	    END IF
	    CALL TUNE(ITWO,'TRI_BAND')
!
	ELSE
	  WRITE(LUER,*)'Error - invalid solution method in SOLVEBA'
	  WRITE(LUER,*)'Solution method is ',METH_SOL
	  STOP
	END IF
!
! This is done to prevent corrections at depth 2 setting the nature of future iterations.
!
	IF(SET_POPS_D2_EQ_D1)THEN                               ! .AND. .NOT. LAMBDA_IT)THEN
	  IF(DST .LE. 2 .AND. DEND .GE. 2)SOL_MAT(:,2)=0.0_LDP
	END IF
!
! 
!
! Lambda iterations (for Ne fixed) should generally have corrections < unity.
! Due to instabilities, large -ve corrections can sometimes arrise. This
! limits these corrections, and potentilly wil help facilitate convergence.
!
	IF(LAMBDA_IT .AND. LAM_SCALE_OPT(1:5) .EQ. 'LIMIT')THEN
	  COUNT(1)=0
	  DO I=DST,DEND
	    DO J=1,NT-1
	      IF(SOL_MAT(J,I) .GT. 1.1_LDP)THEN
	        SOL_MAT(J,I)=0.999_LDP
	        COUNT(1)=COUNT(1)+1
	      END IF
	    END DO
	  END DO
	  IF(COUNT(1) .NE. 0)WRITE(LUER,*)'Warning -- using  LIMIT option for LAMBDA iteration in SOLVEBA_MPI_V1.f'
	END IF
!
! Determine maximum corrections to the 'population parameters', and output summary file:
!
	CALL CREATE_CORRECTION_SUM_MPI_V1(SOL_MAT,DECREASE,INCREASE,ACTUAL_MAX_dT_COR,IDEC,IINC,DST,DEND,ND,NT)
!
	INCREASE_SAVE=INCREASE; DECREASE_SAVE=DECREASE
	DECREASE=100.0_LDP*DECREASE
	INCREASE=100.0_LDP*INCREASE
	IF(MYPE .EQ. 0)THEN
	  IF(LAMBDA_IT)THEN
 	    WRITE(LUER,9000)IINC,ABS(INCREASE),'  (LAMBDA)','          --- iteration ',MAIN_COUNTER
	    WRITE(LUER,9200)IDEC,DECREASE,     '  (LAMBDA)','          --- iteration ',MAIN_COUNTER
	  ELSE IF(BA_COMPUTED)THEN
 	    WRITE(LUER,9000)IINC,ABS(INCREASE),'  (BA computed)','     --- iteration ',MAIN_COUNTER
	    WRITE(LUER,9200)IDEC,DECREASE,     '  (BA computed)','     --- iteration ',MAIN_COUNTER
	  ELSE
 	    WRITE(LUER,9000)IINC,ABS(INCREASE),'  (BA not computed)',' --- iteration ',MAIN_COUNTER
	    WRITE(LUER,9200)IDEC,DECREASE,     '  (BA not computed)',' --- iteration ',MAIN_COUNTER
	  END IF
9000	  FORMAT(' Maximum % increase at depth ',I4,' is',ES10.2,A,A,I4)
9200	  FORMAT(' Maximum % decrease at depth ',I4,' is',ES10.2,A,A,I4)
	END IF
!
! Convert decrease into more useful form. If DECREASE=100, the
! new population value is zero and change is infinite. i.e we
! convert decrease to form 100(old-new)/new. Currently in the form
! 100(old-new)/old [same as INCREASE]. We allow for big decreases -
! limit MAXCH to 1.0D+07 - This value does not halt program execution.
!
	MAXCH=-INCREASE
	IF(DECREASE .LT. 99.999_LDP)THEN
	  DECREASE=100.0_LDP*DECREASE/(100.0_LDP-DECREASE)
	ELSE
	  DECREASE=1.0E+07_LDP
	END IF
	MAXCH=MAX(MAXCH,DECREASE)
!	CALL DO_LEVEL_CHECK_MPI_V1(SOL_MAT,DST,DEND,ND,NT)
!
!**********************************************************************************
!**********************************************************************************
!
! Adjust scale parameter so that the biggest decrease in any variable
! (except T) is a factor of 20, the biggest increase is a factor of CHANGE_LIM.
! T is limited to a maximum change of 20% . Three options are used.
! In case I (LOCAL), the changes are scaled local according to the largest
! single change. In case II (NONE) no scaling is done - the changes for EACH
! variable are limited to the values given above.  Case III (MAJOR) is a
! combination of LOCAL and NONE scaling. We use NONE scaling for species
! whose population is significantly less then the Electron density. In case IV
! (GLOBAL), the biggest change at any depth is used to scale all changes.
! This option is probaly obsolete.
!
!**********************************************************************************
!**********************************************************************************

	IF(CHANGE_LIM .LE. 1.0_LDP)THEN
          WRITE(LUER,'(A,1PE12.4)')' Error in SOLVEBA_MPI_V1'
          WRITE(LUER,'(A,1PE12.4)')' Maximum change for normal iteration must be > 1.'
	  STOP
	END IF
	BIG_LIM=(CHANGE_LIM-1.0_LDP)/CHANGE_LIM		!Prevents division by zero and insures
	LIT_LIM=1.0_LDP-CHANGE_LIM                      !scale=1 for small changes
	MINSCALE=1.0_LDP
!
	IF(SCALE_OPT(1:5) .EQ. 'LOCAL')THEN
	    T1=MAX(BIG_LIM,DECREASE_SAVE)   	!Note + means decrease
	    T2=MIN(LIT_LIM,INCREASE_SAVE)     !Note - means increase
!
! Limit the change in T to a maximum of MAX_dT_COR, and ensure T > T_MIN.
!
	    T3=MAX( MAX_dT_COR,ABS(SOL_MAT(NT,I)) )
	    SCALE=MIN( MAX_dT_COR/T3,SCALE )
	    IF(SOL_MAT(NT,I) .NE. 0 .AND. POPS(NT,I) .GT. T_MIN .AND.
	1                      POPS(NT,I)*(1.0_LDP-SOL_MAT(NT,I)*SCALE) .LT. T_MIN)THEN
	      SCALE=(1.0_LDP-T_MIN/POPS(NT,I))/SOL_MAT(NT,I)
	    END IF
!
	  DO I=DST,DEND
	    DO J=1,NT
	      POPS(J,I)=POPS(J,I)*(1.0_LDP-SOL_MAT(J,I)*SCALE)
	    END DO
	    MINSCALE=MIN(SCALE,MINSCALE)
	  END DO
	  WRITE(LUER,'(A,1PE12.4)')' The local minimum value of scale is:',MINSCALE
!
!
	ELSE IF(SCALE_OPT(1:4) .EQ. 'NONE')THEN
	  DO I=DST,DEND
	    DO J=1,NT-1
	      IF(SOL_MAT(J,I) .GT. BIG_LIM)THEN
	        POPS(J,I)=POPS(J,I)*(1.0_LDP-BIG_LIM)
	        MINSCALE=MIN( BIG_LIM/SOL_MAT(J,I),MINSCALE )
	      ELSE IF(SOL_MAT(J,I) .LT. LIT_LIM)THEN
	        POPS(J,I)=POPS(J,I)*(1.0_LDP-LIT_LIM)
	        MINSCALE=MIN( LIT_LIM/SOL_MAT(J,I),MINSCALE )
	      ELSE
	        POPS(J,I)=POPS(J,I)*(1.0_LDP-SOL_MAT(J,I))
	      END IF
	    END DO
	    IF(POPS(NT,I) .LT. T_MIN)POPS(NT,I)=T_MIN
!
! Limit T to a 20% change if it is a variable.
!
	    SCALE=0.2_LDP/MAX( 0.2_LDP,ABS(SOL_MAT(NT,I)) )
	    POPS(NT,I)=POPS(NT,I)*(1.0_LDP-SOL_MAT(NT,I)*SCALE)
	    MINSCALE=MIN(SCALE,MINSCALE)
	  END DO
	  IF(MYPE .EQ. 0)WRITE(LUER,'(A,1PE12.4)')' The minimum value of scale for all species is:',MINSCALE
!
!
	ELSE IF(SCALE_OPT(1:5) .EQ. 'MAJOR')THEN
	  CALL FIDDLE_POP_CORRECTIONS_MPI_V1(POPS,SOL_MAT,T_MIN,CHANGE_LIM,MAX_dT_COR,
	1        SCALE_OPT,LAMBDA_IT,NT,DST,DEND,ND)
!
!
	ELSE			!Global Scaling !
	  T1=BIG_LIM		!Prevents division by zero and insures
	  T2=LIT_LIM   		!SCALE=1 if small changes.
	  DO I=DST,DEND
	    DO J=1,NT-1
	      T1=MAX(SOL_MAT(J,I),T1)			!Note + means decrease
	      T2=MIN(SOL_MAT(J,I),T2) 			!Note - means increase
	    END DO
	  END DO
	  SCALE=MIN( BIG_LIM/T1, LIT_LIM/T2)
	  DO I=DST,DEND
	    IF(SOL_MAT(NT,I) .NE. 0 .AND. POPS(NT,I) .GT. T_MIN .AND.
	1                      POPS(NT,I)*(1.0_LDP-SOL_MAT(NT,I)*SCALE) .LT. T_MIN)THEN
	        SCALE=(1.0_LDP-T_MIN/POPS(NT,I))/SOL_MAT(NT,I)
	    END IF
	  END DO
!
! Limit the change in T to a maximum of 20%.
!
	  SCALE=MIN( 0.2_LDP/MAX_dT_COR,SCALE )
	  IF(MYPE .EQ. 0)WRITE(LUER,'(A,1PE12.4)')' The value of scale is:',SCALE
	  T1=SCALE
	  CALL MPI_REDUCE(T1,SCALE,IONE,MPI_DOUBLE_PRECISION,MPI_MAX,IZERO,MPI_COMM_WORLD,IERR)
!
! Update the population levels (and the temperature) .
!
	  DO I=DST,DEND
	    DO J=1,NT
	      POPS(J,I)=POPS(J,I)*(1.0_LDP-SOL_MAT(J,I)*SCALE)
	    END DO
	  END DO
	END IF
!
	IF(SET_POPS_D2_EQ_D1)THEN                      ! .AND. .NOT. LAMBDA_IT)THEN
	  POPS(:,2)=POPS(:,1)*POP_ATOM(2)/POP_ATOM(1)
	END IF
!
! This adjust populations to their LTE values. This is only done if LTE_MODEL
! is set to TRUE, which is ascertained in SET_POPS_TO_LTE.
!
!	IF(LTE_MODEL)CALL SET_POPS_TO_PURE_LTE(POPS,NT,ND)
!
	RETURN
	END
