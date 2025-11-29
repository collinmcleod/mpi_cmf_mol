!
! This routine can be compiled instead of set_phot_cross-sections.f. Make
! sure other versions are not beeing compiled, and have been deleted from
! libsubs.a
!
        MODULE MOD_PHOT_THREAD_INFO
        USE SET_KIND_MODULE
        IMPLICIT NONE
        REAL(KIND=LDP), ALLOCATABLE ::  COMP_FREQ(:)
	SAVE
	END MODULE  MOD_PHOT_THREAD_INFO
!
! Subroutine to compute the photoionization cross-sections for all species.
! In this routine,  the cross-sections are computed by each processor.
! There is no parallel processing.
!
	SUBROUTINE SET_PHOT_CROSS_SECTIONS_V1(FREQ,NU,NU_EVAL_CONT,CUR_ML,NCF)
	USE SET_KIND_MODULE
	USE PHOT_DATA_MOD
	USE MOD_PHOT_THREAD_INFO
	USE MOD_CMFGEN
	IMPLICIT NONE
!
	INTEGER NCF
	INTEGER CUR_ML
	REAL(KIND=LDP) FREQ
	REAL(KIND=LDP) NU(NCF)
	REAL(KIND=LDP) NU_EVAL_CONT(NCF)
!
	REAL(KIND=LDP), SAVE, ALLOCATABLE :: ALPHA_VEC(:)
	LOGICAL, SAVE :: FIRST=.TRUE.
!
	REAL(KIND=LDP) LOC_FREQ
	INTEGER ID
	INTEGER PHOT_ID		!Photoionization ID (path)
!
	INTEGER I,K,ML
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	IF(FIRST)THEN
	  ALLOCATE(COMP_FREQ(0:NTHREAD-1))
	  COMP_FREQ=0.0_LDP; K=0
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)K=MAX(K,ATM(ID)%NXZV_F)
	  END DO
	  ALLOCATE(ALPHA_VEC(K))
	  FIRST = .FALSE.
	  IF(MYPE .EQ. 0)THEN
	    WRITE(6,*)'Using single processore versions to compute photoionization cross-sections.'
	    WRITE(6,*)'Each processor computes their own photoionization cross-sections.'
            WRITE(6,*)'Called file is set_phot_cross_sections_v1.f'
	    WRITE(6,*)'Allocated COMP_FREQ,ALPHA_VEC in SET_PHOT_CROSS_SECTIONS_V1'
	  END IF
	  FLUSH(UNIT=6)
	END IF
!
	IF(FREQ .NE. 0)THEN
          LOC_FREQ=FREQ
	  IF(LOC_FREQ .NE. 0.0_LDP)THEN
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        DO I=1,ATM(ID)%N_XzV_PHOT
	          PHOT_ID=I
	          CALL SUB_PHOT_GEN_MPI_V2(ID,PD(ID)%CUR_CROSS(:,PHOT_ID),LOC_FREQ,ATM(ID)%EDGEXZV_F,ATM(ID)%NXzV_F,PHOT_ID,L_FALSE)
	        END DO
	      END IF
	    END DO
	  END IF
!
! 
	ELSE
!
! Compute the photoionization cross-sections for all levels at their edge frequencies.
! Only done for PHOT_ID .EQ. 1 at present.
!
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO I=1,ATM(ID)%N_XzV_PHOT
	        PHOT_ID=I
	        IF(FIRST .AND. PHOT_ID .EQ. 1)THEN
	          CALL SUB_PHOT_GEN_MPI_V2(ID,ALPHA_VEC,FREQ,ATM(ID)%EDGEXZV_F,ATM(ID)%NXzV_F,PHOT_ID,L_TRUE)
	        END IF
	      END DO
	    END IF
	  END DO
	END IF
!
	RETURN
	END
!
	SUBROUTINE GET_PHOT_CROSS_SECTIONS_V1(PHOT,ID,PHOT_ID,NLEVS,FREQ,ADD_EDGE)
	USE SET_KIND_MODULE
	USE PHOT_DATA_MOD
	USE MOD_PHOT_THREAD_INFO
	USE MPI
	IMPLICIT NONE
!
	INTEGER ID
	INTEGER PHOT_ID
	INTEGER NLEVS
	REAL(KIND=LDP) PHOT(NLEVS)
	REAL(KIND=LDP) FREQ
	LOGICAL ADD_EDGE
!
	REAL(KIND=LDP) T1
	INTEGER I,K
	INTEGER L(2)
	INTEGER IERR
	INTEGER NU_THRD
!	
!
! Set the photoionization cross-sections on all nodes for the current
! frequency.
!
	NU_THRD=MYPE
	PHOT(1:NLEVS)=PD(ID)%CUR_CROSS(1:NLEVS,PHOT_ID)
!
	IF(ADD_EDGE .AND. PHOT_ID .EQ. 1)THEN
	  DO I=1,NLEVS
	    T1=PD(ID)%EDGE_FREQ(I,PHOT_ID)
	    IF(FREQ .GE. 0.5_LDP*T1 .AND. FREQ .LT. T1)THEN
	      PHOT(I)=PD(ID)%EDGE_CROSS(I,PHOT_ID)
	    END IF
	  END DO
	END IF
!
!
	RETURN
	END
!
	SUBROUTINE GET_EDGE_CROSS_SECTIONS_V1(PHOT,ID,PHOT_ID,NLEVS)
	USE SET_KIND_MODULE
	USE PHOT_DATA_MOD
	USE MOD_PHOT_THREAD_INFO
!
	INTEGER ID
	INTEGER PHOT_ID
	INTEGER NLEVS
	REAL(KIND=LDP) PHOT(NLEVS)
!
	PHOT(1:NLEVS)=PD(ID)%EDGE_CROSS(1:NLEVS,PHOT_ID)
!
	RETURN
	END
