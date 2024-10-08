!
! This subroutine (which replaces EVAL_LTE_V4.INC) evaluates the LTE populations with respect to
! the ground state for
!
!  (a) The full atoms.
!  (B) The Super-Level model atom.
!
	SUBROUTINE EVAL_ROOT_LTE_MPI_V1(DO_LEV_DISSOLUTION,ND)
	USE SET_KIND_MODULE
	USE MPI
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created 26-Aug-2024
!
	INTEGER ND
	LOGICAL DO_LEV_DISSOLUTION
!
! Local variables.
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER I
	INTEGER ID
!	INCLUDE 'mpif.h'
!
! Revise vector constants for evaluating the level dissolution. These
! constants are the same for all species. These are stored in a common block,
! and are required by SUP_TO_FULL and LTE_POP_WLD.
!
	CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
! The final statements set the population of the ground state of the next
! ionizations stages. These must be set since they are used in determining
! opacities and photoinization terms.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    CALL LTEPOP_WLD_V2(
	1          ROOT(ID)%XzVLTE_F, ROOT(ID)%LOG_XzVLTE_F, ROOT(ID)%W_XzV_F,
	1          ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,  ATM(ID)%ZXzV,
	1          ATM(ID)%GIONXzV_F, ATM(ID)%NXzV_F,
	1          ROOT(ID)%DXzV_F,    ED,T,            IONE, ND, ND)
	    CALL LTE_POP_SL_V2(
	1          ROOT(ID)%XzVLTE,         ROOT(ID)%LOG_XzVLTE,    ROOT(ID)%dlnXzVLTE_dlnT,
	1          ATM(ID)%NXzV,            ROOT(ID)%XzVLTE_F,      ROOT(ID)%LOG_XzVLTE_F,
	1          ROOT(ID)%XzVLTE_F_ON_S,  ATM(ID)%EDGEXzV_F,      ATM(ID)%F_TO_S_XzV,
	1          ATM(ID)%NXzV_F,         ATM(ID)%XzV_PRES, T,     IONE, ND, ND)
	    IF(.NOT. ATM(ID+1)%XzV_PRES)THEN
	      DO I=1,ND
	        ROOT(ID+1)%XzV(1,I)=ROOT(ID)%DXzV_F(I)			!True if not present.
	        ROOT(ID+1)%XzVLTE(1,I)=ROOT(ID)%DXzV_F(I)			!True if not present.
	        ROOT(ID+1)%LOG_XzVLTE(1,I)=LOG(ROOT(ID)%DXzV_F(I))	!True if not present.
	        ROOT(ID+1)%dlnXzVLTE_dlnT(1,I)=0.0_LDP
	      END DO
	    END IF
	  END IF
	END DO
!
	RETURN
	END
