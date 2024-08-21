!
! Subroutine to store the populations back into their individual storage
! locations. Replaces SUP_TO_FULL_V4.INC.
!
	SUBROUTINE SUP_TO_FULL_V4(POPS,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered 05-Apr-2011 : Now call SUP_TO_FULL_V3b (28-Nov-2010).
!                         Changes done to give a wider dynamic rangs in pops.
!
	INTEGER ND
	INTEGER NT
!
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) Z_POP(NT)
	LOGICAL DO_LEV_DISSOLUTION
!
	INTEGER I
	INTEGER J
	INTEGER ID
	INTEGER, PARAMETER :: IONE=1
	include 'mpif.h'
!
! We have now store the revised popuations back in their individual
! storage locations. For some species we have 2 atomic models. For these
! species we need to take the super-level populations and compute:
!
! 1. The LTE population off all leve ls in the FULL atom.
! 2. The population off all levels in the FULL atom.
! 3. The LTE population off all super-levels.
!
! This is done by the routine SUP_TO_FULL.
!
! Compute vector constants for evaluating the level disolution. These
! constants are the same for all species. These are stored in a module,
! and are required by SUP_TO_FULL. We first need to revise POPION.
!
	  DO J=1,ND
	    POPION(J)=0.0_LDP
	    DO I=1,NT
	      IF(Z_POP(I) .GT. 0.01_LDP)POPION(J)=POPION(J)+POPS(I,J)
	    END DO
	  END DO
	  WRITE(130+MYPE,*)ED(1),ED(ND),POPION(1),POPION(ND),T(1),T(ND)
	  FLUSH(UNIT=130+MYPE)
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
	  CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
	  WRITE(6,*)'Starting ATM -- now for root',DST, DEND; FLUSH(UNIT=6)
	  DO ID=1,NUM_IONS-1
	    WRITE(145+MYPE,'(3I5,4E14.4)')ID,DST,DEND,ATM(ID)%XzV(1,DST),ATM(ID)%XzV(1,DEND),
	1                                ATM(ID)%DXzV(DST),ATM(ID)%DXzV(DEND)
	    FLUSH(UNIT=145+MYPE)
	    IF(ATM(ID)%XzV_PRES)THEN
	      I=DEND-DST+1
	      WRITE(145+MYPE,*)ID
	      CALL WR2D(ATM(ID)%XzV(1,DST),ATM(ID)%NXzV,I,'ATM(ID)%XzV',145+MYPE)
	      CALL WR2D(ATM(ID)%DXzV(DST),IONE,I,'ATM(ID)%DXzV',145+MYPE)
	    END IF
	    CALL SUP_TO_FULL_V3b(
	1       ATM(ID)%XzV_F,        ATM(ID)%XzVLTE_F,    ATM(ID)%LOG_XzVLTE_F, ATM(ID)%XzVLTE_F_ON_S,
	1       ATM(ID)%W_XzV_F,      ATM(ID)%EDGEXzV_F,   ATM(ID)%GXzV_F,       ATM(ID)%F_TO_S_XzV,
	1       ATM(ID)%INT_SEQ_XzV , ATM(ID)%NXzV_F,      ATM(ID)%DXzV_F,
	1       ATM(ID)%GIONXzV_F,    ATM(ID)%XzV,         ATM(ID)%XzVLTE,       ATM(ID)%LOG_XzVLTE,
	1       ATM(ID)%NXzV,         ATM(ID)%DXzV,        ATM(ID)%ZXzV,         ATM(ID)%XzV_PRES,
	1       ATM(ID+1)%EDGEXzV_F,  ATM(ID+1)%GXzV_F,
	1       ATM(ID+1)%F_TO_S_XzV, ATM(ID+1)%NXzV_F,
	1       ATM(ID+1)%XzV_PRES,   T,ED, DST, DEND, ND)
	  END DO
50	  CONTINUE
	  WRITE(6,*)'Done ATM -- now for root',DST, DEND; FLUSH(UNIT=6)
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)!
	  IF(DST .EQ. 1)THEN
	    DO ID=1,NUM_IONS-1
	      CALL SUP_TO_FULL_V3b(
	1         ROOT(ID)%XzV_F,       ROOT(ID)%XzVLTE_F,    ROOT(ID)%LOG_XzVLTE_F, ROOT(ID)%XzVLTE_F_ON_S,
	1         ROOT(ID)%W_XzV_F,     ATM(ID)%EDGEXzV_F,   ATM(ID)%GXzV_F,       ATM(ID)%F_TO_S_XzV,
	1         ATM(ID)%INT_SEQ_XzV,  ATM(ID)%NXzV_F,     ROOT(ID)%DXzV_F,
	1         ATM(ID)%GIONXzV_F,   ROOT(ID)%XzV,        ROOT(ID)%XzVLTE,       ROOT(ID)%LOG_XzVLTE,
	1         ATM(ID)%NXzV,        ROOT(ID)%DXzV,        ATM(ID)%ZXzV,         ATM(ID)%XzV_PRES,
	1         ATM(ID+1)%EDGEXzV_F,  ATM(ID+1)%GXzV_F,
	1         ATM(ID+1)%F_TO_S_XzV, ATM(ID+1)%NXzV_F,
	1         ATM(ID+1)%XzV_PRES,   T,ED, IONE, ND, ND)
	    END DO
	  END IF
!
	RETURN
	END
