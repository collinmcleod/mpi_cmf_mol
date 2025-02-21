!
! Subroutine to reset the populations after the velocity has been adjusted.
! It is assumed the POP_SPEC and POP_ATOM have been reset.
!
! This routine was written assuming the radius grid (but not the boundary R values) have changed.
!
        SUBROUTINE WIND_SCALE_POPS_MPI_V1(POPS,R_OLD,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
        USE SET_KIND_MODULE
        USE MOD_CMFGEN
        USE MPI
	IMPLICIT NONE
!
! Under devlopment.
!
! Altered: 21-FEb-2025 -- Fixex compile-check bug (incosistent calls to IONTOPOP)
! Created: 10-Feb-2025 -- Based on wind_scale_pops_v1.f
!
	INTEGER ND
	INTEGER NT
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) R_OLD(ND)
	REAL(KIND=LDP) Z_POP(NT)
	LOGICAL DO_LEV_DISSOLUTION
!
! Local variables and vectors
!
	REAL(KIND=LDP) LOG_R(ND)
	REAL(KIND=LDP) LOG_R_OLD(ND)
!
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) GAM(ND)
	REAL(KIND=LDP) OLD_POPATOM(ND)
	REAL(KIND=LDP) T1
!
	INTEGER ID
	INTEGER ISPEC
	INTEGER I,J,K
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
!
	LOGICAL FIRST
!
	LOG_R=LOG(R)
	LOG_R_OLD=LOG(R_OLD)
!

	WRITE(6,'(/,A)')'***************************'
	WRITE(6,'(A,/)')'***************************'
	WRITE(6,*)' Warning -- wind_scale_pops_mpi_v1 has not been tested since its conversion to mpi.'
	WRITE(6,'(/,A)')'***************************'
	WRITE(6,'(A,/)')'***************************'
!
! Compute revised electron density. We assume the fractioanl ioization is constant.
!
! We only sum to NT-2 as NT-1 is Ne and NT is T.
!
	DO I=1,ND
	  T1=0.0_LDP
	  DO J=1,NT-2
	    T1=T1+POPS(J,I)
	  END DO
	  OLD_POPATOM(I)=T1
	  GAM(I)=ED(I)/T1
	END DO
!
	CALL MON_INTERP(ED,ND,IONE,LOG_R,ND,GAM,ND,LOG_R_OLD,ND)
	ED=ED*POP_ATOM
!
	TA=T
	CALL MON_INTERP(T,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
!
! Compute vector constants for evaluating the level dissolution. These
! constants are the same for all species. These are stored in a common
! block, and are required by SUP_TO_FULL and LTE_POP_WLD.
!
! As a first estimate of POPION, we simply assume the fractional ioization is
! constant.
!
	DO J=1,ND
	  TA(J)=0.0_LDP
	  DO I=1,NT-2
	    IF(Z_POP(I) .GT. 0.01_LDP)TA(J)=TA(J)+POPS(I,J)
	  END DO
	END DO
        DO I=1,ND
          POPION(I)=POP_ATOM(I)*TA(I)/OLD_POPATOM(I)
        END DO
        CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
! Since we need to interpolate on the grid, we need to work with the ROOT populations.
!
! Compute LOG(DC)
!
	IF(MYPE .EQ. 0)THEN
!
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
              IF(ATM(ID)%XzV_PRES)ROOT(ID)%XZV_F=LOG(ROOT(ID)%XzV_F)-ROOT(ID)%LOG_XzVLTE_F
	    END DO
	  END DO
!
! Interpolate onto the new R grid.
!
	  DO ISPEC=1,NUM_SPECIES
!
	    IF(SPECIES_BEG_ID(ISPEC) .GT. 0)THEN
	      ID=SPECIES_END_ID(ISPEC)-1
	      TA=LOG(ROOT(ID)%DxZV_F)
	      CALL MON_INTERP(TB,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	      ROOT(ID)%DXzV_F=EXP(TB)
	    END IF
!
	    DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	      IF(ATM(ID)%XzV_PRES)THEN
	        DO I=1,ATM(ID)%NXzV_F
	           TA(:)=ROOT(ID)%XzV_F(I,:)
	           CALL MON_INTERP(TB,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	        END DO
	      END IF
	    END DO
!
	    FIRST=.TRUE.
	    DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	      IF(ATM(ID)%XzV_PRES)THEN
	        CALL LTEPOP_WLD_V2(ROOT(ID)%XzVLTE_F, ROOT(ID)%LOG_XzVLTE_F,ROOT(ID)%W_XzV_F,
	1              ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,
	1              ATM(ID)%ZXzV,      ATM(ID)%GIONXzV_F,
	1              ATM(ID)%NXzV_F,    ROOT(ID)%DXzV_F,     ED,T,ND)
	        CALL CNVT_FR_DC_V2(ROOT(ID)%XzV_F, ROOT(ID)%LOG_XzVLTE_F,
	1              ROOT(ID)%DXzV_F,    ATM(ID)%NXzV_F,
	1              TB,                TA,ND,
	1              FIRST,             ATM(ID+1)%XzV_PRES)
                IF(ID .NE.  SPECIES_BEG_ID(ISPEC))ROOT(ID-1)%DXzV_F(1:ND)=TB(1:ND)
   	      END IF
	    END DO
!
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      IF(ATM(ID)%XzV_PRES)THEN
	         CALL SCALE_POPS(ROOT(ID)%XzV_F,ROOT(ID)%DXzV_F,
	1               POP_SPECIES(1,ISPEC),TA,ATM(ID)%NXzV_F,ND)
	      END IF
	    END DO
	  END DO
!
	  DO ID=NUM_IONS-1,1,-1
            CALL FULL_TO_SUP_MPI_V1(
	1        ROOT(ID)%XzV,   ROOT(ID)%NXzV,       ROOT(ID)%DXzV,  ATM(ID)%XzV_PRES,
	1        ROOT(ID)%XzV_F, ROOT(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ROOT(ID)%DXzV_F,
	1        ROOT(ID+1)%XzV, ATM(ID+1)%NXzV,      ATM(ID+1)%XzV_PRES, IONE, ND)
          END DO
!
! Store all quantities in POPS array. This is done here as it enables POPION
! to be readily computed. It also ensures that POPS is correct if we  don't
! iterate on T.
!
	  DO ID=1,NUM_IONS-1
	    CALL IONTOPOP(POPS,  ROOT(ID)%XzV, ROOT(ID)%DXzV, ED,T,
	1           ATM(ID)%EQXzV, ATM(ID)%NXzV, NT, IONE, ND, ND, ATM(ID)%XzV_PRES)
	  END DO
	END IF
!
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
! Scatter departure coefficients and ions to all processors.
!
        CALL SCATTER_XzV_F_AND_IONS(ND)
        IF(MYPE .EQ. 0)WRITE(6,*)'Scattered pops in wind_scale_pops_mpi_v1.f';FLUSH(UNIT=6)
!
!
! This section is idential to that in mod_subs/set_new_model_estimates.f.
!
! We now need to compute the populations for the model atom with Super-levels.
! We do this in reverse order (i.e. highest ionization stage first) in order
! that we the ion density for the lower ionization stage is available for
! the next call.
!
! For 1st call to FULL_TO_SUP, Last line contains FeX etc as FeXI not installed.
!
	DO ID=NUM_IONS-1,1,-1
	   CALL FULL_TO_SUP_MPI_V1(
	1      ATM(ID)%XzV,   ATM(ID)%NXzV,       ATM(ID)%DXzV,   ATM(ID)%XzV_PRES,
	1      ATM(ID)%XzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,
	1      ATM(ID+1)%XzV, ATM(ID+1)%NXzV,     ATM(ID+1)%XzV_PRES, DST, DEND)
	END DO
!
	IF(MYPE .EQ. 0)THEN
	  DO ID=NUM_IONS-1,1,-1
	     CALL FULL_TO_SUP_MPI_V1(
	1        ROOT(ID)%XzV,   ATM(ID)%NXzV,      ROOT(ID)%DXzV,   ATM(ID)%XzV_PRES,
	1        ROOT(ID)%XzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ROOT(ID)%DXzV_F,
	1       ROOT(ID+1)%XzV, ATM(ID+1)%NXzV,     ATM(ID+1)%XzV_PRES, IONE, ND)
	  END DO
	END IF
!
! Store all quantities in POPS array. This is done here as it enables POPION
! to be readily computed. It also ensures that POS is correct if we don't
! iterate on T.
!
	IF(MYPE .EQ. 0)THEN
	  DO ID=1,NUM_IONS-1
	    CALL IONTOPOP(POPS,  ROOT(ID)%XzV, ROOT(ID)%DXzV, ED,T,
	1          ATM(ID)%EQXzV, ATM(ID)%NXzV, NT, IONE, ND, ND, ATM(ID)%XzV_PRES)
 	  END DO
	END IF
	K=NT*ND
	CALL MPI_BCAST(POPS,K,MPI_DOUBLE_PRECISION,IZERO,MPI_COMM_WORLD,IERR)
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	IF(MYPE .EQ. 0)WRITE(6,*)'Done ION TO POP'
!
! Compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	DO J=DST,DEND
	  POPION(J)=0.0_LDP
	  DO I=1,NT
	     IF(Z_POP(I) .GT. 0.01_LDP)POPION(J)=POPION(J)+POPS(I,J)
	  END DO
	END DO
!
! Evaluates LTE populations for both the FULL atom, and super levels.
!
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	CALL EVAL_LTE_V5(DO_LEV_DISSOLUTION,ND)
	IF(MYPE .EQ. 0)CALL EVAL_ROOT_LTE_MPI_V1(DO_LEV_DISSOLUTION,ND)
!
	RETURN
	END
