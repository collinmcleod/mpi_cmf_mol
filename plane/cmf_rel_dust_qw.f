!
! Subroutined designed to compute the quadrature weights for anisotropic
! dust scattering. ad Henney-Greensteing phase function is assumed.
!
! The weights are retuned via mod_spac_grid_v2
!
	SUBROUTINE CMF_REL_DUST_QW(G_HG,ND,NP)
	USE SET_KIND_MODULE
	USE MOD_SPACE_GRID_V2
	IMPLICIT NONE
!
! Finalized:  09-Mar-2024
!
	INTEGER NP
	INTEGER ND
!
	REAL(KIND=LDP), ALLOCATABLE :: MU_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: HG_FUN(:)
	REAL(KIND=LDP), ALLOCATABLE :: COEF(:,:)
!
	REAL(KIND=LDP) G_HG
	REAL(KIND=LDP) dMU
	REAL(KIND=LDP), ALLOCATABLE :: CON_SUM(:)
	REAL(KIND=LDP), ALLOCATABLE :: LIN_SUM(:)
!
	INTEGER MAX_NANG
	INTEGER L,K
	INTEGER ID
	INTEGER IP
	INTEGER NANG
	INTEGER NMU
	INTEGER LM,LP
	INTEGER KM,KP
	INTEGER LUOUT
	LOGICAL, PARAMETER :: CHECK_QUAD_ACCURACY=.TRUE.
!
! We get the maximumcray length so the that we don't need to continually allocate/
! deallocate the vectors MU and HG_FUN.
!
	MAX_NANG=2*NP
	ALLOCATE(MU_VEC(MAX_NANG))		!2* as two directions -- _P amd _M
	ALLOCATE(HG_FUN(MAX_NANG))
	ALLOCATE(CON_SUM(NP))
	ALLOCATE(LIN_SUM(NP))
!
! These vectore are defied as follows:
!
!   DW(IP)%WGT_P_TO_M(K,ID) gives the scattering emissivty in the M direction
!      for ray IP arisng from radiaton in the P directions at angles specified by K.
!
	DO L=1,NP
	  ALLOCATE(DW(L)%WGT_P_TO_P(NP,ND));     DW(L)%WGT_P_TO_P=0.0_LDP
	  ALLOCATE(DW(L)%WGT_P_TO_M(NP,ND));     DW(L)%WGT_P_TO_M=0.0_LDP
	  ALLOCATE(DW(L)%WGT_M_TO_P(NP,ND));     DW(L)%WGT_M_TO_P=0.0_LDP
	  ALLOCATE(DW(L)%WGT_M_TO_M(NP,ND));     DW(L)%WGT_M_TO_M=0.0_LDP
	END DO
!
! As we are using a sphrical atmosphere, the mu vary with depth.
!
	DO ID=1,ND
!
! Get all MU at a given depth.
!
          NMU=NP+1-ID
	  NANG=2*NMU-1
          DO IP=1,NMU
            MU_VEC(IP)=RAY(IP)%MU_P(RAY(IP)%LNK(ID))
          END DO
          DO IP=1,NMU
            MU_VEC(NANG-IP+1)=RAY(IP)%MU_M(RAY(IP)%LNK(ID))
	   WRITE(188,'(2I5,2ES18.8)')ID,IP,RAY(IP)%MU_P(RAY(IP)%LNK(ID)),RAY(IP)%MU_M(RAY(IP)%LNK(ID))
	  END DO
	  DO IP=1,NANG
	    WRITE(187,'(2I5,ES18.8)')ID,IP,MU_VEC(IP)
	  END DO
!
! The HG fun is a function of MU, MU'. For each value of MU we need to compute the
! HD function as a function of MU'. To get the scattered contribution for we will sum up
! I(mu')HG(mu,mu') over mu.
!
! To compute the weight we approximate the HG function as a monotonic cubic, and then assume
! that I is a linear function of my.
!
	  IF(ALLOCATED(COEF))DEALLOCATE(COEF)
	  ALLOCATE(COEF(NANG,4)); COEF=0.0_LDP
	  DO L=1,NMU
	    CALL GET_HG_FUN(HG_FUN,MU_VEC,NANG,G_HG,MU_VEC(L))
	    CALL MON_INT_FUNS_V2(COEF,HG_FUN,MU_VEC,NANG)
	    DO K=1,NMU-1
	      dMU=(MU_VEC(K+1)-MU_VEC(K))
	      DW(L)%WGT_P_TO_P(K,ID)  =DW(L)%WGT_P_TO_P(K,ID)  - dMU*(COEF(K,4)/2+dMU*(COEF(K,3)/3 +dMU*(COEF(K,2)/4 + dMU*COEF(K,1)/5)))
	      DW(L)%WGT_P_TO_P(K+1,ID)=DW(L)%WGT_P_TO_P(K+1,ID)- dMU*(COEF(K,4)/2+dMU*(COEF(K,3)/6 +dMU*(COEF(K,2)/12 +dMU*COEF(K,1)/20)))
	    END DO
	    DO K=NMU,NANG-1
	      KM=NANG+1-K
	      dMU=(MU_VEC(K+1)-MU_VEC(K))
	      DW(L)%WGT_M_TO_P(KM,ID)  =DW(L)%WGT_M_TO_P(KM,ID)  -dMU*(COEF(K,4)/2+dMU*(COEF(K,3)/3 +dMU*(COEF(K,2)/4 + dMU*COEF(K,1)/5)))
	      DW(L)%WGT_M_TO_P(KM-1,ID)=DW(L)%WGT_M_TO_P(KM-1,ID)-dMU*(COEF(K,4)/2+dMU*(COEF(K,3)/6 +dMU*(COEF(K,2)/12 +dMU*COEF(K,1)/20)))
	    END DO
	    LM=NMU+1-L
	    DO K=1,NMU-1
	      dMU=(MU_VEC(K+1)-MU_VEC(K))
	      DW(LM)%WGT_P_TO_M(K,ID)  =DW(LM)%WGT_P_TO_M(K,ID)  -dMU*(COEF(K,4)/2+dMU*(COEF(K,3)/3 +dMU*(COEF(K,2)/4 + dMU*COEF(K,1)/5)))
	      DW(LM)%WGT_P_TO_M(K+1,ID)=DW(LM)%WGT_P_TO_M(K+1,ID)-dMU*(COEF(K,4)/2+dMU*(COEF(K,3)/6 +dMU*(COEF(K,2)/12 +dMU*COEF(K,1)/20)))
	    END DO
	    DO K=NMU,NANG-1
	      KM=NANG-K+1
	      dMU=(MU_VEC(K+1)-MU_VEC(K))
	      DW(LM)%WGT_M_TO_M(KM,ID)  =DW(LM)%WGT_M_TO_M(KM,ID)  -dMU*(COEF(K,4)/2+dMU*(COEF(K,3)/3 +dMU*(COEF(K,2)/4 + dMU*COEF(K,1)/5)))
	      DW(LM)%WGT_M_TO_M(KM-1,ID)=DW(LM)%WGT_M_TO_M(KM-1,ID)-dMU*(COEF(K,4)/2+dMU*(COEF(K,3)/6 +dMU*(COEF(K,2)/12 +dMU*COEF(K,1)/20)))
	    END DO
	  END DO
	END DO
	DEALLOCATE(COEF)
!
	IF(CHECK_QUAD_ACCURACY)THEN
!
	  CALL GET_LU(LUOUT,'CMF_REL_DUST_QW')
	  OPEN(UNIT=LUOUT,STATUS='UNKNOWN',ACTION='WRITE',FILE='DUST_QW_CHECK')
!
! Check the accuracy of the quadrature.
!
	  DO ID=1,ND
            NMU=NP+1-ID
	    DO L=1,NMU
	      CON_SUM(L)=0.0_LDP
	      LIN_SUM(L)=0.0_LDP
	      DO K=1,NMU
	        CON_SUM(L)=CON_SUM(L) +
	1                    DW(L)%WGT_P_TO_P(K,ID) +
	1                    DW(L)%WGT_M_TO_P(K,ID)
	        LIN_SUM(L)=LIN_SUM(L) + 
	1                    DW(L)%WGT_P_TO_P(K,ID)*RAY(K)%MU_P(RAY(K)%LNK(ID)) +
	1                    DW(L)%WGT_M_TO_P(K,ID)*RAY(K)%MU_M(RAY(K)%LNK(ID))
	      END DO
	      WRITE(LUOUT,'(2I5,3ES16.8)')ID,L,RAY(L)%MU_P(RAY(L)%LNK(ID)),CON_SUM(L),LIN_SUM(L)
	      IF(CON_SUM(L) .GT. 1)DW(L)%WGT_P_TO_P(1:NMU,ID)=DW(L)%WGT_P_TO_P(1:NMU,ID)/CON_SUM(L)
	      IF(CON_SUM(L) .GT. 1)DW(L)%WGT_P_TO_M(1:NMU,ID)=DW(L)%WGT_P_TO_M(1:NMU,ID)/CON_SUM(L)
	    END DO
!
	    DO L=1,NMU-1
	      CON_SUM(L)=0.0_LDP
	      LIN_SUM(L)=0.0_LDP
	      DO K=1,NMU
	        CON_SUM(L)=CON_SUM(L) +
	1                     DW(L)%WGT_P_TO_M(K,ID) +
	1                     DW(L)%WGT_M_TO_M(K,ID)
	        LIN_SUM(L)=LIN_SUM(L) + 
	1                     DW(L)%WGT_P_TO_M(K,ID)*RAY(K)%MU_P(RAY(K)%LNK(ID)) +
	1                     DW(L)%WGT_M_TO_M(K,ID)*RAY(K)%MU_M(RAY(K)%LNK(ID))
	      END DO
	      WRITE(LUOUT,'(2I5,3ES16.8)')ID,L,RAY(L)%MU_M(RAY(L)%LNK(ID)),CON_SUM(L),LIN_SUM(L)
	      IF(CON_SUM(L) .GT. 1)DW(L)%WGT_P_TO_M(1:NMU,ID)=DW(L)%WGT_P_TO_M(1:NMU,ID)/CON_SUM(L)
	      IF(CON_SUM(L) .GT. 1)DW(L)%WGT_M_TO_M(1:NMU,ID)=DW(L)%WGT_M_TO_M(1:NMU,ID)/CON_SUM(L)
	    END DO
	    FLUSH(LUOUT)
!
	  END DO
	  CLOSE(LUOUT)
	END IF
	WRITE(6,*)'Called cmf_rel_dust_qw.f'; FLUSH(UNIT=6)
!
	RETURN
	END
