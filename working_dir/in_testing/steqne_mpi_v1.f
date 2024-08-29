!
! Routine to increment the charge conservation equation, and the
! variation charge equation. This was originally done in
! steqheii.
!
	SUBROUTINE STEQNE_MPI_V1(ED,NT,DIAG_INDX,ND,COMPUTE_BA,DST,DEND)
	USE SET_KIND_MODULE
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Created 16-Aug-2024 : Based on ! /ihome/dhillier/hillier/cur_cmf/new_main/subs/steqne_v4.f
!
	INTEGER NT
	INTEGER DIAG_INDX
	INTEGER ND
	INTEGER DST,DEND
!
	REAL(KIND=LDP) ED(ND)
	LOGICAL COMPUTE_BA
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local varaiables.
!
	INTEGER K,M,JJ
!
	DO K=DST,DEND
	  STEQ_ED(K)=STEQ_ED(K)-ED(K)
	END DO
!
	IF(COMPUTE_BA)THEN
 	  DO K=DST,DEND
	    BA_ED(NT-1,DIAG_INDX,K)=BA_ED(NT-1,DIAG_INDX,K)-1.0_LDP
	  END DO
	END IF
!
	RETURN
	END
