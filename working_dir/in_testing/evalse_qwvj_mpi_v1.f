!
! Subroutine to increment the statistical equilibrium equations for each
! depth point given the value of the mean intensity at each depth point.
!
! Subroutine also increments the QFV matrix that describe the  variation
! of the SE quations with respect to RJ.
!
! Routine also increments the ionization equilibrium equations.
!
	SUBROUTINE EVALSE_QWVJ_MPI_V1(ID,WSE,WCR,
	1                    HN,HNST,NLEV,ION_LEV,
	1                    DI,LOG_DIST,N_DI,NPHOT,
	1                    JREC,JPHOT,JREC_CR,JPHOT_CR,
	1                    DST,DEND, NT,ND)
	USE SET_KIND_MODULE
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Created 16-Aug-2024 - Based on EVALSE_QWVJ_V8.
!
	INTEGER ION_LEV(NPHOT)		!Levl ID of target in DI (i.e. the ION).
	INTEGER ID            !Ionic species identifier.
	INTEGER NLEV          !Number of atomic levels
	INTEGER N_DI          !Number of atomic levels in next ioization satge.
	INTEGER NT		!Total number of levels
	INTEGER ND		!Number of depth points
	INTEGER DST
	INTEGER DEND
	INTEGER ATM(ID)%N_XzV_PHOT
!
! NB --- NION is the total number of ionic species i.e. for
! HI,HII,CI,CII,CIII,CIV,CV would have NION=5 (dont count HII and CV).
!
	REAL(KIND=LDP) HN(NLEV,DST:DEND)
	REAL(KIND=LDP) HNST(NLEV,DST:DEND)
	REAL(KIND=LDP) WSE(NLEV,DST:DEND)
	REAL(KIND=LDP) WCR(NLEV,DST:DEND)
	REAL(KIND=LDP) DI(N_DI,DST:DEND)
	REAL(KIND=LDP) LOG_DIST(N_DI,DST:DEND)
!
	REAL(KIND=LDP) JREC(DST:DEND)
	REAL(KIND=LDP) JPHOT(DST:DEND)
	REAL(KIND=LDP) JREC_CR(DST:DEND)
	REAL(KIND=LDP) JPHOT_CR(DST:DEND)
!
! Constants for opacity etc.
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local variables.
!
	INTEGER I,J
	REAL(KIND=LDP) NETR
	REAL(KIND=LDP) T1,T2
!
! REV_HNST referes to the LTE population  of the level defined with respect
! to the actual destination (target) level.
!
	REAL(KIND=LDP) REV_HNST
!
	REAL(KIND=LDP) SUM_SE,SUM_VJ_R,SUM_VJ_P
	REAL(KIND=LDP) B_RAT
	REAL(KIND=LDP) PLANCKS_CONSTANT,PC
	EXTERNAL PLANCKS_CONSTANT
!
	INTEGER IPROC
	INTEGER, SAVE :: FIRST=.TRUE.
	INTEGER, SAVE :: MYPE
	INTEGER, SAVE :: NTHREAD
	INTEGER, SAVE :: NUM_DEPTHS_PER_THREAD
	INCLUDE 'mpif.h'
!
	INTEGER GET_DEPTH_INDX
	GET_DEPTH_INDX(IPROC,MYPE,NUM_DEPTHS_PER_THREAD)=IPROC+MYPE*NUM_DEPTHS_PER_THREAD
! 
!^L
	IF(FIRST)THEN
	  FIRST=.FALSE.
	  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYPE,IERR)
	  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NTHREAD,IERR)
	  NUM_DEPTHS_PER_THREAD=(ND-1)/NTHREAD+1
	END IF
!
	PC=1.0E+15_LDP*PLANCKS_CONSTANT()            !1.0D+15 due to units of NU.
!
	DO IL=1,NPHOT
	  IF(ION_LEV(IL) .EQ. 0)RETURN
!
! The net ionization (collisional and radaitive) to the last ionization stage
! must be zero from the sum of the previous equilibrum equations. Hence
! there is no need for a rate equation for the final species - it is
! preserved for the abundance equation.
!
! If there only ionizations to the ground state, the net ionization
! term could be neglected from the rate equation for that level.
! However, ionizations to excited levels, and Auger ionization (to
! another ionization stage) mean the terms have to be explicitly included.
!
! TMP_HST is the LTE population relative to the target level in the ion.
! REV_HNST= HNST * B(ION_LEV)/B(1) where b is the deparure coefficient.
!

	  DO IPROC=1,NUM_DEPTHS_PER_THREAD
	    J=GET_DEPTH_INDX(IPROC,MYPE,NUM_DEPTHS_PER_THREAD)
	    SUM_SE=0.0_LDP
	    SUM_VJ_R=0.0_LDP
	    SUM_VJ_P=0.0_LDP
	    B_RAT=EXP(LOG_DIST(1,J)-LOG_DIST(ION_LEV,J))*(DI(ION_LEV,J)/DI(1,J))
	    DO I=1,NLEV
	      IF(WSE(I,J) .NE. 0.0_LDP)THEN
	        REV_HNST=HNST(I,J)*B_RAT
	        NETR=WSE(I,J)*( REV_HNST*JREC(J)-HN(I,J)*JPHOT(J) )
	        SE(ID)%STEQ(I,J)=SE(ID)%STEQ(I,J)+NETR
	        SUM_SE=SUM_SE+NETR
	        SE(ID)%QFV_R(I,J)=SE(ID)%QFV_R(I,J)+WSE(I,J)*REV_HNST
	        SUM_VJ_R=SUM_VJ_R+WSE(I,J)*REV_HNST
	        SE(ID)%QFV_P(I,J)=SE(ID)%QFV_P(I,J)+WSE(I,J)*HN(I,J)
	        SUM_VJ_P=SUM_VJ_P+WSE(I,J)*HN(I,J)
!
! The next 5 lines update the EHB equation, and compute two matrices
! needed for its variation.
!
	        T1=WSE(I,J)*(HN(I,J)*JPHOT_CR(J)-REV_HNST*JREC_CR(J))
	        T2=WCR(I,J)*(HN(I,J)*JPHOT(J)-REV_HNST*JREC(J))
	        SE(ID)%T_EHB(J)=SE(ID)%T_EHB(J)+PC*(T1+T2)
	        SE(ID)%QFV_R_EHB(I,J)=SE(ID)%QFV_R_EHB(I,J)+WCR(I,J)*REV_HNST
	        SE(ID)%QFV_P_EHB(I,J)=SE(ID)%QFV_P_EHB(I,J)+WCR(I,J)*HN(I,J)
	      END IF
	    END DO
!
	    I=SE(ID)%ION_LEV_TO_EQ_PNT(ION_LEV(IL))
	    SE(ID)%STEQ(I,J)=SE(ID)%STEQ(I,J)-SUM_SE
	    SE(ID)%QFV_R(I,J)=SE(ID)%QFV_R(I,J)-SUM_VJ_R
	    SE(ID)%QFV_P(I,J)=SE(ID)%QFV_P(I,J)-SUM_VJ_P
	  END DO
	END DO
!
	RETURN
	END
