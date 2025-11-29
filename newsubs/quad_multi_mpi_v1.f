!
! Subroutine to compute the quadrature weights for the statistical
! equilibrium equations. These quadrature weight now have to be multplied
! by FQW/NU before use. This change was made to allow for a fixed continuum
! photioization cross-section.
!
! This routine is specifically designed for the handling of super levels.
! That is, we treat the process in a large atom but assume that the populations
! can be described by a smaller set of levels.
!
! Routine also handles level dissolution.
!
! Notation:
!
!         We use _F to denote populations and variables for the FULL atom,
!            with all terms and levels treated separately.
!	  We use _S to denote populations and variables for the SMALL model
!            atom, with many terms and levels treated as one (i.e using
!            SUPER levels).
!
	SUBROUTINE QUAD_MULTI_MPI_V1(WSE_S,dWSE_SdT,WCR,dWCRdT,
	1                       HNST_S,dlnHNST_S_dlnT,N_S,
	1                       HNST_F_ON_S,EDGE_F,N_F,
	1                       F_TO_S_MAPPING,NU_CONT,T,
	1                       DST,DEND,ND,
	1                       COMPUTE_BA,FIXED_T,LAST_ITERATION,
	1                       DESC,ZION,NPHOT,ID)
	USE SET_KIND_MODULE
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
	EXTERNAL SUB_PHOT_GEN
!
! Altered 01-Oct-2025 -- Removed call to SUB_PHOT_GEN_MPI_V2 
!                     -- Fixed bugs that have no effect.
!
	INTEGER ID,NPHOT
	INTEGER N_S,N_F,ND
	INTEGER DST, DEND
	REAL(KIND=LDP) WSE_S(N_S,DST:DEND,NPHOT)
	REAL(KIND=LDP) dWSE_SdT(N_S,DST:DEND,NPHOT)
	REAL(KIND=LDP) WCR(N_S,DST:DEND,NPHOT)
	REAL(KIND=LDP) dWCRdT(N_S,DST:DEND,NPHOT)
!
	REAL(KIND=LDP) HNST_S(N_S,DST:DEND)
	REAL(KIND=LDP) dlnHNST_S_dlnT(N_S,DST:DEND)
	REAL(KIND=LDP) HNST_F_ON_S(N_F,DST:DEND)
!
	REAL(KIND=LDP) EDGE_F(N_F)			!In 10^15 Hz
	INTEGER F_TO_S_MAPPING(N_F)
	REAL(KIND=LDP) T(ND)
!
	REAL(KIND=LDP) NU_CONT
	REAL(KIND=LDP) ZION
	CHARACTER*(*) DESC
	INTEGER IP,PHOT_ID
!
	REAL(KIND=LDP) YDIS(ND)		!Constant for computing level dissolution/
	REAL(KIND=LDP) XDIS(ND)		!Constant for computing level dissolution/
	REAL(KIND=LDP) DIS_CONST(N_F)	!Constant appearing in dissolution formula.
	REAL(KIND=LDP) ALPHA_VEC(N_F)
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) OPLIN,EMLIN
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL COMPUTE_BA,FIXED_T,LAST_ITERATION
!
! Local Variables,
!
	INTEGER I_S,I_F,J
	INTEGER IPROC
	INTEGER DPTH_INDX
	REAL(KIND=LDP) T1,T2,T3,ZION_CUBED,NEFF,FOUR_PI_D_H
	LOGICAL DO_ALL
!
!        INCLUDE 'mpif.h'
! 
!^L
! NB: WSE_OLD=WSE*FQW/NU
!     dWSEdT_OLD=dWSEdT*FQW/NU
!     WCR_OLD=(NU*WSE+WCR)*FQW/NU
!
! The factor of DEX(-10) in FOUR_PI_D_H is due to the definition of the
! cross-section in SUB_GEN_PHOT
! which is DEX(10) times the photoionization cross section so that
! CHI*R is constant.
!
! Note FOUR_PI_D_H differs by 10^-15 from original constant in QUADGEN because
! FQW has C units of Hz, not 10^15 Hz.
!
	FOUR_PI_D_H=1.0_LDP/5.27296E-03_LDP                    !1.8965D+02		!4*PI/H*DEX(-10)*DEX(-15)
!	CALL TUNE(1,'QUAD_ZERO')
	WSE_S=0.0_LDP
	WCR=0.0_LDP
	dWSE_SdT=0.0_LDP
	dWCRdT=0.0_LDP
!	CALL TUNE(2,'QUAD_ZERO')
!
! Get edge frequencies.
!
	DO IP=1,NPHOT
!
	  PHOT_ID=IP
!
! Get photoionization cross-sections for all levels. The first call returns
! the threshold cross-section when NU < EDGE.
!
!	  CALL TUNE(1,'SUB_PH2')
	  IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	    CALL GET_PHOT_CROSS_SECTIONS_V1(ALPHA_VEC,ID,PHOT_ID,N_F,NU_CONT,L_TRUE)
	  ELSE
	    CALL GET_PHOT_CROSS_SECTIONS_V1(ALPHA_VEC,ID,PHOT_ID,N_F,NU_CONT,L_FALSE)
	  END IF
!	  CALL TUNE(2,'SUB_PH2')
!
! DIS_CONST is the constant K appearing in the expression for level dissolution.
! A negative value f_CONST implies that the cross-section is zero.
!
! NB Edge frequencies ordered. Thus dont need to check remaining levels
! once NU_CONT > EDGE_FREQ.
!
!	  CALL TUNE(1,'SUB_DIS_ZERO')
	  DIS_CONST(1:N_F)=-1.0_LDP
!	  CALL TUNE(2,'SUB_DIS_ZERO')
!
!	  CALL TUNE(1,'SUB_DIS')
	  IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	    ZION_CUBED=ZION*ZION*ZION
	    DO I_F=1,N_F
	      IF(NU_CONT .GE. EDGE_F(I_F))EXIT
	      IF(ALPHA_VEC(I_F) .NE. 0 .AND. NU_CONT .GT. 0.8_LDP*EDGE_F(I_F))THEN
	        NEFF=SQRT(3.289395_LDP*ZION*ZION/(EDGE_F(I_F)-NU_CONT))
	        IF(NEFF .GT. 2*ZION)THEN
	          T1=MIN(1.0_LDP,16.0_LDP*NEFF/(1+NEFF)/(1+NEFF)/3.0_LDP)
	          DIS_CONST(I_F)=( T1*ZION_CUBED/(NEFF**4) )**1.5_LDP
	        END IF
	      END IF
	    END DO
	  END IF
!	  CALL TUNE(2,'SUB_DIS')
!
	  DO_ALL=.FALSE.
	  IF(COMPUTE_BA)DO_ALL=.TRUE.
	  IF(FIXED_T)DO_ALL=.FALSE.
	  IF(LAST_ITERATION)DO_ALL=.TRUE.
!
	  DO DPTH_INDX=DST,DEND
!
! Compute dissolution vectors that are independent of level.
!
	    IF(MOD_DO_LEV_DIS)THEN
	      J=DPTH_INDX
	      YDIS(J)=1.091_LDP*(X_LEV_DIS(J)+4.0_LDP*(ZION-1)*A_LEV_DIS(J))*
	1               B_LEV_DIS(J)*B_LEV_DIS(J)
	      XDIS(J)=B_LEV_DIS(J)*X_LEV_DIS(J)
	    END IF
!
! We have to loop over depth (rather than frequency) because of the
! FULL to SUPER level mapping.
!
! Note that WSE_S and WCR have alternate signs.
!
	    IF(DO_ALL)THEN

	        J=DPTH_INDX
	        DO I_F=1,N_F
	          I_S=F_TO_S_MAPPING(I_F)
	          IF(ALPHA_VEC(I_F) .LE. 0.0_LDP)THEN
	          ELSE IF(NU_CONT .GE. EDGE_F(I_F))THEN
	            T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	            WSE_S(I_S,J,IP)=WSE_S(I_S,J,IP) + T1*HNST_F_ON_S(I_F,J)
	            WCR(I_S,J,IP)=WCR(I_S,J,IP) - EDGE_F(I_F)*T1*HNST_F_ON_S(I_F,J)
	            T2=T1*HNST_F_ON_S(I_F,J)*(dlnHNST_S_dlnT(I_S,J)+1.5_LDP+HDKT*EDGE_F(I_F)/T(J))/T(J)
	            dWSE_SdT(I_S,J,IP)=dWSE_SdT(I_S,J,IP) - T2
	            dWCRdT(I_S,J,IP)=dWCRdT(I_S,J,IP) + EDGE_F(I_F)*T2
!
! We only allow for level dissolutions when the ionizations are occurring to
! the ground state.
!
	          ELSE IF(DIS_CONST(I_F) .GE. 0.0_LDP)THEN
	            T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	            T2=7.782_LDP+XDIS(J)*DIS_CONST(I_F)
	            T3=T2/(T2+YDIS(J)*DIS_CONST(I_F)*DIS_CONST(I_F))
	            IF(T3 .GT. PHOT_DIS_PARAMETER)THEN
	              T3=T1*T3
	              WSE_S(I_S,J,IP)=WSE_S(I_S,J,IP) + T3*HNST_F_ON_S(I_F,J)
	              WCR(I_S,J,IP)=WCR(I_S,J,IP) - EDGE_F(I_F)*T3*HNST_F_ON_S(I_F,J)
	              T2=T3*HNST_F_ON_S(I_F,J)*(dlnHNST_S_dlnT(I_S,J)+1.5_LDP+HDKT*EDGE_F(I_F)/T(J))/T(J)
	              dWSE_SdT(I_S,J,IP)=dWSE_SdT(I_S,J,IP) - T2
	              dWCRdT(I_S,J,IP)=dWCRdT(I_S,J,IP) + EDGE_F(I_F)*T2
	            END IF
	          END IF
	        END DO
	      ELSE
!
	        J=DPTH_INDX
	        DO I_F=1,N_F
	          I_S=F_TO_S_MAPPING(I_F)
	          IF(ALPHA_VEC(I_F) .LE. 0.0_LDP)THEN
	          ELSE IF(NU_CONT .GE. EDGE_F(I_F))THEN
	            T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	            WSE_S(I_S,J,IP)=WSE_S(I_S,J,IP) + T1*HNST_F_ON_S(I_F,J)
	            WCR(I_S,J,IP)=WCR(I_S,J,IP) - EDGE_F(I_F)*T1*HNST_F_ON_S(I_F,J)
!
! We only allow for level dissolutions when the ionizations are occurring to
! the ground state.
!
	          ELSE IF(DIS_CONST(I_F) .GE. 0.0_LDP)THEN
	            T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	            T2=7.782_LDP+XDIS(J)*DIS_CONST(I_F)
	            T3=T2/(T2+YDIS(J)*DIS_CONST(I_F)*DIS_CONST(I_F))
	            IF(T3 .GT. PHOT_DIS_PARAMETER)THEN
	              T3=T1*T3
	              WSE_S(I_S,J,IP)=WSE_S(I_S,J,IP) + T3*HNST_F_ON_S(I_F,J)
	              WCR(I_S,J,IP)=WCR(I_S,J,IP) - EDGE_F(I_F)*T3*HNST_F_ON_S(I_F,J)
	            END IF
	          END IF
	        END DO
	      END IF
	  END DO		!Depth index
	END DO                  !Photon route
!
	RETURN
	END
