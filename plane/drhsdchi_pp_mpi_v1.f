!
! Subroutine to compute the coefficient matrix (dimension ND*ND) for
! the variation in opacity for a plane-parallel atmosphere. The
! matrix is tridiagonal.  No velocity filed is present.
!
	SUBROUTINE dRHSdCHI_PP_MPI_V1(SOURCE,CHI,DTAU,COH_VEC,RJ,F,R,DIFF,DBB,DST,DEND,ND)
	USE SET_KIND_MODULE
	USE MOD_TRAP_DERIVATIVES
	USE MOD_VAR_OPAC_J, ONLY : TX, VDST, VDEND
	USE MPI 
	IMPLICIT NONE
!
! Created 19-MAr-2006. Based on VKIEFEAU_IBC & EDD_J_VAR_V5.F (spherical routine).
!                      This routine differs from the former routines in that
!                      Q and R^2 terms have been omitted.
!
	INTEGER DST,DEND,ND
	REAL(KIND=LDP) SOURCE(ND)
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) DTAU(ND)
	REAL(KIND=LDP) COH_VEC(ND)
	REAL(KIND=LDP) RJ(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) F(ND)
	REAL(KIND=LDP) DBB
	LOGICAL DIFF
!
! Local varaibles.
!
	REAL(KIND=LDP) LOC_W(ND,ND)
	REAL(KIND=LDP) dTAU_dCHI(ND,ND)
	REAL(KIND=LDP) ALPHA,BETA,T1
	REAL(KIND=LDP) UIJ,UII
	INTEGER I,J,K
	INTEGER IERR
!
! Compute dDTAU(I)/dCHI
!
	dTAU_dCHI(:,:)=0.0_LDP
	DO I=1,ND-1
	  K=I+1
	  ALPHA=0.5_LDP*(R(I)-R(K))
	  BETA=ALPHA*(R(I)-R(K))/6.0_LDP
	  IF(I .NE. 1)dTAU_dCHI(I,I-1)=-BETA*A(I)
	  dTAU_dCHI(I,I)= ALPHA + BETA*(A(K)-B(I))
	  dTAU_dCHI(I,K)= ALPHA + BETA*(B(K)-C(I))
	  IF(K .NE. ND)dTAU_dCHI(I,K+1)=BETA*C(K)
	END DO
!
! Compute W(I,J)=dRHS(I)/dCHI(J)
!
	TX(:,:,1)=0.0_LDP
	LOC_W=0.0_LDP
	DO I=MAX(2,DST),MIN(DEND,ND-1)
	  K=I+1
	  J=I-1
!
! We use UIJ for the sum of RHS terms containing 1/DTAU(I-1)
! We use UII for the sum of RHS terms containing 1/DTAU(I)
!
	  UIJ=-F(J)*RJ(J)/DTAU(J)/DTAU(J)
	  T1=-0.5_LDP*(1.0_LDP-COH_VEC(I))+F(I)/DTAU(J)/DTAU(J)
	  UIJ=UIJ+T1*RJ(I)+0.5_LDP*SOURCE(I)
	  UII=-F(K)*RJ(K)/DTAU(I)/DTAU(I)
	  T1=-0.5_LDP*(1.0_LDP-COH_VEC(I))+F(I)/DTAU(I)/DTAU(I)
	  UII=UII+T1*RJ(I)+0.5_LDP*SOURCE(I)
!
	  IF(J .NE. 1)LOC_W(I,J-1)=LOC_W(I,J-1)+UIJ*dTAU_dCHI(J,J-1)
	  LOC_W(I,J)=LOC_W(I,J)+UIJ*dTAU_dCHI(J,J)+UII*dTAU_dCHI(I,J)
	  LOC_W(I,I)=LOC_W(I,I)+UIJ*dTAU_dCHI(J,I)+UII*dTAU_dCHI(I,I)
	  LOC_W(I,K)=LOC_W(I,K)+UII*dTAU_dCHI(I,K)
!
	  LOC_W(I,I)=LOC_W(I,I)-0.5_LDP*(DTAU(J)+DTAU(I))*
	1                       (SOURCE(I)+COH_VEC(I)*RJ(I))/CHI(I)
!
	END DO
!
! outer boundary condition.
!
	IF(DST .EQ. 1)THEN
	  UII=-0.5_LDP*(1.0_LDP-COH_VEC(1))*RJ(1)+(RJ(1)*F(1)-RJ(2)*F(2))/DTAU(1)/DTAU(1)
	  UII=UII+0.5_LDP*SOURCE(1)
	  LOC_W(1,1)=LOC_W(1,1)+UII*dTAU_dCHI(1,1)
	  LOC_W(1,2)=LOC_W(1,2)+UII*dTAU_dCHI(1,2)
	  LOC_W(1,3)=LOC_W(1,3)+UII*dTAU_dCHI(1,3)
	  LOC_W(1,1)=LOC_W(1,1)-0.50_LDP*DTAU(1)*(RJ(1)*COH_VEC(1)+SOURCE(1))/CHI(1)
	END IF
!
! Inner boundary condition.
!
	IF(DEND .EQ. ND)THEN
	  UIJ=-0.5_LDP*(1.0_LDP-COH_VEC(ND))*RJ(ND)+(RJ(ND)*F(ND)-F(ND-1)*RJ(ND-1))/DTAU(ND-1)/DTAU(ND-1)
	  UIJ=UIJ+0.5_LDP*SOURCE(ND)
	  LOC_W(ND,ND-2)=LOC_W(ND,ND-2)+UIJ*dTAU_dCHI(ND-1,ND-2)
	  LOC_W(ND,ND-1)=LOC_W(ND,ND-1)+UIJ*dTAU_dCHI(ND-1,ND-1)
	  LOC_W(ND,ND)=LOC_W(ND,ND)+UIJ*dTAU_dCHI(ND-1,ND)
	  IF(DIFF)THEN
	    LOC_W(ND,ND)=LOC_W(ND,ND)-DBB/CHI(ND)/CHI(ND)/3.0_LDP
	  END IF
	  LOC_W(ND,ND)=LOC_W(ND,ND)-0.50_LDP*DTAU(ND-1)*(RJ(ND)*COH_VEC(ND)+SOURCE(ND))/CHI(ND)
	END IF
!
	I=ND*ND
	CALL MPI_ALLREDUCE(MPI_IN_PLACE,LOC_W,I,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	DO J=VDST,VDEND
	  TX(:,J,1)=LOC_W(:,J)
	END DO
!
	RETURN
	END
