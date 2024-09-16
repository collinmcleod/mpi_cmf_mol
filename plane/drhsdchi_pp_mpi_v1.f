!
! Subroutine to compute the coefficient matrix (dimension ND*ND) for
! the variation in opacity for a plane-parallel atmosphere. The
! matrix is tridiagonal.  No velocity filed is present.
!
	SUBROUTINE dRHSdCHI_PP_MPI_V1(W,SOURCE,CHI,DTAU,COH_VEC,RJ,F,R,DIFF,DBB,DST,DEND,ND)
	USE SET_KIND_MODULE
	USE MOD_TRAP_DERIVATIVES
	IMPLICIT NONE
!
! Created 19-MAr-2006. Based on VKIEFEAU_IBC & EDD_J_VAR_V5.F (spherical routine).
!                      This routine differs from the former routines in that
!                      Q and R^2 terms have been omitted.
!
	INTEGER ND
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
! Output:
!
	REAL(KIND=LDP) W(ND,DST:DEND)
!
! Local varaibles.
!
	REAL(KIND=LDP) dTAU_dCHI(ND,DST:DEND)
	REAL(KIND=LDP) ALPHA,BETA,T1
	REAL(KIND=LDP) UIJ,UII
	INTEGER I,J,K
!
! Compute dDTAU(I)/dCHI
!
	dTAU_dCHI(:,:)=0.0_LDP
	DO I=MAX(1,DST-1),MIN(DEND+1,ND-1)
	  K=I+1
	  ALPHA=0.5_LDP*(R(I)-R(K))
	  BETA=ALPHA*(R(I)-R(K))/6.0_LDP
	  IF(I .GT. MAX(DST,2))dTAU_dCHI(I,I-1)=-BETA*A(I)
	  dTAU_dCHI(I,I)= ALPHA + BETA*(A(K)-B(I))
	  IF(K .LE. DEND)dTAU_dCHI(I,K)= ALPHA + BETA*(B(K)-C(I))
	  IF(K .LT. MIN(ND,DEND+1)dTAU_dCHI(I,K+1)=BETA*C(K)
	END DO
!
! We use UIJ for the sum of RHS terms containing 1/DTAU(I-1)
! We use UII for the sum of RHS terms containing 1/DTAU(I)
!
	W(:,:)=0.0_LDP
	DO I=MAX(2,DST-1),MIN(DEND+1,ND-1)
	  K=I+1
	  J=I-1
!
	  UIJ=-F(J)*RJ(J)/DTAU(J)/DTAU(J)
	  T1=-0.5_LDP*(1.0_LDP-COH_VEC(I))+F(I)/DTAU(J)/DTAU(J)
	  UIJ=UIJ+T1*RJ(I)+0.5_LDP*SOURCE(I)
	  UII=-F(K)*RJ(K)/DTAU(I)/DTAU(I)
	  T1=-0.5_LDP*(1.0_LDP-COH_VEC(I))+F(I)/DTAU(I)/DTAU(I)
	  UII=UII+T1*RJ(I)+0.5_LDP*SOURCE(I)
!
	  IF(J .GT. DST)THEN
	    IF(J .NE. 1)W(I,J-1)=W(I,J-1)+UIJ*dTAU_dCHI(J,J-1)
	    W(I,J)=W(I,J)+UIJ*dTAU_dCHI(J,J)+UII*dTAU_dCHI(I,J)
	  END IF
	  W(I,I)=W(I,I)+UIJ*dTAU_dCHI(J,I)+UII*dTAU_dCHI(I,I)
	  W(I,I)=W(I,I)-0.5_LDP*(DTAU(J)+DTAU(I))*
	1                       (SOURCE(I)+COH_VEC(I)*RJ(I))/CHI(I)
	  IF(K .LE. DEND)W(I,K)=W(I,K)+UII*dTAU_dCHI(I,K)
!
	END DO
!
! Outer boundary condition.
!
	IF(DST .LE. 3)THEN
	  UII=-0.5_LDP*(1.0_LDP-COH_VEC(1))*RJ(1)+(RJ(1)*F(1)-RJ(2)*F(2))/DTAU(1)/DTAU(1)
	  UII=UII+0.5_LDP*SOURCE(1)
	  IF(DST .EQ. 1)THEN
	    W(1,1)=W(1,1)+UII*dTAU_dCHI(1,1)
	    W(1,1)=W(1,1)-0.50_LDP*DTAU(1)*(RJ(1)*COH_VEC(1)+SOURCE(1))/CHI(1)
	  END IF
	  IF((DEND-2)*(2-DST) .GE. 0)W(1,2)=W(1,2)+UII*dTAU_dCHI(1,2)
	  IF((DEND-3)*(3-DST) .GE. 0)W(1,3)=W(1,3)+UII*dTAU_dCHI(1,3)
	END IF
!
! Inner boundary condition.
!
	IF(DEND .GE. NE-2)THEN
	  UIJ=-0.5_LDP*(1.0_LDP-COH_VEC(ND))*RJ(ND)+(RJ(ND)*F(ND)-F(ND-1)*RJ(ND-1))/DTAU(ND-1)/DTAU(ND-1)
	  UIJ=UIJ+0.5_LDP*SOURCE(ND)
	  IF((DEND-ND+2)*(ND-2-ST) .GE. 0)W(ND,ND-2)=W(ND,ND-2)+UIJ*dTAU_dCHI(ND-1,ND-2)
	  IF((DEND-ND+1)*(ND-1-ST) .GE. 0)W(ND,ND-1)=W(ND,ND-1)+UIJ*dTAU_dCHI(ND-1,ND-1)
	  IF(DEND .EQ. ND)THEN
	    W(ND,ND)=W(ND,ND)+UIJ*dTAU_dCHI(ND-1,ND)
	    W(ND,ND)=W(ND,ND)-0.50_LDP*DTAU(ND-1)*(RJ(ND)*COH_VEC(ND)+SOURCE(ND))/CHI(ND)
	    IF(DIFF)THEN
	      W(ND,ND)=W(ND,ND)-DBB/CHI(ND)/CHI(ND)/3.0_LDP
	    END IF
	  END IF
	END IF
!
	RETURN
	END
