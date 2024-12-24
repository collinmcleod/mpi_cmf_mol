!
! Routine to compute the opacity & emissivity variation matrices for
! the case with lines. Also computes the matrix dRHS_dCHI which
! multiply's %KI in the V equation. X is the line profile. It is assumed
! that :-
!				KI( , ,1)=dCHI
!				KI( , ,2)=dETA
!
!  				dRHS_dCHI( , ,)=dCHI
!
	SUBROUTINE EDD_JREL_VAR_MPI_V1(
	1                  R,SIGMA,CHI,ESEC,FEDD,dIBCHI_A,dIBCHI_B,
	1                  RHS_JNU,IB_STAB_FACTOR,DBB,
	1                  INNER_BND_METH,METHOD,
	1                  DST,DEND,ND,NM)
	USE SET_KIND_MODULE
	USE MOD_VAR_JREL_MPI_V1
	USE MOD_VAR_OPAC_J, ONLY : KI, RHS_dHdCHI, VDST, VDEND
	USE MPI
	IMPLICIT NONE
!
! Created 13-Dec-2024 : Based on edd_jrel_var_v3.f 
!
	INTEGER DST,DEND,ND,NM
	REAL(KIND=LDP) ESEC(ND)
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) SIGMA(ND)
	REAL(KIND=LDP) FEDD(ND)
	REAL(KIND=LDP) DBB
	REAL(KIND=LDP) RHS_JNU
	REAL(KIND=LDP) IB_STAB_FACTOR
	REAL(KIND=LDP) dIBCHI_A,dIBCHI_B
	CHARACTER(LEN=*) METHOD
	CHARACTER(LEN=*) INNER_BND_METH
!
! Local variables.
!
        REAL(KIND=LDP) WRK_RHS(ND-1,ND)
        REAL(KIND=LDP) WRK_MAT(ND,ND)
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	REAL(KIND=LDP) dUdCHI
	REAL(KIND=LDP) dTAdCHI_J,dTAdCHI_I
	REAL(KIND=LDP) dTCdCHI_I,dTCdCHI_K
	REAL(KIND=LDP) dTBdCHI_J,dTBdCHI_I,dTBdCHI_K,dTBdCHI
	REAL(KIND=LDP) dXM_EPS_J,dXM_EPS_I,dXM_EPS_K
	REAL(KIND=LDP) T1,T2
	INTEGER MDST,MDEND,IERR
	INTEGER I,J,K,L
!
! 
!
	KI=0.0_LDP;  RHS_dHdCHI=0.0_LDP
	MDST=MAX(1,DST-1)
	MDEND=MIN(DEND+1,ND)
!
! NB: The dTAUdCHI_J & dTAUdCHI_H matrices were computed in the calling routine..
!
! The following derivatives are valid for all ML.
!
	DO I=MDST,MIN(MDEND,ND-1)
	  T1=(P_H(I)+W(I))*(CHI_H(I)+CHI_H(I+1))
	  dHUdCHI(I)=HU(I)*W(I)/T1
	  dHUdTAU(I)=-HU(I)/DTAU_H(I)
	  dHLdCHI(I)=HL(I)*W(I)/T1
	  dHLdTAU(I)=-HL(I)/DTAU_H(I)
	  dHSdCHI(I)=-P_H(I)*HS(I)/T1
	  EPS_FAC(I)=-1.0_LDP/T1
	END DO
!
! 
!
! Firstly we compute the variation of the elements with respect to
! DTAU. We then multiply by dTAUdCHI matrix.
!
! DTAU_H terms
!
	DO I=MDST,MIN(MDEND,ND-1)
	  J=I-1
	  K=I+1
	  dTAdCHI_J=-dHLdTAU(J)
	  dTCdCHI_I=-dHUdTAU(I)
	  dTBdCHI_I=dHLdTAU(I)
	  dTBdCHI_J=dHUdTAU(J)
	  dRHSdJ(I)=  - dTAdCHI_J*JNU(J) - dTBdCHI_J*JNU(I)
	  dRHSdI(I)=  - dTCdCHI_I*JNU(K) - dTBdCHI_I*JNU(I)
	END DO
!
	DO I=MAX(2,DST),MIN(DEND,ND-1)
	  J=I-1
	  DO L=1,ND
	    WRK_MAT(I,L)=WRK_MAT(I,L)+dRHSdJ(I)*dTAUdCHI_H(J,L)
	    WRK_MAT(I,L)=WRK_MAT(I,L)+dRHSdI(I)*dTAUdCHI_H(I,L)
	  END DO
	END DO
!
! DTAU_J terms
!
	DO I=MAX(2,DST),MIN(DEND,ND-1)
	  J=I-1
	  T1=0.5_LDP*R(I)*R(I)/Q(I)
	  dTBdCHI=PSI(I)/(DTAU_J(J)+DTAU_J(I)) +
	1          T1*(P_J(I)-COH_VEC(I))*GAM_REL(I)
!
! dDELUB is use as correction because UB(I)=-TB(I)-PSI(I)-PSIPREV(I)
!
	  dUdCHI=PSIPREV(I)/(DTAU_J(J)+DTAU_J(I))
	  dRHSdJ(I)= T1*SOURCE(I) - dTBdCHI*JNU(I)  + dUdCHI*JNU_PREV(I)
	  dRHSdI(I)= T1*SOURCE(I) + dUdCHI*JNU_PREV(I) - dTBdCHI*JNU(I)
!
	END DO
!
	DO I=MAX(2,DST),MIN(DEND,ND-1)
	  J=I-1
	  DO L=1,ND
	    WRK_MAT(I,L)=WRK_MAT(I,L)+dRHSdJ(I)*dTAUdCHI_J(J,L)
	    WRK_MAT(I,L)=WRK_MAT(I,L)+dRHSdI(I)*dTAUdCHI_J(I,L)
	  END DO
	END DO
!
! Can now update KI for direct opacity variation due to CHI_H.
!
	DO I=MAX(2,DST),MIN(DEND,ND-1)
	  J=I-1
	  K=I+1
!
	  dTAdCHI_J=-dHLdCHI(J)
	  dTAdCHI_I=-dHLdCHI(J)
	  dTCdCHI_I=-dHUdCHI(I)
	  dTCdCHI_K=-dHUdCHI(I)
!
	  dTBdCHI_J=dHUdCHI(J)
	  dTBdCHI_I=dHLdCHI(I)+dHUdCHI(J)
	  dTBdCHI_K=dHLdCHI(I)
!
	  dXM_EPS_J=( EPS_A(J)*JNU(J)-EPS_PREV_A(J)*JNU_PREV(J) +
	1             EPS_B(J)*JNU(I)-EPS_PREV_B(J)*JNU_PREV(I) )*EPS_FAC(J)
	  dXM_EPS_K=( EPS_PREV_A(I)*JNU_PREV(I)-EPS_A(I)*JNU(I) +
	1             EPS_PREV_B(I)*JNU_PREV(I+1)-EPS_B(I)*JNU(I+1) )*EPS_FAC(I)
	  dXM_EPS_I=dXM_EPS_J+dXM_EPS_K
!
! NB  :  VB(I)=-HS(J) and VC(I)=HS(I)
!
	  WRK_MAT(I,J)=WRK_MAT(I,J)
	1             - (dTAdCHI_J*JNU(J)+dTBdCHI_J*JNU(I))
	1             - DERIV_SCL_FAC(I)*dHSDCHI(J)*GAM_RSQHNU_PREV(J)
	1             + dXM_EPS_J
!
	  WRK_MAT(I,K)=WRK_MAT(I,K)
	1             - (dTCdCHI_K*JNU(K)+dTBdCHI_K*JNU(I))
	1             + DERIV_SCL_FAC(I)*dHSDCHI(I)*GAM_RSQHNU_PREV(I)
	1             + dXM_EPS_K
!
	  WRK_MAT(I,I)=WRK_MAT(I,I)
	1             - (dTAdCHI_I*JNU(J)+dTCdCHI_I*JNU(K)+dTBdCHI_I*JNU(I))
	1             + DERIV_SCL_FAC(I)*(dHSDCHI(I)*GAM_RSQHNU_PREV(I)-dHSDCHI(J)*GAM_RSQHNU_PREV(J))
	1             + dXM_EPS_I
!
	END DO
!
! Can now update KI for direct opacity variation due to CHI_J.
!
	DO I=MAX(2,DST),MIN(DEND,ND-1)
	  J=I-1
	  T1=0.5_LDP*R(I)*R(I)/Q(I)
!
	  dTBdCHI= -PSI(I)/CHI_J(I)+GAM_RSQ_DTAUONQ(I)*(COH_VEC(I)-VdJdR_TERM(I))/CHI_J(I)
	  dUdCHI=-PSIPREV(I)/CHI_J(I)
!
	  T1=T1*(DTAU_J(J)+DTAU_J(I))/CHI_J(I)
	  WRK_MAT(I,I)=WRK_MAT(I,I)
	1             + (dUdCHI*JNU_PREV(I)-dTBdCHI*JNU(I))
	1             - T1*SOURCE(I)
!
	  KI(I,I,2)=T1
	END DO
!
! Now do the boundary conditions.
!
	IF(DST .EQ. 1)THEN
	  T1=  ( (FEDD(1)+VdHdR_TERM(1))*Q(1)*JNU(1)*GAM_RSQ(1) -
	1      (FEDD(2)+VdHdR_TERM(2))*Q(2)*JNU(2)*GAM_RSQ(2) )/DTAU_H(1)/DTAU_H(1)
	  DO L=1,ND
	    KI(1,L,1)=KI(1,L,1)+T1*dTAUdCHI_H(1,L)
	  END DO
	  KI(1,1,1)=KI(1,1,1)+ (PSI(1)*JNU(1)- PSIPREV(1)*JNU_PREV(1))/CHI_H(1)
	END IF
!!
! NB: We multiply DBB by GAM_REL, as we will divid KI by GAM_REL later.
! The DBB derivative is with  respect to CHI, where as the other derivatives
! were donie with respect to CHI_J or CHI_H.
!
	IF(DEND .EQ. ND)THEN
	  IF(INNER_BND_METH .EQ. 'JEQB')THEN
	
	  ELSE IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	    T1= ( GAM_RSQ(ND)*(FEDD(ND)+VdHdR_TERM(ND))*JNU(ND) -
	1          GAM_RSQ(ND-1)*(FEDD(ND-1)+VdHdR_TERM(ND-1))*Q(ND-1)*JNU(ND-1) )
	1           / DTAU_H(ND-1)/DTAU_H(ND-1)
	    DO L=1,ND
	      WRK_MAT(ND,L)=WRK_MAT(ND,L)+T1*dTAUdCHI_H(ND-1,L)
	    END  DO
	    WRK_MAT(ND,ND)=WRK_MAT(ND,ND)-GAM_RSQ(ND)*GAM_REL(ND)*DBB/3.0_LDP/CHI(ND)/CHI(ND)
!
	  ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	    T1 = GAM_RSQ(ND-1)*(FEDD(ND-1)+VdHdR_TERM(ND-1))*Q(ND-1)
	    T2 = GAM_RSQ(ND)*(FEDD(ND)+VdHdR_TERM(ND))
	    T2 = (T2*JNU(ND)-T1*JNU(ND-1)) + IB_STAB_FACTOR*T2*(JNU(ND)-RHS_JNU)
	    T2 = T2/ DTAU_H(ND-1)/DTAU_H(ND-1)
	    DO L=1,ND
	      WRK_MAT(ND,L)=WRK_MAT(ND,L)+T2*dTAUdCHI_H(ND-1,L)
	    END DO
!
	  ELSE IF(INNER_BND_METH .EQ. 'HOLLOW')THEN
	      WRK_MAT(ND,ND)=WRK_MAT(ND,ND)+dIBCHI_A
	      WRK_MAT(ND,ND-1)=WRK_MAT(ND,ND-1)+dIBCHI_B
!
	  ELSE
	    T1= ( GAM_RSQ(ND)*(FEDD(ND)+VdHdR_TERM(ND))*JNU(ND) -
	1           GAM_RSQ(ND-1)*(FEDD(ND-1)+VdHdR_TERM(ND-1))*Q(ND-1)*JNU(ND-1) )
	1           / DTAU_H(ND-1)/DTAU_H(ND-1)
	    DO L=1,ND
	      KI(ND,L,1)=KI(ND,L,1)+T1*dTAUdCHI_H(ND-1,L)
	    END DO
	  END IF
	END IF
!
! Divide KI(:,:,1) by GAM_REL since dCHI_J/dCHI=dCHI_H/dCHI=1/GAM_REL
!
	DO J=1,ND
	  DO I=DST,DEND
	    WRK_MAT(I,J)=WRK_MAT(I,J)/GAM_REL(J)
	  END DO
	END DO
!
	I=ND*ND
	CALL MPI_ALLREDUCE(MPI_IN_PLACE,WRK_MAT,I,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	DO J=VDST,VDEND
	  KI(:,J,1)=WRK_MAT(:,J)
	END DO
!
! 
!
! Now we compute the variation of the equation which updates the
! flux variation. We don't correct for the line profile, since
! we would the require two matrices. Note that HU(I), HL(I) and
! HS(I) depend directly on CHI(I) and CHI(I+1).
!
	DO I=DST,MIN(DEND,ND-1)
	  T1=dHUdTAU(I)*JNU(I+1)-dHLdTAU(I)*JNU(I)
	  DO L=1,ND
	    WRK_RHS(I,L)=WRK_RHS(I,L)+T1*dTAUdCHI_H(I,L)
	  END DO
	  T1=dHUdCHI(I)*JNU(I+1) - dHLdCHI(I)*JNU(I)
	1                       + dHSdCHI(I)*GAM_RSQHNU_PREV(I) +
	1     EPS_FAC(I)*( EPS_PREV_A(I)*JNU_PREV(I)-EPS_A(I)*JNU(I) +
	1                  EPS_PREV_B(I)*JNU_PREV(I+1)-EPS_B(I)*JNU(I+1) )
	  WRK_RHS(I,I)=WRK_RHS(I,I) + T1
	  WRK_RHS(I,I+1)=WRK_RHS(I,I+1) + T1
	END DO
!
! Recall dCHI_H/dCHI=1/GAM_REL
!
	DO L=VDST,VDEND
	  DO I=DST,MIN(DEND,ND-1)
	    RHS_dHdCHI(I,L)=RHS_dHdCHI(I,L)/GAM_REL(L)
	  END DO
	END DO
!
	I=(ND-1)*ND
	CALL MPI_ALLREDUCE(MPI_IN_PLACE,WRK_RHS,I,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	DO J=VDST,VDEND
	  RHS_dHdCHI(:,J)=WRK_RHS(:,J)
	END DO
!
	RETURN
	END
