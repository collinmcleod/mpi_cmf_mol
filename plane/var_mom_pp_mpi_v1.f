!
! Subroutine to solve the plane-parallel moment equations
! for the variation of J with CHI and ETA. This routine
! should be called in conjunction with: fcomp_pp.f & mom_j_pp_v1.f
!
	SUBROUTINE VAR_MOM_PP_MPI_V1(R,ETA,CHI,ESEC,F,
	1               TX,dJ_DIFF_dT,dJ_DIFF_ddTdR,DO_THIS_MATRIX,
	1               HBC_J,HBC_S,IN_HBC,
	1               DIFF,DBB,dDBBdT,dTdR,
	1               IC,METHOD,COHERENT,DST,DEND,ND,NM)
	USE SET_KIND_MODULE
	USE MPI
	IMPLICIT NONE
!
! Created 4-Jan-2005
!
	INTEGER ND
	INTEGER DST,DEND
	INTEGER NM
!
! Atmospheric variables: these are supplied.
!
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) ETA(ND)
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) ESEC(ND)
	REAL(KIND=LDP) F(ND)		!F (=J/K) must be supplied.
!
! Radiation field variables.These are computed.
!
	REAL(KIND=LDP) TX(ND,DST:DEND,NM)
	REAL(KIND=LDP) dJ_DIFF_dT(ND)
	REAL(KIND=LDP) dJ_DIFF_ddTdR(ND)
	LOGICAL DO_THIS_MATRIX(NM)
!
! Boundary conditions: required input.
!
	REAL(KIND=LDP) HBC_J		!H = HBC_J*J(1) + HBC_J*S(1)
	REAL(KIND=LDP) HBC_S
	REAL(KIND=LDP) IN_HBC		!H = 0.25Ic - hJ -hIc/2 (inner boundary)
!
	REAL(KIND=LDP) DBB		!|dB/dR| -- for diffusion approximation.
	REAL(KIND=LDP) dDBBdT
	REAL(KIND=LDP) dTdR
	REAL(KIND=LDP) IC		!Intensity at inner boundary (used if DIFF=.FALSE.)
!
	LOGICAL DIFF		!Use diffusion approximation?
	LOGICAL COHERENT	!Assume coherent scattering?
	CHARACTER*6 METHOD	!Method to compute dCHI/dR
!
! Local vectors.
!
	REAL(KIND=LDP) dJ_dDBB(ND)
	REAL(KIND=LDP) JNU(ND)
!
	REAL(KIND=LDP), TARGET  :: TRI_VECS(ND,4)
	REAL(KIND=LDP), POINTER :: TA(:)
	REAL(KIND=LDP), POINTER :: DD(:)
	REAL(KIND=LDP), POINTER :: TC(:)
	REAL(KIND=LDP), POINTER ::  RHS(:)
!
	REAL(KIND=LDP) DTAU(ND)
	REAL(KIND=LDP) dCHIdR(ND)
	REAL(KIND=LDP) SOURCE(ND)
	REAL(KIND=LDP) COH_VEC(ND)
	REAL(KIND=LDP) TOR,E2TOR,EXPN
	REAL(KIND=LDP) dR
	EXTERNAL EXPN
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER I,LUER,ERROR_LU
	INTEGER M_DST,M_DEND
	INTEGER IERR
	EXTERNAL ERROR_LU
!
	M_DST=MAX(1,DST-1)
        M_DEND=MIN(DEND+1,ND)
!
	TA(1:ND)=>TRI_VECS(1:ND,1)
	DD(1:ND)=>TRI_VECS(1:ND,2)
	TC(1:ND)=>TRI_VECS(1:ND,3)
	RHS(1:ND)=>TRI_VECS(1:ND,4)
	TRI_VECS=0.0_LDP
!
! As no coupling with earlier frequencies
!
	CALL TUNE(1,'ZER_TX')
	DO I=1,NM
	  IF(DO_THIS_MATRIX(I))TX(:,:,I)=0.0_LDP
        END DO
	CALL TUNE(2,'ZER_TX')
!
! Compute optical depth scale.
!
	dCHIdR=0.0_LDP      					!(to be checked if needed).
	CALL DERIVCHI_MPI_V1(dCHIdR,CHI,R,M_DST,M_DEND,ND,METHOD)
	CALL d_DERIVCHI_dCHI_MPI_V1(dCHIdR,CHI,R,M_DST,M_DEND,ND,METHOD)
        DO I=M_DST,MIN(ND-1,DEND)
          dR=R(I)-R(I+1)
          DTAU(I)=0.5_LDP*dR*(CHI(I)+CHI(I+1)+dR*(dCHIdR(I+1)-dCHIdR(I))/6.0_LDP)
        END DO
!
! If the scattering is incoherent, ETA must contain the full
! emissivity.
!
	DO I=DST,DEND
	  SOURCE(I)=ETA(I)/CHI(I)
	END DO
	IF(COHERENT)THEN
	  DO I=DST,DEND
	    COH_VEC(I)=ESEC(I)/CHI(I)
	  END DO
	ELSE
	  DO I=DST,DEND
	    COH_VEC(I)=0.0_LDP
	  END DO
	END IF
!
! Compute the TRIDIAGONAL operators, and the RHS source vector.
! Note that the diagonal component, TB, is given by
!
!     TB= -DD(I)-TA(I)-TC(I)
!
	DO I=MAX(2,DST),MIN(DEND,ND-1)
	  TA(I)=-F(I-1)/DTAU(I-1)
	  TC(I)=-F(I+1)/DTAU(I)
	  DD(I)=-0.5_LDP*(DTAU(I-1)+DTAU(I))*(1.0_LDP-COH_VEC(I)) -
	1          (F(I)-F(I-1))/DTAU(I-1) - (F(I)-F(I+1))/DTAU(I)
	  RHS(I)=0.5_LDP*(DTAU(I-1)+DTAU(I))*SOURCE(I)
	END DO
!
! Evaluate TA,TB,TC for boundary conditions. These are second order.
!
	IF(DST .EQ. 1)THEN
	  TC(1)=-F(2)/DTAU(1)
	  DD(1)= -0.5_LDP*DTAU(1)*(1.0_LDP-COH_VEC(1))-HBC_J-(F(1)-F(2))/DTAU(1)
	1             +HBC_S*COH_VEC(1)
	  RHS(1)=0.5_LDP*DTAU(1)*SOURCE(1)+HBC_S*SOURCE(1)
	  TA(1)=0.0_LDP
	END IF
!
	IF(DEND .EQ. ND)THEN
	  TA(ND)=-F(ND-1)/DTAU(ND-1)
	  IF(DIFF)THEN
	    DD(ND)=-(F(ND)-F(ND-1))/DTAU(ND-1)-0.5_LDP*DTAU(ND-1)*(1.0_LDP-COH_VEC(ND))
	    RHS(ND)=DBB/3.0_LDP/CHI(ND)+0.5_LDP*DTAU(ND-1)*SOURCE(ND)
	  ELSE
	    DD(ND)=-(F(ND)-F(ND-1))/DTAU(ND-1)-0.5_LDP*DTAU(ND-1)*(1.0_LDP-COH_VEC(ND))-IN_HBC(1)
	    RHS(ND)=IC*(0.25_LDP+0.5_LDP*IN_HBC(1))+0.5_LDP*DTAU(ND-1)*SOURCE(ND)
	  END IF
	  TC(ND)=0.0_LDP
	END IF
!
! Solve for the radiation field along ray for this frequency.
!
	I=4*ND
	CALL MPI_ALLREDUCE(MPI_IN_PLACE,TRI_VECS,I,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	CALL THOMAS_RH(TA,DD,TC,RHS,ND,1)
!
! Save J
!
	JNU(1:ND)=RHS(1:ND)
!
! NB: As there in no velocity field, there is no coupling with
! earlier frequencies.
!
! After this call:
!         TX(I,J,1)= dRHS(I)/dCHI(J)
!
! Evaluate optical depth at the outer boundary. This expression should be the same as in
! procedure that evaluates the Eddington factors (formal solution). After the call to
! dRHSdCII_PP, we update TX(1,1,1) for the dependence of the puer boundary condition
! onthe source function.
!
	TOR=CHI(1)*(R(1)-R(3))/LOG(CHI(3)/CHI(1))
        IF(TOR .LT. 0.0_LDP .OR. HBC_S .LE. 0.0_LDP)TOR=0.0_LDP
	E2TOR=TOR*EXPN(ITWO,TOR)
!
	CALL TUNE(1,'SOL_TX')
 	CALL dRHSdCHI_PP_MPI_V1(TX(:,:,1),SOURCE,CHI,DTAU,COH_VEC,JNU,F,R,DIFF,DBB,DST,DEND,ND)
	IF(DST .EQ. 1)TX(1,1,1)=TX(1,1,1)-(SOURCE(1)+COH_VEC(1)*JNU(1))*(HBC_S-E2TOR)/CHI(1)
	I=DEND-DST+1
	CALL SIMPTH_RH(TA,DD,TC,TX(:,:,1),ND,I)
	CALL TUNE(2,'SOL_TX')
!
! Compute dJ/dETA matrix (=TX(:,:,2); variation of eta).
!
	CALL TUNE(1,'ETA_TX')
        TX(:,:,2)=0.0_LDP
        IF(DST .EQ. 1)TX(1,1,2)=0.5_LDP*DTAU(1)/CHI(1)+HBC_S/CHI(1)
        DO I=MAX(2,DST),MIN(DEND,ND-1)
          TX(I,I,2)=0.5_LDP*(DTAU(I-1)+DTAU(I))/CHI(I)
        END DO
        IF(DEND .EQ. ND)TX(ND,ND,2)=0.5_LDP*DTAU(ND-1)/CHI(ND)		!Diff and non diff
	I=DEND-DST+1
	CALL SIMPTH_RH(TA,DD,TC,TX(:,:,2),ND,I)
	CALL TUNE(2,'ETA_TX')
!
! Allow for the variation of the diffusion approximation boundary condition.
! because their is no frequency coupling, we only need one vector. Both
! dJ__DIF_ddTdR and dJ_DIF_dT can be computed from this vector.
!
	dJ_dDBB(1:ND)=0.0_LDP
	IF(DIFF)THEN
	  dJ_dDBB(ND)=1.0_LDP/3.0_LDP/CHI(ND)
	  CALL SIMPTH_RH(TA,DD,TC,dJ_dDBB,ND,IONE)
	  dJ_DIFF_dT(1:ND)=dJ_dDBB(1:ND)*dDBBdT
	  dJ_DIFF_ddTDR(1:ND)=dJ_dDBB(1:ND)*DBB/dTdR
	END IF
!
	RETURN
	END
