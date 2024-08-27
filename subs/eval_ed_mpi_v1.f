C
C Routine to compute the contibution to the electron density by
C a set of ionization species. U and PHI are non-LTE partition
C functions as computed by PAR_FUN. Also returned are dCHARGEdNE
C (so that a Newton-Rapshon technique can be used to solve for Ne
C as equations are non-linear) and the POPULATION of the highest
C ionization species.
C
C Created 29-Aug-1990.
C
	SUBROUTINE EVAL_ED_MPI_V1(CHARGE,dCHARGEdNE,
	1                  U, PHI, DION_MAT, Z, ED, POPOXY,
	1                  LOC_CH, dDIONdNe, dLOC_CH_dNe,
	1                  SPEC_ID, NSPEC, ION_ID_ST, ION_ID_END, NION,
	1                  DST, DEND, ND, FIRST)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
	LOGICAL FIRST
	INTEGER ION_ID_ST,ION_ID_END,NION
	INTEGER ND
	INTEGER SPEC_ID, NSPEC
	INTEGER DST,DEND
	REAL(KIND=LDP) CHARGE(ND),dCHARGEdNE(ND)

	REAL(KIND=LDP) U(DST:DEND,NION)
	REAL(KIND=LDP) PHI(DST:DEND,NION)
	REAL(KIND=LDP) DION_MAT(DST:DEND,NSPEC)
	REAL(KIND=LDP) Z(NSPEC),ED(ND),POPOXY(ND)
	REAL(KIND=LDP) LOC_CH(ND),dDIONdNe(ND),dLOC_CH_dNe(ND)
C
	REAL(KIND=LDP) T1
	INTEGER I,J,K
C
	IF(FIRST)THEN
	  DO I=DST,DEND
	    CHARGE(I)=0.0_LDP
	    dCHARGEdNe(I)=0.0_LDP
	  END DO
	  FIRST=.FALSE.
	END IF
C
	DO I=DST,DEND
	  DION_MAT(I,SPEC_ID)=0.0_LDP
	  dDIONdNe(I)=0.0_LDP
	  LOC_CH(I)=0.0_LDP
	  dLOC_CH_dNe(I)=0.0_LDP
	END DO
!
	K=ION_ID_END-ION_ID_ST
	DO I=DST,DEND
	  DO J=ION_ID_ST,ION_ID_END-1
	    T1=PHI(I,J)*( U(I,J)/U(I,J+1) )
            DION_MAT(I,SPEC_ID)=T1*ED(I)*( 1.0_LDP+DION_MAT(I,SPEC_ID) )
	    dDIONdNe(I)=T1*( K  + ED(I)*dDIONdNe(I) )
	    LOC_CH(I)=T1*ED(I)*( Z(J)+LOC_CH(I) )
	    dLOC_CH_dNe(I)=T1*( Z(J)*K + ED(I)*dLOC_CH_dNe(I) )
	  END DO
	  DION_MAT(I,SPEC_ID)=DION_MAT(I,SPEC_ID)+1.0_LDP
	  LOC_CH(I)=LOC_CH(I)+Z(ION_ID_END)
	END DO
C
	DO I=DST,DEND
	  T1=POPOXY(I)/DION_MAT(I,SPEC_ID)
	  CHARGE(I)=CHARGE(I)+T1*LOC_CH(I)
	  dCHARGEdNE(I)=dCHARGEdNE(I) + T1*(
	1                  dLOC_CH_dNe(I) - dDIONdNe(I)*(LOC_CH(I)/DION_MAT(I,SPEC_ID)) )
	  DION_MAT(I,SPEC_ID)=T1
	END DO
C
	RETURN
	END
