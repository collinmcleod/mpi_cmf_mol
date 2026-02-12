!
! Subroutine to compute the value of the statistical equilibrium
! equations and the variation of the statistical equilibrium matrix for
! terms which are radiation field independent.
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
!
! The COLLISION routine that is called has a special FORM, which is distinct
! from that in STEQGEN_V2.
!
! NB - ZION is the charge on the ion - thus ZHYD=1.0D0
!
! Routine also increments the ionization equilibrium equations. Routine no
! longer works for NUM_BNDS=ND.
!
! At present only collisional ionizations to ground state are considered.
!
! NION is the the first dimension of STEQ[ION]. In general we
! NION would be the total number of ionic species.
!
	SUBROUTINE STEQ_MULTI_MPI_V1(ED,T,
	1       HN_S,HNST_S,dlnHNST_S_dlnT,AVE_ENERGY,N_S,DI_S,
	1       HN_F,HNST_F_ON_S,W_F,A_F,FEDGE_F,G_F,LEVNAME_F,N_F,
	1       F_TO_S_MAPPING,POP,NEXT_PRES,ZION,
	1       ID,COL_FILE,OMEGA_GEN,
	1       EQGS,NUM_BNDS,DST,DEND,ND,NION,NT,
	1       COMPUTE_BA,FIXED_T,LAST_ITERATION,IS_A_MOLECULE)
	USE SET_KIND_MODULE
	USE MPI
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
	REAL(KIND=LDP) PLANCKS_CONSTANT
	EXTERNAL PLANCKS_CONSTANT
	EXTERNAL OMEGA_GEN
!
	INTEGER NT
	INTEGER NUM_BNDS
	INTEGER ND
	INTEGER NION
	INTEGER EQGS
	INTEGER DST, DEND
!
	INTEGER ID
	INTEGER N_S,N_F
!
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) ED(ND)
!
	REAL(KIND=LDP) DI_S(DST:DEND)
	REAL(KIND=LDP) HN_S(N_S,DST:DEND)
	REAL(KIND=LDP) HNST_S(N_S,DST:DEND)
	REAL(KIND=LDP) dlnHNST_S_dlnT(N_S,DST:DEND)
	REAL(KIND=LDP) AVE_ENERGY(N_S)
!
	REAL(KIND=LDP) HN_F(N_F,DST:DEND)
	REAL(KIND=LDP) HNST_F_ON_S(N_F,DST:DEND)
	REAL(KIND=LDP) W_F(N_F,DST:DEND)
	REAL(KIND=LDP) A_F(N_F,N_F)
	REAL(KIND=LDP) FEDGE_F(N_F)
	REAL(KIND=LDP) G_F(N_F)
	CHARACTER*(*) LEVNAME_F(N_F),COL_FILE
	INTEGER F_TO_S_MAPPING(N_F)
	REAL(KIND=LDP) ZION
!
	REAL(KIND=LDP) POP(ND)		!Population of species.
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	LOGICAL NEXT_PRES
	LOGICAL COMPUTE_BA,FIXED_T,LAST_ITERATION,IS_A_MOLECULE
!
! Local variables.
!
	INTEGER EQION,IONE
	INTEGER I,J,K,L,M
	INTEGER DPTH_INDX
	INTEGER UP,JUP
	INTEGER JJ,KK,LION,JL
	INTEGER VION,VED,VT
	INTEGER EQ_ION_BAL
	INTEGER EQ_NUM_CONV
	INTEGER IERR
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) TMP_VEC_ED(1)
	REAL(KIND=LDP) TMP_VEC_COOL(1)
	PARAMETER (IONE=1)
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
! Work arrays.  DCNM refers to dCNMdT
!
	REAL(KIND=LDP) OMEGA_F(N_F,N_F)
	REAL(KIND=LDP) dln_OMEGA_F_dlnT(N_F,N_F)
	REAL(KIND=LDP) CNM(N_S,N_S)
	REAL(KIND=LDP) DCNM(N_S,N_S)
!
	INTEGER, SAVE :: NUM_DEPTHS_PER_THREAD
!
!	INCLUDE 'mpif.h'
! 
!
	NUM_DEPTHS_PER_THREAD=(ND-1)/NTHREAD+1
        EQION=N_S+1			!Ion equation
        VION=N_S+1
        VED=SE(ID)%N_IV-1
        VT=SE(ID)%N_IV
	EQ_NUM_CONV=SE(ID)%NUMBER_BAL_EQ
	M=(NUM_BNDS/2)+1
!
	DO DPTH_INDX=DST,DEND
!
! Compute collisional cross-sections (and their T derivatives)
! We call this routine ND times so the CNM and DCM arrays can be
! smaller (i.e. no ND dimension).
!
! OMEGA_F,dln_OMEGA_dlnT are work arrays only.
! T1 is returned with the total cooling rate. Not used in this routine.
! We use arrays (even though of length 1) so that some F90 compilers
! don't give an error message because a scaler is passed a vector.
!
	  TMP_VEC_ED(1)=1.0_LDP		!Electron density
	  TMP_VEC_COOL(1)=0.0_LDP		!Initialize cooling rate even
!                                                        though not used here.
          I=DPTH_INDX
	  CALL SUBCOL_MULTI_V7(
	1         OMEGA_F,dln_OMEGA_F_dlnT,
	1         CNM,DCNM,
	1         HN_S(1,I),HNST_S(1,I),dlnHNST_S_dlnT(1,I),N_S,
	1         HN_F(1,I),HNST_F_ON_S(1,I),W_F(1,I),FEDGE_F,
	1         A_F,G_F,LEVNAME_F,N_F,
	1         ZION,ID,COL_FILE,OMEGA_GEN,
	1         F_TO_S_MAPPING,TMP_VEC_COOL,T(I),TMP_VEC_ED,IONE,
	1         COMPUTE_BA,FIXED_T,LAST_ITERATION,IS_A_MOLECULE)
!
!	WRITE(6,*)MYPE,DST,DEND
!	WRITE(6,*)LEVNAME_F(1),LEVNAME_F(N_F)
!	FLUSH(UNIT=6)
!	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!	CALL MPI_FINALIZE()
!	STOP
!
!
! 
!
! *********************** Code for Electron Energy Balance **********************************
! ***********************       Heating - Cooling           **********************************
!
! We use L  for the lower level (JL  refers to the same state in the full vector).
! We use UP for the upper level (JUP refers to the same state in the full vector).
!
	  T1=1.0E+15_LDP*PLANCKS_CONSTANT()*ED(I)              !1.0D+15 due to units of NU.
	  I=DPTH_INDX
	  DO L=1,N_S
	    STEQ_T_EHB(I)=STEQ_T_EHB(I)+T1*CNM(L,L)*(HNST_S(L,I)-HN_S(L,I))*AVE_ENERGY(L)
	    DO UP=L+1,N_S
	      STEQ_T_EHB(I)=STEQ_T_EHB(I)+T1*(HN_S(UP,I)*CNM(UP,L)-HN_S(L,I)*CNM(L,UP))*(AVE_ENERGY(L)-AVE_ENERGY(UP))
	    END DO
	  END DO
!
	  IF(COMPUTE_BA)THEN
	    I=DPTH_INDX
	    LION=EQGS+N_S
	    DO L=1,N_S
	      JL=EQGS+L-1
	      BA_T_PAR_EHB(JL,I)=BA_T_PAR_EHB(JL,I)-T1*CNM(L,L)*AVE_ENERGY(L)
	      T2=T1*HNST_S(L,I)*CNM(L,L)*AVE_ENERGY(L)
	      BA_T_PAR_EHB(NT-1,I)=BA_T_PAR_EHB(NT-1,I)+2*T2/ED(I)
	      BA_T_PAR_EHB(LION,I)=BA_T_PAR_EHB(LION,I)+T2/DI_S(I)
	      BA_T_PAR_EHB(NT,I)=BA_T_PAR_EHB(NT,I)+
	1            T1*DCNM(L,L)*(HNST_S(L,I)-HN_S(L,I))*AVE_ENERGY(L) +
	1            T2*dlnHNST_S_dlnT(L,I)/T(I)
	      DO UP=L+1,N_S
	        JUP=EQGS+UP-1
	        T2=T1*(AVE_ENERGY(L)-AVE_ENERGY(UP))
	        BA_T_PAR_EHB(JL,I)=BA_T_PAR_EHB(JL,I)-T2*CNM(L,UP)
	        BA_T_PAR_EHB(JUP,I)=BA_T_PAR_EHB(JUP,I)+T2*CNM(UP,L)
	        BA_T_PAR_EHB(NT-1,I)=BA_T_PAR_EHB(NT-1,I)+T2*(CNM(UP,L)*HN_S(UP,I)-CNM(L,UP)*HN_S(L,I))/ED(I)
	        BA_T_PAR_EHB(NT,I)=BA_T_PAR_EHB(NT,I)+T2*(DCNM(UP,L)*HN_S(UP,I)-DCNM(L,UP)*HN_S(L,I))
	      END DO
	    END DO
	  END IF
!
! *********************** End Code for Electron Energy Balance **********************************
!
	  I=DPTH_INDX
	  DO J=1,N_S			!Which S.E. equation
	    T1=0.0_LDP
	    T2=0.0_LDP
	    DO L=1,N_S
	      T1=T1+( HN_S(L,I)*CNM(L,J)-HN_S(J,I)*CNM(J,L) )
	      T2=T2+( HN_S(L,I)*DCNM(L,J)-HN_S(J,I)*DCNM(J,L) )
	    END DO
	    SE(ID)%STEQ(J,I)=SE(ID)%STEQ(J,I)+(T1+(HNST_S(J,I)-HN_S(J,I))*CNM(J,J))*ED(I)
	  END DO
!
	  IF(COMPUTE_BA)THEN
	    I=DPTH_INDX
	    DO J=1,N_S			!Which S.E. equation
	      DO K=1,N_S			!Which variable
	        IF(K.EQ.J)THEN
	          T1=0.0_LDP
	          DO L=1,N_S
		    T1=T1+CNM(J,L)
	          END DO
		  SE(ID)%BA(J,K,M,I)=SE(ID)%BA(J,K,M,I)-T1*ED(I)
	        ELSE
	          SE(ID)%BA(J,K,M,I)=SE(ID)%BA(J,K,M,I)+ED(I)*CNM(K,J)
	        END IF
	      END DO
!
	      T1=0.0_LDP
	      T2=0.0_LDP
	      DO L=1,N_S
	        T1=T1+( HN_S(L,I)*CNM(L,J)-HN_S(J,I)*CNM(J,L) )
	        T2=T2+( HN_S(L,I)*DCNM(L,J)-HN_S(J,I)*DCNM(J,L) )
	      END DO
!
	      SE(ID)%BA(J,VION,M,I)=SE(ID)%BA(J,VION,M,I) + (HNST_S(J,I)*CNM(J,J))*ED(I)/DI_S(I)
	      SE(ID)%BA(J,VED,M,I) =SE(ID)%BA(J,VED,M,I) + T1+CNM(J,J)*(2*HNST_S(J,I)-HN_S(J,I))
	      SE(ID)%BA(J,VT,M,I)  =SE(ID)%BA(J,VT,M,I) +
	1                              ED(I)*( T2+(HNST_S(J,I)-HN_S(J,I))*DCNM(J,J)+
	1                              CNM(J,J)*HNST_S(J,I)*dlnHNST_S_dlnT(J,I)/T(I) )
	    END DO
	  END IF
!
! EQION is the ion equation
!
	  T1=0.0_LDP
	  T2=0.0_LDP
	  I=DPTH_INDX
	  DO J=1,N_S
	    T1=T1+(HNST_S(J,I)-HN_S(J,I))*CNM(J,J)
	    T2=T2+HNST_S(J,I)*CNM(J,J)
	  END DO
	  SE(ID)%STEQ(EQION,I)=SE(ID)%STEQ(EQION,I)-T1*ED(I)
!
	  IF(COMPUTE_BA)THEN
	    I=DPTH_INDX
	    DO J=1,N_S
	      SE(ID)%BA(EQION,J,M,I)=SE(ID)%BA(EQION,J,M,I)+CNM(J,J)*ED(I)
	      SE(ID)%BA(EQION,VT,M,I)=SE(ID)%BA(EQION,VT,M,I) -
	1          ED(I)*(  (HNST_S(J,I)-HN_S(J,I))*DCNM(J,J) +
	1          CNM(J,J)*HNST_S(J,I)*dlnHNST_S_dlnT(J,I)/T(I)  )
	    END DO
	    SE(ID)%BA(EQION,VION,M,I)=SE(ID)%BA(EQION,VION,M,I) - T2*ED(I)/DI_S(I)
	    SE(ID)%BA(EQION,VED,M,I)=SE(ID)%BA(EQION,VED,M,I) - T1 - T2
	  END IF
!
!
          IF(COMPUTE_BA)THEN
	    I=DPTH_INDX
	    DO J=1,N_S
              JJ=EQGS+J-1;
	      BA_ED(JJ,M,I)=BA_ED(JJ,M,I)+(ZION-1.0_LDP)
	      SE(ID)%BA(EQ_NUM_CONV,J,M,I)=SE(ID)%BA(EQ_NUM_CONV,J,M,I)+1.0_LDP
	    END DO
	  END IF
!
	  T1=0.0_LDP
	  I=DPTH_INDX
	  DO L=1,N_S
	    T1=T1+HN_S(L,I)
	  END DO
	  SE(ID)%STEQ(EQ_NUM_CONV,I)=SE(ID)%STEQ(EQ_NUM_CONV,I)+T1
	  STEQ_ED(I)=STEQ_ED(I)+(ZION-1.0_LDP)*T1
!
! We only include DI in the population and charge conservation equations
! if the higher ionization species is not present. Necessary to do this as
! DI is the ground state of the next species. We also correct the
! conservation equation for POP if the higher ionization species is not
! present.
!
! Note the charge on DI is ZION. Since we will ALWAYS declare the
! ground state to be an important variable, JJ should never be zero.
!
	  I=DPTH_INDX
	  IF(.NOT. NEXT_PRES)THEN
	     IF (IS_A_MOLECULE) THEN
		SE(ID)%STEQ(EQ_NUM_CONV,I)=SE(ID)%STEQ(EQ_NUM_CONV,I)+DI_S(I)-POP(I)
	     ELSE
		SE(ID)%STEQ(EQ_NUM_CONV,I)=SE(ID)%STEQ(EQ_NUM_CONV,I)+DI_S(I)-POP(I)
	     END IF
	    STEQ_ED(I)=STEQ_ED(I)+DI_S(I)*ZION
	  END IF
	  IF(COMPUTE_BA .AND. .NOT. NEXT_PRES)THEN
	    JJ=EQ_NUM_CONV
	    SE(ID)%BA(JJ,VION,M,I)=SE(ID)%BA(JJ,VION,M,I)+1.0_LDP
	    JJ=EQGS+N_S
	    BA_ED(JJ,M,I)=BA_ED(JJ,M,I)+ZION
	  END IF
!
	END DO		!Loop over depth -- IPROC
!
	RETURN
	END
