!
! Routine to increment the photoionization and recombination rates
! for an arbitrary ion. The bound-free cooling rate (in ergs/cm**3/s)
! is also computed. The Free-Free cooling rate is computed under the
! assumption that is is hydrogenic, and the ion has charge ZHYD.
!
! Separate quadrature weights are passed to compute the cooling rate.
!
!
	SUBROUTINE PRRR_SL_MPI_V1(PR,RR,BFCR,FF,WSE,WCR,
	1                     HN,HNST,NLEV,ZHYD,
	1                     DI,LOG_DIST,N_DI,
	1                     PHOT_ID,ION_LEV,ED,T,
	1                     JREC,JPHOT,JREC_CR,JPHOT_CR,BPHOT_CR,
	1                     NU,NU_CONT,INIT_ARRAYS,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 01-Dec-2023 : Added OMP paralleization.
! Altered 24-Sep-2023 : Adjusted constants for consistency. Some data ouput for checking between
!                            RE and EHB equations (LONG ver -- 15-Oct-2023).
! Altered 04-May-2022 : Changed HMI to H0
! Altered 23-Jun-2015 - Added H- free-free cooling.
! Altered 20-Oct-2011 - Now sum up all in ion levels for FF. Only do this when PHOT_ID=1
! Altered 05-Apr-2011 - Changed to V6.
!                       LOG_DIST (rather than dwHNST_F) is passed in call.
!                         Modifications done to allow lower temperaturs.
!                         Primary editing done 25-Jan-2011
! Altered 25-Jan-2010 : Bug fixed with previous alteration.
! Created 29-Nov-2010 : Based on PRRR_SL_V5
!                         LOG_DIST installed to prevent crashing caused by high ionization
!                         stages at low temperatures.
! Altered 03-Mar-2004 : NU installed in CALL. Computation of FF cooling revised.
! Altered 14-May-2001 : Bug fixed. Arrays were not being initialized
!                         correctly. Using continuum bands, ML may not be
!                         one on first call. Replaced ML by INIT_ARRAYS.
!                         Changed to V4.
! Altered 25-May-1996 : DIM_LIM removed (now use dynamic memory allocation for GFF_VAL)
! Altered 29-Sep-1995 : DI,DIST inserted to allow treatment of ionizations
!                         to multiple final states without the need of
!                         separate LTE population for each target level in
!                         the final ion.
!                       Call changed extensively. Now version V2
!                       FLAG deleted as testing of whether to initialize
!                         arrays can be done using PHOT_ID.
!
! Created 23-Sep-87 - Based on PRRRCOOLGEN_V3
!
!
! CONSTANTS FOR OPACITY ETC.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER NLEV			!Number of levels in species
	INTEGER N_DI			!Number of levels in final ion
	INTEGER ND			!Number of depth points.
	INTEGER PHOT_ID		!Photoionization ID
	INTEGER ION_LEV		!Final (destination) level in ion.
	INTEGER DST,DEND
!
	REAL(KIND=LDP) PR(NLEV,DST:DEND),RR(NLEV,ND),BFCR(NLEV,ND),FF(DST:DEND)
	REAL(KIND=LDP) HN(NLEV,DST:DEND),HNST(NLEV,ND),WSE(NLEV,ND),WCR(NLEV,ND)
	REAL(KIND=LDP) DI(N_DI,DST:DEND),LOG_DIST(N_DI,ND)
!
	REAL(KIND=LDP) ED(ND),T(ND)
	REAL(KIND=LDP) JREC(ND)
	REAL(KIND=LDP) JPHOT(ND)
	REAL(KIND=LDP) JREC_CR(ND)
	REAL(KIND=LDP) JPHOT_CR(ND)
	REAL(KIND=LDP) BPHOT_CR(ND)
	REAL(KIND=LDP) NU
	REAL(KIND=LDP) NU_CONT
!
	LOGICAL INIT_ARRAYS	        !Used to signify initialization
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER I,K
	REAL(KIND=LDP) POP_SUM,T1,T2,A1,TMP_HNST
	REAL(KIND=LDP) JB_RAT, JC_RAT
	REAL(KIND=LDP) H,ZHYD,CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Dynamic memory allocation for free-free gaunt factor as a function of depth.
!
	REAL(KIND=LDP) GFF;    EXTERNAL GFF
	REAL(KIND=LDP) GFF_VAL
	REAL(KIND=LDP) PLANCKS_CONSTANT
	EXTERNAL PLANCKS_CONSTANT
!
! 4PI*1.0E-10 (R scaling) Note that ordering is important or get underflow.
! FQW is approximately 10^15.
!
	H=PLANCKS_CONSTANT()*1.0E+15_LDP   		 !ergs/s (*1.0E+15 due to *nu)
!
! If ML=1 and and PHOT_ID .EQ. 1 then initialize all arrays. This routine
! should be called first for ionizations to the ground state.
!
	IF(INIT_ARRAYS .AND. PHOT_ID .EQ. 1)THEN
	  PR(:,:)=0.0_LDP
	  RR(:,:)=0.0_LDP
	  BFCR(:,:)=0.0_LDP
	  FF(:)=0.0_LDP
	END IF
!
! Note that JREC     = Int [ (2hv^3/c^2 +J) exp(-hv/kT)/v dv ]
!           JREC_CR  = Int [ (2hv^3/c^2 +J) exp(-hv/kT)   dv ]
!           JPHOT    = Int [ J/v dv]
!           JPHOT_CR = Int [ J dv]
!
! Since BFCR = Int (nu-edge)/nu, J?_CR is associated with WSE in the expression
! for BFCR.
!
	DO K=DST,DEND
	  IF(JREC(K) .GT. 0.0_LDP)THEN
	    JB_RAT=LOG(DI(ION_LEV,K)/DI(1,K))-LOG_DIST(ION_LEV,K)+LOG_DIST(1,K)
	    JB_RAT=EXP(LOG(JREC(K))+JB_RAT)
	  ELSE
	    JB_RAT=0.0_LDP
	  END IF
	  IF(JREC_CR(K) .GT. 0.0_LDP)THEN
	    JC_RAT=LOG(DI(ION_LEV,K)/DI(1,K))-LOG_DIST(ION_LEV,K)+LOG_DIST(1,K)
	    JC_RAT=EXP(LOG(JREC_CR(K))+JC_RAT)
	  ELSE
	    JC_RAT=0.0_LDP
	  END IF
	  DO I=1,NLEV
	    IF(WSE(I,K) .NE. 0)THEN
	      PR(I,K)=PR(I,K)+WSE(I,K)*HN(I,K)*JPHOT(K)
	      RR(I,K)=RR(I,K)+WSE(I,K)*HNST(I,K)*JB_RAT
	      BFCR(I,K)=BFCR(I,K)+
	1          ( HNST(I,K)*(WCR(I,K)*JB_RAT+WSE(I,K)*JC_RAT)
	1              -HN(I,K)*(WCR(I,K)*JPHOT(K)+WSE(I,K)*JPHOT_CR(K)) )*H
	    END IF
	  END DO
	END DO
!
! Compute Free-Free cooling.
!
	IF(ZHYD .EQ. 0.0_LDP)THEN
!
! This is for H-. We use T1 as a temporary storage for the ground state
! population of neutral hydrogen.
!
	  DO K=DST,DEND
	    T1=DI(1,K)
	    CALL DO_H0_FF_COOL(FF(K),T1,ED(K),T(K),BPHOT_CR(K),JPHOT_CR(K),NU_CONT,IONE)
	  END DO
!
	ELSE IF(ION_LEV .EQ. 1)THEN
!  
! The opacity is evaluated at NU_CONT. However, the stimulated emission term
! should be evaluated at NU, since we correct CHI and ETA for the change
! in freqency but assuming a constant cross-section.
!
! The constant in T2 is 4PI x 1.0E-10.
!
! Compute free-free gaunt factors.
!
	  T2=1.256637061E-09_LDP*ZHYD*ZHYD*CHIFF/(NU_CONT**3)
	  DO K=DST,DEND
            GFF_VAL=GFF(NU,T(K),ZHYD)
!            CALL FF_RES_GAUNT(GFF_VAL,NU,T(K),ID,GION,ZHYD,IONE)
!
	    POP_SUM=SUM(DI(:,K))
	    A1=EXP(-HDKT*NU/T(K))
	    FF(K) =FF(K)+T2*ED(K)*POP_SUM/SQRT(T(K))*(1.0_LDP-A1)
	1           *GFF_VAL*( BPHOT_CR(K)-JPHOT_CR(K) )
	  END DO
	END IF
!
	RETURN
	END
