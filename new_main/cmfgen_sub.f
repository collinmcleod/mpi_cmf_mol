!
! Main Subroutine to self consistently solve the equation of transfer and the
! equations of statistical equilibrium for a spherically extended atmosphere
! in the presence of outflows.
!
! At present the atmosphere considered to consist of
!      H, He, plus other species such as C, N, O and Fe.
! Eithe H or He must be present.
!
! Abundances are un-normalized relative fractional abundances (i.e. specified
! specified with respect to some arbitrary species X. Generally we have used
! He as X. If negative they are interpreted as mass-fractions.
!
! Several options are available for handeling both the lines and continuum.
! Some options have not been recently tested.
!
! Inclusion of new species is generally straightforward. Only the
! calling routine, CMFGEN, needs to modified. The main work is in
! collating the atomic data. Dynamical memory allocation is used.
!
! 
!
	SUBROUTINE CMFGEN_SUB(ND,NC,NP,NT,
	1                     NUM_BNDS,NION,DIAG_INDX,
	1                     NDMAX,NPMAX,NCF_MAX,NLINE_MAX,
	1                     TX_OFFSET,MAX_SIM,NM,NM_KI,NLF)
	USE SET_KIND_MODULE
	USE MPI
	USE MOD_CMFGEN
	USE ANG_QW_MOD
!	USE CMF_SOB_MOD
	USE CONTROL_VARIABLE_MOD
	USE MOD_VAR_OPAC_J
	USE OPAC_MOD
	USE STEQ_DATA_MOD
	USE MOD_LEV_DIS_BLK
	USE LINE_VEC_MOD
	USE LINE_MOD
	USE LUMINOSITY_MOD 
	USE RADIATION_MOD
	USE SHOCK_POWER_MOD
	USE UPDATE_KEYWORD_INTERFACE
	USE VAR_RAD_MOD_MPI_V1
	USE EDDFAC_REC_DEFS_MOD
	USE MOD_J_TWO_PHOT_MPI_V1, ONLY : SET_POP_FOR_TWOJ, INIT_GET_J_FOR_TWO_PHOT
	IMPLICIT NONE
!
! Altered 10-Feb-2025 : SCRTEMP writes only done by process 0.
! Altered 02-Aug-2024 : Minor cleaning. Changed minimum valuer for RLUMST & warning output.
! Altered 17-Aug-2019 : Default for DO_T_AUTO is now TRUE. Incorporated into IBIS versions.
! Altered 16-Aug-2019 : Now output Planck mean to RVTJ.
! Altered Jul/Aug 2019: Extensive changes to treat electron energy balance equation (done on OSIRIS),
!                          and PNT srce option.
!                          Call SET_ANG_QW_V2 insted of SET_ANG_QW_V2
!                          Call STEQ_MULTI_V10 insted of STEQ_MULTI_V9
!                          Calls to EHB routines added.
! Aleterd 06-Sep-2016 : Output more R digits to MEANOPAC.
! Altered 29-Sep-2015 : Changed back to COMP_OPAC (Different TWO_PHOT in routine)
! Altered 24-Jul-2015 : Added HMI, Changed to COMP_OPAC_V2 (cur_hmi, 19-Aug-2015)
! Altered 27-Mar-2013 : dE_WORK and RAD_DECAY_LUM updated for clumping.
!                       LUWARN inserted (done earlier)
!                       LIN_PROF_SIM is now depth dependent -- code cannow handle depth dependent
!                          line profiles.
!                       Calls to AUTO_CLUMP_REV and SPECIFY_IT_CYCLE added.
! Altered  6-Dec-2013 : When removing lines from the variation set we assume a minimum terminal
!                          velocity of 300 km/s. Important for plane-parallel models with Vinf=0.
! Altered 29-Nov-2011 : Call to STEQ_CO_MOV_DERIV changed to _V3. Done to facilitae
!                         hadling of additional levels with time dependence.
! Altered 07-Nov-2011 : Now output MNL_F, MNUP_F to NETRATE and TOTRATE.
! Altered       -2011 : Extensive alterations to handle low temperatures.
!                       Extensive alterations to include non-thermal ionizations
! Altered 19-Jan-2009 : SL otion inserted; rd_f_to_s_ids_v2.f now used.
! Altered 16-Feb-2006 : Changed and modified over 2 month period. Section solving for
!                         populations and performing NG acceleration etc removed to
!                         subroutine (SOLVE_FOR_POPS). Routines added to allow time
!                         variability of the statistical equilibrium equations.
!                         Currently these are only for a Hubble law flow. Relativistic
!                         terms added to COMP_OBS (now COMP_OBS_V2).
! Altered 20-Feb-2005 : Changed to use FLUX_MEAN & ROSS_MEAN which are defined in
!                         MOD_CMFGEN. Previusly used FLUXMEAN & ROSSMEAN defined in
!                         RADIATION_MOD.
!
	INTEGER ND,NC,NP,NT
	INTEGER NUM_BNDS,NION,DIAG_INDX
	INTEGER NDMAX,NPMAX
	INTEGER NCF_MAX,NLINE_MAX
	INTEGER NM,NM_KI,MAX_SIM,NLF
	INTEGER TX_OFFSET
!
	INTEGER NCF
	LOGICAL, PARAMETER :: IMPURITY_CODE=.FALSE.
!
	CHARACTER(LEN=12), PARAMETER :: PRODATE='12-Feb-2026'		!Must be changed after alterations
!
! 
!
	REAL(KIND=LDP) SOL(NT,ND)		!Temp. stor. area for ST. EQ.
!
! Constants for opacity etc. These are set in CMFGEN.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) OPLIN,EMLIN
!
! Internally used variables
!
	REAL(KIND=LDP) S1,REPA
	REAL(KIND=LDP) MAXCH,MAXCH_SUM
	REAL(KIND=LDP) T1,T2,T3,T4,SRAT
	REAL(KIND=LDP) FL,AMASS,FL_OLD
	REAL(KIND=LDP) FG_COUNT
	REAL(KIND=LDP) SCL_FAC
	REAL(KIND=LDP) SUM_BA
	REAL(KIND=LDP) RLUMST_BND
	REAL(KIND=LDP) LUM_FOR_TEFF
	REAL(KIND=LDP) LUM_SCL_FAC
!
	LOGICAL LST_DEPTH_ONLY
!
! REC_SIZE     is the (maximum) record length in bytes.
! UNIT_SIZE    is the number of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the number of bytes used to represent the number.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
 	INTEGER REC_SIZE
	INTEGER UNIT_SIZE
	INTEGER WORD_SIZE
	INTEGER N_PER_REC
!
! 
!
! Logical Unit assignments. Those indicated with a # after the ! are open in
!  large sections of the code. Other units generally used temprarily.
!
	INTEGER               LUER      	!Output/Error file.
	INTEGER               LUWARN	        !
	INTEGER, PARAMETER :: LUIN=7            !General input unit (closed after accesses).
	INTEGER, PARAMETER :: LUMOD=8           !Model Description file.
!
	INTEGER, PARAMETER :: LU_DC=9      	!Departure coefficient Output.
	INTEGER, PARAMETER :: LU_FLUX=10   	!Flux/Luminosity Data (OBSFLUX)
	INTEGER, PARAMETER :: LU_SE=16     	!Statistical equilibrium and solution arrays.
	INTEGER, PARAMETER :: LU_NET=17    	!# Line Netrate data.
	INTEGER, PARAMETER :: LU_OPAC=18   	!Rosseland mean opacity etc.
	INTEGER, PARAMETER :: LU_DR=19     	!# Downward rate (Nu. Z. A).
	INTEGER, PARAMETER :: LU_EW=20     	!# EW data.
	INTEGER, PARAMETER :: LU_REC_CHK=21	!# EW data.
	INTEGER LU_T_EHB
!
! For writing scratch file (SCRTEMP). Also used in reading in  MODEL data.
!
	INTEGER, PARAMETER :: LUSCR=26
!
	INTEGER, PARAMETER :: LU_HT=27     !#LINEHEAT (i.e Line heating term in R.E. equation)
!
! Used for RVTJ file and POPCARB, POPNIT etc.
!
	INTEGER, PARAMETER :: LU_POP=30
!
	INTEGER, PARAMETER :: LU_IMP=34       !J and CHI for impurity calculation.
	INTEGER, PARAMETER :: LU_EDD=35       !Continuum Eddington factors.
	INTEGER, PARAMETER :: LU_JEW=36       !J for EW computation (JEW)
	INTEGER, PARAMETER :: LU_JCOMP=37     !J_COMP
	INTEGER, PARAMETER :: LU_ES=38        !ES_J_CONV
!
! Following is used output the BA matrix, and its associated
! pointer file.
!
	INTEGER, PARAMETER :: LU_BA=40
!
! For listing of transitions with TOTAL negative opacity values at some depths.
!
	INTEGER, PARAMETER :: LU_NEG=75
!
! 
!
	INTEGER NNM				!Include cont. var in line var.
	INTEGER NL,NUP
	INTEGER MNL,MNUP
	INTEGER MNL_F,MNUP_F
	INTEGER PHOT_ID
	INTEGER DPTH_INDX
	INTEGER ROSS_PHOT_DPTH_INDX
	INTEGER VAR_INDX
	INTEGER I,J,K,L,ML,LS,LINE_INDX,NEXT_LOC
	INTEGER IREC,MATELIM
!
	CHARACTER*80 TMP_STRING
	CHARACTER*20 TMP_KEY
	INTEGER ID,LOC_ID,ID_SAV,JJ
        INTEGER  L1,L2,U1,U2
	INTEGER IT,MNT,NIV
	INTEGER ISPEC
	INTEGER GAMMA_ADD_SLOWLY_COUNTER
!
! Main iteration loop variables.
!
	INTEGER MAIN_COUNTER,NITSF,NUM_ITS_RD,NUM_ITS_TO_DO
	LOGICAL LST_ITERATION
!
! Functions called
!
	INTEGER ICHRLEN,ERROR_LU,WARNING_LU
	REAL(KIND=LDP) DOP_PRO
	REAL(KIND=LDP) S15ADF
	REAL(KIND=LDP) LAMVACAIR
	REAL(KIND=LDP) ATOMIC_MASS_UNIT
	REAL(KIND=LDP) SPEED_OF_LIGHT
	REAL(KIND=LDP) GRAVITATIONAL_CONSTANT
	REAL(KIND=LDP) RAD_SUN
	REAL(KIND=LDP) TEFF_SUN
	REAL(KIND=LDP) MASS_SUN
	LOGICAL EQUAL
	EXTERNAL ICHRLEN,ERROR_LU,WARNING_LU,SPEED_OF_LIGHT,GRAVITATIONAL_CONSTANT
	EXTERNAL MASS_SUN,RAD_SUN,TEFF_SUN
!
	INTEGER GET_DIAG
	INTEGER BNDST
	INTEGER BNDEND
	INTEGER BND_TO_FULL
!
! Photoionization cross-section routines.
!
! Collisional routines.
!
	EXTERNAL OMEGA_GEN_V3
!
! Wind variablity arrays.
!
	REAL(KIND=LDP) POPS(NT,ND)		!Population for all species.
	REAL(KIND=LDP) MEAN_ATOMIC_WEIGHT	!Mean atomic weight of atoms  (neutrals
!                          		! and ions) in atomic mass units.
	REAL(KIND=LDP) ABUND_SUM
!
!
! Arrays for improving on the initial T structure --- partition functions.
! Need one for each atomic species.
!
	REAL(KIND=LDP), ALLOCATABLE :: U_PAR_FN(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: PHI_PAR_FN(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: Z_PAR_FN(:)
!
	REAL(KIND=LDP) TGREY(ND)
	REAL(KIND=LDP) T_SAVE(ND)
	REAL(KIND=LDP) ZNET(ND)
!
! Variables for scaling the line cooling rates in oder that the radiative
! equilibrium equation is more consistent with the electron heating/cooling
! equation. The scaling is done when the line frequency is with a fraction
! of SCL_LINE_HT_FAC of the tmean frequency for the super-level under
! consideration. 0.5 is presently the prefered value.
!
!	REAL(KIND=LDP) AVE_ENERGY(NT)		!Average energy of each super level
	REAL(KIND=LDP) STEQ_T_SCL(DST:DEND)
	REAL(KIND=LDP) STEQ_T_NO_SCL(DST:DEND)
! 
!
! Dielectronic recombination variables and arrays.
!
	INTEGER NMAXDIE
	PARAMETER (NMAXDIE=500)
!
	REAL(KIND=LDP) EDGEDIE(NMAXDIE)		!Ionization frequency (negative)
	REAL(KIND=LDP) EINADIE(NMAXDIE)		!Einstein A coefficient
	REAL(KIND=LDP) GUPDIE(NMAXDIE)		!Stat. weight of autoionizing level.
!
	INTEGER LEVDIE(NMAXDIE)  	!Indicates MNL of low state
	INTEGER INDXDIE(NMAXDIE)
!
! Used for species identification as INDXDIE is not unique (specied not
! present can have same index as a species thats present.)
!
	CHARACTER*10 SPECDIE(NMAXDIE)
	CHARACTER*35 DIENAME(NMAXDIE)
!
	INTEGER NDIETOT
!
! Arrays and variables used for both Dielectronic recombination, and
! the implicit recombination.
!
	INTEGER EQION,EQSPEC
	REAL(KIND=LDP) GLOW,GION
	REAL(KIND=LDP) NUST(ND)			!LTE autoionizing population.
	REAL(KIND=LDP) DION(ND)			!Ion population
	REAL(KIND=LDP), ALLOCATABLE :: DIECOOL(:,:)    !Dielec. cooling check for all spec.
! 
!
! Opacity/emissivity
!
	REAL(KIND=LDP) CHIL(ND)                 !Line opacity (without prof.)
	REAL(KIND=LDP) ETAL(ND)                 !Line emissivity (without prof.)
!
! Quadrature weights.
!
	REAL(KIND=LDP) FQW(NCF_MAX)		!Frequency weights
!
! Transfer equation vectors
	REAL(KIND=LDP) R_OLD(NDMAX)		!Used to store previous R grid in SN models.
!
! Line vectors
	REAL(KIND=LDP) AV(ND)
	REAL(KIND=LDP) VB(NDMAX)		!Used for error calculations
	REAL(KIND=LDP) VC(NDMAX)		!Used for error calculations
	REAL(KIND=LDP) H(ND)
	REAL(KIND=LDP) Q(ND)			!FREQ DEPENDENT.
	REAL(KIND=LDP) QH(ND)			!  "      "
	REAL(KIND=LDP) GAM(ND)			!FREQ INDEPENDENT
	REAL(KIND=LDP) GAMH(ND)			!  "      "
! 
!
! Arrays and variables for computation of the continuum intensity
! using Eddington factors. This is separate to the "inclusion of
! additional points".
!
	LOGICAL EDDINGTON
!
! Variables for EW's and LINE blanketing.
!
	REAL(KIND=LDP) CONT_INT,EW
	INTEGER ACCESS_JEW
	LOGICAL COMPUTE_EW,COMPUTE_JEW,COMPUTE_LAM,MID,FULL_ES
!
! ACESS_F is the current record we are writing in EDDFACTOR.
!
	INTEGER ACCESS_F
	INTEGER NDEXT,NCEXT,NPEXT
!
	REAL(KIND=LDP) CNM(NDMAX,NDMAX)		!For collisions cross-section in
	REAL(KIND=LDP) DCNM(NDMAX,NDMAX)	!STEQGEN
!
	INTEGER, PARAMETER :: N_FLUX_MEAN_BANDS=12
	REAL(KIND=LDP)     LAM_FLUX_MEAN_BAND_END(N_FLUX_MEAN_BANDS)
	REAL(KIND=LDP)     BAND_FLUX_MEAN(ND,N_FLUX_MEAN_BANDS)
	REAL(KIND=LDP)     BAND_FLUX(ND,N_FLUX_MEAN_BANDS)
	DATA LAM_FLUX_MEAN_BAND_END/100.0_LDP,150.0_LDP,200.0_LDP,227.83_LDP,258.90_LDP,300.0_LDP,504.25_LDP,911.75_LDP,
	1                         1200.0_LDP,1500.0_LDP,2000.0_LDP,1.0E+08_LDP/
!
! Continuum frequency variables and arrays.
!
	REAL(KIND=LDP) NU(NCF_MAX)		!Continuum and line frequencies
	REAL(KIND=LDP) NU_EVAL_CONT(NCF_MAX)	!Frequencies to evaluate continuum
	REAL(KIND=LDP) OBS(NCF_MAX)		!Observers spectrum
!
! Vectors and arrays used for the observed flux.
!
	INTEGER N_OBS
	REAL(KIND=LDP) OBS_FREQ(NCF_MAX)		!Since N_OBS < NCF =< NCF_MAX
	REAL(KIND=LDP) OBS_FLUX(NCF_MAX)
	LOGICAL FIRST_OBS_COMP
!
	CHARACTER TIME*20
	CHARACTER FMT*120
	CHARACTER*20 SECTION,FORMAT_DATE*20
	CHARACTER STRING*132
	CHARACTER EW_STRING*132
	CHARACTER TEMP_CHAR*132
	CHARACTER*2 FORMFEED
!
! Global vectors:
!
	REAL(KIND=LDP) AMASS_ALL(NT)
	INTEGER N_LINE_FREQ
!
	INTEGER LINES_THIS_FREQ(NCF_MAX)
!
	REAL(KIND=LDP) NU_DOP
	REAL(KIND=LDP) NU_MAX_OBS
	REAL(KIND=LDP) NU_MIN_OBS
!
	INTEGER FREQ_INDX
	INTEGER X_INDX
	INTEGER FIRST_LINE
	INTEGER LAST_LINE
!
! Variables to limit the computation of the continuum opacities and
! emissivities.
!
	REAL(KIND=LDP) JREC(DST:DEND)
	REAL(KIND=LDP) dJRECdT(DST:DEND)
	REAL(KIND=LDP) JPHOT(DST:DEND)
	REAL(KIND=LDP) JREC_CR(DST:DEND)
	REAL(KIND=LDP) dJREC_CRdT(DST:DEND)
	REAL(KIND=LDP) JPHOT_CR(DST:DEND)
	REAL(KIND=LDP) BPHOT_CR(DST:DEND)
!
	REAL(KIND=LDP) CONT_FREQ
	LOGICAL FINAL_CONSTANT_CROSS
!
! Indicates whether APRXzV, FFXzZ etc should be zeroed.
!
	LOGICAL ZERO_REC_COOL_ARRAYS
!
! 
!
	REAL(KIND=LDP) Z_POP(NT)		!Ionic charge for each species
!
! Variables etc for computation of continuum in comoving frame.
!
	LOGICAL FIRST_FREQ
	LOGICAL RAT_TOO_BIG
	LOGICAL NEW_FREQ
!
! 
!
! X-ray variables.
! We dimension from 0 so that we can access a Null vector for the 1st included
! ioinization stage of each species.
!
	REAL(KIND=LDP) XRAY_HEATING(ND)
!
	REAL(KIND=LDP) OBS_XRAY_LUM_0P1
	REAL(KIND=LDP) OBS_XRAY_LUM_1KEV
	REAL(KIND=LDP) GFF,XCROSS_V2
	EXTERNAL GFF,XCROSS_V2
!
	REAL(KIND=LDP) SPEC_DEN(ND,NUM_SPECIES)		!Used by ELEC_PREP
	REAL(KIND=LDP) AT_NO_VEC(ND,NUM_SPECIES)
!
	REAL(KIND=LDP) AD_COOL_V(DST:DEND)
	REAL(KIND=LDP) AD_COOL_DT(DST:DEND)
!
	REAL(KIND=LDP) ARTIFICIAL_HEAT_TERM(DST:DEND)
	REAL(KIND=LDP) dE_RAD_DECAY(DST:DEND)
	REAL(KIND=LDP) dE_WORK(DST:DEND)
	REAL(KIND=LDP) dE_SHOCK_POWER(DST:DEND)
!
	LOGICAL FIRST
	LOGICAL CHK,SUCCESS
        LOGICAL VAR_SOB_JC
	LOGICAL NEG_OPACITY(ND),FIRST_NEG
	LOGICAL AT_LEAST_ONE_NEG_OPAC
	LOGICAL FILE_OPEN
	LOGICAL VERBOSE
	LOGICAL F_TO_S_RD_ERROR
	LOGICAL TMP_LOGICAL
!
! Inidicates approximate frequencies for which TAU at outer boundary is written
! to OUTGEN on the last iteration.
!
! They are the He2 ege, NIII/CIII egde, HeI, HI, HI(N=2).
!
	INTEGER, PARAMETER :: N_TAU_EDGE=5
	REAL(KIND=LDP) TAU_EDGE(N_TAU_EDGE)
	DATA TAU_EDGE/13.16_LDP,11.60_LDP,5.95_LDP,3.29_LDP,0.83_LDP/
!
!***********************************************************************
!
!*******************FUNCTION DEFINITIONS********************************
!
! This function takes a band-index and converts it the equivalent index
! in the full matrix. L=BND_TO_FULL(J,K) is equivalent to the statements:
!     IF(NUM_BNDS .EQ. ND)THEN L=J ELSE L=K+J-DIAG_INDX END IF
! The second indice is the equation depth.
!
	BND_TO_FULL(J,K)=(NUM_BNDS/ND)*(DIAG_INDX-K)+K+J-DIAG_INDX
!
! This function computes the index L on BA( , ,?,K) corresponding
! to the local depth variable (i.e that at K). It is equivalent
! to IF (NUM_BNDS .EQ. ND)THEN L=K ELSE L=DIAG END IF
!
	GET_DIAG(K)=(NUM_BNDS/ND)*(K-DIAG_INDX)+DIAG_INDX
!
! These two functions compute the start and end indices when updating
! VJ. eg. we do not wish to update VJ( ,1, ) if we are using the banded
! matrix since this refers to a variable beyond the outer atmosphere.
!
	BNDST(K)=MAX( (NUM_BNDS/ND)*(K-DIAG_INDX)+1+DIAG_INDX-K, 1 )
	BNDEND(K)=MIN( (NUM_BNDS/ND)*(K-DIAG_INDX)+ND+DIAG_INDX-K,
	1                 NUM_BNDS )
! 
!
!	INCLUDE 'mpif.h'
!
! Call MPI initialization routines. These must be the first executable
! statements in the code. 
! MPYE will be use to label each process, and NTHREAD is the number of
! processors. To avoid
! multiple output, we will use processor 0 for most output.
!       
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NTHREAD,IERR)
!
!****************************************************************************
!
! Initialization section
!
	DO I=0,NTHREAD-1
	 IF(MYPE .EQ. 0 .AND. MYPE .EQ. I)THEN
	    WRITE(6,'(/,A,/,A)')' Top of CMFGEN_SUB',' Processor depth allocations'
	    WRITE(6,'(T20,3(3X,A))')'MYPE',' DST','DEND'
	    WRITE(6,'(T20,3I7)')MYPE,DST,DEND; FLUSH(UNIT=6)
	    FLUSH(UNIT=6)
	  ELSE IF(MYPE .EQ. I)THEN
	    WRITE(6,'(T20,3I7)')MYPE,DST,DEND; FLUSH(UNIT=6)
	  END IF
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	END DO
	IF(MYPE .EQ. 0)WRITE(6,*)' '; FLUSH(UNIT=6)
!
	LUER=ERROR_LU()
	LUWARN=WARNING_LU()
	COMPUTE_LAM=.FALSE.
	COMPUTE_EW=.TRUE.
	FULL_ES=.TRUE.
	SN_MODEL=.FALSE.
	VINF=0.0_LDP				!Will be reset later
        TREAT_NON_THERMAL_ELECTRONS=.FALSE.
	INCL_RADIOACTIVE_DECAY=.FALSE.
	INC_SHOCK_POWER = .FALSE.
	ZERO_REC_COOL_ARRAYS=.TRUE.
	GAMMA_ADD_SLOWLY_COUNTER=-1
	I=12
	FORMFEED=' '//CHAR(I)
	CNT_FIX_BA=0
	MAXCH_SUM=0.0_LDP
	LST_ITERATION=.FALSE.
	LUM_SCL_FAC=4.1274E-12_LDP              !(4*PI)**2*Dex(+20)/L(sun)
!
	DPTH_INDX=3
	DPTH_INDX=MIN(DPTH_INDX,ND)		!Thus no problem if 84 > ND
	VAR_INDX=366
	VAR_INDX=MIN(VAR_INDX,NT)
	CALL GET_VERBOSE_INFO(VERBOSE)
!
! When TRUE, FIXED_T indicated that T is to be heled fixed (at least at some
! depths) in the linearization. This variable is set automatically by the
! code depending on the magnitude of the corrections.
!
	FIXED_T=.FALSE.
!
! Set NDIETOT to zero in case no dielectronic lines are include in
! the dielectronic section.
!
	NDIETOT=0
!
! Set so that it is defined for the TEST whether to do an accurate flux
! calculation.
!
	MAXCH=100._LDP
!
! A value of 1000 is used to indicate that the last change was greater them
! VAL_DO_NG, or a NEW_MODEL.
!
! A value of 2000 indicates a continuing model. In this case LAST_NG
! take precedence. NB: NEXT_NG is reset to 1000 for a new model in the
! new model section.
!
	NEXT_NG=2000			!Initial value Indicate model just bega
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
!
! 
!
	CALL GEN_ASCI_OPEN(LU_SE,'STEQ_VALS','UNKNOWN','APPEND',' ',IZERO,IOS)
!
! Open a scratch file to record model parameters. This file will eventually
! be renamed MODEL.
!
	WRITE(STRING,'(I3.3)')MYPE; STRING='MODEL_SCR_'//STRING	
	CALL GEN_ASCI_OPEN(LUSCR,STRING,'UNKNOWN',' ',' ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening MODEL_SCR in CMFGEN, IOS=',IOS
	  STOP
	END IF
	CALL SET_LINE_BUFFERING(LUSCR)
	WRITE(LUSCR,'()')
!
! Read in parameters which can change during a single model run. These
! parameters contol the number of iterations, and whether we wish to perform
! a LAMBDA iteration.
!
	  CALL GEN_ASCI_OPEN(LUIN,'IN_ITS','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(LUER,*)'Error opening IN_ITS in CMFGEN, IOS=',IOS
	     STOP
	  END IF
	  CALL RD_OPTIONS_INTO_STORE(LUIN,LUSCR)
	  CALL RD_STORE_INT(NUM_ITS_TO_DO,'NUM_ITS',L_TRUE,'Number of iterations to perform')
	  CALL RD_STORE_LOG(RD_LAMBDA,'DO_LAM_IT',L_TRUE,'Do LAMBDA iterations ?')
	  DO_LAMBDA_AUTO=.TRUE.
	  CALL RD_STORE_LOG(DO_LAMBDA_AUTO,'DO_LAM_AUTO',L_FALSE,
	1                      'Start non-lambda iterations automatically?')
	  DO_GREY_T_AUTO=.TRUE.
	  CALL RD_STORE_LOG(DO_GREY_T_AUTO,'DO_GT_AUTO',L_FALSE,
	1                      'Do a grey temperature iteration after revising USE_FIXED_J?')
	  DO_T_AUTO=.FALSE.
	  CALL RD_STORE_LOG(DO_T_AUTO,'DO_T_AUTO',L_FALSE,
	1                      'Allow temperature to vary when sufficent convergence has been obtained?')
	  SET_POPS_D2_EQ_D1=.FALSE.
	  CALL RD_STORE_LOG(SET_POPS_D2_EQ_D1,'D2_EQ_D1',L_FALSE,
	1                      'Replace pops at depth 2 with those at depth 1 for non-LAMBDA it?')
	  CALL CLEAN_RD_STORE()
	  CLOSE(UNIT=LUIN)
	NUM_ITS_RD=NUM_ITS_TO_DO
!
! These two parameters were originally read in from IN_ITS, however they
! are no longer required by program when using the BAND solution. They
! are still passed, however, to SOLVEBA.
!
	REPA=1.2_LDP
	MATELIM=1

!
! This section does the following :-
!	1) If new model read/determine the new radius scale  and
! populations.
!	2) If old model, populations are read in from scratch file
! previous radius scale etc are used. Ponit1 and point2
! point to the input data record (Note : Single Record)
!
	CALL RD_CONTROL_VARIABLES(LUIN,LUSCR,LUER,NUM_BNDS)
!
! Get the Solar abundance scale that will be used when MOD_SUM is
! created.
!
	CALL RD_SOL_ABUND_SCAL_MPI_V1(SOL_MASS_FRAC,SOL_ABUND_HSCL,AT_MASS,AT_NO,
	1            SOL_ABUND_REF_SET,L_FALSE,NUM_SPECIES)
!
! RMDOT is the density at R=10dex10 cm and V=1km/s (atomic mass units)
!
	RMDOT=RMDOT*3.02286E+23_LDP
!
! LAMBDA_ITERATION controls whether a LAMBDA iteration is performed.
! A Lambda iteration is forced if RD_LAMDA is true. FIXED_T is set true
! as the Radiative equilibrium equation is not linearized if we are
! performing a LAMBDA iteration (could be changed with effort).
! The FIX_IMPURITY option is only used in non-LAMBDA mode.
!
	LAMBDA_ITERATION=RD_LAMBDA
	IF(LAMBDA_ITERATION)THEN
	  FIX_IMPURITY=.FALSE.
          FIXED_T=.TRUE.
	ELSE
	  FIX_IMPURITY=RD_FIX_IMP
	END IF
!
! This ensures that the ITS_DONE keyword is in HYDRO_DEFAULTS.
!
	IF(DO_HYDRO .AND. .NOT. SN_MODEL)THEN
	  CALL CHECK_HYDRO_DEF(STRING,LUIN,LUER)
	END IF
!
! 
!
	IF(ACCURATE)THEN
!
! We first verify that the interpolation range is valid.
!
	  IF(END_INTERP_INDX .GT. ND)END_INTERP_INDX=ND
	  IF(DEEP .GT. ND)DEEP=MIN(5,ND)
	  NDEXT=(END_INTERP_INDX-ST_INTERP_INDX)*NPINS+ND
	  IF(NDEXT .GT. NDMAX)THEN
	    WRITE(LUER,*)' Error - NDEXT larger than NDMAX in CMFGEN'
	    WRITE(LUER,*)' Need to increase NDMAX in CMFGEN'
	    STOP
	  END IF
	  NCEXT=NC
!
! NB: The following expression guarentees that NPEXT has the same relationship
! to NDEXT and NCEXT as does NP to ND and NC.
!
	  NPEXT=NDEXT+NCEXT+(NP-ND-NC)
	  IF(NPEXT .GT. NPMAX)THEN
	    WRITE(LUER,*)' Error - NPEXT larger than NPMAX in CMFGEN'
	    WRITE(LUER,*)' Need to increase NPMAX in CMFGEN'
	    STOP
	  END IF
	ELSE
	  NDEXT=ND; NCEXT=NC; NPEXT=NP
	END IF
!
	CALL SET_LUMINOSITY_MOD(ND)
	CALL SET_RADIATION_MOD(DST,DEND,ND,NDMAX,NPMAX)
	CALL SET_LINE_MOD(DST,DEND,ND,NT,MAX_SIM,NM)
        CALL SET_VAR_RAD_MOD_MPI_V1(DST,DEND,ND,NDEXT,
	1        NT,NUM_BNDS,NM,MAX_SIM,NM_KI,ACCURATE,L_TRUE)
	CALL SET_CMF_SOB_MOD(ND,NUM_BNDS,NT,NM_KI,NLF,LUER)
	CALL SET_VAR_OPAC_J(ACCURATE,DST,DEND,ND,NM_KI,MAX_SIM,NM,NT,NUM_BNDS)
!
!	T1=10.0_LDP/(NLF-1)
!	DO ML=1,NLF
!	  PF(ML)=5.0_LDP-T1*(ML-1)
!	END DO
! 
!
! Read in bound-free gaunt factors for individual n states of hydrogen,
! and hydrogenic cross-sections for individual l states (n =0 to 30,
! l=0 to n-1)
!
	CALL RD_HYD_BF_DATA(LUIN,LUSCR,LUER)
!
! Read in atomic data for 2-photon transitions.
!
	CALL RD_TWO_PHOT(LUIN,INCL_TWO_PHOT)
!
! Read in data for charge exchange reactions.
!
	CALL RD_CHG_EXCH_V3(LUIN,INCL_CHG_EXCH)
!
! Read in X-ray photoionization cross-sections.
!
	CALL RD_XRAY_FITS(LUIN)
!
	IF(XRAYS .AND. .NOT. FF_XRAYS)THEN
	  CALL RD_XRAY_SPEC(T_SHOCK_1,T_SHOCK_2,LUIN)
	END IF
!
	IF(XRAYS .AND. ADD_XRAYS_SLOWLY .AND. RD_LAMBDA)THEN
	   FILL_X1_SAV=FILL_FAC_XRAYS_1
	   FILL_X2_SAV=FILL_FAC_XRAYS_2
	   FILL_FAC_XRAYS_1=FILL_FAC_X1_BEG
	   FILL_FAC_XRAYS_2=FILL_FAC_X2_BEG
	END IF
!
	IF(TREAT_NON_THERMAL_ELECTRONS .AND. (SCL_NT_CROSEC .OR. SCL_NT_ION_CROSEC))THEN
	  CALL RD_NT_CROSEC_SCLFAC_V2(LUIN,LUER)
	END IF
!
	IF(TREAT_NON_THERMAL_ELECTRONS .AND. ADD_DEC_NRG_SLOWLY .AND. RD_LAMBDA)THEN
	  DEC_NRG_SCL_FAC=DEC_NRG_SCL_FAC_BEG
	ELSE
	  DEC_NRG_SCL_FAC=1.0_LDP
	END IF
!
	CALL RD_NUC_DECAY_DATA_V2(INCL_RADIOACTIVE_DECAY,GAMRAY_TRANS,ND,LUIN)
!
	IF (INC_SHOCK_POWER .AND. ADD_SHOCK_POWER_SLOWLY .AND.  RD_LAMBDA) THEN
	  SHOCK_POWER_FAC = SHOCK_POWER_FAC_BEG
	ELSE
	  SHOCK_POWER_FAC = 1.0_LDP
	ENDIF
!
! 
!
! Read in oscillator strengths, the photoionization cross section data,
! dielectronic data, and implicit recombination data for carbon.
! Individual species are grouped together (rather than grouping all the
! oscillator reads) so that the headers of the INPUT files are grouped
! in the MODEL output file.
!
! We do this in reverse order so that GIONXzV can be correctly set.
!
! Note Well - in GENOSICL  - T1 is returned with the ionization energy.
!                            T2 is returned with screened nuclear charge.
!                             I is returned with the number of transitions.
!
! RDGENDIE returns the dielectronic transitions as a line list. These lines
! are treated as individual lines.
!
! RD_XzV_PHOT_DIE associates the dielectronic transitions with the
! photoionization cross-sections. They are then handeled as part of
! the continuum cross-sections.
!
! NB: The passed GF_CUT is set to zero if the Atomic NO. of the species
!      under consideration is less than AT_NO_GF_CUT.
!
	F_TO_S_RD_ERROR=.FALSE.
	WRITE(LUWARN,'(/,A,/)')' Reading in atomic data'
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	    IF( ATM(ID)%XzV_PRES)THEN
	      IF( NINT(AT_NO(SPECIES_LNK(ID))) .LT. NINT(AT_NO_GF_CUT) )THEN
	        T2=0.0_LDP
	      ELSE
	        T2=GF_CUT
	      END IF
	      TMP_STRING=TRIM(ION_ID(ID))//'_F_OSCDAT'
	      CALL GENOSC_V9( ATM(ID)%AXzV_F, ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,ATM(ID)%XzVLEVNAME_F,
	1                 ATM(ID)%ARAD,ATM(ID)%GAM2,ATM(ID)%GAM4,ATM(ID)%OBSERVED_LEVEL,
	1                 T1, ATM(ID)%ZXzV,
	1                 ATM(ID)%XzV_OSCDATE, ATM(ID)%NXzV_F,I,
	1                 'SET_ZERO',T2,GF_LEV_CUT,MIN_NUM_TRANS,L_FALSE,L_FALSE,
	1                 LUIN,LUSCR,TMP_STRING)
	      TMP_STRING=TRIM(ION_ID(ID))//'_F_TO_S'
	      CALL RD_F_TO_S_IDS_V4( ATM(ID)%F_TO_S_XzV, ATM(ID)%INT_SEQ_XzV,
	1           ATM(ID)%XzVLEVNAME_F, ATM(ID)%NXzV_F, ATM(ID)%NXzV,
	1           LUIN,TMP_STRING,SL_OPTION,dE_OPTION,F_TO_S_RD_ERROR)
	      CALL RDPHOT_GEN_V2( ATM(ID)%EDGEXzV_F, ATM(ID)%XzVLEVNAME_F,
	1           ATM(ID)%GIONXzV_F,AT_NO(SPECIES_LNK(ID)),
	1           ATM(ID)%ZXzV, ATM(ID)%NXzV_F,
	1           ATM(ID)%XzV_ION_LEV_ID, ATM(ID)%N_XzV_PHOT,  NPHOT_MAX,
	1           ATM(ID+1)%XzV_PRES,     ATM(ID+1)%EDGEXzV_F, ATM(ID+1)%GXzV_F,
	1           ATM(ID+1)%F_TO_S_XzV,   ATM(ID+1)%XzVLEVNAME_F, ATM(ID+1)%NXzV_F,
	1           SIG_GAU_KMS,FRAC_SIG_GAU,CUT_ACCURACY,ABOVE_EDGE,
	1           XRAYS,ID,ION_ID(ID),LUIN,LUSCR)
              IF(ATM(ID+1)%XzV_PRES) ATM(ID)%GIONXzV_F= ATM(ID+1)%GXzV_F(1)
 	      IF(DIE_AS_LINE .AND. (ATM(ID)%DIE_AUTO_XzV .OR.  ATM(ID)%DIE_WI_XzV) )THEN
	        TMP_STRING='DIE'//TRIM(ION_ID(ID))
	        CALL RDGENDIE_V4( ATM(ID)%XzVLEVNAME_F, ATM(ID)%INDX_XzV,
	1             ATM(ID)%NXzV_F,
	1             EDGEDIE,EINADIE,GUPDIE,
	1             LEVDIE,INDXDIE,SPECDIE,DIENAME, ATM(ID)%GIONXzV_F,
	1             ATM(ID)%DIE_AUTO_XzV, ATM(ID)%DIE_WI_XzV,
	1             ION_ID(ID),LUIN,LUSCR,L_TRUE,TMP_STRING,NMAXDIE,NDIETOT)
	      ELSE IF( ATM(ID)%DIE_AUTO_XzV .OR.  ATM(ID)%DIE_WI_XzV)THEN
	        TMP_STRING='DIE'//TRIM(ION_ID(ID))
	        CALL RD_PHOT_DIE_V1(ID,
	1             ATM(ID)%EDGEXzV_F, ATM(ID)%XzVLEVNAME_F,
	1             ATM(ID)%NXzV_F,    ATM(ID)%GIONXzV_F,
	1             VSM_DIE_KMS, ATM(ID)%DIE_AUTO_XzV, ATM(ID)%DIE_WI_XzV,
	1             ION_ID(ID),LUIN,LUSCR,TMP_STRING)
	      END IF
!
! This is simply to get any data references that are in the collisional data file.
!
	      TMP_STRING=TRIM(ION_ID(ID))//'_COL_DATA'
	      IF(MYPE .EQ. 0)CALL GET_COL_REF(TMP_STRING,LUIN,LUSCR)
	   END IF
	  END DO
	END DO
	IF(F_TO_S_RD_ERROR)THEN
	  WRITE(6,*)' '
	  WRITE(6,*)'There are errors reading in the super level links.'
	  WRITE(6,*)'These need to be fixed before the code can run.'
	  WRITE(6,*)'See F_TO_S_RD_ERRORS for details.'
	  WRITE(6,*)' '
	  STOP
	END IF
!
	CALL ALLOCATE_WSE_ARRAYS_MPI_V1(ND)
!
	IF(MYPE .EQ. 0)THEN
	  DO ID=1,NUM_IONS
	    ROOT(ID)%F_TO_S_XzV=ATM(ID)%F_TO_S_XzV
	    ROOT(ID)%INT_SEQ_XzV=ATM(ID)%INT_SEQ_XzV
	  END DO
	END IF
!
! 
!
! We open a new MODEL file so that the information is at the head of the
! file.
!
	IF(MYPE .EQ. 0)THEN
	  CALL GEN_ASCI_OPEN(LUMOD,'MODEL','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening MODEL in CMFGEN, IOS=',IOS
	    STOP
	  END IF
	  WRITE(6,'(/,A)')' Successfully opened file MODEL'
!
! Output description of model. This is done after the reading of most data
! since some of the model information is read in.
!
	  CALL DATE_TIME(TIME)
	  WRITE(LUMOD,'(//,'' Model Started on:'',20X,(A))')TIME
	  WRITE(LUMOD,
	1       '('' Main program (MPI) last changed on:'',2X,(A))')PRODATE
	  WRITE(LUMOD,'()')
	  FMT='(5X,I8,5X,''!Number of depth points'')'
	  WRITE(LUMOD,FMT)ND
	  FMT='(5X,I8,5X,''!Number of core rays'')'
	  WRITE(LUMOD,FMT)NC
	  FMT='(5X,I8,5X,''!Total number of rays'')'
	  WRITE(LUMOD,FMT)NP
	  FMT='(5X,I8,5X,''!Total number of variables'')'
	  WRITE(LUMOD,FMT)NT
	  FMT='(5X,I8,5X,''!Maximum number of frequencies'')'
	  WRITE(LUMOD,FMT)NCF_MAX
	  FMT='(5X,I8,5X,''!Number of bands'')'
	  WRITE(LUMOD,FMT)NUM_BNDS
!	
	  WRITE(LUMOD,'()')
	  CALL RITE_ATMHD_V4(LUMOD)
!
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      ISPEC=SPECIES_LNK(ID)
	      CALL RITE_ATMDES_V4( ATM(ID)%XzV_PRES, ATM(ID)%NXzV,
	1          ATM(ID)%ZXzV, ATM(ID)%EQXzV, ATM(ID)%XzVLEVNAME_F,
	1          ATM(ID)%NXzV_F, ATM(ID)%GIONXzV_F, ATM(ID)%N_XzV_PHOT,
	1          AT_NO(ISPEC),AT_MASS(ISPEC),ID,ISPEC,LUMOD,ION_ID(ID))
	    END IF
	  END DO
!
! Append VADAT information and atomic data headers to model file.
! Thus only have 1 model file output.
!
	  REWIND(LUSCR)
	  IOS=0
	  TEMP_CHAR=FORMFEED
	  DO WHILE(IOS .EQ. 0)
	    I=ICHRLEN(TEMP_CHAR)
	    IF(I .GT. 0)THEN
	      WRITE(LUMOD,'(A)')TEMP_CHAR(1:I)
	    ELSE
	      WRITE(LUMOD,'()')
	    END IF
	    READ(LUSCR,'(A)',IOSTAT=IOS)TEMP_CHAR
	  END DO
	  CLOSE(UNIT=LUMOD)
	END IF
!
! Finished all data read, so can close LUMOD (output descriptor).
!
	CLOSE(UNIT=LUSCR,STATUS='DELETE') 
! 
!
! Set the vector Z_POP to contain the ionic charge for each species.
!
	DO I=1,NT
	  Z_POP(I)=0.0_LDP
	END DO
!
	DO ID=1,NUM_IONS-1
	  CALL SET_Z_POP(Z_POP, ATM(ID)%ZXzV, ATM(ID)%EQXzV,
	1              ATM(ID)%NXzV, NT, ATM(ID)%XzV_PRES)
	END DO
!
! Store atomic masses in vector of LENGTH NT for later use by line
! calculations. G_ALL and LEVEL_ID  are no longer used due to the use
! of super levels.
!
	AMASS_ALL(1:NT)=0.0_LDP
	DO ID=1,NUM_IONS-1
	  IF(ATM(ID)%XzV_PRES)AMASS_ALL( ATM(ID)%EQXzV: ATM(ID)%EQXzV+ATM(ID)%NXzV-1)=
	1         AT_MASS(SPECIES_LNK(ID))
	END DO
!
! 
!
! Define the average energy of each super level. At present this is
! depth independent, which should be adequate for most models.
! This average energy is used to scale the line cooling rates in
! the radiative equilibrium equation so that is more consistent
! with the electron cooling rate. The need for this scaling
! arises when levels within a super level have a 'relatively large'
! energy separation, and the dominat rates are scattering.
!
	AVE_ENERGY(:)=0.0_LDP
	DO ID=1,NUM_IONS-1
	   CALL AVE_LEVEL_ENERGY(AVE_ENERGY, ATM(ID)%EDGEXzV_F,
	1         ATM(ID)%GXzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%EQXzV,
	1         ATM(ID)%NXzV,   ATM(ID)%NXzV_F, NT, ATM(ID)%XzV_PRES)
	END DO
!
! 
!
! Check to see if old model. If so, read in R,V, SIGMA and POPS arrays.
! If not, set NEWMOD to .TRUE. We also check the format of the file,
! in case we are revising the R grid.
!
	NLBEGIN=0		! Initialize for lines.
	IREC=0                  ! Get last iteration
	CALL SCR_READ_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,LAST_NG,
	1                 WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
        IF( (REVISE_R_GRID .OR. DO_HYDRO) .AND. NEWMOD)THEN
	  WRITE_RVSIG=.TRUE.
	ELSE IF(NEWMOD)THEN
	  WRITE_RVSIG=.FALSE.
	ELSE IF(REVISE_R_GRID .OR. DO_HYDRO)THEN
	   IF(.NOT. WRITE_RVSIG)THEN
	     WRITE(LUER,*)'Error in CMFGEN_SUB with SCRTEMP'
	     WRITE(LUER,*)'Inconsistent format request: RVSIG must be written for each iteration'
	     WRITE(LUER,*)'Restart a fresh model or use REWRITE_SCR to correct file format'
	     STOP
	   END IF
	END IF
	IF(NEWMOD)THEN
	  IF(MYPE .EQ. 0)THEN
	    WRITE(LUER,*)'Starting a new model.'
	    WRITE(LUER,*)'*_IN files will be used to start model'
	    WRITE(LUER,*)'Setting ITS_DONE keyword in HYDRO_DEFAULTS to 0'
	    IF(DO_HYDRO .AND. .NOT. SN_MODEL)THEN
              CALL UPDATE_KEYWORD(IZERO,'[ITS_DONE]','HYDRO_DEFAULTS',L_TRUE,L_TRUE,LUIN)
	    END IF
	  END IF
	ELSE
!
! Generally RP and RMAX will be consistent with R(ND) and R(1). However they
! will need to be reset if we have rewound a model (i.e., changed POINT1 to
! use older iterations) when computing the hyrostatic structure.
!
	  IF(.NOT. SN_HYDRO_MODEL .AND. (RP .NE. R(ND) .OR. R(1) .NE. RMAX))THEN
	    IF(MYPE .EQ. 0)THEN
	      WRITE(LUER,*)'Warning: RP and RMAX in CMFGEN are inconsistent'
	      WRITE(LUER,*)'with values in SCRTEMP. This inconsistency mayhave occured'
              WRITE(LUER,*)'if you have rewound (changd POINT1) a model with DO_HYDRO=T.'
              WRITE(LUER,*)'Alternatively, it could be a rounding error.'
              WRITE(LUER,*)'  RP=',RP,  ' R(ND)=',R(ND)
              WRITE(LUER,*)'RMAX=',RMAX,'  R(1)=',R(1)
	      WRITE(LUER,*)'Please revise VADAT (or rewind SCRTEMP) to set consistency'
	      T1=ABS(R(ND)/RP-1.0_LDP)/(R(ND-1)-R(ND))
	      T2=ABS(R(1)/RMAX-1.0_LDP)/(R(1)-R(2))
	      IF(T1 .GT. 1.0E-03_LDP .OR. T2 .GT. 1.0E-03_LDP)THEN
	        WRITE(6,*)'ABS(R(ND)/RP-1.0_LDP)/(R(ND-1)-R(ND))',T1
	        WRITE(6,*)'ABS(R(1)/RMAX-1.0_LDP)/(R(1)-R(2))',T2
	        WRITE(6,*)'RMAX/RP=',R(1)/R(ND)
	        WRITE(6,*)'Stopping code so error can be fixed'
	        CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	        STOP
	      END IF
	      RP=R(ND); RMAX=R(1)
	    END IF
	  END IF
	END IF
!
! Now does accurate flux calculation for a single iteration provided not a new model.
!
	IF(.NOT. NEWMOD)MAXCH=0.0_LDP
!
!		' OLD MODEL '
!
	IF(.NOT. NEWMOD)THEN
!
! Convert back from POPS array to individual matrices.
!
	  DO ID=1,NUM_IONS-1
	    CALL POPTOION(POPS, ATM(ID)%XzV, ATM(ID)%DXzV,ED,T,
	1            ATM(ID)%EQXzV, ATM(ID)%NXzV,
	1            NT,DST, DEND, ND, ATM(ID)%XzV_PRES)
	  END DO
	  IF(MYPE .EQ. 0)THEN
	    DO ID=1,NUM_IONS-1
	      CALL POPTOION(POPS, ROOT(ID)%XzV, ROOT(ID)%DXzV,ED,T,
	1            ATM(ID)%EQXzV, ATM(ID)%NXzV,
	1            NT, IONE, ND, ND, ATM(ID)%XzV_PRES)
	    END DO
	  END IF
	  ED(1:ND)=POPS(NT-1,1:ND); T(1:ND)=POPS(NT,1:ND)
!
! We have now stored the revised populations back in their individual
! storage locations. For some species we have 2 atomic models. For these
! species we need to take the super-level populations and compute:
!
! 1. The LTE population off all level ls in the FULL atom.
! 2. The population off all levels in the FULL atom.
! 3. The LTE population off all super-levels.
!
! This is done by the following include statement, which  is a sequence of
! calls to the routine SUP_TO_FULL.
!
! Compute the ion population at each depth.
! These are required when evaluation the occupation probabilities, and hence
! the LTE populations. Not that Z_POP is effectivy integer, thus we
! use 0.01 as check.
!
	  DO J=1,ND
	    POPION(J)=0.0_LDP
	    DO I=1,NT
	      IF(Z_POP(I) .GT. 0.01_LDP)POPION(J)=POPION(J)+POPS(I,J)
	    END DO
	  END DO
!	  J=155+MYPE; CALL WR2D(POPION,ND,IONE,'POPION',J); FLUSH(UNIT=J)
!
!	  J=160+MYPE; CALL WR2D(POPS,NT,ND,'POPS',J); FLUSH(UNIT=J)
	  CALL SUP_TO_FULL_V4(POPS,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
!
	ELSE
!
! Compute R, V, and SIGMA separately from DC read so that can evaluate
! POPHE etc, and angle quadrature weights. These are required in NEWMODEL
! section if iterating on the initial temperature structure.
!
	  IF(VELTYPE .EQ. 1)THEN
	    CALL STARNEW(R,V,SIGMA,RMAX,RP,RN,VRP,VINF,EPPS1,GAMMA1
	1    ,RP2,RN2,VRP2,VINF2,EPPS2,GAMMA2,ND,TA,TB,TC)
	  ELSE IF(VELTYPE .EQ. 2)THEN
	    CALL STARFIN(R,V,SIGMA,RMAX,RP,RN,VRP,VINF,EPPS1,GAMMA1
	1   ,RP2,RN2,VRP2,VINF2,EPPS2,GAMMA2,ND,TA,TB,TC)
	  ELSE IF(VELTYPE .EQ. 3 .OR. VELTYPE .EQ. 6)THEN
	    CALL STARPCYG_V3(R,V,SIGMA,RMAX,RP,
	1             SCL_HT,VCORE,VPHOT,VINF1,V_BETA1,V_EPPS1,
	1             VINF2,V_BETA2,V_EPPS2,
	1             N_OB_INS,CONS_FOR_R_GRID,EXP_FOR_R_GRID,
	1             ND,TA,TB,TC,RDINR,LUIN)
          ELSE IF(VELTYPE .EQ. 4)THEN
            CALL STARRAVE(R,V,SIGMA,ND,LUIN,RMAX,RP)
         ELSE IF(VELTYPE .EQ. 7)THEN
	    CALL RD_RV_FILE_V2(R,V,SIGMA,RMAX,RP,VINF,LUIN,ND,VEL_OPTION,NUM_V_OPTS)
         ELSE IF(VELTYPE .EQ. 10)THEN
	    CALL RV_SN_MODEL_V2(R,V,SIGMA,RMAX,RP,VCORE,V_BETA1,RDINR,LUIN,ND)
	 ELSE IF(VELTYPE .EQ. 11)THEN
	   CALL SET_RV_HYDRO_MODEL_V3(R,V,SIGMA,RMAX,RP,RMAX_ON_RCORE,SN_AGE_DAYS,
	1                      PURE_HUBBLE_FLOW,N_IB_INS,N_OB_INS,RDINR,ND,LUIN)
           VINF=V(1)
	 ELSE IF(VELTYPE .EQ. 12)THEN
	   CALL RV_SN_MODEL_SNIIN(R,V,SIGMA,RMAX,RP,VCORE,V_BETA1,RDINR,LUIN,ND)
	   VINF=V(1)
	 ELSE
	   WRITE(LUER,*)'Invalid Velocity Law'
	   STOP
	 END IF
!
	  IF(SN_HYDRO_MODEL .AND. MYPE .EQ. 0)THEN
	    T1=RMAX/RP
	    CALL UPDATE_KEYWORD(RP,'[RSTAR]','VADAT',L_TRUE,L_FALSE,LUIN)
	    CALL UPDATE_KEYWORD(T1,'[RMAX]','VADAT',L_FALSE,L_TRUE,LUIN)
	    WRITE(LUER,*)'Updated RP and RMAX in VADAT as new SN hydro model'
	  END IF
	END IF
!
	IF(VINF .EQ. 0.0_LDP)VINF=V(1)
!
! 
!
! Compute CLUMP_FAC(1:ND) which allow for the possibility that the wind is
! clumped. At the sime time, we compute the vectors which give the density,
! the atom density, and the species density at each depth.
!
! This routine also computes VTURB_VEC.
!
	CALL SET_ABUND_CLUMP(MEAN_ATOMIC_WEIGHT,ABUND_SUM,LUER,ND)
!
! Compute profile frequencies such that for the adopted doppler
! velocity the profile ranges from 5 to -5 doppler widths.
! This section needs to be rewritten if we want the profile to
! vary with depth.
!
! ERF is used in computing the Sobolev incident intensity at the
! outer boundary. ERF = int from "x" to "inf" of -e(-x^2)/sqrt(pi).
! Note that ERF is not the error function. ERF is related to the
! complementary error function by ERF =-0.5D0 . erfc(X).
! S15ADF is a NAG routine which returns erfc(x).
!
! The incident Sobolev intensity is S[ 1.0-exp(tau(sob)*ERF) ]
! NB -from the definition, -1<erf<0 .
!
	T1=4.286299E-05_LDP*SQRT( TDOP/AMASS_DOP + (VTURB/12.85_LDP)**2 )
	J=0
!	DO I=1,NLF
!	  ERF(I)=-0.5_LDP*S15ADF(PF(I),J)
!	  PF(I)=PF(I)*T1
!	END DO
        VDOP_VEC(1:ND)=12.85_LDP*SQRT( TDOP/AMASS_DOP + (VTURB/12.85_LDP)**2 )
!
	IF(GLOBAL_LINE_PROF(1:4) .EQ. 'LIST')THEN
	  CALL RD_STRK_LIST(LUIN)
	END IF
!
! Compute the frequency grid for CMFGEN. Routine also allocates the vectors
! needed for the line data, sets the line data, and puts the line data into
! numerical order.
!
	IF(MYPE .EQ. 0)THEN
	  WRITE(6,*)'Calling routine to set the frequency grid.'; FLUSH(UNIT=6)
	END IF
	CALL SET_FREQUENCY_GRID_V2(NU,FQW,LINES_THIS_FREQ,NU_EVAL_CONT,
	1               NCF,NCF_MAX,N_LINE_FREQ,ND,
	1               OBS_FREQ,OBS,N_OBS,LUIN,IMPURITY_CODE)
	IF(MYPE .EQ. 0)THEN
	  WRITE(6,*)'Frequency grid set'; FLUSH(UNIT=6)
	END IF
!
! 
!
	IF(ACCURATE)THEN
	  I=ND-DEEP
	  IF(INTERP_TYPE .NE. 'LOG')THEN
	    WRITE(LUER,*)'Error in CMFGEN_SUB'
	    WRITE(LUER,*)'The INTERP_TYPE currently implemented is LOG'
	    STOP
	  END IF
	  CALL REXT_COEF_V2(REXT,COEF,INDX,NDEXT,R,POS_IN_NEW_GRID,
	1         ND,NPINS,L_TRUE,I,ST_INTERP_INDX,END_INTERP_INDX)
	  TA(1:ND)=1.0_LDP	!TEXT not required, T currently zero
	  CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,
	1        V,TA,SIGMA,ND)
!
          VDOP_VEC_EXT(1:NDEXT)=12.85_LDP*SQRT( TDOP/AMASS_DOP + (VTURB/12.85_LDP)**2 )
	  CALL SET_POP_FOR_TWOJ(POS_IN_NEW_GRID,EDD_CONT_REC,LU_EDD,NDEXT)
	ELSE
	  CALL SET_POP_FOR_TWOJ(POS_IN_NEW_GRID,EDD_CONT_REC,LU_EDD,ND)
	END IF
	
!
! Need to calculate impact parameters, and angular quadrature weights here
! as these may be required when setting up the initial temperature
! distribution of the atmosphere (i.e. required by JGREY).
!
	CALL SET_ANG_QW_V2(R,NC,ND,NP,REXT,NCEXT,NDEXT,NPEXT,
	1                  R_PNT_SRCE,NC_PNT_SRCE,TRAPFORJ,ACCURATE)
!
! Allocate memory for opacities.
!
        CALL INIT_OPAC_MOD(DST,DEND,ND,NT,L_TRUE)
! 
!
!		'NEW MODEL'
!
! Read in old estimates for the departure coefficents for the new model.
! T and ED are also estimated.
!
	IF(NEWMOD)THEN
	  NITSF=0
	  IREC=0
	  LAST_NG=-1000  			!Must be -1000
	  NEXT_NG=1000				!Must be initialized to 1000
	  CALL SET_NEW_MODEL_ESTIMATES(POPS,Z_POP,NU,NU_EVAL_CONT,FQW,
	1            LUER,LUIN,NC,ND,NP,NT,NCF,N_LINE_FREQ,MAX_SIM)
	END IF
!
! VEXT and SIGMAEXT have already been computed. We need TEXT for
! convolving J with the electron scattering redistribution function.
!
	IF(ACCURATE)THEN
	  CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,
	1        V,T,SIGMA,ND)
	ELSE
	  TEXT(1:ND)=T(1:ND)
	END IF
!
! 
!
! Set TSTAR which is required by TOTOPA_JILA. Its precise value is
! irrelevant (provided >0) as we always adopt the diffusion
! approximation.
!
	TSTAR=T(ND)
!
! 
! Section allows the user to read in a modified solution vector.
! Useful for assisting convergence when the model is experiencing
! difficulty converging.
!
	IF(NUM_ITS_TO_DO .EQ. 0)THEN
	  IF(RDINSOL)THEN
!
	    CALL LOCSOLUT(POPS,SOL,TA,30,NT,ND)
	    DO ID=1,NUM_IONS-1
	      CALL POPTOION(POPS, ATM(ID)%XzV, ATM(ID)%DXzV, ED,T,
	1            ATM(ID)%EQXzV, ATM(ID)%NXzV, NT, DST, DEND, ND,
	1            ATM(ID)%XzV_PRES)
	    END DO
!
! This include block also compute the LTE populations of the FULL model atom,
! and the SUPER level model atom.
!
	    CALL SUP_TO_FULL_V4(POPS,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
!	    INCLUDE 'SUP_TO_FULL_V4.INC'
	  END IF
!
! Write pointer file and output data necessary to begin a new
! iteration.
!
	  IF(RDINSOL .OR. NEWMOD)THEN
	    MAIN_COUNTER=NITSF+1
	    IF(MYPE .EQ. 0)THEN
	      CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
	1                   LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	      CALL MPI_BCAST(IREC,IONE,MPI_INTEGER,IZERO,MPI_COMM_WORLD,IERR)
	    END IF
	  END IF
	  LST_ITERATION=.TRUE.
	  CALL TUNE(IONE,'GIT')
	  GOTO 9999			!End (write out POPS.)
	END IF
!
	ARTIFICIAL_HEAT_TERM=0.0_LDP
	dE_WORK=0.0_LDP
	dE_RAD_DECAY=0.0_LDP 			!For non-SN models.
	dE_SHOCK_POWER= 0.0_LDP
	dE_XRAY_0P1=0.0_LDP 
        dE_XRAY_1kev=0.0_LDP
        dE_XRAY_TOT=0.0_LDP
	DEP_RAD_EQ=0.0_LDP
!
	IF(MYPE .EQ. 0)WRITE(6,'(/,A)')' About to start main iteration section'
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	INCLUDE 'main_it_section.f'
! 
9999	CONTINUE
!
! NB - Need GE because of NG acceleration.
!
	IF(LST_ITERATION)THEN
!
	  CALL TUNE(ITWO,'GIT')
	  IF(MYPE .EQ. 0)CALL TUNE(ITHREE,' ')
!
!	  IF(MYPE .EQ. 0)THEN
!	    CALL CHECK_IONS_PRESENT(ND,NUM_IONS)
!	  END IF
!
! Make sure root has all the needed information.
!
!	  CALL SCATTER_XZV_F_AND_IONS(ROOT)
!
! This file is a direct access file and contains the models
! output (i.e. T,density,population levels etc). No longer
! assumes that RVTJ exists (28-Sep-1990). Altered (1-Dec-1991) to be
! an asci file for CRAY/VAX compatibility. NB. Obs may be ZERO if program
! was inadvertently started during the last iteration.
!
	  CALL DATE_TIME(TIME)
!
	  TA=0.0_LDP; TB=0.0_LDP; TB(DST:DEND)=dE_RAD_DECAY(DST:DEND)
	  CALL MPI_REDUCE(TB,TA,ND,MY_MPI_DP,MPI_SUM,IZERO,MPI_COMM_WORLD,IERR)
	  CALL GATHER_SELF_VEC_MPI_V1(J_INT,DST,DEND,ND)
	  CALL GATHER_SELF_VEC_MPI_V1(H_INT,DST,DEND,ND)
	  CALL GATHER_SELF_VEC_MPI_V1(K_INT,DST,DEND,ND)
	  IF(MYPE .EQ. 0)THEN
	    CALL GEN_ASCI_OPEN(LU_POP,'RVTJ','UNKNOWN',' ',' ',IZERO,IOS)
	    FORMAT_DATE='15-Aug-2019'
	    WRITE(LU_POP,'(1X,A,T30,A)')'Output format date:',FORMAT_DATE
	    WRITE(LU_POP,'(1X,A,T30,A)')'Completion of Model:',TIME
	    WRITE(LU_POP,'(1X,A,T30,A)')'Program Date:',PRODATE
!
	    WRITE(LU_POP,'(1X,A,T30,I5)')'ND:',ND
	    WRITE(LU_POP,'(1X,A,T30,I5)')'NC:',NC
	    WRITE(LU_POP,'(1X,A,T30,I5)')'NP:',NP
	    WRITE(LU_POP,'(1X,A,T30,I6)')'NCF:',N_OBS
!
	    WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')'Mdot(Msun/yr):',RMDOT/3.02286D+23
	    WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')'L(Lsun):',LUM
	    WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')'H/He abundance:',AT_ABUND(1)
	    WRITE(LU_POP,'(1X,A,T30,L1)')'Was T fixed?:',RD_FIX_T
	    WRITE(LU_POP,'(1X,A,T30,A)')'Species naming convention:',NAME_CONVENTION
!
	    WRITE(LU_POP,'(A)')' Radius (10^10 cm)'
	    WRITE(LU_POP,'(1X,8ES18.10)')R
	    WRITE(LU_POP,'(A)')' Velocity (km/s)'
	    WRITE(LU_POP,'(1X,8ES18.10)')V
	    WRITE(LU_POP,'(A)')' dlnV/dlnr-1'
	    WRITE(LU_POP,'(1X,1P8E16.7)')SIGMA
	    WRITE(LU_POP,'(A)')' Electron density'

	    WRITE(LU_POP,'(1X,1P8E16.7)')ED
	    WRITE(LU_POP,'(A)')' Temperature (10^4K)'
	    WRITE(LU_POP,'(1X,1P8E16.7)')T
	    WRITE(LU_POP,'(A)')' Grey temperature (10^4K)'
	    WRITE(LU_POP,'(1X,1P8E16.7)')TGREY
!
	    WRITE(LU_POP,'(A)')' Heating: radioactive decay (ergs/cm^3/s) '
	    WRITE(LU_POP,'(1X,1P8E16.7)')TA(1:ND)
!
	    WRITE(LU_POP,'(A)')' Rosseland Mean Opacity'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(ROSS_MEAN(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Flux Mean Opacity'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(FLUX_MEAN(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Planck Mean Opacity'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(PLANCK_MEAN(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Absorption Mean Opacity'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(ABS_MEAN(I),I=1,ND); FLUSH(LU_POP)
!
	    WRITE(LU_POP,'(A)')' J moment of radiation field'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(J_INT(I)/LUM_SCL_FAC,I=1,ND)
	    WRITE(LU_POP,'(A)')' H moment of radiation field'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(H_INT(I)/LUM_SCL_FAC,I=1,ND)
	    WRITE(LU_POP,'(A)')' K moment of radiation field'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(K_INT(I)/LUM_SCL_FAC,I=1,ND); FLUSH(LU_POP)
!
! Compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	    DO J=1,ND
	      POPION(J)=0.0_LDP
	      DO I=1,NT
	        IF(Z_POP(I) .GT. 0.01_LDP)POPION(J)=POPION(J)+POPS(I,J)
	      END DO
	    END DO
	    WRITE(LU_POP,'(A)')' Atom Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(POP_ATOM(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Ion Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(POPION(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Mass Density (gm/cm^3)'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(DENSITY(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Clumping Factor'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(CLUMP_FAC(I),I=1,ND)
!
	    WRITE(LU_POP,'(A)')' Hydrogen Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')POP_SPECIES(1:ND,1)
	    WRITE(LU_POP,'(A)')' Helium Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')POP_SPECIES(1:ND,2)
!
	    IF(ROOT(1)%XzV_PRES)THEN
	      WRITE(LU_POP,'(A)')TRIM(ION_ID(1))//' populations'
	      WRITE(LU_POP,'(A,T30,I3)')' Number of Hydrogen levels:',ATM(1)%NXzV_F
	      WRITE(LU_POP,'(A,T30,A)')TRIM(ION_ID(1))//
	1                ' oscillator date:',ATM(1)%XzV_OSCDATE
	      WRITE(LU_POP,'(1X,1P8E16.7)')ROOT(1)%XzV_F
	      WRITE(LU_POP,'(A)')' DHYD population'
	      WRITE(LU_POP,'(1X,1P8E16.7)')ROOT(1)%DXzV_F
	    END IF
!
	    IF(ROOT(3)%XzV_PRES)THEN
	      WRITE(LU_POP,'(A)')' HeI populations'
	      WRITE(LU_POP,'(A,T30,I3)')' Number of Helium I levels:',ATM(3)%NXzV_F
	      WRITE(LU_POP,'(A,T30,A)')TRIM(ION_ID(3))//
	1                ' oscillator date:',ATM(3)%XzV_OSCDATE
	      WRITE(LU_POP,'(1X,1P8E16.7)')ROOT(3)%XzV_F
	      WRITE(LU_POP,'(A)')' DHeI population'
	      WRITE(LU_POP,'(1X,1P8E16.7)')ROOT(3)%DXzV_F
	    END IF
!
	    IF(ATM(4)%XzV_PRES)THEN
	      WRITE(LU_POP,'(A)')' He2 populations'
	      WRITE(LU_POP,'(A,T30,I3)')' Number of Helium II levels:',ATM(4)%NXzV_F
	      WRITE(LU_POP,'(A,T30,A)')TRIM(ION_ID(4))//
	1                ' oscillator date:',ATM(1)%XzV_OSCDATE
	      WRITE(LU_POP,'(1X,1P8E16.7)')ROOT(4)%XzV_F
	      WRITE(LU_POP,'(A)')' DHe2 population'
	      WRITE(LU_POP,'(1X,1P8E16.7)')ROOT(4)%DXzV_F
	    END IF
!
	    CLOSE(UNIT=LU_POP)
	  END IF
!
	  IF(TREAT_NON_THERMAL_ELECTRONS)THEN
	    CALL GEN_ASCI_OPEN(LU_POP,'NON_THERM_COOL','UNKNOWN',' ',' ',IZERO,IOS)
	    WRITE(LU_POP,'(1X,A,T30,I5)')'ND:',ND
	    WRITE(LU_POP,*)''
	    DO ID=1,NUM_IONS
	      IF(ROOT(ID)%XzV_PRES)THEN
	        WRITE(LU_POP,'(A)')TRIM(ION_ID(ID))//' ionization cooling (eV)'
	        WRITE(LU_POP,'(1X,1P8E16.7)')ROOT(ID)%NT_ION_CXzV
	        WRITE(LU_POP,*)''
	        WRITE(LU_POP,'(A)')TRIM(ION_ID(ID))//' excitation cooling (eV)'
	        WRITE(LU_POP,'(1X,1P8E16.7)')ROOT(ID)%NT_EXC_CXzV
	        WRITE(LU_POP,*)''
	      END IF
	    END DO
	  CLOSE(UNIT=LU_POP)
	  END IF
!
	  IF(SN_HYDRO_MODEL)THEN
	    CALL OUT_SN_POPS_V3('SN_HYDRO_FOR_NEXT_MODEL',SN_AGE_DAYS,USE_OLD_MF_OUTPUT,ND,LUMOD)
	  END IF
! 
!
	  IF(MYPE .EQ. 0)THEN
	    FORMAT_DATE='27-JAN-1992'
	    DO ISPEC=1,NUM_SPECIES
	      IF(POP_SPECIES(ND,ISPEC) .NE. 0)THEN
	        TMP_STRING='POP'//TRIM(SPECIES(ISPEC))
	        CALL GEN_ASCI_OPEN(LU_POP,TMP_STRING,'UNKNOWN',' ',' ',IZERO,IOS)
	        WRITE(LU_POP,'(1X,A,T30,A)')'Output format date:',FORMAT_DATE
	        WRITE(LU_POP,'(1X,A,T30,A)')'Completion of Model:',TIME
	        WRITE(LU_POP,'(1X,A,T30,I5)')'ND:',ND
	        WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')
	1              TRIM(SPECIES(ISPEC))//'/He abundance:',AT_ABUND(ISPEC)
     	        WRITE(LU_POP,'(1X,1P8E16.7)')POP_SPECIES(1:ND,ISPEC)
	        DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	          CALL RITE_ASC( ATM(ID)%XzV_PRES, ROOT(ID)%XzV_F,ROOT(ID)%DXzV_F,
	1                ATM(ID)%NXzV_F, ND,
	1                ATM(ID)%XzV_OSCDATE,TRIM(ION_ID(ID)),LU_POP)
	        END DO
	      END IF
	    END DO
	  END IF
! 
!
! Write out departure coefficients to ASCI file.
! NB - 1 refers to dimension of DHYD (i.e. DHYD(1,nd)
!      1 refers to format for output.
!      1,NHY - For use with HeI.
!
	IF(MYPE .EQ. 0)THEN
	  CALL EVAL_ROOT_LTE_MPI_V1(DO_LEV_DISSOLUTION,ND)
	END IF
!
! GAM_SPECIES refers to the number of electrons arising from each species (eg
! carbon).
!
	  IF(MYPE .EQ. 0)THEN
	    DO ISPEC=1,NUM_SPECIES
	      GAM_SPECIES(1:ND,ISPEC)=0.0_LDP
	      FIRST=.TRUE.
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        IF( ROOT(ID)%XzV_PRES)THEN
	          CALL UPDATE_GAM( GAM_SPECIES(1,ISPEC),
	1              ROOT(ID)%XzV_F, ROOT(ID)%DXzV_F, ATM(ID)%ZXzV,
	1              ATM(ID)%NXzV_F,ND,
	1              ATM(ID+1)%XzV_PRES,FIRST)
	          TMP_STRING=TRIM(ION_ID(ID))//'OUT'
	          CALL WRITEDC_V3( ROOT(ID)%XzV_F, ROOT(ID)%LOG_XzVLTE_F,
	1              ATM(ID)%NXzV_F, ROOT(ID)%DXzV_F,IONE,
	1              R,T,ED,V,CLUMP_FAC,LUM,ND,
	1              TRIM(TMP_STRING),'DC',IONE)
	        END IF
	      END  DO
	    END DO
!
! We only output GAM when it is not zero.
!
	    CALL RITE_GAM_HEAD(R,ED,T,ND,LUIN,'GAMMAS')
	    DO ISPEC=1,NUM_SPECIES
	      CALL RITE_GAM_V2(POP_SPECIES(1,ISPEC),GAM_SPECIES(1,ISPEC),
	1                       AT_NO(ISPEC),SPECIES(ISPEC),ND,LUIN)
	    END DO
	    CLOSE(LUIN)
!
! If we have a time variability model, we output the model so it can be used by the
! next model in the time sequence. IREC is the ouput record, and is thus simply
! the current TIME_SEQ_NO.
!
	    IF(SN_MODEL)THEN
	      CALL WRITE_SEQ_TIME_FILE_MPI_V1(SN_AGE_DAYS,ND,LUSCR)
	      IF(MYPE .EQ. 0)WRITE(6,*)'Successfully output OLD_MODEL_DATA'
	    END IF
!
	    CLOSE(UNIT=LUER)
	    CLOSE(UNIT=LU_SE)
!
	    IF(TREAT_NON_THERMAL_ELECTRONS .AND. WRITE_RATES)THEN
	      CALL WRITE_NON_THERM_MPI_V1(dE_RAD_DECAY,NT,ND,DEC_NRG_SCL_FAC)
	    END IF
	  END IF
!
	  RETURN
	ELSE
!
	INCLUDE 'set_next_iteration.f'
!
	END IF
!
! 
	END
