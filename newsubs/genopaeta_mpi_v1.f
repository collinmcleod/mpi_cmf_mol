!
! Subroutine to compute the contribution to the opacity AND emissivity
! by FREE-FREE and BOUND-FREE processes for a general ion. The
! contribution is added directly to the opacity CHI and emissivity
! ETA.
!
	SUBROUTINE GENOPAETA_MPI_V1(ID,CHI,ETA,NU,
	1              HN,HNST,LOG_HNST,EDGE,GION,ZION,N,
	1              DI,LOG_DIST,N_DI,PHOT_ID,ION_LEV,
	1              ED,T,EMHNUKT,IONFF,DST,DEND,ND,LST_DEPTH_ONLY)
	USE SET_KIND_MODULE
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
!
	INTEGER ID,N,N_DI,ND
	INTEGER DST,DEND
	INTEGER DPTH_INDX
	LOGICAL IONFF,LST_DEPTH_ONLY
	LOGICAL KEEP_PHOT
!
! Constants for opacity etc.
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	REAL(KIND=LDP) CHI(ND)			!Opacity
	REAL(KIND=LDP) ETA(ND)			!Emissivity
!
! Large Model Atom Populations.
!
	REAL(KIND=LDP) HN(N,DST:DEND)
	REAL(KIND=LDP) HNST(N,DST:DEND)
	REAL(KIND=LDP) LOG_HNST(N,DST:DEND)
	REAL(KIND=LDP) EDGE(N)
!
! Ion populations. These populations should refer to the small model atoms.
! (i.e. the model atom with super levels_
!
	REAL(KIND=LDP) DI(N_DI,DST:DEND)
	REAL(KIND=LDP) LOG_DIST(N_DI,DST:DEND)
!
	REAL(KIND=LDP) T(ND)			!Temperature (K)
	REAL(KIND=LDP) ED(ND)			!Electron density
	REAL(KIND=LDP) EMHNUKT(ND)		!EXP(-hv/kT)
	REAL(KIND=LDP) NU			!Frequency (10^15 Hz)
	REAL(KIND=LDP) ZION			!Charge on resultiong ion.
	REAL(KIND=LDP) GION			!Charge on resultiong ion.
!
	INTEGER PHOT_ID		!Photoionization ID (path)
	INTEGER ION_LEV		!Target level for ionizations in ion.
!
! Vectors to save computational effort.
!
	REAL(KIND=LDP) YDIS(DST:DEND)			!Constant for computing level dissolution/
	REAL(KIND=LDP) XDIS(DST:DEND)			!Constant for computing level dissolution/
	REAL(KIND=LDP) DIS_CONST(N)		!Constant appearing in dissolution formula.
	REAL(KIND=LDP) ALPHA_VEC(N)		!Photionization cross-section
	REAL(KIND=LDP) TMP_CHI(N)		!Photionization cross-section
	REAL(KIND=LDP) TMP_ETA(N)		!Photionization cross-section
!
	REAL(KIND=LDP) GFF_VAL		!g(ff) as a function of depth
	REAL(KIND=LDP) COR_FAC		!Factor to convert HNST for ION_LEV
	REAL(KIND=LDP) LOG_COR_FAC		!Factor to convert HNST for ION_LEV
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
! Local constants.
!
	INTEGER LOC_DST,LOC_DEND
	INTEGER I,K,K_ST,ND_LOC,NO_NON_ZERO_PHOT
	REAL(KIND=LDP) ALPHA,TCHI1,TETA1,TETA2
	REAL(KIND=LDP) T1,T2,ZION_CUBED,NEFF
	REAL(KIND=LDP) GFF
	EXTERNAL GFF
	INTEGER, PARAMETER :: IONE=1
!
!^L
!
	DO DPTH_INDX=DST,DEND
!
! Compute the photo-ionization cross-sections for all levels.
!
	  IF(DPTH_INDX .EQ. DST)THEN
	    IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	      CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU,EDGE,N,PHOT_ID,L_TRUE)
	    ELSE
	      CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU,EDGE,N,PHOT_ID,L_FALSE)
	    END IF
	    NO_NON_ZERO_PHOT=COUNT(ALPHA_VEC .GT. 0.0_LDP)
	    IF(NO_NON_ZERO_PHOT .EQ. 0)RETURN
!
! DIS_CONST is the constant K appearing in the expression for level dissolution.
! A negative value for DIS_CONST implies that the cross-section is zero.
!
	    DIS_CONST(1:N)=-1.0_LDP
	    IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	      ZION_CUBED=ZION*ZION*ZION
	      DO I=1,N
	        IF(NU .LT. EDGE(I) .AND. ALPHA_VEC(I) .NE. 0)THEN
	          NEFF=SQRT(3.289395_LDP*ZION*ZION/(EDGE(I)-NU))
	          IF(NEFF .GT. 2*ZION)THEN
	            T1=MIN(1.0_LDP,16.0_LDP*NEFF/(1+NEFF)/(1+NEFF)/3.0_LDP)
	            DIS_CONST(I)=( T1*ZION_CUBED/(NEFF**4) )**1.5_LDP
	          END IF
	        END IF
	      END DO
	    END IF
	  END IF
!
!
! Add in free-free contribution. Because SN can be dominated by elements other
! than H and He, we now sum over all levels. To make sure that we only do this
! one, we only include the FREE-FREE contribution for the ion when PHOT_ID is one.
!
	  IF(ZION .EQ. 0.0_LDP)THEN
	    I=7				!Used for IO
	    K=DPTH_INDX
	    COR_FAC=DI(1,K)
	    CALL DO_H0_FF(ETA(K),CHI(K),COR_FAC,ED(K),T(K),EMHNUKT(K),NU,I,IONE)
	  ELSE IF(IONFF .AND. PHOT_ID .EQ. 1)THEN
!
! Compute free-free gaunt factors. Replaces call to GFF in following DO loop.
!
	    K=DPTH_INDX
	    GFF_VAL=GFF(NU,T(K),ZION)
	    IF(ION_LEV .EQ. 1)THEN
	      CALL FF_RES_GAUNT(GFF_VAL,NU,T(K),ID,GION,ZION,IONE)
	    END IF
	  END IF
! We use COR_FAC as a temporary vector containing the sum of all level populations in
! the ion at each depth.
!
	  K=DPTH_INDX
	  COR_FAC=SUM(DI(:,K),1)
	  TCHI1=CHIFF*ZION*ZION/(NU*NU*NU)
	  TETA1=CHIFF*ZION*ZION*TWOHCSQ
	  ALPHA=ED(K)*COR_FAC*GFF_VAL/SQRT(T(K))
	  CHI(K)=CHI(K)+TCHI1*ALPHA*(1.0_LDP-EMHNUKT(K))
	  ETA(K)=ETA(K)+TETA1*ALPHA*EMHNUKT(K)
!
! 
! Now add in BOUND-FREE contributions. We first compute vectors which can
! decrease the compuation time.
!
! NB: A clearer way of writing the expressions for ETA and CHI is
! where TMP_HNST is the LTE population defined by the actual population of
! the destination (target) level.
!
! TMP_HNST=HNST(I,K)*(DI(ION_LEV,K)/DIST(ION_LEV,K))*(DIST(1,K)/DI(1,K))
! CHI(K)=CHI(K)+ALPHA*(HN(I,K)-TMP_HNST*EMHNUKT(K))
! ETA(K)=ETA(K)+TETA2*TMP_HNST*EMHNUKT(K)
!
! In case some levels are out of order, we get make sure we get the minimum edge
! frequency.
!
	  T1=MINVAL(EDGE(1:N))
	  IF(NU .GE. T1)THEN
	    K=DPTH_INDX
	    LOG_COR_FAC= LOG(DI(ION_LEV,K)/DI(1,K))+LOG_DIST(1,K)-LOG_DIST(ION_LEV,K)-HDKT*NU/T(K)
	    COR_FAC=EXP(LOG_COR_FAC)
	  END IF
!
! Compute dissolution vectors that are independent of level.
!
	 IF(MOD_DO_LEV_DIS)THEN
	    K=DPTH_INDX
	    YDIS(K)=1.091_LDP*(X_LEV_DIS(K)+4.0_LDP*(ZION-1)*A_LEV_DIS(K))*B_LEV_DIS(K)*B_LEV_DIS(K)
	    XDIS(K)=B_LEV_DIS(K)*X_LEV_DIS(K)
	  END IF
!
! Due to the choice of units, ESEC=6.65E-15*ED(K). We will thus
! ignore the opacity if less than 1.0D-05 times the e.s. opacity.
! This test is fairly mild, since the test has only to be violated at
! on depth for that opacity source to be included.
!
!
! 
! Now do the actual Bound-Free computation.
!
! Evaluate the bound-free contributions. If N is small, the inner loop
! is over ND, otherwise the inner loop is over N. For some species
! (e.g. FeIV) N can be large (e.g. 300). The optimal switching point is
! unclear because some photo-ionization cross-sections can be zero,
! and because of the different overheads. It is proably machine dependent.
!
! Need to refine the 0.5. We can exit since EDGE(1) > EDGE(2) etc
!
	  TETA1=TWOHCSQ*(NU**3)
	  K=DPTH_INDX
	  DO I=N,1,-1
	    IF(NU .LT. 0.5*EDGE(I))THEN
	      EXIT   
	    ELSE IF(NU .GE. EDGE(I) .AND. ALPHA_VEC(I) .GT. 0.0_LDP)THEN
	      TETA2=TETA1*ALPHA_VEC(I)
	      T1=EXP(LOG_COR_FAC+LOG_HNST(I,K))
	      CHI(K)=CHI(K)+ALPHA_VEC(I)*(HN(I,K)-T1)
	      ETA(K)=ETA(K)+TETA2*T1
	    ELSE IF(DIS_CONST(I) .GE. 0.0_LDP)THEN
!
! Add in BOUND-FREE contributions due to level dissolution.
!
	      TETA2=TETA1*ALPHA_VEC(I)
	      T1=7.782_LDP+XDIS(K)*DIS_CONST(I)
	      T2=T1/(T1+YDIS(K)*DIS_CONST(I)*DIS_CONST(I))
	      IF(T2 .GT. PHOT_DIS_PARAMETER)THEN
	        T1=EXP(LOG_HNST(I,K)-HDKT*NU/T(K))
	        CHI(K)=CHI(K)+ALPHA_VEC(I)*T2*(HN(I,K)-T1)
	        ETA(K)=ETA(K)+TETA2*T2*T1
	      END IF		!
	    END IF		!NU > EDGE
	  END DO		!Variable
	END DO			!Depth
!
	RETURN
	END
