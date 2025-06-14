	SUBROUTINE SE_BA_NON_THERM_MPI_V1(dE_RAD_DECAY,dE_SHOCK_POWER,COMPUTE_BA,NT,ND,DEC_NRG_SCL_FAC)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
        USE MOD_NON_THERM
	USE NUC_ISO_MOD
	USE SHOCK_POWER_MOD
	USE CONTROL_VARIABLE_MOD, only : INC_SHOCK_POWER, SHOCK_POWER_FAC
	USE MPI
	IMPLICIT NONE
!
! Altered 14-Aug-2022 : Added SHOCK terms
! Altered 29-Jul-2019 : Added EHB equation.
! Altered 06-Jun-2012 : Now use ATM(ID)%XzVLTE_F_ON_S to avoid floating point errors.
! Altered 12-Feb-2012 : Changed ordering of loops, split ionizations from excitations, and paralleized over depth.
! Created 16-Sep-2010
!
        INTEGER NT
        INTEGER ND
	INTEGER, PARAMETER :: IZERO=0
!
! Output: STEQ_T in MOD_CMFGEN is modified.
!
        REAL(KIND=LDP) dE_RAD_DECAY(DST:DEND)
	REAL(KIND=LDP) dE_SHOCK_POWER(DST:DEND)
!
	REAL(KIND=LDP) LOCAL_ION_HEATING(DST:DEND)
	REAL(KIND=LDP) LOCAL_EXC_HEATING(DST:DEND)
	REAL(KIND=LDP) WRK_EDEP(DST:DEND)
	REAL(KIND=LDP) WRK_EDEP_eV(DST:DEND)
	LOGICAL COMPUTE_BA
!
	REAL(KIND=LDP) MOD_YE(NKT,DST:DEND)
	REAL(KIND=LDP) PI
	REAL(KIND=LDP) ELECTRON_VOLT
	INTEGER GET_INDX_DP
	EXTERNAL GET_INDX_DP, ELECTRON_VOLT
!
	REAL(KIND=LDP) RATE
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) ION_EXC_EN
	REAL(KIND=LDP) SCALE
	REAL(KIND=LDP) GUPPER
	REAL(KIND=LDP) T1_SUM(ND),T2_SUM(ND)
!
	INTEGER ID
	INTEGER ISPEC
	INTEGER IT
	INTEGER IKT
	INTEGER IKT_ST
	INTEGER DPTH_INDX
	INTEGER SE_ION_LEV
	INTEGER I,J
!
	INTEGER NL
	INTEGER NUP
	INTEGER NL_F
	INTEGER NUP_F
!
	INTEGER, PARAMETER :: LU_TH=8
	INTEGER, PARAMETER :: LU_ER=6
	REAL(KIND=LDP) DEC_NRG_SCL_FAC
	REAL(KIND=LDP), PARAMETER :: Hz_to_eV=13.60_LDP/3.2897_LDP
	REAL(KIND=LDP), PARAMETER :: Hz_to_erg=6.6261965E-12_LDP
	LOGICAL ONE_ELEC
!
	LOCAL_ION_HEATING=0.0_LDP
	LOCAL_EXC_HEATING=0.0_LDP
!
! For historical reasons STEQ contains Int[chi.J - eta]dv. Rather than multiply
! this term everywhere by 4pi, we divide the radiactive heating by 4pi.
! Also, since chi, and eta are a factor of 10^10 too large, we need to
! scale by 10^10. The BA matrix does not need to be altered.
!
        PI=ACOS(-1.0_LDP)
        SCALE=1.0E+10_LDP/4.0_LDP/PI

	DO DPTH_INDX=DST,DEND
	  I=DPTH_INDX
          dE_RAD_DECAY(I)=RADIOACTIVE_DECAY_ENERGY(I)
	  dE_SHOCK_POWER(I)=0.0_LDP
          WRK_EDEP (I)= RADIOACTIVE_DECAY_ENERGY(I)
	  IF(INC_SHOCK_POWER)THEN
	    WRK_EDEP(I)= WRK_EDEP(I) + SHOCK_POWER(I) * SHOCK_POWER_FAC
	    dE_SHOCK_POWER(I)= SHOCK_POWER(I) * SHOCK_POWER_FAC
	  END IF
!
! Below, we replace RADIOACTIVE_DECAY_ENERGY by WRK_EDEP,  which can therefore contain decay and shock power
!
          STEQ_T(I)=STEQ_T(I)+SCALE*WRK_EDEP(I)                      !*FRAC_ELEC_HEATING
	  WRK_EDEP_eV(I)=DEC_NRG_SCL_FAC*RADIOACTIVE_DECAY_ENERGY(I)/ELECTRON_VOLT()
	END DO
!
	DO I=1,NUM_IONS
	  ATM(I)%NTCXzV(:)=0.0_LDP
	  ATM(I)%NTIXzV(:)=0.0_LDP
	  ATM(I)%NTIXzV_2E(:)=0.0_LDP
	  ATM(I)%NT_ION_CXzV(:)=0.0_LDP
	  ATM(I)%NT_EXC_CXzV(:)=0.0_LDP
	END DO
!
! Compute the ionization contribution
!
	IF (INCLUDE_NON_THERM_IONIZATION) THEN
!
! Loop over different ionization routes.
!
	  DO IT=1,NUM_THD
	    ID=THD(IT)%LNK_TO_ION
	    ISPEC=THD(IT)%LNK_TO_SPECIES
	    IKT_ST=GET_INDX_DP(THD(IT)%ION_POT,XKT,NKT)
!
	    DO DPTH_INDX=DST,DEND
!
! Get the total excitation route. We ignore slight differences in IP.
!
	      RATE = 0.0_LDP
	      DO IKT=IKT_ST,NKT
	        RATE = RATE + YE(IKT,DPTH_INDX)*dXKT(IKT)*THD(IT)%CROSS_SEC(IKT)
	      END DO
	      RATE=RATE*WRK_EDEP_eV(DPTH_INDX)
!
	      DO J=1,THD(IT)%N_ION_ROUTES
	        IF(ID .EQ. SPECIES_END_ID(ISPEC)-1)THEN
                  SE_ION_LEV=ATM(ID)%NXzV+1
	          GUPPER=THD(IT)%SUM_GION
	          ION_EXC_EN=0.0_LDP
	          ONE_ELEC=.TRUE.
	        ELSE IF(THD(IT)%ION_LEV(J) .EQ. SE(ID)%XRAY_EQ)THEN
                  SE_ION_LEV=SE(ID)%XRAY_EQ
	          GUPPER=THD(IT)%SUM_GION
	          ION_EXC_EN=ATM(ID+1)%EDGEXzV_F(1)
	          ONE_ELEC=.FALSE.
	        ELSE
	          NUP_F=THD(IT)%ION_LEV(J)
!	          NUP=ATM(ID+1)%F_TO_S_XzV(NUP_F)
!                 SE_ION_LEV=SE(ID)%ION_LEV_TO_EQ_PNT(NUP)
	          GUPPER=ATM(ID+1)%GXzV_F(1)                    !NUP_F)
                  SE_ION_LEV=ATM(ID)%NXzV+1                     !NUP -- assume all to ground state at  present.
	          ION_EXC_EN=0.0_LDP                             !=ATM(ID+1)%EDGEXzV_F(1)-ATM(ID+1)%EDGEXzV_F(NUP_F)
	          ONE_ELEC=.TRUE.
	        END IF
	        DO I=1,THD(IT)%N_STATES
	          NL_F=THD(IT)%ATOM_STATES(I); NL=ATM(ID)%F_TO_S_XzV(NL_F)
	          T1=RATE*ATM(ID)%XzV_F(NL_F,DPTH_INDX)*GUPPER/THD(IT)%SUM_GION
	          IF(ONE_ELEC)THEN
	            ATM(ID)%NTIXzV(DPTH_INDX)=ATM(ID)%NTIXzV(DPTH_INDX)+T1             ! non-thermal rates?
	          ELSE
	            ATM(ID)%NTIXzV_2E(DPTH_INDX)=ATM(ID)%NTIXzV_2E(DPTH_INDX)+T1             ! non-thermal rates?
	          END IF
	          ATM(ID)%NTCXzV(DPTH_INDX)=ATM(ID)%NTCXzV(DPTH_INDX)  &
	                          +Hz_to_erg*T1*(ATM(ID)%EDGEXzV_F(NL_F)+ION_EXC_EN)  ! non-thermal cooling?
!
	          SE(ID)%STEQ(NL,DPTH_INDX)=SE(ID)%STEQ(NL,DPTH_INDX)-T1
	          SE(ID)%STEQ(SE_ION_LEV,DPTH_INDX)=SE(ID)%STEQ(SE_ION_LEV,DPTH_INDX)+T1
	          LOCAL_ION_HEATING(DPTH_INDX)=LOCAL_ION_HEATING(DPTH_INDX)+ T1*(ATM(ID)%EDGEXzV_F(NL_F)+ION_EXC_EN)
	          ATM(ID)%NT_ION_CXzV(DPTH_INDX)=ATM(ID)%NT_ION_CXzV(DPTH_INDX) &
	                          +Hz_to_eV*T1*(ATM(ID)%EDGEXzV_F(NL_F)+ION_EXC_EN)
	        END DO
	      END DO
!
! NB:  XzV_F= XzV  (XzVLTE_F/XzVLTE)
!
	      IF(COMPUTE_BA)THEN
	        DO J=1,THD(IT)%N_ION_ROUTES
	          IF(ID .EQ. SPECIES_END_ID(ISPEC)-1)THEN
                    SE_ION_LEV=ATM(ID)%NXzV+1
	            GUPPER=THD(IT)%SUM_GION
	          ELSE IF(THD(IT)%ION_LEV(J) .EQ. SE(ID)%XRAY_EQ)THEN
                    SE_ION_LEV=SE(ID)%XRAY_EQ
	            GUPPER=THD(IT)%SUM_GION
	          ELSE
	            NUP_F=THD(IT)%ION_LEV(J)
! 	            NUP=ATM(ID+1)%F_TO_S_XzV(NUP)
!	            SE_ION_LEV=SE(ID)%ION_LEV_TO_EQ_PNT(NUP)
	            GUPPER=ATM(ID+1)%GXzV_F(1)                    !NUP_F)
	          END IF
                  SE_ION_LEV=ATM(ID)%NXzV+1               !NUP -- assume all to ground state at  present.
	          DO I=1,THD(IT)%N_STATES
	            NL_F=THD(IT)%ATOM_STATES(I); NL=ATM(ID)%F_TO_S_XzV(NL_F)
	            T1=RATE*GUPPER/THD(IT)%SUM_GION * ATM(ID)%XzVLTE_F_ON_S(NL_F,DPTH_INDX)
!&	                         (ATM(ID)%XzVLTE_F(NL_F,DPTH_INDX)/ATM(ID)%XzVLTE(NL,DPTH_INDX))
	            SE(ID)%BA_PAR(NL,NL,DPTH_INDX)=SE(ID)%BA_PAR(NL,NL,DPTH_INDX)-T1
	            SE(ID)%BA_PAR(SE_ION_LEV,NL,DPTH_INDX)=SE(ID)%BA_PAR(SE_ION_LEV,NL,DPTH_INDX)+T1
	            BA_T_PAR_EHB(NL,DPTH_INDX)=BA_T_PAR_EHB(NL,DPTH_INDX)-Hz_to_erg*T1*(ATM(ID)%EDGEXzV_F(NL_F)+ION_EXC_EN)
	          END DO
	        END DO
	      END IF
!
	    END DO   		!Loop over depth
	  END DO		!Ionization species/route
	END IF			!Include ionizations?
!
	DO DPTH_INDX=DST,DEND
	  DO I=1,NKT
	    MOD_YE(I,DPTH_INDX)=YE(I,DPTH_INDX)*dXKT(I)/XKT(I)
	  END DO
	END DO
!
	IF(INCLUDE_NON_THERM_EXCITATION)THEN
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO DPTH_INDX=DST,DEND
	        DO J=2,ATM(ID)%NXzV_F
	          DO I=1,MIN(20,J-1)
	            NL_F=I; NUP_F=J
	            CALL TOTAL_BETHE_RATE_MOD_YE_MPI_V1(RATE,NL_F,NUP_F,MOD_YE(:,DPTH_INDX),XKT,NKT,ID)
	            NUP=ATM(ID)%F_TO_S_XzV(NUP_F)
!
! NB:  XzV_F= XzV  (XzVLTE_F/XzVLTE)
!
	            NL=ATM(ID)%F_TO_S_XzV(NL_F)
	            T1=RATE*WRK_EDEP_eV(DPTH_INDX)
	            IF(COMPUTE_BA)THEN
	              T2=T1*ATM(ID)%XzVLTE_F_ON_S(NL_F,DPTH_INDX)
	              SE(ID)%BA_PAR(NL,NL,DPTH_INDX)=SE(ID)%BA_PAR(NL,NL,DPTH_INDX)-T2
	              SE(ID)%BA_PAR(NUP,NL,DPTH_INDX)=SE(ID)%BA_PAR(NUP,NL,DPTH_INDX)+T2
	              BA_T_PAR_EHB(NL,DPTH_INDX)=BA_T_PAR_EHB(NL,DPTH_INDX) - &
	                             HZ_TO_ERG*T2*(ATM(ID)%EDGEXzV_F(NL_F)-ATM(ID)%EDGEXzV_F(NUP_F))
	            END IF
	            T1=T1*ATM(ID)%XzV_F(NL_F,DPTH_INDX)
	            SE(ID)%STEQ(NL,DPTH_INDX)=SE(ID)%STEQ(NL,DPTH_INDX) - T1
	            SE(ID)%STEQ(NUP,DPTH_INDX)=SE(ID)%STEQ(NUP,DPTH_INDX) + T1
	            T1=T1*(ATM(ID)%EDGEXzV_F(NL_F)-ATM(ID)%EDGEXzV_F(NUP_F))
	            LOCAL_EXC_HEATING(DPTH_INDX)=LOCAL_EXC_HEATING(DPTH_INDX) + T1
	            T1=T1*Hz_to_erg
	            ATM(ID)%NTCXzV(DPTH_INDX)=ATM(ID)%NTCXzV(DPTH_INDX) + T1 			 ! non-thermal cooling?
	            ATM(ID)%NT_EXC_CXzV(DPTH_INDX)=ATM(ID)%NT_EXC_CXzV(DPTH_INDX)  + T1
	          END DO
	        END DO
	      END DO 		!Loop over depth
!	      I=SE(ID)%N_SE*SE(ID)%N_IV*(DEND-DST+1)
!	      CALL CHECK_VEC_NAN(SE(ID)%BA_PAR,I,'SE_NON_THERM(STOP)',J)
	    END IF
	  END DO		!Loop over ionization stage
	END IF			!Include excitations?
!
! It is a minus sign since ATM(ID)%NTCXzV contains the cooling rate.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    STEQ_T_EHB=STEQ_T_EHB-ATM(ID)%NTCXzV
	  END IF
	END DO
!
! Estimate number of non-thermal electrons.
!
	T1_SUM=0.0_LDP
	T2_SUM=0.0_LDP
	DO DPTH_INDX=DST,DEND
	  T1=0.0_LDP
	  T2=0.0_LDP
	  DO IKT=1,NKT
	    T1=T1+YE(IKT,DPTH_INDX)*dXKT(IKT)/SQRT(XKT(IKT))
	    T2=T2+XKT(IKT)*YE(IKT,DPTH_INDX)*dXKT(IKT)/SQRT(XKT(IKT))
	  END DO
	  T1_SUM(DPTH_INDX)=T1*SQRT(0.5_LDP*9.109389E-28_LDP/1.602177E-12_LDP)*WRK_EDEP_eV(DPTH_INDX)
	  T2_SUM(DPTH_INDX)=T2*SQRT(0.5_LDP*9.109389E-28_LDP/1.602177E-12_LDP)*WRK_EDEP_eV(DPTH_INDX)
	END DO
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	CALL MPI_ALLREDUCE(MPI_IN_PLACE,T1_SUM,ND,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,IERR)
	CALL MPI_ALLREDUCE(MPI_IN_PLACE,T2_SUM,ND,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,IERR)
!
	IF(MYPE .EQ. 0)THEN
          OPEN(UNIT=LU_TH,FILE='NON_THERM_SPEC_INFO',STATUS='UNKNOWN',POSITION='APPEND')
	    WRITE(LU_TH,*)''
	    WRITE(LU_TH,'(X,A5,6(A12))')'Depth','Ne(NT)','Ne','Ne(NT)/Ne','Ek(NT)','Ek','Ek(NT)/Ek'
	    DO DPTH_INDX=DST,DEND
	      I=DPTH_INDX
	      WRITE(LU_TH,'(X,I5,6ES12.4)')I,T1_SUM(I),ED(I),T1_SUM(I)/ED(I),                  &
	                             T2_SUM(I),1.5_LDP*8.617343E-5_LDP*ED(I)*T(I)*1.0E4_LDP,   &
	                             T2_SUM(I)/(1.5_LDP*8.617343E-5_LDP*ED(I)*T(I)*1.0E4_LDP)
	    END DO
!
	    WRITE(LU_TH,'(//,A)')'Comparison of heating fractions (SE as evaluated in SE_BA_NON_THERM)'
	    WRITE(LU_TH,'(//,A,9(3X,A))')       &
	            '  Fe_nuc(eV)','    Felec','Felec(SE)','     Fion',' Fion(SE)','     Fexc',' Fexc(SE)', &
	            '  E_ion  ','  E_exc  ','  E_cool '
	    DO I=DST,DEND
	    IF(WRK_EDEP_eV(I) .NE. 0.0_LDP)THEN
	      T2=Hz_TO_eV*LOCAL_ION_HEATING(I)/WRK_EDEP_eV(I)
	      T3=Hz_TO_eV*LOCAL_EXC_HEATING(I)/WRK_EDEP_eV(I)
	    ELSE
	      T2=0.0_LDP
	      T3=0.0_LDP
	    END IF
	    T1=(1.0_LDP-T2-T3)
	    WRITE(LU_TH,'(ES12.4,9(ES12.4))')WRK_EDEP_eV(I),        &
	                 FRAC_ELEC_HEATING(I),T1,         &
	                 FRAC_ION_HEATING(I),T2,          &
	                 FRAC_EXCITE_HEATING(I),T3, &
	                 Hz_TO_eV*LOCAL_ION_HEATING(I),Hz_TO_eV*LOCAL_EXC_HEATING(I), &
	                 Hz_TO_eV*LOCAL_ION_HEATING(I)+Hz_TO_eV*LOCAL_EXC_HEATING(I)
	  END DO
	  CLOSE(LU_TH)
	  IF(MYPE .EQ. 0)WRITE(6,*)'Scale factor is ',DEC_NRG_SCL_FAC
	END IF
!
	RETURN
	END
