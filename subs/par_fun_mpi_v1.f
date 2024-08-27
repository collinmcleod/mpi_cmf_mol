!
! Routine to evaluate the population of the highest ionization stage
! assuming fixed departure coefficients. U AND PHI are defined in
! Mihalas (1978), page 110. U and PHI are also evaluated for
! next ionization stage. This will be overwritten if an additional
! ionization stage is present.
!
! On exit, He2 contains the log of the departure coefficients.
!
	SUBROUTINE PAR_FUN_MPI_V1(U,PHI,ZPFN,HIGH_POP,
	1             HE2,LOG_HE2LTE,W_HE2,DHE2,EDGE,GHE2,GION,ZION,
	1             T,OLD_T,OLD_ED,N,DST,DEND,ND,
	1             ISPEC,NSPEC,ID,NION_MAX,PRES,ION_ID,MODE)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 05-Apr-2011 - Update from PAR_FUN_V3.
!                         LOG_HE2LTE insted of HE2LTE passed in call (18-Dec-2010).
! Altered 22-Feb-2010 - Altered computation of PHI to help problems with dynamic range.
! Altered 15-May-2009 - VERBOSE option installed.
! Altered 03-Feb-2008 - Can now use excitation temperatures (MODE=TX), as well
!                         as departure coefficients (MODE=DC), to evaluate
!                         partition functions.
! Altered 07-Mar-2006 - 2.07D-22 pulled into exponential.
! Altered 25-May-1996 - ERROR_LU inserted.
! Altered 16-Jan-1995 - Occupation probabilities included in calculation of
!                         partition function. CALL altered, but still V2.
! Altered 15-Jan-1995 - HIGH_POP and DHE2 inserted. Allowest population of
!                         the highest ionization state to be set through
!                         successive calls. CALL altered.
! Altered 27-Sep-1990 - HE2 is now converted to departure coefficients
!                       in routine.
!
	INTEGER N
	INTEGER DST,DEND,ND
	INTEGER ISPEC,NSPEC
	INTEGER ID,NION_MAX
!
	REAL(KIND=LDP) HE2(N,DST:DEND)
	REAL(KIND=LDP) LOG_HE2LTE(N,DST:DEND)
	REAL(KIND=LDP) W_HE2(N,DST:DEND)
	REAL(KIND=LDP) DHE2(DST:DEND)
!
	REAL(KIND=LDP) EDGE(N)
	REAL(KIND=LDP) GHE2(N)
	REAL(KIND=LDP) GION,ZION
!
	REAL(KIND=LDP) U(DST:DEND,NION_MAX)
	REAL(KIND=LDP) PHI(DST:DEND,NION_MAX)
	REAL(KIND=LDP) HIGH_POP(DST:DEND,NSPEC)
	REAL(KIND=LDP) ZPFN(NION_MAX)
!
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) OLD_T(ND)
	REAL(KIND=LDP) OLD_ED(ND)
	CHARACTER(LEN=*) MODE
	CHARACTER(LEN=*) ION_ID
	LOGICAL PRES
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
!
	REAL(KIND=LDP) T1,T2,T3,LOG_RGU
	REAL(KIND=LDP) TEXCITE(N)
	REAL(KIND=LDP) DELTA_T
	REAL(KIND=LDP) CONST
!
	INTEGER I,J
	INTEGER COUNT
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
	LOGICAL VERBOSE
!
	IF(.NOT. PRES)RETURN
	LUER=ERROR_LU()
	CALL GET_VERBOSE_INFO(VERBOSE)
!
	IF(ID .GT. NION_MAX)THEN
	  WRITE(LUER,*)'Error in PAR_FUN ---- ID too small'
	  STOP
	END IF
!
	DO I=DST,DEND
	  HIGH_POP(I,ISPEC)=DHE2(I)
	  IF(W_HE2(1,I) .NE. 1.0_LDP)THEN
	    WRITE(LUER,*)'Error in PAR_FUN --- occupation probability for ',
	1               'ground state must be unity.'
	    WRITE(LUER,*)'N=',N
	    WRITE(LUER,*)'EDGE(1)=',EDGE(1)
	    STOP
	  END IF
	END DO
!
! Ouput departure coefficents.
!
	IF(VERBOSE)THEN
	  OPEN(UNIT=7,FILE=TRIM(ION_ID)//'_ORIG_DC',STATUS='UNKNOWN')
	    WRITE(7,'(A)')' '
	    WRITE(7,*)N,ND
	    DO I=DST,DEND
	      WRITE(7,'(A)')' '
	      WRITE(7,'(4ES14.6)')DHE2(I),OLD_ED(I),T(I),OLD_T(I)
	      WRITE(7,'(5ES14.6)')(LOG(HE2(J,I))-LOG_HE2LTE(J,I),J=1,N)
	    END DO
	  CLOSE(UNIT=7)
	END IF
!
! Convert HE2 array to departure coefficients.
!
	IF(MODE .EQ. 'DC')THEN
	  DO I=DST,DEND
	    DO J=1,N
	      HE2(J,I)=LOG(HE2(J,I))-LOG_HE2LTE(J,I)
	    END DO
	  END DO
	ELSE IF(MODE .EQ. 'TX')THEN
	  DO I=DST,DEND
!
! We first get the excitation temperature for each level, and then compute
! the revised departure coefficients. Remeber the He2 contains log b.
!
	    DO J=1,N
	      T2=HDKT*EDGE(J)
	      CONST=HE2(J,I)+1.5_LDP*LOG(T2)+LOG( GION/GHE2(J)/W_HE2(J,I)/2.07078E-22_LDP/OLD_ED(I)/DHE2(I) )
	      IF(CONST .LT. 1.0_LDP)THEN
	        TEXCITE(J)=EXP(CONST*0.67_LDP)
	      ELSE
	        TEXCITE(J)=CONST
	      END IF
              COUNT=0
	      DELTA_T=1.0E+10_LDP
	      DO WHILE( ABS(DELTA_T/TEXCITE(J)) .GT. 1.0E-10_LDP)
	        COUNT=COUNT+1
	        T1=SQRT(TEXCITE(J))
	        T3=EXP(CONST-TEXCITE(J))
	        DELTA_T=(T1*TEXCITE(J) - CONST)/T1/(1.5_LDP+TEXCITE(J))
	        IF(DELTA_T .GT. 0.8_LDP*TEXCITE(J))DELTA_T=0.8_LDP*TEXCITE(J)
	        IF(DELTA_T .LT. -0.8_LDP*TEXCITE(J))DELTA_T=-0.8_LDP*TEXCITE(J)
	        TEXCITE(J)=TEXCITE(J)-DELTA_T
	        IF(COUNT .GE. 100)THEN
	          WRITE(LUER,*)'Error in PAR_FUN_V3'
	          WRITE(LUER,*)'Too many iterations to get the excitation temperature'
	          STOP
	        END IF
	      END DO
              TEXCITE(J)=T2/TEXCITE(J)
	      TEXCITE(J)=T(I)+(TEXCITE(J)-OLD_T(I))
	      HE2(J,I)=HDKT*EDGE(J)*(1.0_LDP/TEXCITE(J)-1.0_LDP/T(I)) + 1.5_LDP*LOG(T(I)/TEXCITE(J))
	    END DO
	  END DO
	ELSE
	  WRITE(LUER,*)'Error in PAR_FN_V3: interplation mode not recognized'
	  WRITE(LUER,*)'Passed mode is ',TRIM(MODE)
	  STOP
	END IF
!
! We can now compute the non-LTE partition function.
!
	ZPFN(ID)=ZION-1
	ZPFN(ID+1)=ZION
	T1=HDKT*EDGE(1)
	DO I=DST,DEND
	  U(I,ID)=GHE2(1)
	  U(I,ID+1)=GION
	  PHI(I,ID+1)=0.0_LDP
	  LOG_RGU=LOG(2.07078E-22_LDP)+HE2(1,I)
	  PHI(I,ID)=EXP( LOG_RGU+T1/T(I) )/T(I)/SQRT(T(I))
	  DO J=2,N
	    T2=HDKT*(EDGE(J)-EDGE(1))
	    U(I,ID)=U(I,ID) + GHE2(J)*EXP(T2/T(I)+HE2(J,I)-HE2(1,I))*W_HE2(J,I)
	  END DO
	END DO
!
! Ouput departure coefficents. If MODE=DC, these will be the same as in _ORID_DC files.
! If MODE=TX, the files will contain the revised departure coefficients.
!
	IF(VERBOSE)THEN
	  OPEN(UNIT=7,FILE=TRIM(ION_ID)//'_REV_DC',STATUS='UNKNOWN')
	    WRITE(7,'(A)')' '
	    WRITE(7,*)DST,DEND
	    DO I=DST,DEND
	      WRITE(7,'(A)')' '
	      WRITE(7,'(4ES14.6)')DHE2(I),OLD_ED(I),T(I),OLD_T(I)
	      WRITE(7,'(5ES14.6)')(HE2(J,I),J=1,N)
	    END DO
	  CLOSE(UNIT=7)
	END IF
!
	RETURN
	END
