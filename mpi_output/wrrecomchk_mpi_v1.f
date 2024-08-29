!
! Subroutine to write out the recombination and collision rates for
! any ion. The valuse are output to an existing file (overwriting), or
! to a new file.
!
!
	SUBROUTINE WRRECOMCHK_MPI_V1(PR,RR,CPR,CRR,CHG_PR,CHG_RR,ADVEC_RR,
	1                 DIERECOM,ADDRECOM,X_RECOM_1,X_RECOM_2,
	1                 NT_ION_RATE,NT_ION_RATE_2E_1,NT_ION_RATE_2E_2,
	1                 R,T,ED,DHYD,N,ND,LU,FILNAM,STRDESC,ID)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 28-Aug-2024 : Based on WRRECOM_CHK_V5
!
	INTEGER ID
	INTEGER N,ND,LU,ML,MF,I,MS,J,IOS
	REAL(KIND=LDP) PR(N,ND)			!Radiative photioization rate
	REAL(KIND=LDP) RR(N,ND)			!Radiative recombination rate
	REAL(KIND=LDP) CPR(ND)			!Collisional ioization rate
	REAL(KIND=LDP) CRR(ND)			!Collisional recombination rate
	REAL(KIND=LDP) CHG_PR(ND)		!Charge ionization rate
	REAL(KIND=LDP) CHG_RR(ND)		!Charge recombination rate
	REAL(KIND=LDP) ADVEC_RR(ND)		!Advection recombination rate
	REAL(KIND=LDP) DIERECOM(ND)
	REAL(KIND=LDP) ADDRECOM(ND)
	REAL(KIND=LDP) X_RECOM_1(ND),X_RECOM_2(ND)
	REAL(KIND=LDP) NT_ION_RATE(ND)
	REAL(KIND=LDP) NT_ION_RATE_2E_1(ND)
	REAL(KIND=LDP) NT_ION_RATE_2E_2(ND)
	REAL(KIND=LDP) R(ND),T(ND),ED(ND),DHYD(ND)
	CHARACTER*(*) FILNAM
	CHARACTER*(*) STRDESC(ID+2)
!
! Local data.
!
	REAL(KIND=LDP) TOTRR(ND),NETRR(ND)
	REAL(KIND=LDP) ABS_SUM
	REAL(KIND=LDP) ADVEC_SUM
	REAL(KIND=LDP) T1
!
	INTEGER ERROR_LU,LUER
	INTEGER, PARAMETER :: IZERO=0
	EXTERNAL ERROR_LU
!
	NETRR(:)=0.0_LDP                 !ND
	TOTRR(:)=0.0_LDP                 !ND
	ADVEC_SUM=SUM(ADVEC_RR)
!
	CALL GEN_ASCI_OPEN(LU,FILNAM,'UNKNOWN',' ',' ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error opening RECOM file',FILNAM
	  WRITE(LUER,*)'IOSTAT=',IOS
	  RETURN
	END IF
!
	MS=1
	DO 10 ML=0,ND-1,10
	  MF=ML+10
	  IF(MF .GT. ND)MF=ND
	  IF(ML .NE. 0)WRITE(LU,'(1H1)')
!
	  WRITE(LU,'(3X,''Depth index'')')
	  WRITE(LU,'(2X,10(I11,'' ''))')(J,J=MS,MF)
	  WRITE(LU,'(/,3X,''Radius [1.0E+10cm] '')')
	  WRITE(LU,999)(R(J),J=MS,MF)
	  WRITE(LU,'(/,3X,''Temperature [1.0E+4K] '')')
	  WRITE(LU,999)(T(J),J=MS,MF)
	  WRITE(LU,'(/,3X,''Electron Density'')')
	  WRITE(LU,999)(ED(J),J=MS,MF)
	  WRITE(LU,'(/,3X,''Ion Density '')')
	  WRITE(LU,999)(DHYD(J),J=MS,MF)
!
	  WRITE(LU,'(/,3X,(A),'' Photoionization Rates'')')TRIM(STRDESC(ID))
	  DO I=1,N
	    WRITE(LU,999)(PR(I,J),J=MS,MF)
	  END DO
!
	  WRITE(LU,'(/3X,''Colisional Ionization Rate '') ')
	  WRITE(LU,999)(CPR(J),J=MS,MF)
!
	  IF(CHG_PR(MS) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Charge Transfer Ionization Rate '') ')
	    WRITE(LU,999)(CHG_PR(J),J=MS,MF)
	  END IF
!
	  WRITE(LU,'(/,3X,(A),'' Recombination Rates'')')TRIM(STRDESC(ID))
	  DO  I=1,N
	    WRITE(LU,999)(RR(I,J),J=MS,MF)
	  END DO
!
	  WRITE(LU,'(/3X,''Colisional Recombination Rate '') ')
	  WRITE(LU,999)(CRR(J),J=MS,MF)
!
	  IF(CHG_PR(MS) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Charge Transfer Recombination Rate '') ')
	    WRITE(LU,999)(CHG_RR(J),J=MS,MF)
	  END IF
!
	  IF(ADVEC_SUM .NE. 0)THEN
	    WRITE(LU,'(/3X,''Effective Advection Recombination Rate '') ')
	    WRITE(LU,999)(ADVEC_RR(J),J=MS,MF)
	  END IF
!
	  IF(NT_ION_RATE(MS) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Non-Thermal Ionization  Rate (one electron) '') ')
	    WRITE(LU,999)(NT_ION_RATE(J),J=MS,MF)
	  END IF
!
! Write out net X-ray recombinations from i+1 to i-1.  eg If NIV, net recom's from NV to NIII).
!
	  IF(NT_ION_RATE_2E_1(1) .NE. 0 .OR. NT_ION_RATE_2E_1(ND) .NE. 0)THEN
	    WRITE(LU,'(/,3X,4A)')'NT ionization rate (2e): ',TRIM(STRDESC(ID-1)),' to ',TRIM(STRDESC(ID+1))
	    WRITE(LU,999)(NT_ION_RATE_2E_1(J),J=MS,MF)
	  END IF
!
! Write out net X-ray recombinations from i+2 to i. (eg If NIV, net ion's from NIV to NVI).
!
	  IF(NT_ION_RATE_2E_2(1) .NE. 0 .OR. NT_ION_RATE_2E_2(ND) .NE. 0)THEN
	    WRITE(LU,'(/,3X,4A)')'NT ionization rate (2e): ',TRIM(STRDESC(ID)),' to ',TRIM(STRDESC(ID+2))
	    WRITE(LU,999)(NT_ION_RATE_2E_2(J),J=MS,MF)
	  END IF
	  FLUSH(LU)
!
	  IF(DIERECOM(1) .NE. 0 .AND. DIERECOM(ND) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Dielectronic Recombination Rate '') ')
	    WRITE(LU,999)(DIERECOM(J),J=MS,MF)
	  END IF
!
	  IF(ADDRECOM(1) .NE. 0 .AND. ADDRECOM(ND) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Implicit Recombination Rate '') ')
	    WRITE(LU,999)(ADDRECOM(J),J=MS,MF)
	  END IF
!
! Write out net X-ray recombinations from i+1 to i-1.
! (eg If NIV, net recom's from NV to NIII).
!
	  IF(X_RECOM_1(1) .NE. 0 .AND. X_RECOM_1(ND) .NE. 0)THEN
	    WRITE(LU,'(/,3X,4A)')'Net X-ray ionization rate: ',TRIM(STRDESC(ID-1)),' to ',TRIM(STRDESC(ID+1))
	    WRITE(LU,999)(-X_RECOM_1(J),J=MS,MF)
	  END IF
!
! Write out net X-ray recombinations from i+2 to i.
! (eg If NIV, net recom's from NVI to NIV).
!
	  IF(X_RECOM_2(1) .NE. 0 .AND. X_RECOM_2(ND) .NE. 0)THEN
	    WRITE(LU,'(/,3X,4A)')'Net X-ray ionization rate: ',TRIM(STRDESC(ID)),' to ',TRIM(STRDESC(ID+2))
	    WRITE(LU,999)(-X_RECOM_2(J),J=MS,MF)
	  END IF
	  FLUSH(LU)
!
	  DO J=MS,MF
	    ABS_SUM=0.0_LDP
	    DO I=1,N
	      NETRR(J)=NETRR(J)+(RR(I,J)-PR(I,J))
	      TOTRR(J)=TOTRR(J)+RR(I,J)
	      ABS_SUM=ABS_SUM+RR(I,J)+PR(I,J)
	    END DO
	    ABS_SUM=ABS_SUM+CRR(J)+CPR(J)+
	1             CHG_PR(J)+CHG_RR(J)+
	1             ABS(DIERECOM(J))+ABS(ADDRECOM(J))+ABS(ADVEC_RR(J))+
	1             ABS(X_RECOM_1(J))+ABS(X_RECOM_2(J))+
	1             ABS(NT_ION_RATE(J))+ABS(NT_ION_RATE_2E_1(J))+ABS(NT_ION_RATE_2E_2(J))
	    NETRR(J)=200.0_LDP*(NETRR(J)+(CRR(J)-CPR(J))+
	1                    (CHG_RR(J)-CHG_PR(J))+
	1                    DIERECOM(J)+ADDRECOM(J)+ADVEC_RR(J)+
	1                    X_RECOM_1(J)+X_RECOM_2(J)-
	1                    NT_ION_RATE(J)-
	1                    NT_ION_RATE_2E_1(J) -NT_ION_RATE_2E_2(J))/ABS_SUM
	    T1=TOTRR(J)+ADDRECOM(J)
	    IF(T1 .NE. 0)DIERECOM(J)=DIERECOM(J)/T1
	    ADDRECOM(J)=ADDRECOM(J)/ED(J)/DHYD(J)
	    TOTRR(J)=TOTRR(J)/ED(J)/DHYD(J)
	  END DO
!
	  WRITE(LU,'(/3X,''Net Recombination Rate (% of total) '') ')
	  WRITE(LU,999)(NETRR(J),J=MS,MF)
	  WRITE(LU,'(/3X,''Radiative Recombination Coefficient for '',
	1                       ''explicitly treated levels.'') ')
	  WRITE(LU,999)(TOTRR(J),J=MS,MF)
	  IF(ADDRECOM(1) .NE. 0 .AND. ADDRECOM(ND) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Radiative Recombination for'',
	1                         ''implicit levels.'') ')
	    WRITE(LU,999)(ADDRECOM(J),J=MS,MF)
	  END IF
	  IF(DIERECOM(1) .NE. 0 .AND. DIERECOM(ND) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Ratio of Total Dielectronic Recombination'',
	1                       '' to Total Radiatve Recombination'') ')
	    WRITE(LU,999)(DIERECOM(J),J=MS,MF)
	  END IF
!
	  MS=MS+10
10	CONTINUE
!
	CLOSE(UNIT=LU)
	RETURN
999	FORMAT(1X,1P,10E12.4)
	END
