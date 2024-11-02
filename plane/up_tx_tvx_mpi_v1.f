!
! Routine to compute the matrices TX and TVX which depend on the integration
! at the previous frequency.
!
! TX(, , )  describes the variation of J
! TVX(, , ) describes the variation of RSQH
!
! It is assumed that TA,TB and TC (the tridiagonal matrix to be inverted) have
! already been modified by a call to THOMAS.
!
! KI is a 2 dimensional matrix with
!
!       KI( , ,1) variation of transfer equation w.r.t. CHI
!       KI( , ,2) variation of transfer equation w.r.t. ETA
!
! Upon entry TX( , ,1)  and TX( , ,2) should be zero, as these now reflect
! change in J with respect to CHI and ETA at the current frequency.
! The other matrices are used to reflect changes in J at an earlier frequency.
!
! In general K=1 denotes dCHI and K=2 denotes dETA.
!
	SUBROUTINE UP_TX_TVX_MPI_V1(
	1                 TA,TB,TC,PSIPREV_MOD,
	1                 VB,VC,HU,HL,HS,
	1                 EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                 ND,NM_TX,NM_KI,
	1                 DTAU_BND,OUT_BC_TYPE,
	1                 INIT,DO_THIS_TX_MATRIX)
	USE SET_KIND_MODULE
	USE MPI
	USE MOD_VAR_OPAC_J, ONLY : TX, TVX, RHS_dHdCHI, KI, VDST, VDEND
	IMPLICIT NONE
!
! Altered 28-May-1996 : Calls to DP_ZERO removed.
!                       (, at end of subroutine specification deleted).
!
! Created 10-Mar-1995.   Based on UPTX_J_EDD_V3
!                        EPS_A,EP_B etc passed in call.
!                        Argument odering altered.
!                        TX and TVX modified in the same routine, since
!                          TVX may depend on TX at the previus frequency.
!
	INTEGER ND,NM_TX,NM_KI
!
	REAL(KIND=LDP) TA(ND),TB(ND),TC(ND)
	REAL(KIND=LDP) PSIPREV_MOD(ND),VB(ND),VC(ND)
	REAL(KIND=LDP) HU(ND),HL(ND),HS(ND)
!
! NB: _A denotes that EPS(I) multiples RJ(I)
!     _B denotes that EPS(I) multiples RJ(I+1)
!
	REAL(KIND=LDP) EPS_A(ND),EPS_B(ND)
	REAL(KIND=LDP) EPS_PREV_A(ND),EPS_PREV_B(ND)
	REAL(KIND=LDP) DTAU_BND
	LOGICAL INIT,DO_THIS_TX_MATRIX(NM_TX)
	INTEGER OUT_BC_TYPE
!
! Work Array.
!
	REAL(KIND=LDP) OLD_TX(ND,VDST:VDEND)
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local varables.
!
	INTEGER I,J,K
	LOGICAL USE_EPS
!
! Determine whether we are using the G eddington factor to describe N
! (in terms of H) or whether N is also being described interms of J.
!
! If only G is being used, a faster section of code is excuted.
! Need to check both EPS_A and EPS_PREV_A because only 1 frequency may be using
! RSQN_ON_RSQJ. NO need to check _B, since differ by constant (factor) only.
!
	USE_EPS=.FALSE.			!G Eddington factor only
	DO I=1,ND
	  IF(EPS_A(I) .NE. 0.0_LDP .OR. EPS_PREV_A(I) .NE. 0.0_LDP)USE_EPS=.TRUE.
	END DO
!
	IF(NM_TX .LT. 2 .OR. NM_KI .LT. 2)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Invalid NM in UP_TX_TVX'
	  WRITE(I,*)'NM_TX=',NM_TX
	  WRITE(I,*)'NM_KI=',NM_KI
	  STOP
	END IF
!
! INIT will be true for the very first frequency. We initialize all storage
! locations, even those not in use.
!
	IF(INIT)THEN
	  TX=0.0_LDP
	  TVX=0_LDP
	  OLD_TX=0.0_LDP
	END IF
!
! Now modify the matrices, operating on each matrix (labeled by K) separately.
! We use OLD_TX to store TX( , ,K) at the previous frequency. Only necessary
! when N is being (at least partially) specified in terms of J.
!
	DO K=1,NM_TX
	  IF(DO_THIS_TX_MATRIX(K))THEN
	    IF(USE_EPS .AND. .NOT. INIT)THEN
	      OLD_TX(:,:)=TX(:,:,K)
	      IF(OUT_BC_TYPE .LE. 1)THEN
	        DO J=VDST,VDEND
	          TX(1,J,K)=PSIPREV_MOD(1)*OLD_TX(1,J) + VC(1)*TVX(1,J,K)
	        END DO
	      ELSE
	        DO J=VDST,VDEND
	          TX(1,J,K)=-PSIPREV_MOD(1)*OLD_TX(1,J) - HS(1)*TVX(1,J,K)/DTAU_BND -
	1             (EPS_PREV_A(1)*OLD_TX(1,J)+EPS_PREV_B(1)*OLD_TX(2,J))/DTAU_BND
	        END DO
	      END IF
!
	      DO J=VDST,VDEND
	        DO I=2,ND-1
 	          TX(I,J,K)=PSIPREV_MOD(I)*OLD_TX(I,J)
	1                 + VB(I)*TVX(I-1,J,K) + VC(I)*TVX(I,J,K)
	1                 + ( EPS_PREV_B(I)*OLD_TX(I+1,J)
	1                       - EPS_PREV_A(I-1)*OLD_TX(I-1,J) )
	        END DO
 	        TX(ND,J,K)= PSIPREV_MOD(ND)*OLD_TX(ND,J) + VB(ND)*TVX(ND-1,J,K)
	      END DO
	    ELSE IF(.NOT. INIT)THEN
	      DO J=VDST,VDEND
	        TX(1,J,K)=PSIPREV_MOD(1)*TX(1,J,K) + VC(1)*TVX(1,J,K)
	        DO I=2,ND-1
 	          TX(I,J,K)=PSIPREV_MOD(I)*TX(I,J,K) + VB(I)*TVX(I-1,J,K) + VC(I)*TVX(I,J,K)
	        END DO
 	        TX(ND,J,K)= PSIPREV_MOD(ND)*TX(ND,J,K) + VB(ND)*TVX(ND-1,J,K)
	      END DO
	    END IF
!
	    IF(K .EQ. 1 .OR. K .EQ. 2)THEN
	      TX(:,:,K)=TX(:,:,K)+KI(:,:,K)
	    END IF
!
! Solve the simultaneous equations.
!
	    CALL SIMPTH(TA,TB,TC,TX(1,1,K),ND,ND)
!
	    IF(USE_EPS)THEN
	      DO J=VDST,VDEND
	        DO I=1,ND-1
	           TVX(I,J,K)= HU(I)*TX(I+1,J,K) - HL(I)*TX(I,J,K)
	1               + HS(I)*TVX(I,J,K) +
	1               (EPS_PREV_A(I)*OLD_TX(I,J)-EPS_A(I)*TX(I,J,K)) +
	1               (EPS_PREV_B(I)*OLD_TX(I+1,J)-EPS_B(I)*TX(I+1,J,K))
	        END DO
	      END DO
	    ELSE
	      DO J=VDST,VDEND
	        DO I=1,ND-1
	           TVX(I,J,K)= HU(I)*TX(I+1,J,K) - HL(I)*TX(I,J,K)
	1               + HS(I)*TVX(I,J,K)
	        END DO
	      END DO
	    END IF
!
	    IF(K .EQ. 1)THEN
	      DO J=VDST,VDEND
	        DO I=1,ND-1
	          TVX(I,J,K)=TVX(I,J,K)+RHS_dHdCHI(I,J)
	        END DO
	      END DO
	    END IF
!
	  END IF	!DO_THIS_MATRIX
	END DO		!K
!
	RETURN
	END
