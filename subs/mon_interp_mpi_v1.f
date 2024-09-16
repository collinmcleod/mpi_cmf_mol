!
! Subroutine to interpolate an array onto a new grid. The grid vector must be
! either a monotonically decreasing or increasing function. A modified cubic
! polynomial is used to do the interpolation. Instead of using
! the exact cubic estimates for the first derivative at the two nodes,
! we use revised estimates which insure that the interpolating function
! is monotonic in the interpolating interval.
!
! Disadvantages: The interpolating weights can only be defined when the
!                function is known. In principal could use these modified
!                first derivatives to compute an accurate integration
!                formulae. However, the integration weights cannot be defined
!                independently of the function values, as desired in many
!                situations.
!
! Ref: Steffen. M, 1990, A/&A, 239, 443-450
!
	SUBROUTINE MON_INTERP_MPI_V1(INTERPOLATED_MAT,NEW_R,NX,LIN_END,
	1                R_PNT,DST,DEND,IN_MAT,R,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER LIN_END
	INTEGER NX
	INTEGER ND
	INTEGER DST,DEND
!
	REAL(KIND=LDP) INTERPOLATED_MAT(NX,LIN_END)
	REAL(KIND=LDP) NEW_R(NX)
	INTEGER R_PNT(NX)
!
	REAL(KIND=LDP) IN_MAT(ND,LIN_END)
	REAL(KIND=LDP) R(ND)
!
	REAL(KIND=LDP) S(ND)		!Slopes
	REAL(KIND=LDP) H(ND)
!
	REAL(KIND=LDP), PARAMETER :: ONE=1.0_LDP
!
	INTEGER I,J,ML
	INTEGER IMIN,IMAX
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) A(ND)
	REAL(KIND=LDP) B(ND)
	REAL(KIND=LDP) C(ND)
	REAL(KIND=LDP) D(ND)		!Used for derivative at I.
	REAL(KIND=LDP) E(ND)
	REAL(KIND=LDP) SGN
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Determine intervals and slopes to minimize computational effort.
!
	DO I=1,ND-1
	  H(I)=R(I+1)-R(I)
	END DO
!
! Loop over frequency space.
!
	DO ML=1,LIN_END
!
! Compute the slopes.
!
	  IMIN=MAX(1,R_PNT(DST)-2);  IMAX=MIN(ND,R_PNT(DEND)+1)
	  DO I=IMIN,MIN(IMAX,ND-1)
	    S(I)=(IN_MAT(I+1,ML)-IN_MAT(I,ML))/H(I)
	  END DO
!
! Compute the first derivatives at node I.
!
          IF(IMIN .EQ. 1)D(1)=S(1) +(S(1)-S(2))*H(1)/(H(1)+H(2))
	  DO I=MAX(2,IMIN),MIN(IMAX,ND-1)
            D(I)=(S(I-1)*H(I)+S(I)*H(I-1))/(H(I-1)+H(I))
	  END DO
!
! Adjust first derivatives so that function is monotonic  in each interval.
!
	  IF(IMIN .EQ. 1)D(1)=( SIGN(ONE,S(1))+SIGN(ONE,D(1)) )*MIN(ABS(S(1)),0.5_LDP*ABS(D(1)))
	  DO I=MAX(2,IMIN),MIN(IMAX,ND-1)
	    D(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5_LDP*ABS(D(I)))
	  END DO
          D(ND)=S(ND-1)+(S(ND-1)-S(ND-2))*H(ND-1)/(H(ND-2)+H(ND-1))
	  D(ND)=( SIGN(ONE,S(ND-1))+SIGN(ONE,D(ND)) )*
	1      MIN(ABS(S(ND-1)),0.5_LDP*ABS(D(ND)))
!
! Determine the coefficients of the monotonic cubic polynomial.
!
! If T1=X-R(I) then
!             Y=A(I)*T1^3 + B(I)*T1^3 + C(I)*T1 + E(I)
!
	  DO I=IMIN,IMAX
            A(I)=(D(I)+D(I+1)-2.0_LDP*S(I))/H(I)/H(I)
	    B(I)=(3.0_LDP*S(I)-2.0_LDP*D(I)-D(I+1))/H(I)
	    C(I)=D(I)
	    E(I)=IN_MAT(I,ML)
	  END DO
!
! Perform the interpolations.
!
	  DO J=R_PNT(DST),R_PNT(DEND)
	    I=R_PNT(J)
	    T1=(NEW_R(J)-R(I))
            INTERPOLATED_MAT(J,ML)=((A(I)*T1+B(I))*T1+C(I))*T1+E(I)
	  END DO
!
	END DO
!
	RETURN
	END
