!
! Subroutine to solve a "Block-Banded" system of simultaneous equations
!                 BA. X. = STEQ
! The solution X is returned in STEQ. BA is corrupted.
!
! The "BANDED" matrix can either be 'Diagonal' or 'Tridiagonal'
!
! Routine is currently designed to operate on the SMALL variation matrix
! which is of dimension N.NIV.NUM_BNDS.ND where NUM_BNDS refers to the number
! of bands.
!
! A integer matrix, LNK_IV_TO_F indicates how BA is expanded into BA_BIG,
! which has dimension N,N,NUM_BNDS,ND.
!
! This routine was desined specifically to handle BA_SM, while conserving
! memory.
!
	SUBROUTINE CMF_DIAG_BAND_MPI_V1(STEQ,POPS,SOL_TYPE,FLAG,
	1                DIAG_INDX,N,NION,NUM_BNDS,DST,DEND,ND,
	1                BA_COMPUTED,WR_BA_INV,WR_PRT_INV)
	USE SET_KIND_MODULE
	USE MPI
	IMPLICIT NONE
!
! Altered 13-Aug-2025: Now writes out which g.s. rate equations are replaced (2-Sep-2025).
! Altered 12-Aug-2025: ION equations written to STEQ_VALS.
! Created 24-Sep-2024: Based on CMF_BLKBND_MPI_V1 (Also see cmf_blkband_v3.f for earlier alterations).
!
! 
! The description here is from BLKBAND, which also allowed for a PENTDIAGONAL
! matrix. FOr simplicity,w e retained the smae notation.
!
! Let A[k], B[k], C[k], D[k], and E[k] be the sub-matrices (dimension N*N) of
! the large Block-Pentadiagonal matrix. Then
!                               A[k]=BA( , ,DI-2,K)
!                               B[k]=BA( , ,DI-1,K)
!                               C[k]=BA( , ,DI,K)
!                               D[k]=BA( , ,DI+1,K)
!                               E[k]=BA( , ,DI+2,K)
!
! where DI=K if ND=NUM_BNDS, else DI=DIAG_INDX.
!
! Likewise we have              L[k] =STEQ( , K)
!
! Elimination defined by:
!                  phi[1]= C^{-1) . L[1]
!                  alpha[1]= C^{-1) . D[1]
!                  gamma[1]= C^{-1) . E[1]
!                  SOL[1] = phi[1] - alpha[1].SOL[2] - gamma[1].SOL[3]
! Then
!                  m[k]=C[k]-B[k].alpha[k-1] - A[k].( gamma[k-2]
!                                            - alpha[k-2].alpha[k-1] )
!                  phi[k]=m[k]^{-1}(L[k] - B[k].phi[k-1]- A[k].( phi[k-2]
!                                            - alpha[k-2].phi(k-1) )
!                  alpha[k]=m[k]^{-1}(D[k]-B[k].gamma[k-1] +
!                                            + A[k].alpha[k-2].gamma[k-1]
!                  gamma[k]=m[k]^{-1} . E[k]
! with
!                  SOL[k]= phi[k] - alpha[k]. SOL[k+1] - gamma[k].SOL[k+2]
!
! Note that for k=2, we take A[k]==0.
!           for k=ND-1, we hace E[ND-1]==0 and hence gamma[ND-1]==0
!           for k=ND, we have D[ND]=E[ND]=0 and hence alpha[ND]=gamma[ND]==0.
! Thus
!                  m[ND]=C[ND]-B[ND].alpha[ND-1] - A[ND].( gamma[ND-2]
!                                            - alpha[ND-2].alpha[ND-1] )
!                  phi[ND]=m[ND]^{-1}(L[ND] - B[ND].phi[ND-1]
!                                            - A[ND].( phi[ND-2]
!                                            - alpha[ND-2].phi(ND-1) )
! and              SOL(ND)=phi[ND]
!
! The back substitution is simly obtained from the exprssion for SOL[k]
! since SOL[ND] is now known.
!
! The tridiagonal case can be obtained from these expressinons by setting
! A[k]=E[k]=0, and A[1]=D[ND]=0. When the tridiagonal option is indicated,
! these matrices are ASSUMED to be zero.
! 
!**************************************************************************
! SUBROUTINES CALLED
! ******************
!
!          DGEMV(TRANS,M,N,ALPHA,A,IDA,X,INCX,BETA,Y,INCY)
!
! BLAS2 routine.
!
! Compute Y <--- alpha A.X + beta Y
!
! Where
!
!  	A = Double preciosion A(IDA,N)
!       M = Number of rows in A
!	N = Number of columns in A
!       IDA = First dimension of A [ > max(1,m) ]
!
!	ALPHA, BETA are REAL
!
!       X=REAL   X(1+(N-1}INCX) if TRANS='N'
!       X=REAL   X(1+{M-1}INCY) if TRANS='T'
!	INCX = Stride for X (>0, generally 1)
!
!       Y=REAL   Y(1+(M-1}INCX) if TRANS='N'
!       Y=REAL   Y(1+(N-1}INCX) if TRANS='T'
!	INCY = Stride for Y (>0, generally 1)
! 
! ****************************************************************************
!
!          DGETRF(M,N,A,LDA,IPRIV,INFO)
!
! LAPACK routine.
!
! Perform the LU decomposition of the general matrix A where:
!
!		M   = Number of rows in A
!		N   = Number of columns in A
!       	IDA = First dimension of A [ > max(1,m) ]
!  		A   = Double precision A(IDA,N)
!
!	      IPIV  = Integer work vector with pivot
!                   = IPIV(min[m,n])
!
!             Integer INFO = 0 (successful exit)
!                          = -i (i th argument has illeagal value)
!                          = i (Pivot A(i,i) is exactly zero)
!
! LU decomposition is stored in A.
!
!******************
!
!                   DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
!
! LAPACK routine.
!
! Solves for the solution of A^{-1} . B . DGETRF must have
! previously been called, and A and IPIV must be passed to DGETRS unchanged.
! The solution is returned in B.
!
!		N      = Order of A
!		NRHS   = Number of RHS
!       	IDA    = First dimension of A [ > max(1,m) ]
!  		A      = Double precision A(IDA,N)
!
!	      IPIV  = Integer work vector with pivot
!                   = IPIV(min[m,n])
!
!  		B      = Double precision B(IDB,NRHS)
!
!             Integer INFO = 0 (successful exit)
!                          = -i (i th argument has illeagal value)
! 
!*******************
!
!	    DGEEQU(M,N,A,LDA,ROW_SF,COL_SF,ROW_CND,COL_CND,MAX_VAL,IFAIL)
!
! LAPACK routine
!
! Compute the Equibriation factors for the M by N matrix A whose first
! dimension is LDA. These factors help improve the condition of the matrix.
!
!
! The equilibrated matrix is found from
!
!         E(I,J)=ROW_SF(I)*A(I,J)*COL_SF(J)
!
! To solve  a set of simultaneous equtions multiply RHS by ROW_SF before
! solution, and by COL_SF after solution.
!
!******************************************************************************
! 
!
	INTEGER N,NION,NUM_BNDS,DIAG_INDX
	INTEGER DST,DEND,ND
	REAL(KIND=LDP) STEQ_STORE(N,DST:DEND)
        REAL(KIND=LDP) STEQ(N,DST:DEND)
        REAL(KIND=LDP) POPS(N,ND)
	CHARACTER(LEN=*) SOL_TYPE
	LOGICAL FLAG
!
	LOGICAL BA_COMPUTED
	LOGICAL WR_BA_INV
	LOGICAL WR_PRT_INV
!
! Local variables
!
        REAL(KIND=LDP), ALLOCATABLE :: C_MAT(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: ORIG_POPS(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: RUB(:,:)	      	      !Not accessed when passed.
!
        REAL(KIND=LDP) ROW_SF(N)
        REAL(KIND=LDP) COL_SF(N)
        REAL(KIND=LDP) ROW_CND,COL_CND,MAX_VAL
	INTEGER IPIVOT(N)
	INTEGER RUB_VEC(N)
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL REPLACE_EQ(NION)
	LOGICAL ZERO_STEQ(N)
	LOGICAL USE_PASSED_REP
!
	INTEGER DEPTH_INDX
	INTEGER BAND_INDX
!
        INTEGER I,J,K,JJ
        INTEGER IOS,IFAIL
	INTEGER, PARAMETER :: NSNG=1
	CHARACTER(LEN=10) DESC
	LOGICAL FIRST_MATRIX,LAST_MATRIX
!
        CHARACTER(LEN=1), PARAMETER :: NO_TRANS='N'
	CHARACTER(LEN=50) TMP_STR
!
	INTEGER IERR
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	LUER=ERROR_LU()
!
!
!
! If the BA matrix was not computed, it means that its inverse must have been previously
! computed. IF WR_BA_INV is true, we can use the previous inverse, which is stored on disk.
!
	IF(SOL_TYPE(1:4) .EQ. 'DIAG' .AND. .NOT. BA_COMPUTED .AND. WR_BA_INV)THEN
	  FIRST_MATRIX=.TRUE.
	  LAST_MATRIX=.FALSE.
	  USE_PASSED_REP=.TRUE.
          ALLOCATE (C_MAT(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (RUB(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (ORIG_POPS(N,ND),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in CMF_BLKBAND'
            WRITE(LUER,*)'Unable to allocate C_MAT etc'
            WRITE(LUER,*)'STAT=',IOS
            STOP
          END IF
	  C_MAT=0.0_LDP; RUB=0.0_LDP; ORIG_POPS=0.0_LDP
!
	  DO K=DST,DEND
!
! Read in  LU decomposition, and other necessary data.
!
	    DEPTH_INDX=K
	    CALL READ_BCD_MAT(RUB,C_MAT,RUB,ROW_SF,COL_SF,IPIVOT,
	1          ORIG_POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'C')
!
! Map the small BA rray onto the full BA array (one depth at a time).
! Really only need to do this for STEQ. We pass RUB, since we don't
! want to corrupt C_MAT. RUB is not used. NB: In earlier version we used
! D_MAT in place of RUB.
!
	    IF(K .EQ. ND)LAST_MATRIX=.TRUE.
	    CALL GENERATE_FULL_MATRIX_V3(
	1         RUB,STEQ(1,K),POPS,REPLACE_EQ,ZERO_STEQ,
	1         N,ND,NION,NUM_BNDS,
	1         DIAG_INDX,DIAG_INDX,DEPTH_INDX,
	1         FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	    FIRST_MATRIX=.FALSE.
!
	    STEQ_STORE(:,K)=STEQ(:,K)
!
! Correct STEQ for the eqilibrization of C
!
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
	    END DO
!
! Now perform the solution.
!
	    CALL DGETRS(NO_TRANS,N,NSNG,C_MAT,N,IPIVOT,STEQ(1,K),N,IFAIL)
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	    END DO
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRS_DIAG'
	      GOTO 9999
	    END IF
!
! Correct for the use of old POPS populations in doing the scaling.
!
	   DO I=1,N
	     STEQ(I,K)=STEQ(I,K)*(ORIG_POPS(I,K)/POPS(I,K))
	   END DO
!
	  END DO
          FLAG=.TRUE.
          DEALLOCATE (C_MAT,ORIG_POPS,RUB)
	  CALL WRITE_STEQ_ION(NION,DST,DEND,ND)
	  CALL WRITE_REPLACE(NION,DST,DEND,ND)
	  CALL WR2D_GATH_MPI_V1(STEQ_STORE,N,DST,DEND,ND,'STEQ_ARRAY','*',L_TRUE,16)
!
	  RETURN
	END IF
!
!
! Check to see if we only require the diagonal solution. For this special
! case, the solutions at each depth are independent. We must compute the
! inverse, which will be saved on DISK if WR_BA_INV is true.
!
	IF(SOL_TYPE(1:4) .EQ. 'DIAG')THEN
!
	  FIRST_MATRIX=.TRUE.
	  LAST_MATRIX=.FALSE.
	  USE_PASSED_REP=.FALSE.
          ALLOCATE (C_MAT(N,N),STAT=IOS)
          ALLOCATE (RUB(N,N),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in CMF_BLKBAND'
            WRITE(LUER,*)'Unable to allocate C_MAT etc'
            WRITE(LUER,*)'STAT=',IOS
            STOP
	  END IF
	  C_MAT=0.0_LDP
!
          DO K=DST,DEND
!
! Map the small BA rray onto the full BA array (one depth at a time).
!
	    DEPTH_INDX=K
	    IF(K .EQ. ND)LAST_MATRIX=.TRUE.
	    CALL GENERATE_FULL_MATRIX_V3(
	1         C_MAT,STEQ(1,K),POPS,REPLACE_EQ,ZERO_STEQ,
	1         N,ND,NION,NUM_BNDS,
	1         DIAG_INDX,DIAG_INDX,DEPTH_INDX,
	1         FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	    FIRST_MATRIX=.FALSE.
	    STEQ_STORE(:,K)=STEQ(:,K)
!
! Perform the LU decomposition using DGETRF. We first equilibrize the matrix
! using DGEEQU sot the the maximum row and column values are approximately
! unity.
!
	    CALL TUNE(1,'DGEEQU')
	    CALL DGEEQU(N,N,C_MAT,N,ROW_SF,COL_SF,
	1               ROW_CND,COL_CND,MAX_VAL,IFAIL)
	    CALL TUNE(2,'DGEEQU')
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
              DO I=1,N
	        C_MAT(I,J)=C_MAT(I,J)*ROW_SF(I)*COL_SF(J)
	      END DO
	    END DO
	    CALL TUNE(1,'DGERTF')
	    CALL DGETRF(N,N,C_MAT,N,IPIVOT,IFAIL)
	    CALL TUNE(2,'DGERTF')
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRF_DIAG'
	      WRITE(LUER,*)'Error in CMF_BLKBAND_V3'
	      WRITE(LUER,*)'Unable to get solution at depth',K
	      WRITE(LUER,*)'Setting fractional corrections to zero'
	      WRITE(LUER,*)'IFAIL=',IFAIL
	      STEQ(:,K)=0.0_LDP
	      GOTO 500			!9999
	    END IF
	    IF(WR_BA_INV)THEN
	      CALL WRITE_BCD_MAT(RUB,C_MAT,RUB,ROW_SF,COL_SF,IPIVOT,
	1              POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'C')
	    END IF
!
! Now perform the solution.
!
	    CALL TUNE(1,'DGETRS')
	    CALL DGETRS(NO_TRANS,N,NSNG,C_MAT,N,IPIVOT,STEQ(1,K),N,IFAIL)
	    CALL TUNE(2,'DGETRS')
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	    END DO
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRS_DIAG'
	      WRITE(LUER,*)'Error in CMF_BLKBAND_V3'
	      WRITE(LUER,*)'Unable to get solution at depth',K
	      WRITE(LUER,*)'Setting fractional corrections to zero'
	      WRITE(LUER,*)'IFAIL=',IFAIL
	      STEQ(:,K)=0.0_LDP
	    END IF
500	    CONTINUE
	  END DO
          FLAG=.TRUE.
          DEALLOCATE (C_MAT,RUB)
	  CALL WRITE_STEQ_ION(NION,DST,DEND,ND)
	  CALL WRITE_REPLACE(NION,DST,DEND,ND)
	  CALL WR2D_GATH_MPI_V1(STEQ_STORE,N,DST,DEND,ND,'STEQ_ARRAY','*',L_TRUE,16)
	  RETURN
	END IF
!
! Error handling routine.
!
9999	CONTINUE
	WRITE(LUER,*)'Error in LINPAC (or BLAS) routine',DESC,' in CMF_BLKBAND'
	WRITE(LUER,100)K,IFAIL
100	FORMAT(1x,'depth=',I3,10x,'IFAIL=',I3)
	FLAG=.FALSE.
	RETURN
!
	END
