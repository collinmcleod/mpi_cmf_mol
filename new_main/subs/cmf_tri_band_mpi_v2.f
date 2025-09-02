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
	SUBROUTINE CMF_TRI_BAND_MPI_V2(STEQ,POPS,SOL_TYPE,FLAG,
	1                  DIAG_INDX,N,NION,NUM_BNDS,DST,DEND,ND,
	1                  BA_COMPUTED,WR_BA_INV,WR_PRT_INV,TRI_SOL_OPTIONS)
	USE SET_KIND_MODULE
	USE MPI
	IMPLICIT NONE
!
! Altered 13-Aug-2025: Now writes out which g.s. rate equations are replaced (2-Sep-2025).
! Altered 12-Aug-2025: ION equations written to STEQ_VALS.
! Altered 10-Jun-2024: Now use diagonal solution when tridiagonal solution fails.
!                         Cleaning done earlier -- 12-May-2025.
! Created 24-Sep-2023: Based on CMF_BLK_BAND_V3
!
! 
! The description here is from TRI_BAMD_MPI_V1, which also allowed for a PENTDIAGONAL
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
	CHARACTER(LEN=*) SOL_TYPE
	CHARACTER(LEN=*) TRI_SOL_OPTIONS
!
	REAL(KIND=LDP) STEQ_STORE(N,DST:DEND)
        REAL(KIND=LDP) STEQ(N,DST:DEND)
        REAL(KIND=LDP) POPS(N,ND)
	LOGICAL FLAG
!
	LOGICAL BA_COMPUTED
	LOGICAL WR_BA_INV
	LOGICAL WR_PRT_INV
!
! Local variables
!
        REAL(KIND=LDP), SAVE, ALLOCATABLE :: OLD_EST(:,:)
!
        REAL(KIND=LDP), ALLOCATABLE :: B_MAT(:,:,:)
        REAL(KIND=LDP), ALLOCATABLE :: C_MAT(:,:,:)
        REAL(KIND=LDP), ALLOCATABLE :: D_MAT(:,:,:)
        REAL(KIND=LDP), ALLOCATABLE :: RUB(:,:)
!
        REAL(KIND=LDP), ALLOCATABLE :: ORIG_POPS(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: NEW_EST(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: ROW_SF(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: COL_SF(:,:)
        REAL(KIND=LDP) ROW_CND,COL_CND,MAX_VAL
	INTEGER, ALLOCATABLE ::  IPIVOT(:,:)
	REAL(KIND=LDP) ERR_EST(ND)
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) RELAX_PARAM
!
	INTEGER, PARAMETER :: MAX_NUM_ITS=100
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL REPLACE_EQ(NION,DST:DEND)
	LOGICAL ZERO_STEQ(N,DST:DEND)
	LOGICAL USE_PASSED_REP
	LOGICAL NAN_PRES
!
	INTEGER DEPTH_INDX
	INTEGER BAND_INDX
!
	INTEGER LU_IT
	INTEGER, SAVE :: ENTRY_COUNTER=0
!
        INTEGER I,J,K,JJ
        INTEGER IOS,IFAIL,IERR
	INTEGER IT_COUNTER
	INTEGER NG_CNT
!
	CHARACTER*10 DESC
	LOGICAL ANONZERO,ENONZERO
	LOGICAL FIRST_MATRIX,LAST_MATRIX
	LOGICAL WR_D_MAT
!
        REAL(KIND=LDP),      PARAMETER :: DP_NEG_ONE=-1.0_LDP
        REAL(KIND=LDP),      PARAMETER :: DP_ONE=1.0_LDP
        INTEGER,   PARAMETER :: INT_ONE=1
        INTEGER,   PARAMETER :: NSNG=1
        CHARACTER*1, PARAMETER :: NO_TRANS='N'
	CHARACTER(LEN=3) OUT_TYPE
	CHARACTER(LEN=100) TMP_STR
	LOGICAL UPDATE_RELAX
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
!	INCLUDE 'mpif.h'
!
	RELAX_PARAM=0.8_LDP
	LUER=ERROR_LU()
	NG_CNT=30
	IF( INDEX(TRI_SOL_OPTIONS,'UPDATE_RELAX') .NE. 0 )UPDATE_RELAX=.TRUE.
	IF(MYPE .EQ. 0)THEN
	  WRITE(6,*)'UPDATE_RELAX=',UPDATE_RELAX,TRIM(TRI_SOL_OPTIONS)
	  FLUSH(UNIT=6)
	END IF
!
! If we get here, we must be performing a TRI diagonal solution.
!
	IF(SOL_TYPE(1:3) .EQ. 'TRI')THEN
	  IF(NUM_BNDS .LT. 3)THEN
	    WRITE(LUER,*)'Error in CMF_TRI_BAMD_MPI_V2 : NUM_BNDS too small for solution type'
	    FLAG=.FALSE.
	    STOP
	  END IF
	ELSE
	  WRITE(LUER,*)'Error in CMF_TRI_BAMD_MPI_V2 - invalid SOL_TYPE - ','SOLTYPE= ',SOL_TYPE
	  STOP
        END IF
!
        IF(.NOT. ALLOCATED(OLD_EST))THEN
	  ALLOCATE (OLD_EST(N,ND),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in CMF_TRI_BAMD_MPI_V2 -- unable to allocate OLD_EST etc'
            WRITE(LUER,*)'STAT=',IOS
            STOP
          END IF
	  OLD_EST=0.0_LDP
	END IF
!
! Perform the TRIDIAGONAL solution. If we reach here, we need to redo
! the LU decomposition of BA.
!
	CALL TUNE(1,'TRI_ALLOCATION')
	  ALLOCATE (B_MAT(N,N,DST:DEND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (C_MAT(N,N,DST:DEND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (D_MAT(N,N,DST:DEND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (RUB(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (COL_SF(N,DST:DEND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (ROW_SF(N,DST:DEND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (IPIVOT(N,DST:DEND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (NEW_EST(N,DST:DEND),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in CMF_TRI_BAMD_MPI_V2'
            WRITE(LUER,*)'Unable to allocate D_MAT etc'
            WRITE(LUER,*)'STAT=',IOS
            STOP
          END IF
	  B_MAT=0.0_LDP; C_MAT=0.0_LDP; D_MAT=0.0_LDP
	CALL TUNE(2,'TRI_ALLOCATION')
!
	IF(.NOT. BA_COMPUTED .AND. WR_BA_INV)THEN
          ALLOCATE (ORIG_POPS(N,DST:DEND),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in CMF_TRI_BAMD_MPI_V2'
            WRITE(LUER,*)'Unable to allocate D_MAT etc: STAT=',IOS
             STOP
          END IF
!
	  USE_PASSED_REP=.TRUE.
	  DO K=DST,DEND
	    IF(K .EQ. ND)LAST_MATRIX=.TRUE.
	    DEPTH_INDX=K
!
! Read in LU decompostion of C, and the original BD matrices. This must be done before the
! call to GENERATE sice we need REPALCE_EQ and ZERO_STEQ.
!
	    CALL TUNE(1,'TRI_RD_BCD')
	    OUT_TYPE='BCD'
            CALL READ_BCD_MAT(B_MAT(:,:,K),C_MAT(:,:,K),D_MAT(:,:,K),ROW_SF(:,K),COL_SF(:,K),
	1          IPIVOT(:,K),ORIG_POPS(:,K),REPLACE_EQ(:,K),ZERO_STEQ(:,K),
	1          N,NION,DEPTH_INDX,OUT_TYPE)
	    CALL TUNE(2,'TRI_RD_BCD')
!
! Don't need D_MAT -- just the STEQ array.
!
	    CALL TUNE(1,'TRI_GEN')
	    CALL GENERATE_FULL_MATRIX_V3(
	1         RUB,STEQ(1,K),POPS,REPLACE_EQ(:,K),ZERO_STEQ(:,K),
	1         N,ND,NION,NUM_BNDS,
	1         DIAG_INDX,DIAG_INDX,DEPTH_INDX,
	1         FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	     FIRST_MATRIX=.FALSE.
	    CALL TUNE(2,'TRI_GEN')
!
	     STEQ_STORE(:,K)=STEQ(:,K)
	  END DO
	ELSE
!
! 
!
! Solve the tridiagonal equations.
!
	  FIRST_MATRIX=.TRUE.
	  LAST_MATRIX=.FALSE.
	  USE_PASSED_REP=.FALSE.
!
	  DO K=DST,DEND
!
! Map the small BA rray onto the full BA array (one depth at a time).
!    To check for NaNs:
!	    I=N*N; CALL CHECK_VEC_NAN(C_MAT(:,:,K),I,'CMAT(STOP)',NAN_PRES)
!
	    DEPTH_INDX=K
	    CALL TUNE(1,'TRI_GEN')
	    CALL GENERATE_FULL_MATRIX_V3(
	1           C_MAT(:,:,K),STEQ(1,K),POPS,REPLACE_EQ(:,K),ZERO_STEQ(:,K),
	1           N,ND,NION,NUM_BNDS,
	1           DIAG_INDX,DIAG_INDX,DEPTH_INDX,
	1           FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	    FIRST_MATRIX=.FALSE.
	    IF(K .NE. 1)THEN
	      IF(K .EQ. ND)LAST_MATRIX=.TRUE.
	      BAND_INDX=DIAG_INDX-1
	      CALL GENERATE_FULL_MATRIX_V3(
	1             B_MAT(:,:,K),STEQ(1,K),POPS,REPLACE_EQ(:,K),ZERO_STEQ(:,K),
	1             N,ND,NION,NUM_BNDS,
	1             BAND_INDX,DIAG_INDX,DEPTH_INDX,
	1             FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	    END IF
	    IF(K .NE. ND)THEN
	      BAND_INDX=DIAG_INDX+1
	      CALL GENERATE_FULL_MATRIX_V3(
	1             D_MAT(:,:,K),STEQ(1,K),POPS,REPLACE_EQ(:,K),ZERO_STEQ(:,K),
	1             N,ND,NION,NUM_BNDS,
	1             BAND_INDX,DIAG_INDX,DEPTH_INDX,
	1             FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	    END IF
	    STEQ_STORE(:,K)=STEQ(:,K)
	    CALL TUNE(2,'TRI_GEN')
!
! Do LU decompostion of m[k]. We first equilibrilze C_MAT. DGEEQU takes
! a minimal amount of time compared to DGETRF.
!
	     CALL DGEEQU(N,N,C_MAT(:,:,K),N,ROW_SF(:,K),COL_SF(:,K),
	1               ROW_CND,COL_CND,MAX_VAL,IFAIL)
!
	     CALL TUNE(1,'TRI_DGETRF')
	     DO J=1,N
               DO I=1,N
	         C_MAT(I,J,K)=C_MAT(I,J,K)*ROW_SF(I,K)*COL_SF(J,K)
	       END DO
	     END DO
	     CALL DGETRF(N,N,C_MAT(:,:,K),N,IPIVOT(:,K),IFAIL)
	     CALL TUNE(2,'TRI_DGETRF')
	     IF(IFAIL .NE. 0)THEN
	       DESC='DGETRF_2'
	       GOTO 9999
	     END IF
!
	     IF(WR_BA_INV)THEN
	       CALL TUNE(1,'WR_BA_INV')
	       OUT_TYPE='BCD'
	       CALL WRITE_BCD_MAT(B_MAT(:,:,K),C_MAT(:,:,K),D_MAT(:,:,K),
	1              ROW_SF(:,K),COL_SF(:,K),
	1              IPIVOT(:,K),POPS(:,K),REPLACE_EQ(:,K),ZERO_STEQ(:,K),
	1              N,NION,DEPTH_INDX,OUT_TYPE)
	       CALL TUNE(2,'WR_BA_INV')
	     END IF
!
	  END DO
	END IF
	CALL TUNE(1,'TRI_GATH')
	CALL WRITE_STEQ_ION(NION,DST,DEND,ND)
	CALL WRITE_REPLACE(NION,DST,DEND,ND)
	CALL WR2D_GATH_MPI_V1(STEQ_STORE,N,DST,DEND,ND,'STEQ_ARRAY','*',L_TRUE,16)
	CALL TUNE(2,'TRI_GATH')
!
! If we have not computed BA, we are probably in a regime where we are converging. Since
! convergence is slow, previous estimates of the solution will generally be very close to
! the current estimates, so we can use them as a starting guess. We check the maximum 
! value of OLD_EST to check whether this model has been restarted (in which case OLD will
! not be available).
! 
	T1=MAXVAL(OLD_EST)
	IF(.NOT. BA_COMPUTED .AND. T1 .NE. 0.0_LDP)THEN
	  DO K=DST,DEND
	    STEQ(:,K)=STEQ_STORE(:,K)
	    IF(K .EQ. 1)THEN
	      DO J=1,N
	        DO I=1,N
	          STEQ(I,K)=STEQ(I,K)-D_MAT(I,J,K)*OLD_EST(J,K+1)
	        END DO
	      END DO  
	    ELSE IF(K .EQ. ND)THEN
	      DO J=1,N
	        DO I=1,N
	          STEQ(I,K)=STEQ(I,K)-B_MAT(I,J,K)*OLD_EST(J,K-1)
	        END DO
	      END DO  
	    ELSE
	      DO J=1,N
	        DO I=1,N
	          STEQ(I,K)=STEQ(I,K)-B_MAT(I,J,K)*OLD_EST(J,K-1)-
	1                                   D_MAT(I,J,K)*OLD_EST(J,K+1)
	        END DO
	      END DO  
	    END IF
	  END DO
	END IF
!
	IF(MYPE .EQ. 0)THEN
	  CALL GET_LU(LU_IT,'LU for iteration information on CMF_TRI_BAND_MPI_V1')
	  IF(ENTRY_COUNTER .EQ. 0)THEN
	    OPEN(LU_IT,FILE='TRI_BAND_IT_INFO',STATUS='UNKNOWN')
	  ELSE
	    OPEN(LU_IT,FILE='TRI_BAND_IT_INFO',STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
	  END IF
	  ENTRY_COUNTER=ENTRY_COUNTER+1
	  WRITE(LU_IT,'(/,A,I4)')   ' ENTRY_COUNTER=',ENTRY_COUNTER
	  WRITE(LU_IT,'(A,3X,L1,/)')'   BA_COMPUTED=',BA_COMPUTED
	END IF
!
	CALL TUNE(1,'IT_COUNTER')
	IT_COUNTER=1
	DO WHILE(IT_COUNTER .LE. MAX_NUM_ITS)
!
	  ERR_EST=0.0_LDP
	  CALL TUNE(1,'TRI_DGETRS')
	  DO K=DST,DEND
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*ROW_SF(J,K)
	    END DO
	    CALL DGETRS(NO_TRANS,N,NSNG,C_MAT(:,:,K),N,IPIVOT(:,K),STEQ(:,K),N,IFAIL)
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRS_4'
	      CALL TUNE(2,'TRI_DGETRS')
	      CALL TUNE(2,'IT_COUNTER')
	      GOTO 9999
	    END IF
!
! NB: The term 1.0E-16 is needed as STEQ(NT,J) will be zero when we hold T fixed.
!
	    DO J=1,N
	      NEW_EST(J,K)=STEQ(J,K)*COL_SF(J,K)
	    END DO
	    IF(RELAX_PARAM .LT. 0.1_LDP)THEN
	       
	    ELSE IF(IT_COUNTER .GE. 2)THEN
	      DO J=1,N
	        NEW_EST(J,K)=OLD_EST(J,K)+RELAX_PARAM*(NEW_EST(J,K)-OLD_EST(J,K))
	        ERR_EST(K)=MAX(ERR_EST(K),ABS(NEW_EST(J,K)-OLD_EST(J,K)) / 
	1                               (ABS(NEW_EST(J,K))+ABS(OLD_EST(J,K))+1.0E-20_LDP))
	      END DO
	    ELSE
	      ERR_EST(K)=1.0_LDP
	    END IF
	  END DO
	  CALL TUNE(2,'TRI_DGETRS')
!
!	  CALL ACCEL_TRI_SOL_IT(NEW_EST,DST,DEND,ND,N,IT_COUNTER,NG_CNT)
!	  IF(IT_COUNTER .GT. NG_CNT)NG_CNT=NG_CNT+20
!
	  OLD_EST=0.0_LDP; OLD_EST(:,DST:DEND)=NEW_EST(:,DST:DEND); I=N*ND
	  CALL MPI_ALLREDUCE(MPI_IN_PLACE,OLD_EST,I,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
!
	  CALL GATHER_SELF_VEC_MPI_V1(ERR_EST,DST,DEND,ND)
	  J=0; CALL MPI_BCAST(ERR_EST,ND,MPI_DOUBLE_PRECISION,J,MPI_COMM_WORLD,IERR)
	  IF(MYPE .EQ. 0)THEN
	    WRITE(TMP_STR,'(I4)')IT_COUNTER
	    TMP_STR='IT_CNT='//TRIM(TMP_STR)
	    CALL WRITE_VEC(ERR_EST,ND,TMP_STR,275)
	    FLUSH(UNIT=275)
	  END IF
!
	  IF(MAXVAL(ERR_EST) .LT. 1.0E-06_LDP)THEN
	    STEQ(:,DST:DEND)=NEW_EST(:,DST:DEND)
	    IF(MYPE .EQ. 0)WRITE(6,*)'Solution of tri-diaginal equations converged'
	    EXIT
!
	  ELSE IF(IT_COUNTER .EQ. MAX_NUM_ITS .AND. MAXVAL(ERR_EST) .LT. 1.0E-02_LDP)THEN
	    WRITE(6,*)'Solution of tri-diaginal equations may not be full converged'
	    STEQ(:,DST:DEND)=NEW_EST(:,DST:DEND)
	    EXIT
!
	  ELSE IF(IT_COUNTER .EQ. MAX_NUM_ITS)THEN
!
! Use diagonal solution.
!
	    WRITE(6,'(/,A)')' Solution of tri-diagonal equations failed to converged'
	    WRITE(6,'(A,/)')' Using diagonal solution'
	    DO K=DST,DEND
	      STEQ(:,K)=STEQ_STORE(:,K)
	      DO J=1,N
	        STEQ(J,K)=STEQ(J,K)*ROW_SF(J,K)
	      END DO
	      CALL DGETRS(NO_TRANS,N,NSNG,C_MAT(:,:,K),N,IPIVOT(:,K),STEQ(:,K),N,IFAIL)
	      IF(IFAIL .NE. 0)THEN
	        DESC='DGETRS_4'
	        CALL TUNE(2,'TRI_DGETRS')
	        CALL TUNE(2,'IT_COUNTER')
	        GOTO 9999
	      END IF
!
	      DO J=1,N
	        NEW_EST(J,K)=STEQ(J,K)*COL_SF(J,K)
	      END DO
	      STEQ(:,DST:DEND)=NEW_EST(:,DST:DEND)
            END DO
	    EXIT
!
! Check if converging. If not converging we restart the iterative procedure, and lower
! the relaxation parameter.
!
	  ELSE IF(IT_COUNTER .GT. 30 .AND. MAXVAL(ERR_EST) .GT. 0.7_LDP .AND. UPDATE_RELAX)THEN
	    RELAX_PARAM=RELAX_PARAM/2.0_LDP
	    OLD_EST=0.0_LDP
	    IT_COUNTER=0
	    IF(MYPE .EQ. 0)THEN
	      WRITE(6,'(A,F5.2)')' Updated relaxation paraemeter in CMF_TRI_BAMD_MPI_V2 to:',RELAX_PARAM
	    END IF
	  END IF
	  IF(MYPE .EQ. 0)THEN
	    WRITE(LU_IT,'(1X,A,I4,A,F10.5,2X,ES14.4,F6.2)')'Maximum error and correction on CM_TRI_BAND iteration',
	1                       IT_COUNTER,' is (in %): ',200.0_LDP*MAXVAL(ERR_EST),MAXVAL(OLD_EST),RELAX_PARAM
	  END IF
!	  IF(MYPE .EQ. 0)THEN
!	    CALL WR2D_V2(OLD_EST,N,ND,'STEQ_ARRAY','*',L_TRUE,450)
!	  END IF
!
!	  PREV_EST(:,:,2:4)=PREV_EST(:,:,1:3)
!	  PREV_EST(:,:,1)=NEW_EST(:,:,1)
!	  CALL ACCEL_IT(OLD_EST,PREV_EST,DST,DEND,ND,NT)
!
! Update RHS for off block diaginal elements and then reiterate.
!
	  DO K=DST,DEND
	    STEQ(:,K)=STEQ_STORE(:,K)
	    IF(K .EQ. 1)THEN
	      DO J=1,N
	        DO I=1,N
	          STEQ(I,K)=STEQ(I,K)-D_MAT(I,J,K)*OLD_EST(J,K+1)
	        END DO
	      END DO  
	    ELSE IF(K .EQ. ND)THEN
	      DO J=1,N
	        DO I=1,N
	          STEQ(I,K)=STEQ(I,K)-B_MAT(I,J,K)*OLD_EST(J,K-1)
	        END DO
	      END DO  
	    ELSE
	      DO J=1,N
	        DO I=1,N
	          STEQ(I,K)=STEQ(I,K)-B_MAT(I,J,K)*OLD_EST(J,K-1)-
	1                                   D_MAT(I,J,K)*OLD_EST(J,K+1)
	        END DO
	      END DO  
	    END IF
	  END DO
	  IT_COUNTER=IT_COUNTER+1
	END DO
	CALL TUNE(2,'IT_COUNTER')
!
	IF(.NOT. BA_COMPUTED .AND. WR_BA_INV)THEN
	  DO K=DST,DEND
	    STEQ(:,K)=STEQ(:,K)*(ORIG_POPS(:,K)/POPS(:,K))
	  END DO
	END IF	
	FLAG=.TRUE.
!
	DEALLOCATE (B_MAT,C_MAT,D_MAT,RUB,NEW_EST)
	DEALLOCATE(IPIVOT,COL_SF,ROW_SF)
	IF(ALLOCATED(ORIG_POPS))DEALLOCATE (ORIG_POPS)
	IF(MYPE .EQ. 0)CLOSE(LU_IT)
!
	RETURN
!
! Error handling setcion
!
9999	CONTINUE
	WRITE(LUER,*)'Error in LINPAC (or BLAS) routine ',DESC,' in CMF_TRI_BAMD_MPI_V2'
	WRITE(LUER,100)K,IFAIL
100	FORMAT(1x,'depth=',I5,10x,'IFAIL=',I5)
	FLAG=.FALSE.
	RETURN
!
	END
