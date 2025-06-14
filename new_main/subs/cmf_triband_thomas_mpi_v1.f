!
! Subroutine to solve a "Block-Banded" system of simultaneous equations
!                 BA. X. = STEQ
! The solution X is returned in STEQ. BA is corrupted.
!
! The "BANDED" matrix is assumed to be 'Tridiagonal'
!
	SUBROUTINE CMF_TRIBAND_THOMAS_MPI_V1(STEQ,POPS,SOL_TYPE,FLAG,
	1                 DIAG_INDX,N,NION,NUM_BNDS,DST,DEND,ND,
	1                 BA_COMPUTED,WR_BA_INV,WR_PRT_INV,TRI_SOL_OPTIONS)
	USE SET_KIND_MODULE
	USE MPI
	IMPLICIT NONE
!
! Altered 1-0Jun-2025 : Comments (partially fixed). Earlier a MPI_BARRIER statement was pu in its
!                          correct location.
!
! 
! The description here is from CMF_TRIBAND_THOMAS_MPI_V1, which also allowed for a PENTDIAGONAL
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
!
!*******************
!                  MAT5PEN(A,B,C,D,E,F,WRKMAT,VEC,N,NS,ENONZERO,DNONZERO)
!
! Evaluates A = A - B. C + D. (F . C - E)
!
! where:
!               A = real A[N,NS] - Retuned with solution
!               B = real B[N,N]  - Unchanged
!               C = real C[N,NS] - Unchanged
!               D = real D[N,N]  - Unchanged
!               E = real E[N,NS] - Unchanged
!               F = real F[N,N] -  Unchanged
!               WRKMAT = real A[N,NS] - Working matrix (corrupted)
!               VEC = real VEC[N] - Working vector (passed but no longer used)
!               ENONZERO - If FALSE, E is assumed to be zero, and E is not
!                             accessed.
!               DNONZERO - If FALSE, D is assumed to be zero, and D, F and E
!                             are not accessed.
!
!
!******************************************************************************
! 
!
	INTEGER N,NION,DST,DEND,ND,NUM_BNDS,DIAG_INDX
	CHARACTER*(*) SOL_TYPE
	CHARACTER*(*) TRI_SOL_OPTIONS
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
        REAL(KIND=LDP), ALLOCATABLE :: B_MAT(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: C_MAT(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: D_MAT(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: ORIG_POPS(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: PREV_D_MAT(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: D_MAT_STORE(:,:,:)
        REAL(KIND=LDP), ALLOCATABLE :: RUB(:,:)	      	      !Not accessed when passed.
        REAL(KIND=LDP), ALLOCATABLE :: STEQ_WRK_VEC(:)
!
        REAL(KIND=LDP), ALLOCATABLE :: ROW_SF(:)
        REAL(KIND=LDP), ALLOCATABLE :: COL_SF(:)
        REAL(KIND=LDP) ROW_CND,COL_CND,MAX_VAL
	INTEGER, ALLOCATABLE ::  IPIVOT(:)
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
        INTEGER IOS,IFAIL,IERR,ITAG
	CHARACTER*10 DESC
	LOGICAL ANONZERO,ENONZERO
	LOGICAL FIRST_MATRIX,LAST_MATRIX
	LOGICAL WR_D_MAT
	INTEGER REC_STATUS(MPI_STATUS_SIZE)
!
        REAL(KIND=LDP),      PARAMETER :: DP_NEG_ONE=-1.0_LDP
        REAL(KIND=LDP),      PARAMETER :: DP_ONE=1.0_LDP
        INTEGER,   PARAMETER :: INT_ONE=1
        INTEGER,   PARAMETER :: NSNG=1
        CHARACTER*1, PARAMETER :: NO_TRANS='N'
	CHARACTER(LEN=50) TMP_STR
!
	INTEGER ITHREAD
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	LUER=ERROR_LU()
!
!
!
! If we get here, we must be performing a TRI diagonal solution.
!
	IF(SOL_TYPE(1:3) .EQ. 'TRI')THEN
	  IF(NUM_BNDS .LT. 3)THEN
	    WRITE(LUER,*)'Error in CMF_TRIBAND_THOMAS_MPI_V1 : NUM_BNDS too small ',
	1             'for solution type'
	    FLAG=.FALSE.
	    STOP
	  END IF
	ELSE
	  WRITE(LUER,*)'Error in CMF_TRIBAND_THOMAS_MPI_V1 - invalid SOL_TYPE - ',
	1          'SOLTYPE= ',SOL_TYPE
	  STOP
        END IF
!
! Solve the TRIDIAGONAL system of equations assuming that we have already performed
! an LU decomposition on the BA matrix.
!
	IF(SOL_TYPE(1:3) .EQ. 'TRI' .AND. .NOT. BA_COMPUTED .AND. WR_BA_INV)THEN
!
! Allocate needed work arrays. These have been saved, 1 depth at a time.
!
          ALLOCATE (B_MAT(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (C_MAT(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (RUB(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (ROW_SF(N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (COL_SF(N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (IPIVOT(N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (ORIG_POPS(N,DST:DEND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (STEQ_WRK_VEC(N),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in CMF_TRIBAND_THOMAS_MPI_V1'
            WRITE(LUER,*)'Unable to allocate B_MAT, C_MAT & D_MAT'
            WRITE(LUER,*)'STAT=',IOS
            STOP
          END IF
	  B_MAT=0.0_LDP; C_MAT=0.0_LDP; RUB=0.0_LDP; ORIG_POPS=0.0_LDP
!
! Needed for MAT5PEN and GENERATE_FULL_MATRIX.
!
          ANONZERO=.FALSE.
          ENONZERO=.FALSE.
	  FIRST_MATRIX=.TRUE.
	  LAST_MATRIX=.FALSE.
	  USE_PASSED_REP=.TRUE.
! 
!
! Do the forward elimination etc.
!
	  DO K=DST,DEND
	    DEPTH_INDX=K
!
	    CALL READ_BCD_MAT(B_MAT,C_MAT,RUB,ROW_SF,COL_SF,IPIVOT,
	1             ORIG_POPS(1,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'BC')
!
! Map the small BA array onto the full BA array (one depth at a time).
! We just need STEQ as a modified C_MAT will be read in. To avoid
! corrupting C_MAT, we use RUB
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
	  END DO
!
! Computes phi[k]' (i.e. ' denotes not yet multiplied by m[k]^{-1} )
! Stored in L[k].
!
	  DO ITHREAD=0,NTHREAD-1
	    IF(MYPE .EQ. ITHREAD)THEN
!
! In this section, STEQ_WRK_VEC is used to STEQ at depth K-1
!
	      DO K=DST,DEND
	        IF(K .EQ. DST .AND. K .NE. 1)THEN
	          ITAG=MYPE
	          CALL MPI_RECV(STEQ_WRK_VEC,N,MY_MPI_DP,ITHREAD-1,ITAG,MPI_COMM_WORLD,REC_STATUS,IERR)
	        ELSE IF(K .NE. 1)THEN
	          STEQ_WRK_VEC=STEQ(:,K-1)
	        END IF
!
	        IF(K .NE. 1)THEN
	          ENONZERO=.FALSE.
	          ANONZERO=.FALSE.
	          CALL MAT5PEN(STEQ(1,K),B_MAT,
	1                 STEQ_WRK_VEC,RUB,RUB,RUB,RUB,RUB,
	1                 N,NSNG,ENONZERO,ANONZERO)
	        END IF
!
! Computes phi[k] (Stored in L[k])
!
	        DO J=1,N
	          STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
	        END DO
	        CALL DGETRS(NO_TRANS,N,NSNG,C_MAT,N,IPIVOT,STEQ,N,IFAIL)
	        IF(IFAIL .NE. 0)THEN
	          DESC='DGETRS_4'
	          GOTO 9999
	        END IF
	        DO J=1,N
	          STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	        END DO
	        IF(K .EQ. DEND .AND. K .LT. ND)THEN
	          ITAG=MYPE+1
	          CALL MPI_SEND(STEQ(1,DEND),N,MY_MPI_DP,ITHREAD+1,ITAG,MPI_COMM_WORLD,IERR)
	        END IF
	      END DO 		!Loop over depth
	    END IF		!Correct processors
	  END DO		!Thread
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
!**************************************************************************
!
! Now we start the back substitution. In this section, STEQ_WRK_VEC is used to store STEQ(:,K+1).
!
	DO ITHREAD=NTHREAD-1,0,-1
	  IF(MYPE .EQ. ITHREAD)THEN
	    DO K=MIN(ND-1,DEND),DST,-1
	      IF(K .EQ. DEND)THEN
	        ITAG=MYPE
	        CALL MPI_RECV(STEQ_WRK_VEC,N,MY_MPI_DP,ITHREAD+1,ITAG,MPI_COMM_WORLD,REC_STATUS,IERR)
	      ELSE
	        STEQ_WRK_VEC=STEQ(:,K+1)
	      END IF
!
! Since we no longer need the C_MAT, we use it to store the D matrix.
!
	      DEPTH_INDX=K
	      CALL READ_BCD_MAT(RUB,RUB,C_MAT,RUB,RUB,RUB_VEC,
	1            ORIG_POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'D')
	      CALL DGEMV(NO_TRANS,N,N,DP_NEG_ONE,C_MAT,N,STEQ_WRK_VEC,
	1                 INT_ONE,DP_ONE,STEQ(1,K),INT_ONE)
	      IF(IFAIL .NE. 0)THEN
	        DESC='DGEMV_2'
	        GOTO 9999
	      END IF
	      IF(K .EQ. DST .AND. DST .NE. 1)THEN
	        ITAG=MYPE-1
	        CALL MPI_SEND(STEQ(1,DST),N,MY_MPI_DP,ITHREAD-1,ITAG,MPI_COMM_WORLD,IERR)
	      END IF
	    END DO
	  END IF
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	END DO
!
! Correct for the use of old POPS populations in doing the scaling.
!
	DO K=DST,DEND
	  DO I=1,N
	    STEQ(I,K)=STEQ(I,K)*(ORIG_POPS(I,K)/POPS(I,K))
	  END DO
	END DO
!
! Successfull solution obtained.
!
	FLAG=.TRUE.
!
	DEALLOCATE (B_MAT,C_MAT,ROW_SF,COL_SF)
	DEALLOCATE (ORIG_POPS,RUB,IPIVOT)
	CALL WR2D_GATH_MPI_V1(STEQ_STORE,N,DST,DEND,ND,'STEQ_ARRAY','*',L_TRUE,16)
!
	RETURN
	END IF
!
! 
!
! Perform the TRIDIAGONAL solution. If we reach here, we need to redo the LU decomposition of BA.
!
	IF(MYPE .EQ. 0)WRITE(6,*)'Beginning TRI solution in CMF_TRIBAND_THOMAS_MPI_V1'
!
	ALLOCATE (B_MAT(N,N),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (C_MAT(N,N),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (D_MAT(N,N),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (PREV_D_MAT(N,N+1),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (RUB(N,N),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (ROW_SF(N),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (COL_SF(N),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (IPIVOT(N),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (STEQ_WRK_VEC(N),STAT=IOS)
        IF(IOS .NE. 0)THEN
          WRITE(LUER,*)'Error in CMF_TRIBAND_THOMAS_MPI_V1'
          WRITE(LUER,*)'Unable to allocate D_MAT etc'
          WRITE(LUER,*)'STAT=',IOS
          STOP
        END IF
	B_MAT=0.0_LDP; C_MAT=0.0_LDP;       D_MAT=0.0_LDP
	RUB=0.0_LDP;   PREV_D_MAT=0.0_LDP
!
! WR_BA_INV is false, we perform the LU decomosition using dynamic memory
! allocation only. If insufficent storage can be allocated to store D_MAT
! (in D_MAT_STORE), we output it to disk, one depth at a time, for later
! use. In this case, D_MAT will be output, independent of WR_BA_INV.
!
	WR_D_MAT=WR_PRT_INV
	IF(.NOT. WR_BA_INV .AND. .NOT. WR_PRT_INV)THEN
          ALLOCATE (D_MAT_STORE(N,N,DST:DEND),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in CMF_TRIBAND_THOMAS_MPI_V1'
            WRITE(LUER,*)'Unable to allocate D_MAT_STORE etc'
            WRITE(LUER,*)'STAT=',IOS
            WR_D_MAT=.TRUE.
          END IF
	END IF
!
! 
!
! Solve the tridiagonal equations.
!
        ANONZERO=.FALSE.
        ENONZERO=.FALSE.
	FIRST_MATRIX=.FALSE.
	IF(MYPE .EQ. 0)FIRST_MATRIX=.TRUE.
	LAST_MATRIX=.FALSE.
	USE_PASSED_REP=.FALSE.
	I=MAX(1,NTHREAD-2)
	CALL OMP_SET_NUM_THREADS(I)
!
	DO ITHREAD=0,NTHREAD-1
	  IF(MYPE .EQ. ITHREAD)THEN
!
	    DO K=DST,DEND
!
! Map the small BA rray onto the full BA array (one depth at a time).
!
	      DEPTH_INDX=K
	      CALL TUNE(1,'TRI_GEN')
	      CALL GENERATE_FULL_MATRIX_V3(
	1            C_MAT,STEQ(1,K),POPS,REPLACE_EQ,ZERO_STEQ,
	1            N,ND,NION,NUM_BNDS,
	1            DIAG_INDX,DIAG_INDX,DEPTH_INDX,
	1            FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	      FIRST_MATRIX=.FALSE.
	      IF(K .NE. 1)THEN
	        IF(K .EQ. ND)LAST_MATRIX=.TRUE.
	        BAND_INDX=DIAG_INDX-1
	        CALL GENERATE_FULL_MATRIX_V3(
	1             B_MAT,STEQ(1,K),POPS,REPLACE_EQ,ZERO_STEQ,
	1             N,ND,NION,NUM_BNDS,
	1             BAND_INDX,DIAG_INDX,DEPTH_INDX,
	1           FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	      END IF
	      IF(K .NE. ND)THEN
	        BAND_INDX=DIAG_INDX+1
	        CALL GENERATE_FULL_MATRIX_V3(
	1             D_MAT,STEQ(1,K),POPS,REPLACE_EQ,ZERO_STEQ,
	1             N,ND,NION,NUM_BNDS,
	1             BAND_INDX,DIAG_INDX,DEPTH_INDX,
	1             FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	      END IF
	      CALL TUNE(2,'TRI_GEN')
	      STEQ_STORE(:,K)=STEQ(:,K)
!
! Computes phi[k]' (i.e. ' denotes not yet multiplied by m[k]^{-1} )
! Stored in L[k].
!
	      CALL TUNE(1,'TRI_MAT5')
	      IF(K .NE. 1)THEN
!
! STEQ_WRK_VEC used to store STEQ(:,K-1)
!
	        IF(K .EQ. DST)THEN
	          I=N*(N+1); ITAG=MYPE
	          CALL MPI_RECV(PREV_D_MAT,I,MY_MPI_DP,MYPE-1,ITAG,MPI_COMM_WORLD,REC_STATUS,IERR)
	          STEQ_WRK_VEC(:)=PREV_D_MAT(:,N+1)
	        ELSE
	          STEQ_WRK_VEC(:)=STEQ(:,K-1)
	        END IF
!
	        CALL MAT5PEN(STEQ(1,K),B_MAT,
	1                  STEQ_WRK_VEC,RUB,RUB,RUB,RUB,RUB,
	1                  N,NSNG,ENONZERO,ANONZERO)
!
! Computes m[k] - Stored in C[k]
!
                CALL MAT5PEN(C_MAT,B_MAT,PREV_D_MAT,
	1                   RUB,RUB,RUB,RUB,RUB,
	1                   N,N,ENONZERO,ANONZERO)
	      END IF
	      CALL TUNE(2,'TRI_MAT5')
!
! Do LU decompostion of m[k]. We first equilibrilze C_MAT.
!
	      CALL TUNE(1,'TRI_DGEEQU')
	      CALL DGEEQU(N,N,C_MAT,N,ROW_SF,COL_SF,
	1                   ROW_CND,COL_CND,MAX_VAL,IFAIL)
	      CALL TUNE(2,'TRI_DGEEQU')
!
	      CALL TUNE(1,'TRI_DGETRF')
!$OMP PARALLEL DO PRIVATE(J,I)
	      DO J=1,N
                DO I=1,N
	          C_MAT(I,J)=C_MAT(I,J)*ROW_SF(I)*COL_SF(J)
	        END DO
	      END DO
	      CALL DGETRF(N,N,C_MAT,N,IPIVOT,IFAIL)
	      CALL TUNE(2,'TRI_DGETRF')
	      IF(IFAIL .NE. 0)THEN
	        DESC='DGETRF_2'
	        GOTO 9999
	      END IF
!
! Computes phi[k] (Stored in L[k])
!
	      CALL TUNE(1,'TRI_DGETRS')
	      DO J=1,N
	        STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
	      END DO
	      CALL DGETRS(NO_TRANS,N,NSNG,C_MAT,N,IPIVOT,STEQ(1,K),N,IFAIL)
	      CALL TUNE(2,'TRI_DGETRS')
	      IF(IFAIL .NE. 0)THEN
	        DESC='DGETRS_4'
	        GOTO 9999
	      END IF
	      DO J=1,N
	        STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	      END DO
!
! Computes alpha[k] (Stored in D[k)).
!
	      IF(K .NE. ND)THEN
!
!$OMP PARALLEL DO PRIVATE(J,I)
	        DO J=1,N
                  DO I=1,N
	             D_MAT(I,J)=D_MAT(I,J)*ROW_SF(I)
	          END DO
	        END DO
	        CALL TUNE(1,'TRI_DGETRS')
	        CALL DGETRS(NO_TRANS,N,N,C_MAT,N,IPIVOT,D_MAT,N,IFAIL)
	        CALL TUNE(2,'TRI_DGETRS')
	        IF(IFAIL .NE. 0)THEN
	          DESC='DGETRS_5'
	          GOTO 9999
	        END IF
!
!$OMP PARALLEL DO PRIVATE(J,I)
	        DO J=1,N
                  DO I=1,N
	            D_MAT(I,J)=D_MAT(I,J)*COL_SF(I)
	          END DO
	        END DO
!
	        PREV_D_MAT(:,1:N)=D_MAT
	        IF(WR_BA_INV .OR. WR_D_MAT)THEN
	          CALL WRITE_BCD_MAT(RUB,RUB,D_MAT,RUB,RUB,RUB_VEC,
	1                POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'D')
	        ELSE
	          D_MAT_STORE(:,:,K)=D_MAT(:,:)
	        END IF
	      END IF
	      IF(K .EQ. DEND .AND. K .NE. ND)THEN
	        I=N*(N+1); ITAG=MYPE+1
	        PREV_D_MAT(:,N+1)=STEQ(:,K)
	        CALL MPI_SEND(PREV_D_MAT,I,MY_MPI_DP,MYPE+1,ITAG,MPI_COMM_WORLD,IERR)
	      END IF
!
	      IF(WR_BA_INV)THEN
	        CALL WRITE_BCD_MAT(B_MAT,C_MAT,RUB,ROW_SF,COL_SF,IPIVOT,
	1              POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'BC')
	      END IF
	    END DO
	  END IF
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	END DO
!
! 
!**************************************************************************
!
! Now we start the back substitution.
!
	DO ITHREAD=NTHREAD-1,0,-1
	  IF(MYPE .EQ. ITHREAD)THEN
	    DO K=MIN(ND-1,DEND),DST,-1
!
! Here STEQ_WRK_VEC is used to store STEQ(:,K+1).
!
	      IF(K .EQ. DEND)THEN
	        ITAG=MYPE
	        CALL MPI_RECV(STEQ_WRK_VEC,N,MY_MPI_DP,ITHREAD+1,ITAG,MPI_COMM_WORLD,REC_STATUS,IERR)
	      ELSE
	        STEQ_WRK_VEC=STEQ(:,K+1)
	      END IF
	      DEPTH_INDX=K
	      IF(WR_BA_INV .OR. WR_D_MAT)THEN
	        CALL READ_BCD_MAT(RUB,RUB,D_MAT,RUB,RUB,RUB_VEC,
	1           POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'D')
	      ELSE
	        D_MAT=D_MAT_STORE(:,:,K)
	      END IF
	      CALL DGEMV(NO_TRANS,N,N,DP_NEG_ONE,D_MAT,N,STEQ_WRK_VEC,
	1                 INT_ONE,DP_ONE,STEQ(1,K),INT_ONE)
	      IF(IFAIL .NE. 0)THEN
	        DESC='DGEMV_2'
	        GOTO 9999
	      END IF
	      IF(K .EQ. DST .AND. DST .NE. 1)THEN
	        ITAG=MYPE-1
	        CALL MPI_SEND(STEQ(1,DST),N,MY_MPI_DP,ITHREAD-1,ITAG,MPI_COMM_WORLD,IERR)
	      END IF
	    END DO
	  END IF
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	END DO
	CALL OMP_SET_NUM_THREADS(1)
	IF(MYPE .EQ. 0)THEN
	  WRITE(6,*)' Obtained solutions of rate equations in CMF_TRIBAND_THOMAS_MPI_V1'
	  FLUSH(UNIT=6)
	END IF
!
! Successfull solution obtained.
!
	FLAG=.TRUE.
	DEALLOCATE (B_MAT,C_MAT,D_MAT,ROW_SF,COL_SF,IPIVOT,PREV_D_MAT,RUB)
	IF(ALLOCATED(D_MAT_STORE))DEALLOCATE (D_MAT_STORE)
!
        CALL TUNE(1,'TRI_GATH')
        CALL WR2D_GATH_MPI_V1(STEQ_STORE,N,DST,DEND,ND,'STEQ_ARRAY','*',L_TRUE,16)
        CALL TUNE(2,'TRI_GATH')
!
	IF(INDEX(TRI_SOL_OPTIONS,'CHECK_SOL') .NE. 0)THEN
	  IF(MYPE .EQ. 0)WRITE(6,*)TRIM(TRI_SOL_OPTIONS)
	  CALL CHECK_BA_SOL(STEQ,POPS,SOL_TYPE,
	1                 DIAG_INDX,N,NION,NUM_BNDS,DST,DEND,ND)
	END IF
!
	RETURN
!
! Error handling routine.
!
9999	CONTINUE
	WRITE(LUER,*)'Error in LINPAC (or BLAS) routine',DESC,' in CMF_TRIBAND_THOMAS_MPI_V1'
	WRITE(LUER,100)K,IFAIL
100	FORMAT(1x,'depth=',I3,10x,'IFAIL=',I3)
	FLAG=.FALSE.
	RETURN
!
	END
