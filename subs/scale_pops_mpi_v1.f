C
C Routine to scale the populations so that the species conservation
C equation is satisfied. Here, SUM, which  must be supplied, is the
C current species population obtained by summing over all sates.
C
C Created 11-Apr-1989
C
	SUBROUTINE SCALE_POPS_MPI_V1(HYD,DHYD,POPHYD,SUM,NHYD,DST,DEND,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
	INTEGER DST,DEND
	INTEGER NHYD,ND
	REAL(KIND=LDP) HYD(NHYD,DST:DEND),DHYD(DST:DEND)
	REAL(KIND=LDP) POPHYD(ND),SUM(ND)
C
	INTEGER I,J
	REAL(KIND=LDP) T1
C
	IF(POPHYD(ND) .EQ. 0)RETURN
C
	DO J=DST,DEND
	  T1=POPHYD(J)/SUM(J)
	  DO I=1,NHYD
	    HYD(I,J)=HYD(I,J)*T1
	  END DO
	  DHYD(J)=DHYD(J)*T1
	END DO
C
	RETURN
	END
