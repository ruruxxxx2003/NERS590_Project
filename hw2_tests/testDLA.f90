PROGRAM testDLA
  USE ISO_FORTRAN_ENV
  USE ISO_C_BINDING
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE solveGaussElim(n,A,b,x,ierr) BIND(C,NAME="solveGaussElim")
      USE ISO_C_BINDING
      INTEGER(C_INT),INTENT(IN) :: n
      REAL(C_DOUBLE),INTENT(IN) :: A(n,*)
      REAL(C_DOUBLE),INTENT(IN) :: b(*)
      REAL(C_DOUBLE),INTENT(INOUT) :: x(*)
      INTEGER(C_INT),INTENT(INOUT) :: ierr
    ENDSUBROUTINE
  ENDINTERFACE

  CALL testIdentity()
  CALL testSmall()
  CALL testLaplacian()
  CALL testBadArgs()
  CALL testBadMatrix()

  WRITE(OUTPUT_UNIT,*) "All unit tests passed"
CONTAINS
!
!-------------------------------------------------------------------------------
  SUBROUTINE testIdentity()
    INTEGER :: i
    REAL(C_DOUBLE) :: a(10,10),b(10),x(10)
    INTEGER(C_INT) :: ierr

    a=0.0d0
    x=0.0d0
    DO i=1,10
      a(i,i)=1.0d0
      b(i)=REAL(i,C_DOUBLE)
    ENDDO 

    CALL solveGaussElim(10_C_INT,a,b,x,ierr)

    IF(ANY(ABS(b-x) > EPSILON(x)) .OR. ierr /= 0) CALL fail('Identity')
  ENDSUBROUTINE testIdentity
!
!-------------------------------------------------------------------------------
  SUBROUTINE testSmall()
    REAL(C_DOUBLE) :: a(4,4),b(4),x(4),refsol(4)
    INTEGER(C_INT) :: ierr

    a(1,1)=1.0d0; a(1,2)=2.0d0; a(1,3)=3.0d0; a(1,4)=4.0d0
    a(2,1)=1.0d0; a(2,2)=3.0d0; a(2,3)=2.0d0; a(2,4)=3.0d0
    a(3,1)=3.0d0; a(3,2)=2.0d0; a(3,3)=3.0d0; a(3,4)=1.0d0
    a(4,1)=1.0d0; a(4,2)=1.0d0; a(4,3)=1.0d0; a(4,4)=1.0d0

    b=(/10.0d0,9.0d0,9.0d0,4.0d0/)

    x=0.0d0
    refsol=1.0d0

    CALL solveGaussElim(4_C_INT,a,b,x,ierr)
    IF(ANY(ABS(refsol-x) > EPSILON(x)) .OR. ierr /=0) CALL fail('Small')
  ENDSUBROUTINE testSmall
!
!-------------------------------------------------------------------------------
  SUBROUTINE testLaplacian()
    INTEGER,PARAMETER :: gridn=4,n=gridn*gridn
    INTEGER :: k
    REAL(C_DOUBLE) :: a(n,n),b(n),x(n),refsol(n)
    INTEGER(C_INT) :: ierr

    a=0.0d0
    DO k=1,n
      a(k,k)=4.0d0
      IF(k > 1) THEN
        a(k-1,k)=-1.0d0
      ENDIF
      IF(k < n) THEN
        a(k+1,k)=-1.0d0
      ENDIF
      IF(k > gridn) THEN
        a(k-gridn,k)=-1.0d0
      ENDIF
      IF(k < n-gridn) THEN
        a(k+gridn,k)=-1.0d0
      ENDIF
      refsol(k)=REAL(k,C_DOUBLE)
    ENDDO

    b=MATMUL(a,refsol)
    x=0.0d0
    
    CALL solveGaussElim(n,a,b,x,ierr)
    IF(ierr /= 0 .OR. ANY(ABS(x-refsol) > 1d-14)) CALL fail('Laplacian')
  ENDSUBROUTINE testLaplacian
!
!-------------------------------------------------------------------------------
  SUBROUTINE testBadMatrix()
    REAL(C_DOUBLE) :: a(4,4),b(4),x(4),refsol(4)
    INTEGER(C_INT) :: ierr

    a=1.0d0
    b=(/10.0d0,9.0d0,9.0d0,4.0d0/)
    x=0.0d0

    CALL solveGaussElim(4_C_INT,a,b,x,ierr)
    IF(ierr == 0) CALL fail('Bad Matrix')
  ENDSUBROUTINE testBadMatrix
!
!-------------------------------------------------------------------------------
  SUBROUTINE testBadArgs()
    INTEGER(C_INT) :: ierr
    REAL(C_DOUBLE) :: a(1,1),b(1),x(1)
    CALL solveGaussElim(0_C_INT,a,b,x,ierr)
    IF(ierr /= -1) CALL fail('Bad input arguments')
  ENDSUBROUTINE testBadArgs
!
!-------------------------------------------------------------------------------
  SUBROUTINE fail(tname)
    CHARACTER(LEN=*),INTENT(IN) :: tname
    WRITE(ERROR_UNIT,*) 'Test "'//TRIM(tname)//'" failed!'
    STOP 1
  ENDSUBROUTINE fail
!
ENDPROGRAM testDLA
