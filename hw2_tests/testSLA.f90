PROGRAM testSLA
  USE ISO_FORTRAN_ENV
  USE ISO_C_BINDING
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE solveGMRES(max_iters,tol,n,ia,ja,aa,b,x,iters,r2,ierr) &
      BIND(C,NAME="solveGMRES")
      USE ISO_C_BINDING
      INTEGER(C_INT),INTENT(IN),VALUE :: max_iters
      REAL(C_DOUBLE),INTENT(IN),VALUE :: tol
      INTEGER(C_INT),INTENT(IN),VALUE :: n
      INTEGER(C_INT),INTENT(IN) :: ia(*)
      INTEGER(C_INT),INTENT(IN) :: ja(*)
      REAL(C_DOUBLE),INTENT(IN) :: aa(*)
      REAL(C_DOUBLE),INTENT(IN) :: b(*)
      REAL(C_DOUBLE),INTENT(INOUT) :: x(*)
      INTEGER(C_INT),INTENT(INOUT) :: iters
      REAL(C_DOUBLE),INTENT(INOUT) :: r2
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
    REAL(C_DOUBLE) :: a(10),b(10),x(10),r2
    INTEGER(C_INT) :: ia(11),ja(11),iters,ierr

    ia(1)=1
    DO i=2,11
      ia(i)=ia(i-1)+1
    ENDDO
    DO i=1,10
     ja(i)=i
     a(i)=1.0d0
     b(i)=REAL(i,C_DOUBLE)
    ENDDO

    x=0.0d0
    CALL solveGMRES(10_C_INT,EPSILON(x),10_C_INT,ia,ja,a,b,x,iters,r2,ierr)
    IF(ANY(ABS(b-x) > EPSILON(x)) .OR. ierr /= 0) CALL fail('Identity')
  ENDSUBROUTINE testIdentity
!
!-------------------------------------------------------------------------------
  SUBROUTINE testSmall()
    REAL(C_DOUBLE) :: a(4,4),b(4),x(4),refsol(4),r2
    REAL(C_DOUBLE),ALLOCATABLE :: aa(:)
    INTEGER(C_INT),ALLOCATABLE :: ia(:),ja(:)
    INTEGER(C_INT) :: ierr,iters
    

    a(1,1)=1.0d0; a(1,2)=2.0d0; a(1,3)=3.0d0; a(1,4)=4.0d0
    a(2,1)=1.0d0; a(2,2)=3.0d0; a(2,3)=2.0d0; a(2,4)=3.0d0
    a(3,1)=3.0d0; a(3,2)=2.0d0; a(3,3)=3.0d0; a(3,4)=1.0d0
    a(4,1)=1.0d0; a(4,2)=1.0d0; a(4,3)=1.0d0; a(4,4)=1.0d0

    CALL full2csr(a,ia,ja,aa)

    b=(/10.0d0,9.0d0,9.0d0,4.0d0/)

    x=0.0d0
    refsol=1.0d0

    CALL solveGMRES(10_C_INT,EPSILON(x),4_C_INT,ia,ja,aa,b,x,iters,r2,ierr)
    IF(ANY(ABS(refsol-x) > r2) .OR. ierr < 0) CALL fail('Small')
  ENDSUBROUTINE testSmall
!
!-------------------------------------------------------------------------------
  SUBROUTINE testLaplacian()
    INTEGER(C_INT),PARAMETER :: gridn=4,n=gridn*gridn
    INTEGER :: k
    REAL(C_DOUBLE) :: a(n,n),b(n),x(n),refsol(n),r2
    INTEGER(C_INT) :: ierr,iters

    REAL(C_DOUBLE),ALLOCATABLE :: aa(:)
    INTEGER(C_INT),ALLOCATABLE :: ia(:),ja(:)

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

    CALL full2csr(a,ia,ja,aa)

    CALL solveGMRES(n,EPSILON(x),n,ia,ja,aa,b,x,iters,r2,ierr)
    IF(ierr < 0 .OR. ANY(ABS(x-refsol) > 1.0d1*r2)) CALL fail('Laplacian')
  ENDSUBROUTINE testLaplacian
!
!-------------------------------------------------------------------------------
  SUBROUTINE testBadMatrix()
    REAL(C_DOUBLE) :: a(4,4),b(4),x(4),refsol(4),r2
    INTEGER(C_INT) :: ierr,iters

    REAL(C_DOUBLE),ALLOCATABLE :: aa(:)
    INTEGER(C_INT),ALLOCATABLE :: ia(:),ja(:)

    a=1.0d0
    b=(/10.0d0,9.0d0,9.0d0,4.0d0/)
    x=0.0d0

    CALL full2csr(a,ia,ja,aa)

    CALL solveGMRES(4_C_INT,EPSILON(x),4_C_INT,ia,ja,aa,b,x,iters,r2,ierr)
    IF(ierr == 0) CALL fail('Bad Matrix')
  ENDSUBROUTINE testBadMatrix
!
!-------------------------------------------------------------------------------
  SUBROUTINE testBadArgs()
    INTEGER(C_INT) :: ia(1),ja(1),iters,ierr
    REAL(C_DOUBLE) :: a(1),b(1),x(1),r2

    CALL solveGMRES(0_C_INT,EPSILON(x),1_C_INT,ia,ja,a,b,x,iters,r2,ierr)
    IF(ierr /= -1) CALL fail('Bad input arguments (max_iters)')
    CALL solveGMRES(10_C_INT,-EPSILON(x),1_C_INT,ia,ja,a,b,x,iters,r2,ierr)
    IF(ierr /= -1) CALL fail('Bad input arguments (tol)')
    CALL solveGMRES(10_C_INT,EPSILON(x),0_C_INT,ia,ja,a,b,x,iters,r2,ierr)
    IF(ierr /= -1) CALL fail('Bad input arguments (n)')
    x=0.0d0
    CALL solveGMRES(10_C_INT,-EPSILON(x),1_C_INT,ia,ja,a,b,x,iters,r2,ierr)
    IF(ierr /= -1) CALL fail('Bad input arguments (x)')
  ENDSUBROUTINE testBadArgs
!
!-------------------------------------------------------------------------------
  SUBROUTINE fail(tname)
    CHARACTER(LEN=*),INTENT(IN) :: tname
    WRITE(ERROR_UNIT,*) 'Test "'//TRIM(tname)//'" failed!'
    STOP 1
  ENDSUBROUTINE fail
!
!-------------------------------------------------------------------------------
  SUBROUTINE full2csr(a,ia,ja,aa)
    REAL(C_DOUBLE),INTENT(IN) :: a(:,:)
    INTEGER(C_INT),ALLOCATABLE,INTENT(INOUT) :: ia(:),ja(:)
    REAL(C_DOUBLE),ALLOCATABLE,INTENT(INOUT) :: aa(:)

    INTEGER(C_INT) :: n,k,i,j,nnz

    n=SIZE(a,DIM=1)

    ALLOCATE(ia(n+1)); ia(1)=1
    k=0
    DO i=1,n
      DO j=1,n
        IF(a(i,j) /= 0.0d0) k=k+1
      ENDDO
      ia(i+1)=k+1
    ENDDO

    nnz=ia(n+1)-1
    ALLOCATE(aa(nnz))
    ALLOCATE(ja(nnz))

    k=0
    DO i=1,n
      DO j=1,n
        IF(a(i,j) /= 0.0d0) THEN
          k=k+1
          aa(k)=a(i,j)
          ja(k)=j
        ENDIF
      ENDDO
    ENDDO

!    DO i=1,n
!      DO j=1,n
!      IF(a(i,j) /= 0.0d0) THEN
!        WRITE(0,*) i,j,a(i,j)
!      ENDIF
!      ENDDO
!    ENDDO
!WRITE(0,*)
!    DO i=1,n
!      DO k=ia(i),ia(i+1)-1
!        WRITE(0,*) i,ja(k),aa(k)
!      ENDDO
!    ENDDO
  ENDSUBROUTINE
!
ENDPROGRAM testSLA
