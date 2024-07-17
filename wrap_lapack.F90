SUBROUTINE diagonalize_matrix(N,A,e)

! Diagonalize a square matrix

  implicit none

! Input variables

  integer,intent(in)            :: N
  double precision,intent(inout):: A(N,N)
  double precision,intent(out)  :: e(N)

! Local variables

  integer                       :: lwork,info
  double precision,allocatable  :: work(:)

! Memory allocation

  allocate(work(3*N))
  lwork = size(work)

  call dsyev('V','U',N,A,N,e,work,lwork,info)
 
  if(info /= 0) then 
    print*,'Problem in diagonalize_matrix (dsyev)!!'
    stop
  endif

END SUBROUTINE diagonalize_matrix

SUBROUTINE svd(N,A,U,D,Vt)

  ! Compute A = U.D.Vt
  ! Dimension of A is NxN

  implicit none

  integer, intent(in)             :: N
  double precision,intent(in)     :: A(N,N)
  double precision,intent(out)    :: U(N,N)
  double precision,intent(out)    :: Vt(N,N)
  double precision,intent(out)    :: D(N)
  double precision,allocatable    :: work(:)
  integer                         :: info,lwork

  double precision,allocatable    :: scr(:,:)

  allocate (scr(N,N))

  scr(:,:) = A(:,:)

  ! Find optimal size for temporary arrays

  allocate(work(1))

  lwork = -1
  call dgesvd('A','A',N,N,scr,N,D,U,N,Vt,N,work,lwork,info)
  lwork = int(work(1))

  deallocate(work)

  allocate(work(lwork))

  call dgesvd('A','A',N,N,scr,N,D,U,N,Vt,N,work,lwork,info)

  deallocate(work,scr)

  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif

END SUBROUTINE svd

SUBROUTINE inverse_matrix(N,A,B)

! Returns the inverse of the square matrix A in B

  implicit none

  integer,intent(in)             :: N
  double precision, intent(in)   :: A(N,N)
  double precision, intent(out)  :: B(N,N)

  integer                        :: info,lwork
  integer, allocatable           :: ipiv(:)
  double precision,allocatable   :: work(:)

  allocate (ipiv(N),work(N*N))
  lwork = size(work)

  B(1:N,1:N) = A(1:N,1:N)

  call dgetrf(N,N,B,N,ipiv,info)

  if (info /= 0) then

    print*,info
    stop 'error in inverse (dgetrf)!!'

  endif

  call dgetri(N,B,N,ipiv,work,lwork,info)

  if (info /= 0) then

    print *,  info
    stop 'error in inverse (dgetri)!!'

  endif

  deallocate(ipiv,work)

END SUBROUTINE inverse_matrix

SUBROUTINE linear_solve(N,A,b,x)

! Solve the linear system A.x = b where A is a NxN matrix
! and x and x are vectors of size N

  implicit none

  integer,intent(in)             :: N
  double precision,intent(in)    :: A(N,N),b(N)
  double precision,intent(out)   :: x(N)

  integer                        :: info,lwork
  integer,allocatable            :: ipiv(:)
  double precision,allocatable   :: work(:)

  allocate(ipiv(N),work(N*N))
  lwork = size(work)

  x = b

  call dsysv('U',N,1,A,N,ipiv,x,N,work,lwork,info)

  if (info /= 0) then

    print *,  info
    stop 'error in linear_solve (dsysv)!!'

  endif

END SUBROUTINE linear_solve

SUBROUTINE eval_determinant (ldm, n, matrix, lu_matrix, ipiv, det)
    !-------------------------------------------------------------!
    ! Given MATRIX(1:N,1:N) return its determinant DET, its LU    !
    ! decomposition LU_MATRIX(1:N,1:N) and its pivot vector IPIV. !
    !-------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER(4)      , INTENT(in)    :: ldm, n
    DOUBLE PRECISION, INTENT(in)    :: matrix(ldm,n)
    DOUBLE PRECISION, INTENT(inout) :: lu_matrix(ldm,n), det
    INTEGER(4)      , INTENT(inout) :: ipiv(n)
    INTEGER(4)                      :: ierr, i
    lu_matrix = matrix
    call dgetrf (n, n, lu_matrix, ldm, ipiv, ierr)
    det = product( (/ (lu_matrix(i,i), i=1,n) /) )
    if ( mod( count( ipiv /= (/(i, i=1,n) /) ), 2) /= 0 ) det = -det
END SUBROUTINE eval_determinant

SUBROUTINE eval_complex_determinant (ldm, n, matrix, lu_matrix, ipiv, det)
    !----------------------------------------------------------------!
    ! Given MATRIX(1:N,1:N) COMPLEX return its determinant DET, its  !
    ! LU decomposition LU_MATRIX(1:N,1:N) and its pivot vector IPIV. !
    !----------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER(4)    , INTENT(in)    :: ldm, n
    DOUBLE COMPLEX, INTENT(in)    :: matrix(ldm,n)
    DOUBLE COMPLEX, INTENT(inout) :: lu_matrix(ldm,n), det
    INTEGER(4)    , INTENT(inout) :: ipiv(n)
    INTEGER(4)                    :: ierr, i,j
    lu_matrix = matrix
    call zgetrf (n, n, lu_matrix, ldm, ipiv, ierr)
    det = product( (/ (lu_matrix(i,i), i=1,n) /) )
    if ( mod( count( ipiv /= (/(i, i=1,n) /) ), 2) /= 0 ) det = -det
    !write(*,*) 'Number of PIVOTS   ', count( ipiv /= (/(i, i=1,n) /) )
    !write(*,*) 'Determinant        ', det
    !write(*,*) 'IERR               ', ierr
    !do i = 1,n,1
    !   do j = 1,n,1
    !      write(*,*) 'matrix[',i,',',j,'] = ', lu_matrix(i,j)
    !   enddo
    !enddo
END SUBROUTINE eval_complex_determinant


