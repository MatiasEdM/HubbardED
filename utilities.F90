!----------------------------------------------------------------------------------------------------------

INTEGER(8) FUNCTION fact(n)

!-------------------------------------------------------------------!        
! FUNCTION to calculate the FACTORIAL of a given integer number n!. !
!-------------------------------------------------------------------!

 implicit none

 integer(8), intent(in) :: n
 
 integer(8)             :: i

 fact = 1

 if (n.eq.0) then
    fact = 1
 elseif (n.eq.1) then
    fact = 1
 else
    do i = 2,n,1
       fact = fact * i
    enddo
 endif    

 return     

END FUNCTION fact

!----------------------------------------------------------------------------------------------------------

INTEGER(8) FUNCTION binom(n,k)

!----------------------------------------------------------------------------------------!        
! FUNCITON to calculate the BINOMIAL COEFFICIENT of a given pair (n,k) = n!/(k!x(n-k)!). !
!----------------------------------------------------------------------------------------!
 
 implicit none

 integer(8), intent(in) :: n, k

 integer(8)             :: fact

 if (k==n) then
    binom = 1
 elseif (k==1) then
    binom = n
 else
    binom = int( fact(n) / ( fact(k) * fact(n-k) )  )
 endif

 return

END FUNCTION binom

!--------------------------------------------------------------------------------------------------------- 

SUBROUTINE cross_product(c,v,w)

!----------------------------------------------------------------------------!
! SUBROUTINE to calculate the cross product between two vectors v, and w as: !
! c = v x w                                                                  !
! in Cartesian coordinates,                                                  !
! c(1) = v(2)*w(3) - v(3)*w(2)                                               !
! c(2) = v(3)*w(1) - v(1)*w(3)                                               !
! c(3) = v(1)*w(2) - v(2)*w(1)                                               !
!----------------------------------------------------------------------------!

  implicit none

  integer(8), parameter         :: dddm = 3

  double precision, intent(in)  :: v(dddm), w(dddm)
  double precision, intent(out) :: c(dddm)

  c(:) = 0.d0

  c(1) = v(2)*w(3) - v(3)*w(2)
  c(2) = v(3)*w(1) - v(1)*w(3)
  c(3) = v(1)*w(2) - v(2)*w(1)

  return

END SUBROUTINE cross_product

!----------------------------------------------------------------------------------------------------------

SUBROUTINE check_MATdiag(O,U,n)

!--------------------------------------------------------------------------------------!        
! SUBROUTINE to check the correct diagonalization of a matrix O(n,n) by performing:    !
! D = Ut * O * U                                                                       !
! where U is the unitary transformation matrix (the matrix containing the coefficients !
! that relate the vectors of the basis in which O is diagonal with the vectors of the  !
! original basis in a column-wise manner).                                             !
!--------------------------------------------------------------------------------------!

 implicit none 

 double precision, parameter  :: thrs = 1E-8

 integer(8), intent(in)       :: n
 double precision, intent(in) :: O(n,n), U(n,n)

 double precision             :: D(n,n)
 integer(8)                   :: i, j, r, s
 logical                      :: check

 D(:,:) = 0.d0

 do j = 1,n,1
    do i = 1,n,1
       do r = 1,n,1
          do s = 1,n,1
             D(i,j) = D(i,j) + U(r,i) * O(r,s) * U(s,j)
          enddo
       enddo
    enddo
 enddo

 check = .FALSE.

 do i = 1,n,1
    do j = 1,i-1,1
       if ( abs(D(i,j)) .gt. thrs ) then
         check = .TRUE.
         exit
       endif
    enddo
 enddo

 if ( check ) then
    write(*,fmt=200) ' - WARNING: Hamiltonian not correctly diagonalized !!!'
    write(*,*)
 endif

 200 format(A54)

 return

END SUBROUTINE check_MATdiag

!----------------------------------------------------------------------------------------------------------

SUBROUTINE normCI_VEC(CI,nDET)

!-------------------------------------------------------------------------------!
! SUBROUTINE that returns the CI Vector for the Ground State and Excited States !
! normalized according to the L2-NORM.                                          !
!-------------------------------------------------------------------------------!

 implicit none

 integer(8)      , intent(in)    :: nDET
 double precision, intent(inout) :: CI(nDET,nDET)

 integer(8)                      :: I
 double precision                :: norm, normCOEFF

 do I = 1,nDET,1
    norm    = normCOEFF( CI(:,I),nDET )
    CI(:,I) = norm * CI(:,I)
 enddo

 return

END SUBROUTINE normCI_VEC

!---------------------------------------------------------------------------------------------------------

DOUBLE PRECISION FUNCTION normCOEFF(v,n)

!-----------------------------------------------------------------------------------!        
! FUNCTION that calculates the normalization coefficient of a given vector v(n) as: !
! N = sqrt[ 1 / SUM_{i}[ |C(i)|^2] ]                                                !
!-----------------------------------------------------------------------------------!

 implicit none

 integer(8)      , intent(in) :: n
 double precision, intent(in) :: v(n)

 double precision             :: summ
 integer(8)                   :: i

 summ = 0.d0

 do i = 1,n,1
    summ = summ + v(i) * v(i)
 enddo

 normCOEFF = 1 / dsqrt( summ )

 return

END FUNCTION

!---------------------------------------------------------------------------------------------------------

DOUBLE PRECISION FUNCTION L2_norm(v,n)

!----------------------------------------------------------------!        
! FUNCTION that calculates de L2 norm of a given vector A(n) as: !
!  ||v|| = sqrt[ SUM_{i}[ |C(i)|^2 ] ]                           !
!----------------------------------------------------------------!

 implicit none

 integer(8)      , intent(in) :: n
 double precision, intent(in) :: v(n)

 integer(8)                   :: i

 L2_norm = 0.d0

  do i = 1,n,1
    L2_norm = L2_norm + v(i) * v(i)
 enddo

 L2_norm = dsqrt( L2_norm )

 return

END FUNCTION L2_norm

!---------------------------------------------------------------------------------------------------------

SUBROUTINE writeINTMAT(A,n,m)

!-----------------------------------------------!        
! SUBROUTINE to print out a MATRIX of INTEGERS. !
!-----------------------------------------------!

 implicit none

 integer(8), intent(in) :: n, m
 integer(8), intent(in) :: A(n,m)

 integer(8)             :: i,j

 do i = 1,n,1
    do j = 1,m-1,1
       write(*,fmt=200,advance='No') A(i,j)
    enddo
    write(*,fmt=200) A(i,m)
 enddo

 200 format(1I2)

 return

END SUBROUTINE writeINTMAT

!----------------------------------------------------------------------------------------------------------

SUBROUTINE writeREALMAT(A,n,m)

!---------------------------------------------------!        
! SUBROUTINE to print out a MATRIX of REAL NUMBERS. !
!---------------------------------------------------!

 implicit none

 integer(8)      , intent(in) :: n, m
 double precision, intent(in) :: A(n,m)

 integer(8)             :: i,j

 do i = 1,n,1
    do j = 1,m-1,1
       write(*,fmt=200,advance='No') A(i,j)
    enddo
    write(*,fmt=200) A(i,m)
 enddo

 200 format(1F6.2)

 return

END SUBROUTINE writeREALMAT

!----------------------------------------------------------------------------------------------------------

SUBROUTINE writeEREALMAT(A,n,m)

!---------------------------------------------------!        
! SUBROUTINE to print out a MATRIX of REAL NUMBERS. !
!---------------------------------------------------!

 implicit none

 double precision, parameter  :: thrs = 1e4

 integer(8)      , intent(in) :: n, m
 double precision, intent(in) :: A(n,m)

 integer(8)             :: i,j

 do i = 1,n,1
    do j = 1,m-1,1
       if ( abs(A(i,j)) .lt. thrs ) then
          write(*,fmt=200,advance='No') A(i,j)
        else
          write(*,fmt=210,advance='No') A(i,j)
        endif
    enddo
    if ( abs(A(i,m)) .lt. thrs ) then
       write(*,fmt=200) A(i,m)
     else
       write(*,fmt=210) A(i,m)
     endif
 enddo
 
 200 format(1F9.2)
 210 format(1E9.2)

 return

END SUBROUTINE writeEREALMAT

!----------------------------------------------------------------------------------------------------------

SUBROUTINE matout(m,n,A,u)

!---------------------------------------------------------------------!        
! Extra SUBROUTINE to print out matrices with row and column indices. !
! Print the MxN array A.                                              !
!---------------------------------------------------------------------!

  implicit none

  integer(8)      , parameter  :: ncol = 5
  double precision, parameter  :: small = 1d-10

  integer(4)      , intent(in) :: u

  integer(8)      , intent(in) :: m,n
  double precision, intent(in) :: A(m,n)
  double precision            :: B(ncol)
  integer(8)                  :: ilower,iupper,num,i,j

  do ilower=1,n,ncol
    iupper = min(ilower + ncol - 1,n)
    num = iupper - ilower + 1
    write(u,'(3X,10(9X,I6))') (j,j=ilower,iupper)
    do i=1,m
      do j=ilower,iupper
        B(j-ilower+1) = A(i,j)
      enddo
      do j=1,num
        if(abs(B(j)) < small) B(j) = 0d0
      enddo
      write(u,'(I7,10F15.8)') i,(B(j),j=1,num)
    enddo
  enddo

END SUBROUTINE matout

!----------------------------------------------------------------------------------------------------------

SUBROUTINE swap(x,y)

!--------------------------------------------!        
! SUBROUTINE that swaps two numbers X <-> Y. !
!--------------------------------------------!

 implicit none

 double precision, intent(out) :: x,y
 double precision              :: tmp

 tmp = x
 x   = y
 y   = tmp

 return

END SUBROUTINE swap

!----------------------------------------------------------------------------------------------------------

DOUBLE PRECISION FUNCTION complex_module(z)
        
!------------------------------------------------------!        
! FUNCTION to calculate the MODULE of a COMPLEX NUMBER !
!  z  = a + i*b                                        !
! |z| = sqrt[z * congj(z)]                             !
!     = sqrt[(a+i*b)(a-i*b)]                           !
!     = a^2 + b^2                                      !
!------------------------------------------------------!

  implicit none

  double complex, intent(in) :: z

  complex_module =  z * conjg(z)

  complex_module = dsqrt(complex_module)

  return

END FUNCTION  complex_module

!----------------------------------------------------------------------------------------------------------

LOGICAL FUNCTION sign_change(x,y)

!-------------------------------------------------------!
! FUNCTION that gives .TRUE. as a result if there is a  !
! sign change between real numbers X and Y, and .FALSE. !
! as result if there is no sign change.                 !
! If X=0.0 or Y=0.0 it is consider to be a sign change, !
! but not if both are zero.                             !
!-------------------------------------------------------!

   implicit none

   double precision, parameter  :: thrs = 1.d-6

   double precision, intent(in) :: x,y

   if ( (abs(x).lt.thrs) .and. (abs(y).lt.thrs) ) then
      sign_change = .FALSE.
   elseif ( (abs(x).lt.thrs) .and. (abs(y).gt.thrs) ) then
      sign_change = .TRUE.
   elseif ( (abs(x).gt.thrs) .and. (abs(y).lt.thrs) ) then
      sign_change = .TRUE.
   elseif ( (abs(x).gt.thrs) .and. (abs(y).gt.thrs) ) then
      if ( (x*y).gt.0 ) then
         sign_change = .FALSE.
      elseif ( (x*y).lt.0 ) then
         sign_change = .TRUE.
      endif
   endif

   return

END FUNCTION sign_change

