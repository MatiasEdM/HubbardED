!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_hopMAP(hopMAP,DMAT,LMAT,nSITES,nDET)

!--------------------------------------------------------------------!
! SUBROUTINE to construct the HOPPING MAP for a given spin subspace. !
! C(I,J) = {-1,0,1}                                                  !
!  0 if determinants I and J are disconnected.                       !
! +1 if determinants I and J are connected by an even permutation.   !
! -1 if determinants I and J are connected by an odd permutation.    !
!--------------------------------------------------------------------!

  implicit none
  
 !--------------------------------------------------!
 ! INOUT VARIABLES                                  !
 !--------------------------------------------------!

  integer(8), intent(in)  :: nSITES, nDET
  integer(8), intent(in)  :: DMAT(nSITES,nDET)
  integer(8), intent(in)  :: LMAT(nSITES,nSITES)
  integer(8), intent(out) :: hopMAP(nDET,nDET)
  
 !--------------------------------------------------!
 ! LOCAL VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)              :: I, J
  integer(8)              :: diff                   !Number of differences in the determinants I and J.
  integer(8)              :: latt(2)                !In case of diff=2, lattice sites (i) and (j) of the
                                                     !hopping of the electron, otherwisre is garbage.
  integer(8)              :: p, parity_IJ           !Parity of connected determinants.

 !--------------------------------------------------!
 ! SUBROUTINE: ALGORITHM FOR HOPPING MAP GENERATION !
 !--------------------------------------------------!

  hopMAP(:,:) = 0

  do I = 1,nDET,1
     do J = 1,I-1,1
        call count_diff( diff,latt(:),DMAT(:,I),DMAT(:,J),nSITES)
        if ( diff .ne. 2 ) then 
           hopMAP(I,J) = 0
           hopMAP(J,I) = 0
        elseif ( (diff .eq. 2) .and. (LMAT(latt(1),latt(2)) .eq. 0) ) then
           hopMAP(I,J) = 0
           hopMAP(J,I) = 0
        elseif ( (diff .eq. 2) .and. (LMAT(latt(1),latt(2)) .eq. 1) ) then
           p           = parity_IJ( DMAT(:,I),DMAT(:,J),latt(:),nSITES )
           hopMAP(I,J) = p
           hopMAP(J,I) = p
        elseif ( (diff .eq. 2) .and. (LMAT(latt(1),latt(2)) .eq. 2) ) then
           p           = parity_IJ( DMAT(:,I),DMAT(:,J),latt(:),nSITES )
           hopMAP(I,J) = 2 * p
           hopMAP(J,I) = 2 * p
        endif    
     enddo
  enddo

  return

END SUBROUTINE generate_hopMAP

!----------------------------------------------------------------------------------------------------------

SUBROUTINE count_diff(diff,latt,DET_I,DET_J,nSITES)

!----------------------------------------------------------------------------!
! SUBROUTINE to count the number of differences in two determinants I and J, !
! and in case of 2-differences, it returns the indices of the lattice sites  !
! involved in the electron hopping process.                                  !
!----------------------------------------------------------------------------!

 implicit none

 integer(8), intent(in)  :: nSITES
 integer(8), intent(in)  :: DET_I(nSITES), DET_J(nSITES)

 integer(8), intent(out) :: diff
 integer(8), intent(out) :: latt(2)

 integer(8)              :: i, k
 logical                 :: tmp(nSITES)

 diff    = 0
 latt(1) = 0
 latt(2) = 0
 k       = 1

 do i = 1,nSITES,1
    if ( DET_I(i) .ne. DET_J(i) ) then
       diff    = diff + 1
       latt(k) = i
       k       = 2
     endif
 enddo

 return

END SUBROUTINE count_diff

!----------------------------------------------------------------------------------------------------------

INTEGER(8) FUNCTION  parity_IJ(DET_I,DET_J,latt,nSITES)


!-----------------------------------------------------------------!
! FUNCTION that returns the PARITY of two connected determinants. !
!-----------------------------------------------------------------!

 implicit none

 integer(8), intent(in) :: nSITES
 integer(8), intent(in) :: DET_I(nSITES), DET_J(nSITES)
 integer(8), intent(in) :: latt(2)

 integer(8)             :: i, ctr

 ctr  = 0

 do i = latt(1),latt(2),1
    if ( (DET_I(i) .eq. 1) .and. (DET_J(i) .eq. 1) )  then
       ctr = ctr + 1
    endif
 enddo

 if ( mod(ctr,2) .eq. 0 ) then
    parity_IJ =  1
 else
    parity_IJ = -1
 endif

 return

END FUNCTION parity_IJ
