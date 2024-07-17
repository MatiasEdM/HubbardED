!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_rsDETS(DMAT,nSITES,nEL,nDET)

!---------------------------------------------------------------------------!
! SUBROUTINE to construct a given set of spin-up or -down determinants      !
! given a set of lattice sites (Ns) and electrons (N) belongin to either    !
! the spin-up or spin-down channel. In practice, it generates all possible  !
! configurations of distributing N electrons in Ns sites, i.e. the binomial !
! coefficient binom(Ns,N) - distributing N 1's in an array of Ns sites.     !
!---------------------------------------------------------------------------!

  implicit none
  
 !-------------------------------------------------!
 ! INOUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8), intent(in)  :: nSITES, nEL, nDET
  integer(8), intent(out) :: DMAT(nSITES,nDET)
  
 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)              :: i

 !-------------------------------------------------!
 ! SUBROUTINE: ALGORITHM FOR DET. GENERATION       !
 !-------------------------------------------------!

  DMAT(:,:) = 0
  
  DMAT(1:nEL,1) = 1

  do i = 2,nDET,1
     call generate_lexicorderDET( DMAT(:,i-1),DMAT(:,i),nSITES,nEL )
  enddo
 
  return

END SUBROUTINE generate_rsDETS

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_lexicorderDET(initCONFIG, nextCONFIG, nSITES, nEL)

!------------------------------------------------------------------------------!
! SUBROUTINE for the generation of configurations by LEXICOGRAPHIC SORTING     !
! Given a certain configuration (Slater Determinanat in the second quantized   !
! formalism) coresponding to  a distribution of nEL electrons in nSITES sites, !
! generates another configuration which corresponds to the following one in    !
! the LEXICOGRAPHIC ORDERING.                                                  !
!------------------------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! INOUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8), intent(in)  :: nSITES, nEL
  integer(8), intent(in)  :: initCONFIG(nSITES)
  integer(8), intent(out) :: nextCONFIG(nSITES)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)              :: i, ctr
  integer(8)              :: pos, mov_pos
  integer(8)              :: count_ones

 !-------------------------------------------------!
 ! SUBROUTINE: ALGORITHM FOR CONFIG. GENERATION    !
 !-------------------------------------------------!
 
  i = pos( initCONFIG(:),nSITES )

  if ( i .eq. nSITES ) then
     i      = mov_pos( initCONFIG(:),nSITES )
     if ( (nSITES - i) .eq. 1) then
        ctr = 1
     else
        ctr = count_ones( initCONFIG(i+1:nSITES),nSITES-i )
     endif
     nextCONFIG(:)               = initCONFIG(:)
     nextCONFIG(i)               = 0
     nextCONFIG(i+1:i+1+ctr)     = 1
     nextCONFIG(i+2+ctr:nSITES)  = 0
  elseif ( i .ne. nSITES ) then
     nextCONFIG(:)   = initCONFIG(:)
     nextCONFIG(i)   = 0
     nextCONFIG(i+1) = 1 
  endif

  return

END SUBROUTINE generate_lexicorderDET

!----------------------------------------------------------------------------------------------------------

INTEGER(8) FUNCTION pos(VECT,n)

!-----------------------------------------------------------!        
! FUNCTION to find the first (1) in a vector of 0's and 1's !
! starting from the last element.                           !
!-----------------------------------------------------------!

 implicit none

 integer(8), intent(in) :: n
 integer(8), intent(in) :: VECT(n)

 integer(8)             :: i

 do i = n,1,-1
    if ( VECT(i) .eq. 1) then
       pos = i
       exit
    endif
  enddo

  return

END FUNCTION pos

!----------------------------------------------------------------------------------------------------------

INTEGER(8) FUNCTION mov_pos(VECT,n)

!-----------------------------------------------------------!
! FUNCTION to find the first movable (1) in a vector of 0's !
! and 1's when the last element is already a (1).           !
!-----------------------------------------------------------!

 implicit none

 integer(8), intent(in) :: n
 integer(8), intent(in) :: VECT(n)

 integer(8)             :: i

 do i = n-1,1,-1
    if ( (VECT(i) .eq. 1) .and. VECT(i+1) .eq. 0 ) then
       mov_pos = i
       exit
    endif
  enddo

  return

END FUNCTION mov_pos

!----------------------------------------------------------------------------------------------------------

INTEGER(8) FUNCTION count_ones(VECT,n)

!------------------------------------------------------!        
! FUNCTION that returns the number of 1's on a vector. !
!------------------------------------------------------!

 implicit none

 integer(8), intent(in) :: n
 integer(8), intent(in) :: VECT(n)

 integer(8)             :: i

 count_ones = 0

 do i = 1,n,1
    if ( VECT(i) .eq. 1) then
       count_ones = count_ones + 1
    endif
 enddo

 return

 END FUNCTION count_ones 
