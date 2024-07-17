INTEGER(8) FUNCTION globalDETindex(Iup,Jdown,nDETup,nDETdown)

!----------------------------------------------------------------------------------------------------------------------------------------!
! FUNCTION that maps the indices of the up-spin determinants (Iup) and the down-spin                                                     !
! determinants (Jdown) to the global numbering of the composed determinants, with the                                                    !
! convention:                                                                                                                            !
! (Iup,Jdown) = (1,1) , (1,2), ..., (1,nDETdown), (2,1), (2,2), ..., (2,nDETdown), ..., (nDETup,1), (nDETup,2), ..., (nDETup,nDETdown)   !
!    (K)      =  (1)  ,  (2) , ...,  (nDETdown) , (nDETdown+1), ...,               ((nDETdown-1)xnDETup)      , ..., (nDETup x nDETdown) !                                                 
! which gives the following mapping                                                                                                      !
!                                  K = ( Iup - 1 ) x nDETdown + Jdown for Iup = {1,...,nDETup} and Jup = {1,...,nDETdown}.               !
!----------------------------------------------------------------------------------------------------------------------------------------!

 implicit none

 integer(8), intent(in) :: nDETUp, nDETdown
 integer(8), intent(in) :: Iup, Jdown

 globalDETindex = ( Iup - 1 ) * nDETdown + Jdown 

 return

END FUNCTION globalDETindex

!----------------------------------------------------------------------------------------------------------

SUBROUTINE spinDETindex(K,Kup,Kdown,nDETup,nDETdown)

!------------------------------------------------------------------------------------------!
! SUBROUTINE inverse of globalDETindex, that for a given determinant K in global notation, !
! returns the indices of the up-spin and down-spin determinants (Iup,Jdown).               !
!------------------------------------------------------------------------------------------!

 implicit none

 integer(8), intent(in)  :: nDETUp, nDETdown
 integer(8), intent(in)  :: K
 integer(8), intent(out) :: Kup, Kdown

 if ( mod(K,nDETdown) .eq. 0 ) then
    Kup   = K/nDETdown
    Kdown = nDETdown 
 else
    Kup   = int( K/nDETdown ) + 1
    Kdown = mod( K,nDETdown )
 endif

 return

END SUBROUTINE spinDETindex

!----------------------------------------------------------------------------------------------------------

INTEGER(8) FUNCTION newDETindex(DET,DMAT,nSITES,nDETs)

!--------------------------------------------------------------------!
! FUNCTION that given an UP-SPIN or DOWN-SPIN DETERMINANT determines !
! its index:                                                         !
!            |K> = |Kup> (x) |Kdown>                                 !
! determines either Kup or Kdown (depending on the arguments.        !
! -------------------------------------------------------------------!

 implicit none

 integer(8), intent(in) :: nSITES
 integer(8), intent(in) :: nDETs
 integer(8), intent(in) :: DMAT(nSITES,nDETs)
 integer(8), intent(in) :: DET(nSITES)

 integer(8)             :: I
 integer(8)             :: site

 logical                :: check

 do I = 1,nDETs,1
    check = .TRUE.
    do site = 1,nSITES,1
       if ( DET(site) .ne. DMAT(site,I) ) then
          check = .FALSE.
          exit
       endif
    enddo
    if ( check ) then
       exit
    endif
 enddo

 newDETindex = I

 return

END FUNCTION newDETindex

!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_OCC_UNOCC_index(OCC,UNOCC,DET,N,nEL)

!--------------------------------------------------------------------!
! SUBROUTINE that given an UP-SPIN or DOWN-SPIN DETERMINANT of dim.  !
! N (nSITES / nKPOINTS) hosting nEL electrons returns a list with    !
! indices (following the adopted lableling of SITES or KPOINTS) of   !
! the OCCUPIED and UNOCCUPIED single-particle states (SITE / KPOINT) !
!            OCC( nEL )                                              !
!            UNOCC( N-nEL )                                          !
! -------------------------------------------------------------------!

  implicit none

  integer(8), intent(in)  :: N
  integer(8), intent(in)  :: nEL
  integer(8), intent(in)  :: DET(N)
  integer(8), intent(out) :: OCC(nEL)
  integer(8), intent(out) :: UNOCC(N-nEL)

  integer(8)              :: I
  integer(8)              :: r, s

  r = 0
  s = 0

  do I = 1,N,1
     if ( DET(I) .eq. 1 ) then
        r        = r + 1
        OCC(r)   = I
     elseif ( DET(I) .eq. 0 ) then
        s        = s + 1
        UNOCC(s) = I
     endif
  enddo

  if ( (r.ne.nEL).or.(s.ne.(N-nEL)) ) then
     write(*,*) 'AN ERROR HAS OCCURRED WITH get_OCC_UNOCC_indx SUBROUTINE'
     stop
  endif

  return

END SUBROUTINE get_OCC_UNOCC_index
