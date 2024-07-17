!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_doubleOCC(dOCC,upDMAT,downDMAT,nSITES,nDET,nDETup,nDETdown)

!--------------------------------------------------------------------------------!
! SUBROUTINE to construct the DOUBLE OCCUPANCY MAP, dOCC(Ns,Ndet), a matrix that !
! contains in its columns the list of doubly occupied sites for each determinant !
! constructed as the tensor product of the up-spin and down-spin determinants,   !
! i.e. |K> = |Iup> x |Jdown>.                                                    !
!--------------------------------------------------------------------------------!

  implicit none
  
 !--------------------------------------------------!
 ! INOUT VARIABLES                                  !
 !--------------------------------------------------!

  integer(8), intent(in)  :: nSITES
  integer(8), intent(in)  :: nDET, nDETup, nDETdown
  integer(8), intent(in)  :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)
  integer(8), intent(out) :: dOCC(nSITES,nDET)
  
 !--------------------------------------------------!
 ! LOCAL VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)              :: Iup, Jdown
  integer(8)              :: Kdet                    !Each (UP)X(DOWN) determinant Kdet is any pair of 
                                                     !(Iup,Jdown), Iup determinant of the up-spin space 
                                                     !and Jdown determinant of the down-spin space.

  integer(8)              :: site

 !--------------------------------------------------!
 ! SUBROUTINE: ALGORITHM FOR DOUBLE OCCUPANCY GEN.  !
 !--------------------------------------------------!

  dOCC(:,:) = 0

  Kdet = 1

  do Iup = 1,nDETup,1
     do Jdown = 1,nDETdown,1
        do site = 1,nSITES,1
           if ( (upDMAT(site,Iup) .eq. 1) .and. (downDMAT(site,Jdown) .eq. 1 ) ) then
              dOCC(site,Kdet) = 1
           endif
        enddo
        Kdet = Kdet + 1
     enddo
  enddo

  return

END SUBROUTINE generate_doubleOCC
