!----------------------------------------------------------------------------------------------------------

SUBROUTINE solve_stoquastized(H,E,nSITES,nDET,nDETup,nDETdown,upDMAT,downDMAT,tot_spin_stoq)

!---------------------------------------------------------------------------------------------!
! SUBROUTINE to solve the STOQUASTIZED PROBLEM by the diagonalization of H_stoq., defined as: !
!  H_stoq(I,J) = -| H(I,J) | for I /= J                                                       !
!  H_stoq(I,J) =    H(I,J)   for I =  J                                                       !        
! It returns the CI VEC and the EIGENVALUES of the stoquastized wave function.                !
!---------------------------------------------------------------------------------------------!

  implicit none

 !--------------------------------------------------!
 ! INOUT VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)      , intent(in)  :: nSITES
  integer(8)      , intent(in)  :: nDET, nDETup, nDETdown
  
  integer(8)      , intent(in)  :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)

  character(3)    , intent(in)  :: tot_spin_stoq

  double precision, intent(in)  :: H(nDET,nDET)
  double precision, intent(in ) :: E(nDET)

 !--------------------------------------------------!
 ! LOCAL VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)                    :: I, J

  double precision, allocatable :: CI_stoq(:,:), E_stoq(:)

  double precision              :: DE_stoq

  double precision              :: S_stoq
  double precision              :: tot_spin_sq

 !--------------------------------------------------!
 ! MEMORY ALLOCATION                                !
 !--------------------------------------------------!

  allocate( CI_stoq(nDET,nDET), E_stoq(nDET) )

 !--------------------------------------------------!
 ! SUBROUTINE: ALGORITHM FOR H_STOQ GENERATION      !
 !--------------------------------------------------!
 
  !-------------------------------------------------!
  ! H_stoq is stored in CI_stoq because it will be  !
  ! overwritten by the EIGENVECTORS during DSYEV    !
  ! diagonalization.                                !
  !-------------------------------------------------!

  CI_stoq(:,:)   = 0.d0

  do I = 1,nDET,1
     CI_stoq(I,I) = H(I,I)
     do J = 1,I-1,1
        CI_stoq(I,J) = - abs( H(I,J) )
        CI_stoq(J,I) = CI_stoq(I,J)
     enddo
  enddo

 !--------------------------------------------------!
 ! SUBROUTINE: DIAGONALIZATION OF H_STOQ            !
 !--------------------------------------------------!

  call diagonalize_matrix( nDET,CI_stoq(:,:),E_stoq(:) )

 !--------------------------------------------------!
 ! SUBROUTINE: TOTAL SPIN S CALCULATION             !
 !--------------------------------------------------!

  if ( tot_spin_stoq .eq. 'yes' ) then
     S_stoq = tot_spin_sq( CI_stoq(:,1),upDMAT(:,:),downDMAT(:,:),nSITES,nDET,nDETup,nDETdown )   !ATTENTION THIS IS <S^2> NOT S. 
  endif

 !--------------------------------------------------!
 ! SUBROUTINE: WRITING OUTPUTS                      !
 !--------------------------------------------------!

  DE_stoq = E(1) - E_stoq(1)
  call write_output_stoq( CI_stoq(:,:),E_stoq(:),DE_stoq,tot_spin_stoq,S_stoq,nSITES,nDET,nDETup,nDETdown, &
          & upDMAT,downDMAT ) 

 !--------------------------------------------------!
 ! MEMORY DEALLOCATION                              !
 !--------------------------------------------------!

  deallocate( CI_stoq, E_stoq )


  return

END SUBROUTINE solve_stoquastized

!----------------------------------------------------------------------------------------------------------

SUBROUTINE solve_stoquastized_tr(H_tr,E_tr,noNODESindx,nSITES,nDET,nDET_tr,nDETup,nDETdown, &
                & upDMAT,downDMAT,tot_spin_stoq)

!---------------------------------------------------------------------------------------------!
! SUBROUTINE to solve the STOQUASTIZED PROBLEM by the diagonalization of H_stoq., defined as: !
!  H_stoq(I,J) = -| H(I,J) | for I /= J                                                       !
!  H_stoq(I,J) =    H(I,J)   for I =  J                                                       !        
! It returns the CI VEC and the EIGENVALUES of the stoquastized wave function.                !
!---------------------------------------------------------------------------------------------!

  implicit none

 !--------------------------------------------------!
 ! INOUT VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)      , intent(in)  :: nSITES
  integer(8)      , intent(in)  :: nDET
  integer(8)      , intent(in)  :: nDET_tr
  integer(8)      , intent(in)  :: nDETup, nDETdown
  
  integer(8)      , intent(in)  :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)

  integer(8)      , intent(in)  :: noNODESindx(nDET_tr)

  character(3)    , intent(in)  :: tot_spin_stoq

  double precision, intent(in)  :: H_tr(nDET_tr,nDET_tr)
  double precision, intent(in ) :: E_tr(nDET_tr)

 !--------------------------------------------------!
 ! LOCAL VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)                    :: I, J

  double precision, allocatable :: CI_stoq(:,:), E_stoq(:)
  double precision, allocatable :: CI_rg(:)

  double precision              :: DE_stoq

  double precision              :: S_stoq
  double precision              :: tot_spin_sq

 !--------------------------------------------------!
 ! MEMORY ALLOCATION                                !
 !--------------------------------------------------!

  allocate( CI_stoq(nDET_tr,nDET_tr), E_stoq(nDET_tr) )
  allocate( CI_rg(nDET) )

 !--------------------------------------------------!
 ! SUBROUTINE: ALGORITHM FOR H_STOQ GENERATION      !
 !--------------------------------------------------!
 
  !-------------------------------------------------!
  ! H_stoq is stored in CI_stoq because it will be  !
  ! overwritten by the EIGENVECTORS during DSYEV    !
  ! diagonalization.                                !
  !-------------------------------------------------!

  CI_stoq(:,:)   = 0.d0

  do I = 1,nDET_tr,1
     CI_stoq(I,I) = H_tr(I,I)
     do J = 1,I-1,1
        CI_stoq(I,J) = - abs( H_tr(I,J) )
        CI_stoq(J,I) = CI_stoq(I,J)
     enddo
  enddo

 !--------------------------------------------------!
 ! SUBROUTINE: DIAGONALIZATION OF H_STOQ            !
 !--------------------------------------------------!

  call diagonalize_matrix( nDET_tr,CI_stoq(:,:),E_stoq(:) )

 !--------------------------------------------------!
 ! SUBROUTINE: REORGANIZATION OF THE CI COEFF.      !
 !--------------------------------------------------!

  call reorganize_CI( CI_rg(:),CI_stoq(:,1),noNODESindx(:),nDET,nDET_tr )

 !--------------------------------------------------!
 ! SUBROUTINE: TOTAL SPIN S CALCULATION             !
 !--------------------------------------------------!

  if ( tot_spin_stoq .eq. 'yes' ) then
     S_stoq = tot_spin_sq( CI_rg(:),upDMAT(:,:),downDMAT(:,:),nSITES,nDET,nDETup,nDETdown )   !ATTENTION THIS IS <S^2> NOT S. 
  endif

 !--------------------------------------------------!
 ! SUBROUTINE: WRITING OUTPUTS                      !
 !--------------------------------------------------!

  DE_stoq = E_tr(1) - E_stoq(1)
  call write_output_stoq_tr( CI_stoq(:,:),E_stoq(:),DE_stoq,tot_spin_stoq,S_stoq,nSITES,nDET,nDET_tr, &
          & nDETup,nDETdown,upDMAT,downDMAT ) 

 !--------------------------------------------------!
 ! MEMORY DEALLOCATION                              !
 !--------------------------------------------------!

  deallocate( CI_stoq, E_stoq )
  deallocate( CI_rg )

  return

END SUBROUTINE solve_stoquastized_tr
