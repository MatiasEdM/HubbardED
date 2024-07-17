!----------------------------------------------------------------------------------------------------------

SUBROUTINE solve_truncated(upDMAT,downDMAT,H,CI,E,nSITES,nDET,nDETup,nDETdown,U,tot_spin_trunc)

!------------------------------------------------------------------------------------!        
! SUBROUTINE that diagonalizes the Hamiltonian MAT in a truncated space generated    !
! only by determinants with non-zero amplitude in the HF determinant. That is, the   !
! Hilbert space in which the node determinants of the HF many-body wave function     !
! have been eliminated. This is equivalent to an exact diagonalization in the        !
! non-interaction limit (U/t = 0.0), but becomes an approximation in the interacting !
! case ( U /= 0.0), when the many-body wave function in reciprocal space has other   !
! contributions other that the HF solution.                                          !
!------------------------------------------------------------------------------------!

  implicit none


 !--------------------------------------------------!
 ! INOUT VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)      , intent(in) :: nSITES
  integer(8)      , intent(in) :: nDET, nDETup, nDETdown
  integer(8)      , intent(in) :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)
  double precision, intent(in) :: H(nDET,nDET), CI(nDET,nDET), E(nDET)
  double precision, intent(in) :: U
  character(3)    , intent(in) :: tot_spin_trunc

 !--------------------------------------------------!
 ! LOCAL VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)                    :: nDET_tr, nNODES
  integer(8)                    :: I, J

  integer(8)      , allocatable :: NODESindx(:), noNODESindx(:)

  double precision, allocatable :: Ho(:,:)  , CIo(:,:)  , Eo(:)
  double precision, allocatable :: H_tr(:,:), CI_tr(:,:), E_tr(:)

  double precision, allocatable :: CI_rg(:)

  double precision              :: DE_error

  double precision              :: S_tr
  double precision              :: tot_spin_sq

 !--------------------------------------------------!
 ! SUBROUTINE: CALCULATION OF NON-INT. LIMIT U=0.0  !
 !--------------------------------------------------!

  allocate( Ho(nDET,nDET),CIo(nDET,nDET),Eo(nDET),CI_rg(nDET) )
  
  if ( U .ne. 0.0 ) then
     Ho(:,:) = H(:,:)
     do I = 1,nDET,1
        Ho(I,I) = 0.d0
     enddo
     CIo(:,:)   = Ho(:,:)
     call diagonalize_matrix( nDET,CIo(:,:),Eo(:) )
   else if ( U .eq. 0.0 ) then
     CIo(:,:)   = CI(:,:)
     Eo(:)      = E(:)
   endif

 !--------------------------------------------------!
 ! SUBROUTINE: GENERATION OF TRUNCATED QUANTITIES   !
 !--------------------------------------------------!

 ! Count of number of nodes, i.e. number of . to
 ! block and remove from the Hilbert space.

  call count_nodes( nNODES,CIo(:,1),nDET ) 

  nDET_tr = nDET - nNODES

  allocate( H_tr(nDET_tr,nDET_tr),CI_tr(nDET_tr,nDET_tr),E_tr(nDET_tr) )

  allocate( NODESindx(nNODES), noNODESindx(nDET_tr) )

  call get_NODESindx( NODESindx(:),noNODESindx(:),CIo(:,1),nDET,nDET_tr,nNODES )

  do I = 1,nDET_tr,1
     do J = 1,I-1,1
        H_tr(I,J) = H( noNODESindx(I), noNODESindx(J) )
        H_tr(J,I) = H_tr(I,J)
     enddo
     H_tr(I,I)    = H( noNODESindx(I), noNODESindx(I) )
   enddo

 !--------------------------------------------------!
 ! SUBROUTINE: DIAGONALIZATION IN TRUNCATED SPACE   !
 !--------------------------------------------------!

   CI_tr(:,:) = H_tr(:,:)

   call diagonalize_matrix( nDET_tr,CI_tr(:,:),E_tr(:) )
 
 !--------------------------------------------------!
 ! SUBROUTINE: REORGANIZATION OF THE CI COEFF.      !
 !--------------------------------------------------!

   call reorganize_CI( CI_rg(:),CI_tr(:,1),noNODESindx(:),nDET,nDET_tr )

 !--------------------------------------------------!
 ! SUBROUTINE: TOTAL SPIN S CALCULATION             !
 !--------------------------------------------------!

  if ( tot_spin_trunc .eq. 'yes' ) then
      S_tr = tot_spin_sq( CI_rg(:),upDMAT(:,:),downDMAT(:,:),nSITES,nDET,nDETup,nDETdown )   !ATTENTION THIS IS <S^2> NOT S. 
  endif

 !--------------------------------------------------!
 ! SUBROUTINE: WRITE OUTPUT                         !
 !--------------------------------------------------!

   DE_error = E_tr(1) - E(1)

   call write_output_trunc( upDMAT(:,:),downDMAT(:,:),H_tr(:,:),CI_tr(:,:),E_tr(:),noNODESindx(:), &
           & DE_error,tot_spin_trunc,S_tr,nSITES,nDET,nDETup,nDETdown,nDET_tr )

 !--------------------------------------------------!
 ! MEMORY DEALLOCATION                              !
 !--------------------------------------------------!
 
  deallocate( Ho,CIo,Eo,CI_rg )
  deallocate( H_tr,CI_tr,E_tr )
  deallocate( NODESindx,noNODESindx )

  return

END SUBROUTINE solve_truncated

!----------------------------------------------------------------------------------------------------------

SUBROUTINE count_nodes(nNODES,CI,nDET)

!-------------------------------------------------------------!        
! SUBROUTINE to count the number of nodes in the CI expansion !
! of the many-body wave function (number of zero-amplitude    !
! determinants).                                              !
!-------------------------------------------------------------!

  implicit none

  double precision, parameter   :: thrs = 10E-8

  integer(8)      , intent(in)  :: nDET
  double precision, intent(in)  :: CI(nDET)

  integer(8)      , intent(out) :: nNODES

  integer(8)                    :: I

  nNODES = 0

  do I = 1,nDET,1
     if ( abs(CI(I)) .le. thrs ) then
        nNODES = nNODES + 1
     endif
  enddo

 return

END SUBROUTINE count_nodes
     
!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_NODESindx( NODESindx,noNODESindx,CI,nDET,nDET_tr,nNODES )

!---------------------------------------------------------------------!        
! SUBROUTINE to get the (global) indices of the determinants that are !
! nodes and those that are not nodes.                                 !
!---------------------------------------------------------------------!

  implicit none

  double precision, parameter   :: thrs = 10E-8

  integer(8)      , intent(in)  :: nDET, nDET_tr, nNODES
  integer(8)      , intent(out) :: NODESindx(nNODES)
  integer(8)      , intent(out) :: noNODESindx(nDET_tr)

  double precision, intent(in)  :: CI(nDET)

  integer(8)                    :: I, J, K

  J = 0
  K = 0

  do I = 1,nDET,1
     if ( abs(CI(I)) .le. thrs ) then
        J              = J + 1
        NODESindx(J)   = I
     else
        K              = K + 1
        noNODESindx(K) = I
     endif
  enddo

  if ( (J + K) .ne. nDET ) then
     write(*,*) 'ERROR DURING TUNCATED HILBERT SPACE CALCULATION'
     stop
  endif
  
  return

END SUBROUTINE get_NODESindx

!----------------------------------------------------------------------------------------------------------

SUBROUTINE reorganize_CI(CI_rg,CI_tr,noNODESindx,nDET,nDET_tr)

!------------------------------------------------------------------------!        
! SUBROUTINE that reorganizes the truncated CI VEC of dimension nDET_tr  !
! to a CI VEC of dimension nDET to keep the correct original labeling of !
! DETs of the complete Hilbert space, placing 0.0 in those DETs that are !
! excluded in the truncated Hilbert Space.                               !
!------------------------------------------------------------------------!

 implicit none

 integer(8)      , intent(in)  :: nDET
 integer(8)      , intent(in)  :: nDET_tr

 integer(8)      , intent(in)  :: noNODESindx(nDET_tr)

 double precision, intent(in)  :: CI_tr(nDET_tr)

 double precision, intent(out) :: CI_rg(nDET)

 integer(8)                    :: I

 CI_rg = 0.d0

 do I = 1,nDET_tr,1
    CI_rg(noNODESindx(I)) = CI_tr(I)
 enddo

 return

END SUBROUTINE reorganize_CI
