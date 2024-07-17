!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_HAMILTONIAN(H,upDMAT,downDMAT,hopMAPup,hopMAPdown,dOCC,nSITES,nDET,nDETup,nDETdown,t,U)

!----------------------------------------------------------------!
! SUBROUTINE to construct the full real-space HAMILTONIAN MATRIX !
! in the many-body determinantal basis.                          !
!----------------------------------------------------------------!

  implicit none
  
 !--------------------------------------------------!
 ! INOUT VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)      , intent(in)  :: nSITES 
  integer(8)      , intent(in)  :: nDET  , nDETup, nDETdown
  integer(8)      , intent(in)  :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)
  integer(8)      , intent(in)  :: hopMAPup(nDETup,nDETup), hopMAPdown(nDETdown,nDETdown) 
  integer(8)      , intent(in)  :: dOCC(nSITES,nDET)
  double precision, intent(in)  :: t, U

  double precision, intent(out) :: H(nDET,nDET)
 
 !--------------------------------------------------!
 ! LOCAL VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)                    :: I, J, Iup, Idown, Jup, Jdown, Kaux
  integer(8)                    :: globalDETindex
  integer(8)                    :: delta

 !--------------------------------------------------!
 ! SUBROUTINE: GENERATION OF HAMILTONIAN MATRIX     !
 !--------------------------------------------------!

 !--------------------------------------------------!
 ! Calculation of the Hamiltonian matrix elements as:
 ! <I|H|J> = < Iup,Idown | H | Jup,Jdown >          !
 !         = < Iup|H|Jup> <Idown|H|Jdown>           !
 !--------------------------------------------------!

  H(:,:) = 0.d0

  do I = 1,nDET,1
     call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
     do J = 1,I-1,1
        call spinDETindex( J,Jup,Jdown,nDETup,nDETdown )
        do Kaux = 1,nDETup,1
           
           H(I,J) = H(I,J) + hopMAPup(Kaux,Jup)*delta(upDMAT(:,Iup),upDMAT(:,Kaux),nSITES)* &
                   & delta(downDMAT(:,Idown),downDMAT(:,Jdown),nSITES)
        enddo
        do Kaux = 1,nDETdown,1
           H(I,J) = H(I,J) + hopMAPdown(Kaux,Jdown)*delta(upDMAT(:,Iup),upDMAT(:,Jup),nSITES)* &
                   & delta(downDMAT(:,Idown),downDMAT(:,Kaux),nSITES)
        enddo
        H(I,J) = -t * H(I,J)
        H(J,I) = H(I,J)
     enddo
     H(I,I) = sum(dOCC(:,I)) * U
  enddo

  return

END SUBROUTINE generate_HAMILTONIAN

!----------------------------------------------------------------------------------------------------------

INTEGER(8) FUNCTION delta(DET_I,DET_J,nSITES)

!-------------------------------------------------------------------------------------!
! FUNCTION to return the value of DELTA_{I,J} on two determinants belonging to either !
! the up-spin or down-spin subspaces.                                                 !
!-------------------------------------------------------------------------------------!

 implicit none

 integer(8), intent(in) :: nSITES
 integer(8), intent(in) :: DET_I(nSITES), DET_J(nSITES)

 integer(8)             :: site

 delta = 1

 do site = 1,nSITES,1
    if ( DET_I(site) .ne. DET_J(site) ) then
       delta = 0
       exit
    endif
 enddo

 return

END FUNCTION delta
