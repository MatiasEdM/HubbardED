DOUBLE PRECISION FUNCTION tot_spin_sq(CI,upDMAT,downDMAT,nSITES,nDET,nDETup,nDETdown)

!---------------------------------------------------!       
! FUNCTION that calculates the EXPECTATION VALUE of !
! TOTAL SPIN (squared) OPERATOR of the MB WFC as:   !
! <Psi| S^2 |Psi>.                                  !
!---------------------------------------------------!

 implicit none

!--------------------------------------------------!
! INOUT VARIABLES                                  !
!--------------------------------------------------!

 integer(8), intent(in)       :: nSITES
 integer(8), intent(in)       :: nDET, nDETup, nDETdown
 integer(8), intent(in)       :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)
 
 double precision, intent(in) :: CI(nDET)

!--------------------------------------------------!
! LOCAL VARIABLES                                  !
!--------------------------------------------------!

 integer(8)                   :: I, J, Iup, Idown, Jup, Jdown
 integer(8)                   :: globalDETindex
 integer(8)                   :: newDETindex

 integer(8)                   :: nOSup, nOSdown               !Number of Open-Shell up- and down-electrons.
 integer(8)                   :: nEXC                         !Number of possible EXCHANGE OPERATIONS a DET admits.
 integer(8), allocatable      :: upINDX(:), downINDX(:)       !INDICES of OS up- and down-electrons in a DET.

 integer(8)                   :: upDETaux(nSITES), downDETaux(nSITES)

 integer(8)                   :: ii, jj

 integer(8)                   :: phase_factor
 integer(8)                   :: sgnIJ

 double precision             :: tot_proj_spin
 double precision             :: Sz
 
!--------------------------------------------------!
! S^2 EXPECTATION VALUE CALCULATION                !
!--------------------------------------------------!

 tot_spin_sq = 0.d0

 do I = 1,nDET,1

    !-------------------------------------------------------------------------------------!
    ! Contribution from the diagonal element of the TOTAL SPIN (squared) OPERATOR in the  !
    ! SD Basis:                                                                           !
    !             < D(I) |S^2| D(I) > = Sz ( Sz - 1 ) + nOSup                             !
    ! where Sz is the TOTAL SPIN PROJECTION and nOSup is the number of open shell up-spin !
    ! electrons in the SD |D(I)>.                                                         !
    ! The total contribution is then weighted  by the CI coefficients:                    !
    !                 CI(I)*CI(I) *  < D(I) |S^2| D(I) >                                  !
    !-------------------------------------------------------------------------------------!

    call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
    
    call open_shell( nOSup,nOSdown,upDMAT(:,Iup),downDMAT(:,Idown),nSITES )

    allocate( upINDX(nOSup), downINDX(nOSdown) )
    
    call get_open_shell_index( upINDX(:), downINDX(:),upDMAT(:,Iup),downDMAT(:,Idown),nSITES,nOSup,nOSdown )

    Sz = tot_proj_spin( upDMAT(:,Iup), downDMAT(:,Idown), nSITES)
    tot_spin_sq = tot_spin_sq + CI(I)*CI(I)*( Sz * ( Sz - 1.d0 ) + dble( nOSUp ) )

    !----------------------------------------------------------------------------------!
    ! Contribution from the off-diagonal elements of the TOTAL SPIN (squared) OPERATOR !
    ! in the SD Basis:                                                                 !
    !                < D(I) |S^2| D(J) > = sgn[D(I),D(J)]                              !
    ! where |D(J)> is obtained as an EXHANGE OPERATION from |D(I)>.                    !
    !----------------------------------------------------------------------------------!

   nEXC = nOSup * nOSdown

   if ( nEXC .ne. 0 ) then
      do ii = 1,nOSup,1
         do jj = 1,nOSdown,1
            
            upDETaux(:)                = upDMAT(:,Iup)
            downDETaux(:)              = downDMAT(:,Idown)
         
            upDETaux(   upINDX(ii)   ) = 0
            upDETaux(   downINDX(jj) ) = 1
            downDETaux( downINDX(jj) ) = 0
            downDETaux( upINDX(ii)   ) = 1

            sgnIJ = phase_factor( upDMAT(:,Iup),downDMAT(:,Idown),upINDX(ii),downINDX(jj),nSITES )
          
            Jup   = newDETindex( upDETaux(:)  ,upDMAT(:,:)  ,nSITES,nDETup   )
            Jdown = newDETindex( downDETaux(:),downDMAT(:,:),nSITES,nDETdown )
            J     = globalDETindex( Jup,Jdown,nDETup,nDETdown )

            tot_spin_sq = tot_spin_sq + CI(I)*CI(J)*sgnIJ
            
         enddo
      enddo

   endif 

   deallocate( upINDX,downINDX )

 enddo

 return

END FUNCTION tot_spin_sq

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_S2MAT( S2,upDMAT,downDMAT,nSITES,nDET,nDETup,nDETdown)

!---------------------------------------------------!       
! FUNCTION that calculates the EXPECTATION VALUE of !
! TOTAL SPIN (squared) OPERATOR of the MB WFC as:   !
! <Psi| S^2 |Psi>.                                  !
!---------------------------------------------------!

 implicit none

!--------------------------------------------------!
! INOUT VARIABLES                                  !
!--------------------------------------------------!

 integer(8)      , intent(in)  :: nSITES
 integer(8)      , intent(in)  :: nDET, nDETup, nDETdown
 integer(8)      , intent(in)  :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)

 double precision, intent(out) :: S2(nDET,nDET)

!--------------------------------------------------!
! LOCAL VARIABLES                                  !
!--------------------------------------------------!

 integer(8)                   :: I, J, Iup, Idown, Jup, Jdown
 integer(8)                   :: globalDETindex
 integer(8)                   :: newDETindex

 integer(8)                   :: nOSup, nOSdown               !Number of Open-Shell up- and down-electrons.
 integer(8)                   :: nEXC                         !Number of possible EXCHANGE OPERATIONS a DET admits.
 integer(8), allocatable      :: upINDX(:), downINDX(:)       !INDICES of OS up- and down-electrons in a DET.

 integer(8)                   :: upDETaux(nSITES), downDETaux(nSITES)

 integer(8)                   :: ii, jj

 integer(8)                   :: phase_factor
 integer(8)                   :: sgnIJ

 double precision             :: tot_proj_spin
 double precision             :: Sz

!--------------------------------------------------!
! S^2 EXPECTATION VALUE CALCULATION                !
!--------------------------------------------------!

 S2(:,:) = 0.d0

 do I = 1,nDET,1

    !-------------------------------------------------------------------------------------!
    ! Contribution from the diagonal element of the TOTAL SPIN (squared) OPERATOR in the  !
    ! SD Basis:                                                                           !
    !             < D(I) |S^2| D(I) > = Sz ( Sz - 1 ) + nOSup                             !
    ! where Sz is the TOTAL SPIN PROJECTION and nOSup is the number of open shell up-spin !
    ! electrons in the SD |D(I)>.                                                         !
    ! The total contribution is then weighted  by the CI coefficients:                    !
    !                 CI(I)*CI(I) *  < D(I) |S^2| D(I) >                                  !
    !-------------------------------------------------------------------------------------!

    call spinDETindex( I,Iup,Idown,nDETup,nDETdown )

    call open_shell( nOSup,nOSdown,upDMAT(:,Iup),downDMAT(:,Idown),nSITES )

    allocate( upINDX(nOSup), downINDX(nOSdown) )

    call get_open_shell_index( upINDX(:), downINDX(:),upDMAT(:,Iup),downDMAT(:,Idown),nSITES,nOSup,nOSdown )

    Sz      = tot_proj_spin( upDMAT(:,Iup), downDMAT(:,Idown), nSITES)
    S2(I,I) =  Sz * ( Sz - 1.d0 ) + dble( nOSUp )

    !------------------------------------------------------------------------------------!
    ! Contribution from the off-diagonal elements of the TOTAL SPIN (squared) OPERATOR   !
    ! in the SD Basis:                                                                   !
    !                < D(I) |S^2| D(J) > = sgn[D(I),D(J)]                                !
    ! where |D(J)> is obtained as an EXHANGE OPERATION from |D(I)>.                      !
    !------------------------------------------------------------------------------------!

    nEXC = nOSup * nOSdown

    if ( nEXC .ne. 0 ) then
       do ii = 1,nOSup,1
          do jj = 1,nOSdown,1

             upDETaux(:)                = upDMAT(:,Iup)
             downDETaux(:)              = downDMAT(:,Idown)

             upDETaux(   upINDX(ii)   ) = 0
             upDETaux(   downINDX(jj) ) = 1
             downDETaux( downINDX(jj) ) = 0
             downDETaux( upINDX(ii)   ) = 1

             sgnIJ   = phase_factor( upDMAT(:,Iup),downDMAT(:,Idown),upINDX(ii),downINDX(jj),nSITES )

             Jup     = newDETindex( upDETaux(:)  ,upDMAT(:,:)  ,nSITES,nDETup   )
             Jdown   = newDETindex( downDETaux(:),downDMAT(:,:),nSITES,nDETdown )
             J       = globalDETindex( Jup,Jdown,nDETup,nDETdown )

             S2(I,J) = sgnIJ

          enddo
       enddo

    endif

 deallocate( upINDX,downINDX )

 enddo

 return

END SUBROUTINE generate_S2MAT

!----------------------------------------------------------------------------------------------------------

DOUBLE PRECISION FUNCTION tot_proj_spin(upDET,downDET,nSITES)

!---------------------------------------------------------------------------------!
! FUNCTION that calculates the EXPECTATION VALUE of the                           !
! TOTAL SPIN PROJECTION of a given determinant |D> composed                       !
! by the tensor product                                                           !
!                      |D> = |D_up> (x) |D_down>                                  !
! as:                                                                             !
! <Sz> = SUM_{i=1,...,nSITES}[ 1/2 ( nUP(i) - nDOWN(i) ) ] = ( nUP - nDOWN ) / 2  !
! where nUP and nDOWN is the total number up spin-up and spin-down electrons.     !
! COMMENT: <Sz> is a constant for this problem because the Hubbard Hamiltonian is !
! diagonalized in a fixed Ms-SUBSPACE.                                            !
!---------------------------------------------------------------------------------!

 implicit none

 integer(8), intent(in) :: nSITES
 integer(8), intent(in) :: upDET(nSITES), downDET(nSITES)

 tot_proj_spin = 0.5d0 * ( dble(sum(upDET(:)))  - dble(sum(downDET(:))) )

 return

END FUNCTION tot_proj_spin

!----------------------------------------------------------------------------------------------------------

SUBROUTINE  open_shell(nOsup,nOSdown,upDET,downDET,nSITES)
    
!------------------------------------------------!        
! SUBROUTINE that calculates the total number of !
! OPEN SHELL SPIN-UPi and SPIN-DOWN electrons    !
! of a given SD |D>.                             !
!------------------------------------------------!

 implicit none

 integer(8), intent(in)  :: nSITES
 integer(8), intent(in)  :: upDET(nSITES), downDET(nSITES)

 integer(8), intent(out) :: nOSUp, nOSdown

 integer(8)             :: site

 nOSup   = 0
 nOSdown = 0

 do site = 1,nSITES,1
    if ( (upDET(site) .eq. 1) .and. (downDET(site) .eq. 0) ) then
       nOSup   = nOSup   + 1
    elseif ( (upDET(site) .eq. 0) .and. (downDET(site) .eq. 1) ) then
       nOSdown = nOSdown + 1
    endif
 enddo

 return

END SUBROUTINE open_shell

!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_open_shell_index( upINDX,downINDX,upDET,downDET,nSITES,nOSup,nOSdown )

!------------------------------------------------------------!        
! SUBROUTINE to get the indices where the OPEN SHELL SPIN-UP !
! and SPIN-DOWN electrons are located.                       !
!------------------------------------------------------------!

 implicit none

 integer(8), intent(in)  :: nSITES
 integer(8), intent(in)  :: nOSup, nOSdown
 integer(8), intent(in)  :: upDET(nSITES), downDET(nSITES)

 integer(8), intent(out) :: upINDX(nOSup), downINDX(nOSdown)

 integer(8)              :: i, j, site

 i = 1
 j = 1

 do site = 1,nSITES,1
     if ( (upDET(site) .eq. 1) .and. (downDET(site) .eq. 0) ) then
       upINDX(i)   = site
       i           = i + 1
    elseif ( (upDET(site) .eq. 0) .and. (downDET(site) .eq. 1) ) then
       downINDX(j) = site
       j           = j + 1
    endif
 enddo


return

END SUBROUTINE get_open_shell_index

!----------------------------------------------------------------------------------------------------------

INTEGER(8) FUNCTION phase_factor(upDET,downDET,site_i,site_j,nSITES)

!-----------------------------------------------------------------!       
! FUNCTION to calculate the PHASE FACTOR of the UP-DOWN-SPIN DET  ! 
! based on how many electrons are jumped over to flip down an     !
! up-electron in site_i and flip up a down-electron in site_j due !
! to FERMIONIC SIGN FLIPS during an EXCHANGE OPERATION.           !
!-----------------------------------------------------------------!

 implicit none

 integer(8), intent(in) :: nSITES
 integer(8), intent(in) :: site_i, site_j
 integer(8), intent(in) :: upDET(nSITES), downDET(nSITES)

 integer(8)             :: DET( nSITES+nSITES )
 integer(8)             :: expn

 !----------------------------------------------------!
 ! [0] Auxiliar DET = [upDET;downDET].                !
 !----------------------------------------------------!

 DET(1:nSITES)               = upDET(:)
 DET(nSITES+1:nSITES+nSITES) = downDET(:)

 !----------------------------------------------------!
 ! [1.1] First Phase Factor Gamma_{up-down} in DET.   !
 !----------------------------------------------------!

 expn = sum( DET(site_i+1:nSITES+site_i-1) )

 phase_factor = (-1) ** expn

 !----------------------------------------------------!
 ! [1.2] DET' resulting from the flipping down of the !
 ! up-electron in site_i.                             !
 !----------------------------------------------------!

 DET(site_i)        = 0
 DET(nSITES+site_i) = 1

 !----------------------------------------------------!
 ! [2.1] Second Phase Factor Gamma_{down-up} in DET'. !
 !----------------------------------------------------!

 expn = sum( DET(site_j+1:nSITES+site_j-1) )

 phase_factor = phase_factor * ( (-1) ** expn )
 
 !----------------------------------------------------!
 ! [2.2] DET'' resulting from the flipping up of the  !
 ! down-electron in site_j (checking purposes).       !
 !----------------------------------------------------!

 DET(site_j)        = 1
 DET(nSITES+site_j) = 0

 return

END FUNCTION phase_factor
