SUBROUTINE get_TWFC_CIEXPANSION( CIindx,CI_K,nDET,nTWFCDET )

!------------------------------------------------------------------!
! SUBROUTINE to get the (global) indices of the first nTWFCDET and !
! most-contributing DET. in K-SPACE of the exact GS WFC[k] for the !
! short CI-EXPANSION of the MC-TWFC[k].                            !
!          TWFC    = SUM_{D=1,...,nTWFCDET} CI(D) |D[K]>           ! 
!------------------------------------------------------------------!

 implicit none

 !---------------------------------------------------!
 ! PARAMETERS                                        !
 !---------------------------------------------------!

  double precision, parameter   :: eps = 1.E-15

 !---------------------------------------------------!
 ! INPUT VARIABLES                                   !
 !---------------------------------------------------!

  integer(8)      , intent(in)  :: nDET
  integer(8)      , intent(in)  :: nTWFCDET

  double precision, intent(in)  :: CI_K(nDET)

  integer(8)      , intent(out) :: CIindx(nTWFCDET)

 !---------------------------------------------------!
 ! LOCAL VARIABLES                                   !
 !---------------------------------------------------!

  double precision, allocatable :: CIaux(:)

  integer(8)                    :: D
  integer(8)                    :: I

 !---------------------------------------------------!
 ! SUBROUTINE: GET THE INDICES                       !
 !---------------------------------------------------!

  allocate( CIaux(nDET) )

  CIaux(:) = CI_K(:)
  do D = 1,nDET,1
     if ( CI_K(D) .eq. 0.d0 ) then
        CIaux(D) = eps
     endif
  enddo

  do D = 1,nTWFCDET,1
     I         = maxloc(abs(CIaux(:)),nDET)
     CIindx(D) = I
     CIaux(I)  = 0.d0
  enddo

  deallocate( CIaux )

  return

END SUBROUTINE get_TWFC_CIEXPANSION

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_kDET(kDET,rot_kDET,KDMAT,KPOINTS,rot_KP,nEL,nKPOINTS)

!-------------------------------------------------------!
! SUBROUTINE to construct the DET. in K-SPACE by        !
! occupying the single-particle states |K> according to !
! the given INDEX Ik to extract information from the    !
! list of DETs. KDMAT.for the up-electrons or the       !
! down-electrons.                                       !
!-------------------------------------------------------!

  implicit none

 !--------------------------------------------------!
 ! PARAMETERS                                       !
 !--------------------------------------------------!

  integer(8)       , parameter  :: dm = 2
  
 !--------------------------------------------------!
 ! INOUT VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)      , intent(in)  :: nEL
  integer(8)      , intent(in)  :: nKPOINTS
  integer(8)      , intent(in)  :: KDMAT(nKPOINTS)
  integer(8)      , intent(in)  :: rot_KP(nKPOINTS)
  integer(8)      , intent(out) :: rot_kDET(nEL)
  
  double precision, intent(in)  :: KPOINTS(dm,nKPOINTS)
  double precision, intent(out) :: kDET(dm,nEL)

 !--------------------------------------------------!
 ! LOCAL VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)                    :: kpoint
  integer(8)                    :: iEL

 !--------------------------------------------------!
 ! SUBROUTINE: GENERATION OF K-SPACE HF DET.        !
 !--------------------------------------------------!

  kDET(:,:)   = 0.d0
  rot_kDET(:) = 0
  iEL         = 0

  do kpoint = 1,nKPOINTS,1
     if ( KDMAT(kpoint) .eq. 1 ) then 
        iEL           = iEL + 1
        kDET(1,iEL)   = KPOINTS(1,kpoint)
        kDET(2,iEL)   = KPOINTS(2,kpoint)
        rot_kDET(iEL) = rot_KP(kpoint)
      endif 
  enddo

  if ( iEL .ne. nEL ) then
     write(*,*) 'AN ERROR HAS OCCURED WITH THE CONSTRUCTION OF THE kDET (generate_kDET)'
     stop
  endif

  return

END SUBROUTINE generate_kDET

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_SGN( SGN,TWFC,H,nDET,Vsf_MAX )

!--------------------------------------------------------------------------------------------!
! SUBROUTINE that calculates the SGN MAT. to categorize the sign-problematic transitions     !
! |R> -> |R'> based on sign changes of the TRIAL WFC (HF) projected on connected |R>-CONFIG. !
! i.e. DETS.                                                                                 !
! The matrix elements <R|SGN|R'> are calculated as:                                          !
!                                                   <R|SGN|R'> = <HF|R'> <R'|H|R> <R|HF>     !
! where the projection of the TRIAL WFC (HF) corresponds to,                                 !
!                                                           <R|HF> = det[M]                  !
!--------------------------------------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  double precision, parameter   :: thrs    = 1.D-6 ! Threshold value below which the projection <R|HF> ~ 0. [CHECK WITH RESPECT TO EXPECTED QUANTITIES

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)      , intent(in)  :: nDET

  double precision, intent(in)  :: H(nDET,nDET)
  double precision, intent(in)  :: Vsf_MAX

  double precision, intent(in)  :: TWFC(nDET)

  double precision, intent(out) :: SGN(nDET,nDET)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                    :: I
  integer(8)                    :: J

!  logical                       :: sign_change
!  logical                       :: auxRE, auxIM

 !-------------------------------------------------!
 ! SUBROUTINE: CONSTRUCT THE SGN MATRIX            !
 !-------------------------------------------------!

  SGN(:,:) = 0.d0

  do I = 1,nDET,1
    do J = 1,nDET,1
       if ( J .eq. I ) cycle
       if ( H(I,J) .ne. 0 ) then
          if ( abs(TWFC(I)) .lt. thrs ) then
             SGN(I,J) = Vsf_MAX
          else
             SGN(I,J) = H(I,J) * ( ( TWFC(J) )/( TWFC(I) ) )
          endif
       endif
    enddo
  enddo

  return

END SUBROUTINE generate_SGN

!----------------------------------------------------------------------------------------------------------

DOUBLE PRECISION FUNCTION SF_POT( I,SGN,H,nDET,Vsf_MAX )

!--------------------------------------------------------------------------------------------!
! FUNCTION that calculates the SIGN FLIP POTENTIAL Vsf(R) contribution to the CONFIG. (DET.) !
! I (|R>) as the diagonal matrix element:                                                    !
!                                        <R|Vsf|R> = SUM_{R'} <R|H|R'> <R'|HF>/<R|HF>        !
!                                                                                            !      
! where the sum runs over all CONFIG. (DET.) for which the transition |R'> -> |R> is SIGN-   !
! VIOLATING, i.e. <R'|SGN|R> > 0                                                             !
! The projection of the TRIAL WFC (HF) corresponds to,                                       !
!                                                           <R|HF> = det[M]                  !
!--------------------------------------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)      , intent(in)  :: I
  integer(8)      , intent(in)  :: nDET

  double precision, intent(in)  :: Vsf_MAX
  double precision, intent(in)  :: H(nDET,nDET)
  double precision, intent(in)  :: SGN(nDET,nDET)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                    :: J

 !-------------------------------------------------!
 ! SUBROUTINE: CALCULATE THE SIGN-FLIP POTENTIAL   !
 !-------------------------------------------------!

  SF_POT = 0.d0
  do J = 1,nDET,1
     if ( SGN(I,J) .eq. Vsf_MAX ) then
        SF_POT = Vsf_MAX
        exit
     elseif ( SGN(I,J) .gt. 0.d0 ) then
        SF_POT = SF_POT + SGN(I,J)
     endif
  enddo

  return

END FUNCTION SF_POT

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_Heff(Heff,Vsf,SGN,H,nDET,Vsf_MAX)

!-----------------------------------------------------------------------!
! SUBROUTINE to generate the EFFECTIVE (FNA) HAMILTONIAN Heff according !
! to the FIXED-NODE APPROX. and the (diagonal) SIGN-FLIP POTENTIAL.     !
! The HAM. MAT. elements <R|Heff|R'> are calculated as:                 !
! - DIAGONAL    :                                                       !
!                <R|Heff|R>  = <R|H|R> + <R|Vsf|R>                      !
! - OFF-DIAGONAL:                                                       !
!                              -                                        !
!                             |  <R|H|R'> if <R|SGN|R'> < 0             !
!                <R|Heff|R'> =|                                         !
!                             |     0     if <R|SGN|R'> > 0 [OR A NODE] !
!                              -                                        !
!-----------------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)      , intent(in)  :: nDET

  double precision, intent(in)  :: H(nDET,nDET)
  double precision, intent(in)  :: SGN(nDET,nDET)
  double precision, intent(in)  :: Vsf_MAX

  double precision, intent(out) :: Heff(nDET,nDET)
  double precision, intent(out) :: Vsf(nDET)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                    :: I, J

  double precision              :: SF_POT

 !-------------------------------------------------!
 ! SUBROUTINE: EFFECTIVE (FNA) HAMILTONIAN MAT.    !
 !-------------------------------------------------!

  Heff(:,:) = 0.d0
  Vsf(:)    = 0.d0

  do I = 1,nDET,1
     Vsf(I)    = SF_POT( I,SGN(:,:),H(:,:),nDET,Vsf_MAX )
     Heff(I,I) = H(I,I) + Vsf(I)
  enddo

  do I = 1,nDET,1
     do J = 1,nDET,1
        if ( J .eq. I ) cycle
        if ( SGN(I,J) .lt. 0 ) then
           Heff(I,J) = H(I,J)
        endif
     enddo
  enddo

  return

END SUBROUTINE generate_Heff

!----------------------------------------------------------------------------------------------------------

SUBROUTINE count_TWFCnodes(nNODES,Heff,nDET,Vsf_MAX)

!-------------------------------------------------------------!
! SUBROUTINE to count the number of nodes in the FNA EFF. HAM.!
! according to the SIGN-FLIP POTENTIAL:                       !
! - If Heff(I,I) = Vmax then the CONFIG. I is considered to be!
!   a NODE of the TRIAL WFC, i.e. <i|Psi_T> = 0.0             !
!-------------------------------------------------------------!

  implicit none

  integer(8)      , intent(in)  :: nDET

  double precision, intent(in)  :: Vsf_MAX
  double precision, intent(in)  :: Heff(nDET,nDET)

  integer(8)      , intent(out) :: nNODES

  integer(8)                    :: I

  nNODES = 0

  do I = 1,nDET,1
     if ( abs(Heff(I,I)) .ge. Vsf_MAX ) then
        nNODES = nNODES + 1
     endif
  enddo

 return

END SUBROUTINE count_TWFCnodes

!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_TWFC_NODESindx( NODESindx,noNODESindx,Heff,nDET,nDET_tr,nNODES,Vsf_MAX )

!---------------------------------------------------------------------!
! SUBROUTINE to get the (global) indices of the determinants that are !
! nodes and those that are not nodes.                                 !
!---------------------------------------------------------------------!

  implicit none

  integer(8)      , intent(in)  :: nDET
  integer(8)      , intent(in)  :: nDET_tr
  integer(8)      , intent(in)  :: nNODES
  integer(8)      , intent(out) :: NODESindx(nNODES)
  integer(8)      , intent(out) :: noNODESindx(nDET_tr)

  double precision, intent(in)  :: Vsf_MAX
  double precision, intent(in)  :: Heff(nDET,nDET)

  integer(8)                    :: I, J, K

  J = 0
  K = 0

  do I = 1,nDET,1
     if ( abs(Heff(I,I)) .ge. Vsf_MAX ) then
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

END SUBROUTINE get_TWFC_NODESindx

!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_Heff_tr(Heff_tr,Heff,noNODESindx,nDET,nDET_tr)


!---------------------------------------------------------------------!
! SUBROUTINE to construct the TRUNCATED FNA EFF. HAM. from the FNA    !
! EFF. HAM. by keeping those CONFIG. (DET.) that are not NODES of the !
! TRIAL WFC, i.e. whose SIGN FLIP POTENTIAL is not Vmax.              !
!---------------------------------------------------------------------!

  implicit none

  integer(8)      , intent(in)  :: nDET
  integer(8)      , intent(in)  :: nDET_tr
  integer(8)      , intent(in)  :: noNODESindx(nDET_tr)

  double precision, intent(in)  :: Heff(nDET,nDET)

  double precision, intent(out) :: Heff_tr(nDET_tr,nDET_tr)

  integer(8)                    :: I, J

  Heff_tr(:,:) = 0.d0

  do I = 1,nDET_tr,1
     do J = 1,I-1,1
        Heff_tr(I,J) = Heff( noNODESindx(I), noNODESindx(J) )
        Heff_tr(J,I) = Heff_tr(I,J)
     enddo
     Heff_tr(I,I)    = Heff( noNODESindx(I), noNODESindx(I) )
  enddo

  return

END SUBROUTINE get_Heff_tr

!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_SUBSPACE_COEFF( CI_TWFC,CIk,CIindx,Hk,nDET,nTWFCDET )

!---------------------------------------------------------------------!
! SUBROUTINE to re-calculate the CI COEFFICIENTS of the (short) CI    !
! EXPANSION of the MC-TWFC by EXACT DIAGONALIZATION in the SUBSPACE   !
! spanned by the DET[k] (K-SPACE DET.) included in the TWFC.          !
!---------------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)      , intent(in)  :: nDET
  integer(8)      , intent(in)  :: nTWFCDET
  integer(8)      , intent(in)  :: CIindx(nTWFCDET)

  double precision, intent(in)  :: Hk(nDET,nDET)
  double precision, intent(in)  :: CIk(nDET)

  double precision, intent(out) :: CI_TWFC(nDET)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  double precision, allocatable :: Hs(:,:)
  double precision, allocatable :: CIs(:,:)
  double precision, allocatable :: Es(:)

  integer(8)                    :: I
  integer(8)                    :: J

 !-------------------------------------------------!
 ! SUBROUTINE: GET THE SUBSPACE HAMILTONIAN        !
 !-------------------------------------------------!

  allocate( Hs(nTWFCDET,nTWFCDET)  )
  allocate( CIs(nTWFCDET,nTWFCDET) )
  allocate( Es(nTWFCDET)           )

  Hs(:,:) = 0.d0

  do I = 1,nTWFCDET,1
     do J = 1,nTWFCDET,1
        Hs(I,J) = Hk( CIindx(I), CIindx(J) )
     enddo
  enddo

 !-------------------------------------------------!
 ! SUBROUTINE: DIAGONALIZE THE SUBSPACE HAMILT.    !
 !-------------------------------------------------!

  CIs(:,:) = Hs(:,:)

  call diagonalize_matrix( nTWFCDET,CIs(:,:),Es(:) )

 !-------------------------------------------------!
 ! SUBROUTINE: GET THE MODIFIED TWFC COEFFICIENTS  !
 !-------------------------------------------------!

  CI_TWFC(:) = CIk(:)

  do I = 1,nTWFCDET,1
     CI_TWFC( CIindx(I) ) = CIs(I,1)
  enddo

  deallocate( Hs  )
  deallocate( CIs )
  deallocate( Es  )

  return

END SUBROUTINE get_SUBSPACE_COEFF

!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_TWF_EN( E_T,TWFC,H,nDET )

!---------------------------------------------------------------------!
! SUBROUTINE to calculate the TRIAL ENERGY of the TRIAL WAVE FUNCTION !
! of the FIXED NODE APPROXIMATION as:                                 !
!                                                                     !
!       E[T] = SUM_{R,R'} <HF|R> <R|H|R'> <R'|HF>                     !
!---------------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)      , intent(in)  :: nDET

  double precision, intent(in)  :: H(nDET,nDET)

  double precision, intent(in)  :: TWFC(nDET)

  double precision, intent(out) :: E_T

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  double precision              :: ET_aux
  double precision              :: normCNST

  integer(8)                    :: I
  integer(8)                    :: J

 !-------------------------------------------------!
 ! SUBROUTINE: CALCULATE THE TRIAL WFC ENERGY      !
 !-------------------------------------------------!

  ET_aux   = 0.d0
  normCNST = 0.d0


  do I = 1,nDET,1
     normCNST  = normCNST + TWFC(I) * TWFC(I)
     do J = 1,nDET,1
        ET_aux = ET_aux   + TWFC(I) * H(I,J) * TWFC(J)
     enddo
  enddo

  E_T = ET_aux / normCNST

  return

END SUBROUTINE get_TWF_EN

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_Vsf_gsWFC( Vsf,H,CI,nDET )

!---------------------------------------------------------------!
! SUBROUTINE to calculate the SIGN FLIP POTENTIAL of the EXACT  !
! FCI WFC as:                                                   !
!                                                               !
!             Vsf = SUM_{R' in S} <R|H|R'> <R'|PSI> / <R|PSI>   !
!                                                               !
! The projections on to DET. |R> are the FCI COEFF.             !
!                                                               !
!                     <R|PSI> = C(R)                            !
!                                                               !
! and the summation runs over:                                  !
! for a given |R> , |R'> in S / S(R,R')= C(R) H(R,R') C(R') > 0 !
!---------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!
  
  double precision, parameter  :: thrs = 10E-8
  double precision, parameter  :: Vmax = 10E+4

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)      , intent(in)  :: nDET

  double precision, intent(in)  :: H(nDET,nDET)
  double precision, intent(in)  :: CI(nDET)

  double precision, intent(out) :: Vsf(nDET)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                    :: I,J
 
  double precision              :: S

 !-------------------------------------------------!
 ! SUBROUTINE: CALCULATION OF Vsf FOR THE GS-WFC   !
 !-------------------------------------------------!

  Vsf(:) = 0.d0

  do I = 1,nDET,1
     if ( abs(CI(I)) .lt. thrs ) then
        Vsf(I) = Vmax
        cycle
     endif
     do J = 1,nDET,1
        if ( J .eq. I ) cycle
        S  = CI(I) * H(I,J) * CI(J)
        if ( S .gt. 0.d0 ) then
           Vsf(I) = Vsf(I) + H(I,J) * ( CI(J)/CI(I) )
        endif
     enddo
  enddo

  return

END SUBROUTINE generate_Vsf_gsWFC

!----------------------------------------------------------------------------------------------------------

SUBROUTINE FNA_ERROR()
  implicit none
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=210) 'IMPLEMENTATION OF FNA FOR OPEN SYSTEM CURRENTLY NOT AVAILABLE'
  write(*,*)
  write(*,fmt=220) 'EXITING CALCULATION ...'
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  stop

  200 format(A57)
  210 format(A61)
  220 format(A23)
END SUBROUTINE FNA_ERROR

!----------------------------------------------------------------------------------------------------------
