SUBROUTINE FNA_LADDER(BoundCond,upDMAT,downDMAT,H,lattVEC,lattPOS,recVEC,A1o,A2o,nSITES,nSITESx,nSITESy, &
                & nEL,nELup,nELdown,nDET,nDETup,nDETdown,t,U,nTWFCDET,tk,Uk,gsVsf,orb_rot,fna_trunc,     &
                &fna_stoq,tot_spin_fna,tot_spin_sq_fna)

!-----------------------------------------------------------------------!
! SUBROUTINE to construct and diagonalize the EFFECTIVE HAMILTONIAN of  !
! the FIXED NODE APPROXIMATION incorporating the SIGN-FLIP POTENTIAL    !
! using the HARTEE-FOCK DET. (U=0.0) as TRIAL WAVE FUNCTION for LADDER  !
! systems with PBC in both one [X|Y] direction.                         !
! The single particle states of the HF TRIAL WFC are mixed PLANE WAVES  !
! of the form:                                                          !
!             <r|k'> = [1/sqrt(N)] exp{i*(k'1 x r1)} sin[(k'2+1)*r2]/N2 !
! with:                                                                 !
!       k'  = (k'1,k'2)                                                 !
!       k'1 = 2*pi/N1 * n1 ; n1 = 0,+-1,...,N1/2 (-N1/2)                !
!       k'2 = 0,1,...,N2                                                !
! (relations are inverted if the PBC are applied to the orthogonal      !
! direction).                                                           !
!-----------------------------------------------------------------------!
 
 implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter    :: dm = 2

  double precision, parameter   :: Vsf_MAX = 1.E+15 ! Threshold MAX value to be asigned to Vsf(R) in case Vsf(R) is a node of the TRIAL (HF) WFC.
  
 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  character(5)    , intent(in)  :: BoundCond

  integer(8)      , intent(in)  :: nSITES, nSITESx, nSITESy
  integer(8)      , intent(in)  :: nEL   , nELup  , nELdown
  integer(8)      , intent(in)  :: nDET  , nDETup , nDETdown

  integer(8)      , intent(in)  :: nTWFCDET
  
  integer(8)      , intent(in)  :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)
 
  double precision, intent(in)  :: A1o
  double precision, intent(in)  :: A2o

  double precision, intent(in)  :: lattVEC(dm,dm)
  double precision, intent(in)  ::  recVEC(dm,dm)
  double precision, intent(in)  :: lattPOS(dm,nSITES)

  double precision, intent(in)  :: H(nDET,nDET)
  double precision, intent(in)  :: gsVsf(nDET)

  double precision, intent(in)  :: t
  double precision, intent(in)  :: U

  double precision, intent(in)  :: tk
  double precision, intent(in)  :: Uk

  character(3)    , intent(in)  :: orb_rot
  character(3)    , intent(in)  :: fna_trunc
  character(3)    , intent(in)  :: fna_stoq
  character(3)    , intent(in)  :: tot_spin_fna
  character(3)    , intent(in)  :: tot_spin_sq_fna

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                    :: mMAXx, mMAXy    !MAX. values of {m1, m2} quantum numbers of K vectors K = m1*B1 + m2*B2
  integer(8)                    :: nKPOINTS        !TOTAL number of K POINTS.
  
  integer(8)                    :: nNODES, nDET_tr !Quantity of NODES and NO-NODES of the TRIAL WFC for the FNA TRUNC. CALC.
  integer(8)      , allocatable :: NODESindx(:)
  integer(8)      , allocatable :: noNODESindx(:)
  
  integer(8)      , allocatable :: KupDMAT(:,:)    !LIST of UP-SPIN   K-DETS.
  integer(8)      , allocatable :: KdownDMAT(:,:)  !LIST of DOWN-SPIN K-DETS.
  
  integer(8)      , allocatable :: rot_KP(:)       !VEC. containing (+) or (-) for the linear combination of KNPOINTS for ORBITAL ROTATION.

  integer(8)      , allocatable :: CIindx(:)       !List of GLOBAL INDICES of the most-contributing DETs[k] to the TWFC CI EXPANSION.

  double precision              :: E_T             !ENERGY of the TRIAL WFC.
  double precision              :: S               !SPIN S.
  double precision              :: tot_spin_sq     !<S^2> FUNCTION.  
  
  double precision, allocatable :: KPOINTS(:,:)    !LIST of K POINTS.
  double precision, allocatable :: KP_ENRG(:)      !LIST of K-POINTS ENERGIES.
  double precision, allocatable :: CI_K(:)         !COEFFICIENTS of CI EXPANSION of GS-WFC in K-SPACE HF-ORBITALS.

  double precision, allocatable :: SGN(:,:)
  double precision, allocatable :: Vsf(:)
  double precision, allocatable :: H_eff(:,:), CI_eff(:,:), E_fn(:)

  double precision, allocatable :: Heff_tr(:,:)
  double precision, allocatable :: CI_rg(:)

  double complex  , allocatable :: TWFC(:)

  INTEGER(8) :: I

 !-------------------------------------------------!
 ! SUBROUTINE: GENERATION OF LIST OF ALL K'-POINTS !
 !-------------------------------------------------!

  if ( BoundCond .eq. 'pbcx' ) then
     if ( mod(nSITESx,2) .ne. 0 ) then
        mMAXx = int( (nSITESx - 1 ) / 2 )
     else
        mMAXx = int( nSITESx / 2 )
     endif
     mMAXy    = nSITESy
     if     ( mod(nSITESx,2).eq.0 ) then
            nKPOINTS = (2*mMAXx) * mMAXy
     elseif ( mod(nSITESx,2).ne.0 ) then
            nKPOINTS = (1 + 2*mMAXx) * mMAXy
     endif
  elseif ( BoundCond .eq. 'pbcy' ) then
      if ( mod(nSITESy,2) .ne. 0 ) then
         mMAXy = int( (nSITESy - 1 ) / 2 )
      else
         mMAXy = int( nSITESy / 2 )
      endif
      mMAXx    = nSITESx
      if     ( mod(nSITESy,2).eq.0 ) then
             nKPOINTS = mMAXx * (2*mMAXy)
      elseif ( mod(nSITESy,2).ne.0 ) then
             nKPOINTS = mMAXx * (1 + 2*mMAXy)
     endif
  endif

  if ( nKPOINTS.ne.nSITES ) then
     write(*,*) 'ERROR IN THE GENERATION OF [FNA] K-POINTS'
     stop
  endif

  allocate( KPOINTS(dm,nKPOINTS) )
  allocate( KP_ENRG(nKPOINTS)    )

  call generate_KPOINTS_LADDER( KPOINTS(:,:),BoundCond,A1o,A2o,nSITESx,nSITESy,mMAXx,mMAXy,nKPOINTS )
  call sort_KPOINTS( KPOINTS(:,:),KP_ENRG(:),A1o,A2o,nKPOINTS,t )

 !-------------------------------------------------!
 ! SUBROUTINE: DIAG. OF K-SPACE HAM. MAT. AND DET. !
 !-------------------------------------------------!

  allocate( KupDMAT(nKPOINTS,nDETup), KdownDMAT(nKPOINTS,nDETdown) )
  allocate( CI_K(nDET)                                             )

  KupDMAT(:,:)   = upDMAT(:,:)
  KdownDMAT(:,:) = downDMAT(:,:)


!  call solve_KSPACE_HUBB(BoundCond,CI_K(:),KupDMAT(:,:),KdownDMAT(:,:),KPOINTS(:,:),KP_ENRG(:),recVEC(:,:), &
!          & nKPOINTS,nEL,nELup,nELdown,nDET,nDETup,nDETdown,tk,Uk)

  call solve_KSPACE_HUBB(BoundCond,CI_K(:),KupDMAT(:,:),KdownDMAT(:,:),KPOINTS(:,:),KP_ENRG(:),recVEC(:,:), &
          & nKPOINTS,nEL,nELup,nELdown,nDET,nDETup,nDETdown,tk,0.d0)

 !-------------------------------------------------!
 ! SUBROUTINE: ROTATION OF KPOINTS ORBITALS        !
 !-------------------------------------------------!

  allocate( rot_KP(nKPOINTS) )

  call rotate_KPOINTS_LADDER( rot_KP(:), KPOINTS(:,:), nKPOINTS, BoundCond )

 !-------------------------------------------------!
 ! SUBROUTINE: GENERATION OF PROJECTED TRIAL WFC   !
 !-------------------------------------------------!

  allocate( TWFC(nDET)       )
  allocate( CIindx(nTWFCDET) )

  call get_PROJ_TWFC_LADDER( BoundCond,TWFC(:),CIindx(:),CI_K(:),KupDMAT(:,:),KdownDMAT(:,:),KPOINTS(:,:),  &
          & rot_KP(:),lattPOS(:,:),upDMAT(:,:),downDMAT(:,:),nKPOINTS,nSITES,nSITESx,nSITESy,nELup,nELdown, &
          & nDET,nDETup,nDETdown,nTWFCDET,orb_rot )

 !--------------------------------------------------!
 ! SUBROUTINE: TRIAL WFC ENERGY CALCULATION         !
 !--------------------------------------------------!

  call get_TWF_EN( E_T,TWFC(:),H(:,:),nDET )

 !--------------------------------------------------!
 ! SUBROUTINE: SIGN-FLIP MATRIX S                   !
 !--------------------------------------------------!

  allocate( SGN(nDET,nDET) )

  call generate_SGN( SGN(:,:),TWFC(:),H(:,:),nDET,Vsf_MAX )

 !--------------------------------------------------!
 ! SUBROUTINE: EFFECTIVE HAMILTONINAN (FNA)         !
 !--------------------------------------------------!

  allocate( H_eff(nDET,nDET), Vsf(nDET) )

  call generate_Heff( H_eff(:,:),Vsf(:),SGN(:,:),TWFC(:),H(:,:),nDET,Vsf_MAX )

  if ( fna_trunc .eq. 'no' ) then

     !--------------------------------------------------!
     ! SUBROUTINE: DIAGONALIZATION OF THE EFF. HAM.     !
     !--------------------------------------------------!

      allocate( CI_eff(nDET,nDET),E_fn(nDET) )

      CI_eff(:,:) = H_eff(:,:)

      call diagonalize_matrix( nDET,CI_eff(:,:),E_fn(:) )

     !-------------------------------------------------!
     ! FNA GROUND STATE TOTAL SPIN S CALCULATION       !
     !-------------------------------------------------!

      if ( tot_spin_fna .eq. 'yes' ) then
         S = tot_spin_sq( CI_eff(:,1),upDMAT(:,:),downDMAT(:,:),nSITES,nDET,nDETup,nDETdown )   !ATTENTION THIS IS <S^2> NOT S.
      endif

     !-------------------------------------------------!
     ! OUTPUT WRITING: FIXED-NODE APP. CALCULATION     !
     !-------------------------------------------------!

      call write_output_fna( upDMAT(:,:),downDMAT(:,:),KupDMAT(:,:),KdownDMAT(:,:),CI_K(:),CIindx(:),  &
              & SGN(:,:),Vsf(:),gsVsf(:),H_eff(:,:),CI_eff(:,:),E_fn(:),nKPOINTS,nSITES,nELup,nELdown, &
              & nDET,nDETup,nDETdown,nTWFCDET,tot_spin_fna,E_T,S )

     !-------------------------------------------------!
     ! STOQ. FNA HAM. MATRIX EXACT DIAGONALIZATION     !
     !-------------------------------------------------!

      if ( fna_stoq .eq. 'yes' ) then
         write(*,*)
         write(*,fmt=200) '---------------------------------------------------------'
         write(*,fmt=200) '----- DIAGONALIZATION OF THE STOQ. FNA-HAM. MATRIX  -----'
         write(*,fmt=200) '---------------------------------------------------------'
         call solve_stoquastized( H_eff(:,:),E_fn(:),nSITES,nDET,nDETup,nDETdown,upDMAT,downDMAT,tot_spin_sq_fna )
      endif

  elseif ( fna_trunc .eq. 'yes' ) then 

     !-------------------------------------------------!
     ! SUBROUTINE: GET NUMBER AND POSITION OF NO-NODES !
     !-------------------------------------------------!
      
      call count_TWFCnodes( nNODES,H_eff(:,:),nDET,Vsf_MAX )
      
      nDET_tr = nDET - nNODES
   
      allocate( NODESindx(nNODES),noNODESindx(nDET_tr),Heff_tr(nDET_tr,nDET_tr), &
              & CI_eff(nDET_tr,nDET_tr),CI_rg(nDET),E_fn(nDET_tr) )

      call get_TWFC_NODESindx( NODESindx(:),noNODESindx(:),H_eff(:,:),nDET,nDET_tr,nNODES,Vsf_MAX )  

     !-------------------------------------------------!
     ! SUBROUTINE: CONSTRUCTION OF THE TRUNC. EFF. HAM.!
     !-------------------------------------------------!

      call get_Heff_tr( Heff_tr(:,:),H_eff(:,:),noNODESindx(:),nDET,nDET_tr)

     !-------------------------------------------------!
     ! SUBROUTINE: DIAGONALIZATION OF THE TR. EFF. HAM.!
     !-------------------------------------------------!

      CI_eff(:,:) = Heff_tr(:,:)

      call diagonalize_matrix( nDET_tr,CI_eff(:,:),E_fn(:) )

     !--------------------------------------------------!
     ! SUBROUTINE: REORGANIZATION OF THE CI COEFF.      !
     !--------------------------------------------------!

      call reorganize_CI( CI_rg(:),CI_eff(:,1),noNODESindx(:),nDET,nDET_tr )

     !-------------------------------------------------!
     ! FNA GROUND STATE TOTAL SPIN S CALCULATION       !
     !-------------------------------------------------!

      if ( tot_spin_fna .eq. 'yes' ) then
         S = tot_spin_sq( CI_rg(:),upDMAT(:,:),downDMAT(:,:),nSITES,nDET,nDETup,nDETdown )   !ATTENTION THIS IS <S^2> NOT S.
      endif

     !-------------------------------------------------!
     ! OUTPUT WRITING: TRUNCATED FIXED-NODE APP. CALC. !
     !-------------------------------------------------!

      call write_output_fna_tr( upDMAT(:,:),downDMAT(:,:),KupDMAT(:,:),KdownDMAT(:,:),CI_K(:),CIindx(:),    &
              & SGN(:,:),Vsf(:),gsVsf(:),Heff_tr(:,:),CI_rg(:),E_fn(:),nKPOINTS,nSITES,nELup,nELdown,nDET, &
              & nDET_tr,nDETup,nDETdown,nTWFCDET,tot_spin_fna,E_T,S )

     !-------------------------------------------------!
     ! STOQ. FNA HAM. MATRIX EXACT DIAGONALIZATION     !
     !-------------------------------------------------!

      if ( fna_stoq .eq. 'yes' ) then
         write(*,*)
         write(*,fmt=200) '---------------------------------------------------------'
         write(*,fmt=200) '--- DIAGONALIZATION OF THE STOQ. TR. FNA-HAM. MATRIX ----'
         write(*,fmt=200) '---------------------------------------------------------'
         call solve_stoquastized_tr( Heff_tr(:,:),E_fn(:),noNODESindx(:),nSITES,nDET,nDET_tr,nDETup, &
                 & nDETdown,upDMAT,downDMAT,tot_spin_sq_fna )
      endif

     deallocate( NODESindx, noNODESindx, Heff_tr, CI_rg )         
  
  endif

  deallocate( KPOINTS, KP_ENRG    )
  deallocate( KupDMAT, KdownDMAT  )
  deallocate( CI_K                )
  deallocate( CIindx              )
  deallocate( rot_KP              )
  deallocate( TWFC                )
  deallocate( SGN, Vsf            )
  deallocate( H_eff, CI_eff, E_fn )

  200 format(A57)

  return

END SUBROUTINE FNA_LADDER

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_KPOINTS_LADDER(KPOINTS,BoundCond,A1o,A2o,nSITESx,nSITESy,mMAXx,mMAXy,nKPOINTS)

!-----------------------------------------------------!
! SUBROUTINE to generate an UNSORTED LIST of KPOINTS  !
! for the LADDER SYSTEM.                              !
!-----------------------------------------------------!

 implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter    :: dm = 2
  double precision, parameter    :: pi = dacos(-1.d0)

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  character(5)    , intent(in)  :: BoundCond

  integer(8)      , intent(in)  :: nSITESx, nSITESy
  integer(8)      , intent(in)  :: mMAXx  , mMAXy
  integer(8)      , intent(in)  :: nKPOINTS

  double precision, intent(in)  :: A1o
  double precision, intent(in)  :: A2o

  double precision, intent(out) :: KPOINTS(dm,nKPOINTS)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                    :: m1, m2
  integer(8)                    :: kpoint, i

 !-------------------------------------------------!
 ! SUBROUTINE: GENERATION OF KPOINTS LIST          !
 !-------------------------------------------------!

  KPOINTS(:,:) = 0.d0

  kpoint       = 1

  if     ( BoundCond .eq. 'pbcx' ) then

         if     ( mod(nSITESx,2).eq.0 ) then

                do m1 = -(mMAXx) + 1,mMAXx,1
                   do m2 = 1,mMAXy,1
                      KPOINTS(1,kpoint) = (2.d0*pi)*( dble(m1)/(dble(nSITESx)*A1o) )
                      KPOINTS(2,kpoint) = pi*( dble(m2)/(dble(nSITESy + 1.d0)*A2o) )
                      kpoint            = kpoint + 1
                   enddo
                enddo
                
         elseif ( mod(nSITESx,2).ne.0 ) then

                do m1 = -(mMAXx),mMAXx,1
                   do m2 = 1,mMAXy,1
                      KPOINTS(1,kpoint) = (2.d0*pi)*( dble(m1)/(dble(nSITESx)*A1o) )
                      KPOINTS(2,kpoint) = pi*( dble(m2)/(dble(nSITESy + 1.d0)*A2o) )
                      kpoint            = kpoint + 1
                   enddo
                enddo
         endif

  elseif ( BoundCond .eq. 'pbcy' ) then

         if     ( mod(nSITESy,2).eq.0 ) then

                do m1 = 1,mMAXx,1
                   do m2 = -(mMAXy) + 1,mMAXy,1
                      KPOINTS(1,kpoint) = pi*( dble(m1)/(dble(nSITESx + 1.d0)*A1o) )
                      KPOINTS(2,kpoint) = (2.d0*pi)*( dble(m2)/(dble(nSITESy)*A2o) )
                      kpoint            = kpoint + 1
                   enddo
                enddo

         elseif ( mod(nSITESy,2).ne.0 ) then

                do m1 = 1,mMAXx,1
                   do m2 = -(mMAXy),mMAXy,1
                      KPOINTS(1,kpoint) = pi*( dble(m1)/(dble(nSITESx + 1.d0)*A1o) )
                      KPOINTS(2,kpoint) = (2.d0*pi)*( dble(m2)/(dble(nSITESy)*A2o) )
                      kpoint            = kpoint + 1
                   enddo
                enddo
         endif
  endif

  return

END SUBROUTINE generate_KPOINTS_LADDER

!----------------------------------------------------------------------------------------------------------

SUBROUTINE rotate_KPOINTS_LADDER(rot_KP,KPOINTS,nKPOINTS,BoundCond)

!--------------------------------------------------------!
! SUBROUTINE to asign the sign (+) or (-) to a given     !
! KPOINT to ROTATE the ORBITALS in Momentum Space as:    !
!                                                        !
!        |K'> = 1/sqrt(2) [ |K> +- |-K> ]                !
!                                                        !
! In the case in which both |K> and |-K> are present in  !
! the Basis Set, the algorithm asign (+) to the first    !
! KPOINT to be occupied, and (-) to the second, in order !
! to preserve the number of Basis Functions and the DIM  !
! of the Hilbert Space.                                  !
!--------------------------------------------------------!

  implicit none

 !--------------------------------------------------!
 ! PARAMETERS                                       !
 !--------------------------------------------------!

  integer(8)       , parameter  :: dm = 2

 !--------------------------------------------------!
 ! INOUT VARIABLES                                  !
 !--------------------------------------------------!

  character(5)    , intent(in)  :: BoundCond

  integer(8)      , intent(in)  :: nKPOINTS
  integer(8)      , intent(out) :: rot_KP(nKPOINTS)

  double precision, intent(in)  :: KPOINTS(dm,nKPOINTS)

 !--------------------------------------------------!
 ! LOCAL VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)                    :: k, kp

 !--------------------------------------------------!
 ! SUBROUTINE: ASSIGNATION OF SIGN FOR ORB. ROT.    !
 !--------------------------------------------------!

  rot_KP(:) = 1


  if ( BoundCond .eq. 'pbcx' ) then
     do k = 2,nKPOINTS,1
        do kp = 2,k-1,1
           if ( (KPOINTS(1,k) .eq. -KPOINTS(1,kp)) .and. (KPOINTS(2,k) .eq. KPOINTS(2,kp)) ) then
              rot_KP(k) = -1
              exit
           endif
        enddo
     enddo
  elseif ( BoundCond .eq. 'pbcy' ) then
     do k = 2,nKPOINTS,1
        do kp = 2,k-1,1
           if ( (KPOINTS(2,k) .eq. -KPOINTS(2,kp)) .and. (KPOINTS(1,k) .eq. KPOINTS(1,kp)) ) then
              rot_KP(k) = -1
              exit
           endif
        enddo
     enddo
  endif

  return

END SUBROUTINE rotate_KPOINTS_LADDER

!----------------------------------------------------------------------------------------------------------


DOUBLE COMPLEX FUNCTION detM_LADDER(BoundCond,R,kDET,rot_kDET,nEL,nSITES,nSITESx,nSITESy,orb_rot)

!-----------------------------------------------------------------------------!
! FUNCTION that calculates the PROJECTION of the TRIAL WAVE FUNCTION,         !
! the HF DET. (U = 0.0) onto a given spatial configuration |R> of the         !
! electrons in the lattice, for a given spin channel (up- or down-electrons). !
! It involves:                                                                !
!             |R> = |Rup> (x) |Rdown>                                         !
!             |K> = |Kup> (x) |Kdown>                                         !
! and:                                                                        !
!           <R|HF> = <Rup|Kup> <Rdown|Kdown>                                  !
!                  = det[M]                                                   !
!                  = det[upM] x det[downM]                                    !
! where M is matrix composed by the following elements:                       !
!           M(i,j) = <r_i|k_j> = mPW (mixed plane wave)                       !
! and is block-diagonal according to spin symmetry.                           !
!-----------------------------------------------------------------------------!

 implicit none

 !--------------------------------------------------!
 ! PARAMETERS                                       !
 !--------------------------------------------------!

 integer(8)      , parameter  :: dm = 2

 double precision, parameter  :: pi = dacos(-1.d0)  !PI PARAMETER.

 double complex  , parameter  :: ic = (0.d0,1.d0)   !IMAGINARY UNIT.

 !--------------------------------------------------!
 ! INOUT VARIABLES                                  !
 !--------------------------------------------------!

 character(5)    , intent(in) :: BoundCond
 character(3)    , intent(in) :: orb_rot

 integer(8)      , intent(in) :: nSITES
 integer(8)      , intent(in) :: nSITESx
 integer(8)      , intent(in) :: nSITESy
 integer(8)      , intent(in) :: nEL
 integer(8)      , intent(in) :: rot_kDET(nEL)

 double precision, intent(in) :: R(dm,nEL)
 double precision, intent(in) :: kDET(dm,nEL)

 !--------------------------------------------------!
 ! LOCAL VARIABLES                                  !
 !--------------------------------------------------!

  double complex              :: M(nEL,nEL)
  double complex, allocatable :: LU(:,:)
  double complex              :: det

  double precision            :: normCNST

  integer(8), allocatable     :: IPIV(:)
  integer(8)                  :: i, j

 !--------------------------------------------------!
 ! SUBROUTINE: GENERATION OF MATRIX M[sigma]        !
 !--------------------------------------------------!

  allocate( LU(nEL,nEL),IPIV(nEL) )

  M(:,:)   = 0.d0

  normCNST = 1.d0

  if ( orb_rot .eq. 'no' ) then
     if     ( BoundCond .eq. 'pbcx' ) then
            normCNST = 1.d0/dsqrt(dble(nSITESx)) * dsqrt(2.d0/dble(nSITESy + 1.d0))
            do i = 1,nEL,1
               do j = 1,nEL,1
                  M(i,j) = normCNST * exp(ic*R(1,i)*kDET(1,j)) * sin( kDET(2,j) * (R(2,i) + 1.d0) )
               enddo
            enddo
     elseif ( BoundCond .eq. 'pbcy' ) then
            normCNST = 1.d0/dsqrt(dble(nSITESy)) * dsqrt(2.d0/dble(nSITESx + 1.d0))
            do i = 1,nEL,1
               do j = 1,nEL,1
                  M(i,j) = normCNST * exp(ic*R(2,i)*kDET(2,j)) * sin( kDET(1,j) * (R(1,i) + 1.d0) )
               enddo
            enddo
     endif
  elseif ( orb_rot .eq. 'yes' ) then
     if     ( BoundCond .eq. 'pbcx' ) then
            normCNST = 1.d0/dsqrt(dble(nSITESx)) * dsqrt(2.d0/dble(nSITESy + 1.d0)) * 1.d0/dsqrt(2.d0)
            do i = 1,nEL,1
               do j = 1,nEL,1
                  if ( rot_kDET(j) .eq. 1 ) then
                     M(i,j) = normCNST * 2.d0 * dcos(R(1,i)*kDET(1,j)) * sin( kDET(2,j) * (R(2,i) + 1.d0) )
                  elseif ( rot_kDET(j) .eq. -1 ) then
                     M(i,j) = normCNST * 2.d0 * dsin(R(1,i)*kDET(1,j)) * sin( kDET(2,j) * (R(2,i) + 1.d0) )
                  endif
               enddo
            enddo
     elseif ( BoundCond .eq. 'pbcy' ) then
            normCNST = 1.d0/dsqrt(dble(nSITESy)) * dsqrt(2.d0/dble(nSITESx + 1.d0)) * 1.d0/dsqrt(2.d0)
            do i = 1,nEL,1
                do j = 1,nEL,1
                   if ( rot_kDET(j) .eq. 1 ) then
                      M(i,j) = normCNST * 2.d0 * dcos(R(2,i)*kDET(2,j)) * sin( kDET(1,j) * (R(1,i) + 1.d0) )
                   elseif ( rot_kDET(j) .eq. -1 ) then
                      M(i,j) = normCNST * 2.d0 * dsin(R(2,i)*kDET(2,j)) * sin( kDET(1,j) * (R(1,i) + 1.d0) )
                   endif
               enddo
            enddo
     endif
  endif

  call eval_complex_determinant( nEL,nEL,M(:,:),LU(:,:),IPIV(:),det )

  detM_LADDER = det

  deallocate( LU, IPIV )

  return

END FUNCTION detM_LADDER

!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_PROJ_TWFC_LADDER( BoundCond,TWFC,CIindx,CI_K,KupDMAT,KdownDMAT,KPOINTS,rot_KP, &
                & lattPOS,upDMAT,downDMAT,nKPOINTS,nSITES,nSITESx,nSITESy,nELup,nELdown,nDET, &
                & nDETup,nDETdown,nTWFCDET,orb_rot )

!--------------------------------------------------------------------------------------------!
! SUBROUTINE that calculates the ONE-DIMENSIONAL VECT. of the TRIAL WFC projected on the     !
! REAL SPACE DET. as:                                                                        !
!                    TWFC[R] =  <R|TWFC> for all R                                           !
! with |R> = |I>                                                                             !
!                    TWFC[I] =  <I|TWFC> for all I = 1,...,nDET                              !
! The TWFC follows a short CI expansion (MC-TWFC)                                            !
!                    TWFC    = SUM_{D=1,...,nTWFCDET} CI(D) |D[K]>                           !
! such that the PROJ. is calculated as:                                                      !
!                    TWFC[I] =  SUM_{D=1,...,nTWFCDET} CI(D) <I|D[K]> for all I = 1,...,nDET !
!--------------------------------------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8), parameter :: dm = 2

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  character(5)    , intent(in)  :: BoundCond
  character(3)    , intent(in)  :: orb_rot

  integer(8)      , intent(in)  :: nKPOINTS
  integer(8)      , intent(in)  :: nSITES
  integer(8)      , intent(in)  :: nSITESx
  integer(8)      , intent(in)  :: nSITESy
  integer(8)      , intent(in)  :: nELup , nELdown
  integer(8)      , intent(in)  :: nDET  , nDETup , nDETdown
  integer(8)      , intent(in)  :: nTWFCDET

  integer(8)      , intent(in)  :: upDMAT(nSITES,nDETup)
  integer(8)      , intent(in)  :: downDMAT(nSITES,nDETdown)
  integer(8)      , intent(in)  :: KupDMAT(nKPOINTS,nDETup)
  integer(8)      , intent(in)  :: KdownDMAT(nKPOINTS,nDETdown)

  integer(8)      , intent(in)  :: rot_KP(nKPOINTS)

  integer(8)      , intent(out) :: CIindx(nTWFCDET)

  double precision, intent(in)  :: CI_K(nDET)
  double precision, intent(in)  :: KPOINTS(dm,nKPOINTS)
  double precision, intent(in)  :: lattPOS(dm,nSITES)

  double complex  , intent(out) :: TWFC(nDET)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)      , allocatable :: rot_Kup(:), rot_Kdown(:)

  double precision, allocatable :: Rup(:,:)  , Rdown(:,:)
  double precision, allocatable :: Kup(:,:)  , Kdown(:,:)

  double complex                :: detM
  double complex                :: detM_LADDER
  double complex                :: detMup  , detMdown

  integer(8)                    :: kDET
  integer(8)                    :: I , Iup , Idown
  integer(8)                    :: Ik, Ikup, Ikdown


 !-------------------------------------------------!
 ! ALLOCATION OF VARIABLES                         !
 !-------------------------------------------------!

  allocate( rot_Kup(nELup), rot_Kdown(nELdown) )
  allocate( Kup(dm,nELup) , Kdown(dm,nELdown)  )
  allocate( Rup(dm,nELup) , Rdown(dm,nELdown)  )

 !-------------------------------------------------!
 ! SUBROUTINE: GET THE FIRST nTWFCDET TERMS OF CI  !
 !-------------------------------------------------!

  call get_TWFC_CIEXPANSION( CIindx(:),CI_K(:),nDET,nTWFCDET )

 !-------------------------------------------------!
 ! SUBROUTINE: CALCULATE THE PROJECTED TRIAL WFC   !
 !-------------------------------------------------!

  TWFC = (0.d0,0.d0)

  if ( (nELup .ne. 0) .and. (nELdown .ne. 0) ) then

     do I = 1,nDET,1
        call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
        call get_DETPOS( Rup(:,:)  ,lattPOS(:,:),upDMAT(:,Iup)    ,nSITES,nDETup  ,nELup   )
        call get_DETPOS( Rdown(:,:),lattPOS(:,:),downDMAT(:,Idown),nSITES,nDETdown,nELdown )
        do kDET = 1,nTWFCDET,1
           Ik       = CIindx(kDET)
           call spinDETindex( Ik,Ikup,Ikdown,nDETup,nDETdown )
           call generate_kDET( Kup(:,:)  ,rot_Kup(:)  ,KupDMAT(:,Ikup)    ,KPOINTS(:,:),rot_KP(:),nELup  ,nKPOINTS)
           call generate_kDET( Kdown(:,:),rot_Kdown(:),KdownDMAT(:,Ikdown),KPOINTS(:,:),rot_KP(:),nELdown,nKPOINTS)
           detMup   = detM_LADDER( BoundCond,Rup(:,:)  ,Kup(:,:)  ,rot_Kup(:)  ,nELup  ,nSITES,nSITESx,nSITESy,orb_rot )
           detMdown = detM_LADDER( BoundCond,Rdown(:,:),Kdown(:,:),rot_Kdown(:),nELdown,nSITES,nSITESx,nSITESy,orb_rot )
           detM     = detMup * detMdown
           TWFC(I)  = TWFC(I) + CI_K(Ik) * detM
        enddo
     enddo

  elseif ( (nELup .ne. 0) .and. (nELdown .eq. 0) ) then

     do I = 1,nDET,1
        call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
        call get_DETPOS( Rup(:,:),lattPOS(:,:),upDMAT(:,Iup),nSITES,nDETup,nELup )
        do kDET = 1,nTWFCDET,1
           Ik       = CIindx(kDET)
           call spinDETindex( Ik,Ikup,Ikdown,nDETup,nDETdown )
           call generate_kDET( Kup(:,:),rot_Kup(:),KupDMAT(:,Ikup),KPOINTS(:,:),rot_KP(:),nELup,nKPOINTS)
           detMup   = detM_LADDER( BoundCond,Rup(:,:),Kup(:,:),rot_Kup(:),nELup,nSITES,nSITESx,nSITESy,orb_rot )
           detM     = detMup
           TWFC(I)  = TWFC(I) + CI_K(Ik) * detM
        enddo
     enddo

  elseif ( (nELup .eq. 0) .and. (nELdown .ne. 0) ) then

     do I = 1,nDET,1
        call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
        call get_DETPOS( Rdown(:,:),lattPOS(:,:),downDMAT(:,Idown),nSITES,nDETdown,nELdown )
        do kDET = 1,nTWFCDET,1
           Ik       = CIindx(kDET)
           call spinDETindex( Ik,Ikup,Ikdown,nDETup,nDETdown )
           call generate_kDET( Kdown(:,:),rot_Kdown(:),KdownDMAT(:,Ikdown),KPOINTS(:,:),rot_KP(:),nELdown,nKPOINTS)
           detMdown = detM_LADDER( BoundCond,Rdown(:,:),Kdown(:,:),rot_Kdown(:),nELdown,nSITES,nSITESx,nSITESy,orb_rot )
           detM     = detMdown
           TWFC(I)  = TWFC(I) + CI_K(Ik) * detM
        enddo
     enddo

  endif

 !-------------------------------------------------!
 ! DEALLOCATION OF VARIABLES                       !
 !-------------------------------------------------!

  deallocate( rot_Kup, rot_Kdown )
  deallocate( Kup    , Kdown     )
  deallocate( Rup    , Rdown     )

  return

END SUBROUTINE get_PROJ_TWFC_LADDER

!----------------------------------------------------------------------------------------------------------
