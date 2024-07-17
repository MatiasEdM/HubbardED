SUBROUTINE FNA_BULK(BoundCond,upDMAT,downDMAT,H,CI,lattVEC,lattPOS,recVEC,A1o,A2o,nSITES,nSITESx,nSITESy, &
                & nEL,nELup,nELdown,nDET,nDETup,nDETdown,t,U,tk,Uk,nTWFCDET,nPDET,gsVsf,fna_trunc,        &
                & subspace_coeff,fna_relax,fna_stoq,tot_spin_fna,tot_spin_sq_fna)

!----------------------------------------------------------------------!
! SUBROUTINE to construct and diagonalize the EFFECTIVE HAMILTONIAN of !
! the FIXED NODE APPROXIMATION incorporating the SIGN-FLIP POTENTIAL   !
! using the HARTEE-FOCK DET. (U=0.0) as TRIAL WAVE FUNCTION for BULK   !
! systems with PBC in both [XY] directions.                            !
! The SINGLE-PARTICLE STATES of the HF TRIAL WFC are PLANE WAVES of    !
! the form:                                                            !
!            <r|k> = [1/sqrt(N)] exp{i*(k x r)}                        !
!----------------------------------------------------------------------!
 
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

  integer(8)      , intent(in)  :: nPDET
  
  integer(8)      , intent(in)  :: upDMAT(nSITES,nDETup)
  integer(8)      , intent(in)  :: downDMAT(nSITES,nDETdown)
 
  double precision, intent(in)  :: A1o
  double precision, intent(in)  :: A2o

  double precision, intent(in)  :: lattVEC(dm,dm), recVEC(dm,dm)
  double precision, intent(in)  :: lattPOS(dm,nSITES)

  double precision, intent(in)  :: H(nDET,nDET)
  double precision, intent(in)  :: CI(nDET)
  double precision, intent(in)  :: gsVsf(nDET)

  double precision, intent(in)  :: t
  double precision, intent(in)  :: U

  double precision, intent(in)  :: tk
  double precision, intent(in)  :: Uk

  character(3)    , intent(in)  :: fna_trunc
  character(3)    , intent(in)  :: subspace_coeff
  character(3)    , intent(in)  :: fna_relax
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

  integer(8)      , allocatable :: class_KP(:)     !VEC. containing the CLASS of KPOINT: (0) [High Symmetry Point GAMMA,X,M], (1) [Edge of BZ] or (2) [Internal Point].

  integer(8)      , allocatable :: CIindx(:)       !List of GLOBAL INDICES of the most-contributing DETs[k] to the TWFC CI EXPANSION.

  double precision              :: E_T             !ENERGY of the TRIAL WFC.
  double precision              :: S               !SPIN S.
  double precision              :: tot_spin_sq     !<S^2> FUNCTION.  

  double precision              :: normDC
  
  double precision, allocatable :: KPOINTS(:,:)    !LIST of K POINTS.
  double precision, allocatable :: KP_ENRG(:)      !LIST of K-POINT ENERGIES.
  double precision, allocatable :: H_K(:,:)        !HUBBARD HAMILTONIAN in K-SPACE.
  double precision, allocatable :: CI_K(:)         !COEFFICIENTS of CI EXPANSION of GS-WFC in K-SPACE HF-ORBITALS.
  double precision, allocatable :: CI_TWFC(:)      !COEFFICIENTS of CI EXPANSION of GS-WFC in K-SPACE HF-ORBITALS for the TWFC.

  double complex  , allocatable :: OVERLAP(:,:)    !MATRIX of single-particle orbitals overlap M[i,j] = <R[i] | K[j]> with ORBITAL ROTATION.

  double precision, allocatable :: SGN(:,:)
  double precision, allocatable :: Vsf(:)
  double precision, allocatable :: H_eff(:,:), CI_eff(:,:), E_fn(:)

  double precision, allocatable :: Heff_tr(:,:)
  double precision, allocatable :: CI_rg(:)

  double precision, allocatable :: TWFC(:)

 !-------------------------------------------------!
 ! SUBROUTINE: GENERATION OF LIST OF ALL K-POINTS  !
 !-------------------------------------------------!

  if ( mod(nSITESx,2) .ne. 0 ) then
     mMAXx = int( (nSITESx - 1 ) / 2 )
  else
     mMAXx = int( nSITESx / 2 )
  endif

  if ( mod(nSITESy,2) .ne. 0 ) then
     mMAXy = int( (nSITESy - 1 ) / 2 )
  else
     mMAXy = int( nSITESy / 2 )
  endif

  if     ( ( mod(nSITESx,2).eq.0 ).and.( mod(nSITESy,2).eq.0 ) ) then
         nKPOINTS = (2*mMAXx) * (2*mMAXy)
  elseif ( ( mod(nSITESx,2).eq.0 ).and.( mod(nSITESy,2).ne.0 ) ) then
         nKPOINTS = (2*mMAXx) * (1 + 2*mMAXy)
  elseif ( ( mod(nSITESx,2).ne.0 ).and.( mod(nSITESy,2).eq.0 ) ) then
         nKPOINTS = (1 + 2*mMAXx) * (2*mMAXy)
  elseif ( ( mod(nSITESx,2).ne.0 ).and.( mod(nSITESy,2).ne.0 ) ) then
         nKPOINTS = (1 + 2*mMAXx) * (1 + 2*mMAXy)
  endif

  if ( nKPOINTS.ne.nSITES ) then
     write(*,*) 'ERROR IN THE GENERATION OF [FNA] K-POINTS'
     stop
  endif

  allocate( KPOINTS(dm,nKPOINTS) )
  allocate( KP_ENRG(nKPOINTS)    )

  call generate_KPOINTS( KPOINTS(:,:),recVEC(:,:),nSITESx,nSITESy,mMAXx,mMAXy,nKPOINTS )
  call sort_KPOINTS(     KPOINTS(:,:),KP_ENRG(:),A1o,A2o,nKPOINTS,t )

 !-------------------------------------------------!
 ! SUBROUTINE: DIAG. OF K-SPACE HAM. MAT. AND DET. !
 !-------------------------------------------------!

  allocate( KupDMAT(nKPOINTS,nDETup), KdownDMAT(nKPOINTS,nDETdown) )
  allocate( H_K(nDET,nDET)                                         )
  allocate( CI_K(nDET)                                             )

  KupDMAT(:,:)   = upDMAT(:,:)
  KdownDMAT(:,:) = downDMAT(:,:)

  call solve_KSPACE_HUBB(BoundCond,H_K(:,:),CI_K(:),KupDMAT(:,:),KdownDMAT(:,:),KPOINTS(:,:),KP_ENRG(:),recVEC(:,:), &
          & nKPOINTS,nEL,nELup,nELdown,nDET,nDETup,nDETdown,tk,Uk)

 !-------------------------------------------------!
 ! SUBROUTINE: CLASSIFICATION OF KPOINTS ORBITALS  !
 !-------------------------------------------------!

  allocate( class_KP(nKPOINTS) )

  call classify_KPOINTS_BULK( class_KP(:), KPOINTS(:,:), A1o, A2o, nKPOINTS )

 !-------------------------------------------------!
 ! SUBROUTINE: CALCULATION OF ORBITALS OVERLAP     !
 !-------------------------------------------------!

  allocate( OVERLAP(nSITES,nKPOINTS) )

  call get_ORB_OVERLAP_BULK( OVERLAP(:,:),lattPOS(:,:),KPOINTS(:,:),KP_ENRG(:),class_KP(:),A1o,A2o,nSITES,nKPOINTS )

 !-------------------------------------------------!
 ! SUBROUTINE: GET THE FIRST nTWFCDET TERMS OF CI  !
 !-------------------------------------------------!

  allocate( CIindx(nTWFCDET)  )
  allocate( CI_TWFC(nDET) )
  
  call get_TWFC_CIEXPANSION( CIindx(:),CI_K(:),nDET,nTWFCDET )

  if ( subspace_coeff.eq.'no' ) then
     CI_TWFC(:) = CI_K(:)
  elseif ( subspace_coeff.eq.'yes' ) then
     call get_SUBSPACE_COEFF( CI_TWFC(:),CI_K(:),CIindx(:),H_K(:,:),nDET,nTWFCDET )
  endif

 !-------------------------------------------------!
 ! SUBROUTINE: GENERATION OF PROJECTED TRIAL WFC   !
 !-------------------------------------------------!

  allocate( TWFC(nDET) )

  call get_PROJ_TWFC_BULK( TWFC(:),CIindx(:),CI_TWFC(:),OVERLAP(:,:),KupDMAT(:,:),KdownDMAT(:,:),  &
          & upDMAT(:,:),downDMAT(:,:),nKPOINTS,nSITES,nELup,nELdown,nDET,nDETup,nDETdown,nTWFCDET )

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

  if ( fna_relax .eq. 'no' ) then
     call generate_Heff( H_eff(:,:),Vsf(:),SGN(:,:),H(:,:),nDET,Vsf_MAX )
  elseif ( fna_relax .eq. 'yes' ) then
     call generate_Hpq( H_eff(:,:),Vsf(:),SGN(:,:),H(:,:),CI(:),nDET,nPDET,Vsf_MAX  )
  endif

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
     ! SING-PROBLEM ESTIMATION: BY L2 NORM ||DC||2     !
     !-------------------------------------------------!

      call estimate_sign_problem( normDC,H_eff(:,:),CI_eff(:,1),nDET )

     !-------------------------------------------------!
     ! OUTPUT WRITING: FIXED-NODE APP. CALCULATION     !
     !-------------------------------------------------!

      call write_output_fna( upDMAT(:,:),downDMAT(:,:),KupDMAT(:,:),KdownDMAT(:,:),CI_TWFC(:),CIindx(:),  &
              & SGN(:,:),Vsf(:),gsVsf(:),H_eff(:,:),CI_eff(:,:),E_fn(:),nKPOINTS,nSITES,nELup,nELdown,    &
              & nDET,nDETup,nDETdown,nTWFCDET,nPDET,subspace_coeff,fna_relax,tot_spin_fna,E_T,S,normDC )

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
     ! SING-PROBLEM ESTIMATION: BY L2 NORM ||DC||2     !
     !-------------------------------------------------!

      call estimate_sign_problem( normDC,Heff_tr(:,:),CI_eff(:,1),nDET_tr )

     !-------------------------------------------------!
     ! OUTPUT WRITING: TRUNCATED FIXED-NODE APP. CALC. !
     !-------------------------------------------------!

      call write_output_fna_tr( upDMAT(:,:),downDMAT(:,:),KupDMAT(:,:),KdownDMAT(:,:),CI_TWFC(:),CIindx(:), &
               & SGN(:,:),Vsf(:),gsVsf(:),Heff_tr(:,:),CI_rg(:),E_fn(:),nKPOINTS,nSITES,nELup,nELdown,nDET, &
               & nDET_tr,nDETup,nDETdown,nTWFCDET,nPDET,subspace_coeff,fna_relax,tot_spin_fna,E_T,S,normDC )

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
  deallocate( H_K                 )
  deallocate( CI_K                )
  deallocate( CIindx              )
  deallocate( CI_TWFC             )
  deallocate( class_KP            )
  deallocate( OVERLAP             )
  deallocate( TWFC                )
  deallocate( SGN, Vsf            )
  deallocate( H_eff, CI_eff, E_fn )

  200 format(A57)

  return

END SUBROUTINE FNA_BULK

!----------------------------------------------------------------------------------------------------------

SUBROUTINE rotate_KPOINTS_BULK(rot_KP,KPOINTS,nKPOINTS)

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

  do k = 2,nKPOINTS,1
     do kp = 2,k-1,1
        if ( (KPOINTS(1,k) .eq. -KPOINTS(1,kp)) .and. (KPOINTS(2,k) .eq. -KPOINTS(2,kp)) ) then
           rot_KP(k) = -1
           exit
        endif
     enddo
  enddo

  return

END SUBROUTINE rotate_KPOINTS_BULK

!----------------------------------------------------------------------------------------------------------

SUBROUTINE classify_KPOINTS_BULK(class_KP,KPOINTS,A1o,A2o,nKPOINTS)

!--------------------------------------------------------!
! SUBROUTINE to classify the set of K-POINTS included in !
! the computational basis according to the following     !
! criteria:                                              !
!          (0) High Symmetry Point {GAMMA,M,X,Y}         !
!          (1) Edge of the BZ (but not of high symmetry) !
!          (2) Internal K-POINT                          !
!--------------------------------------------------------!

  implicit none

 !--------------------------------------------------!
 ! PARAMETERS                                       !
 !--------------------------------------------------!

  integer(8)       , parameter  :: dm = 2

  double precision , parameter  :: pi = dacos(-1.d0)

 !---------------------------------------------------!
 ! INOUT VARIABLES                                   !
 !---------------------------------------------------!

  integer(8)      , intent(in)  :: nKPOINTS
  integer(8)      , intent(out) :: class_KP(nKPOINTS)

  double precision, intent(in)  :: A1o
  double precision, intent(in)  :: A2o
  double precision, intent(in)  :: KPOINTS(dm,nKPOINTS)

 !---------------------------------------------------!
 ! LOCAL VARIABLES                                   !
 !---------------------------------------------------!

  double precision, allocatable :: G(:) 
  double precision, allocatable :: M(:)
  double precision, allocatable :: X(:)
  double precision, allocatable :: Y(:)

  integer(8)                    :: k

 !---------------------------------------------------!
 ! SUBROUTINE: HIGH SYMMETRY POINTS 2D RECT. LATTICE !
 !---------------------------------------------------!

  allocate( G(dm) )
  allocate( M(dm) )
  allocate( X(dm) )
  allocate( Y(dm) )

  G(1) = 0.d0   ; G(2) = 0.d0
  M(1) = pi/A1o ; M(2) = pi/A2o
  X(1) = pi/A1o ; X(2) = 0.d0
  Y(1) = 0.d0   ; Y(2) = pi/A2o

 !---------------------------------------------------!
 ! SUBROUTINE: ASSIGNATION OF CLASS TO K-POINTS      !
 !---------------------------------------------------!

  !(0) Check for HIGH-SYMMETRY K-POINT.
  !(1) Check for BZ EDGE K-POINT.
  !(2) Check for INTERNAL K-POIN.

  do k = 1,nKPOINTS,1
     if ( ((KPOINTS(1,k).eq.G(1)).and.(KPOINTS(2,k).eq.G(2))).or.((KPOINTS(1,k).eq.M(1)).and.(KPOINTS(2,k).eq.M(2))).or. &
             & ((KPOINTS(1,k).eq.X(1)).and.(KPOINTS(2,k).eq.X(2))).or.((KPOINTS(1,k).eq.Y(1)).and.(KPOINTS(2,k).eq.Y(2))) ) then
        class_KP(k) = 0
     elseif ( (KPOINTS(1,k).eq.M(1)).or.(KPOINTS(2,k).eq.M(2)) ) then
        class_KP(k) = 1
     else
        class_KP(k) = 2 
     endif
  enddo

  deallocate( G )
  deallocate( M )
  deallocate( X )
  deallocate( Y )

  return

END SUBROUTINE classify_KPOINTS_BULK

!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_ORB_OVERLAP_BULK( OVERLAP,lattPOS,KPOINTS,KP_ENRG,class_KP,A1o,A2o,nSITES,nKPOINTS )

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

 !-------------------------------------------------------!
 ! PARAMETERS                                            !
 !-------------------------------------------------------!

  integer(8)      , parameter  :: dm   = 2

  double precision, parameter  :: pi   = dacos(-1.d0)
  double precision, parameter  :: thrs = 1.E-5

  double complex  , parameter  :: ic   = (0.d0,1.d0)   !IMAGINARY UNIT.

 !-------------------------------------------------------!
 ! INPUT VARIABLES                                       !
 !-------------------------------------------------------!

  integer(8)      , intent(in)  :: nSITES
  integer(8)      , intent(in)  :: nKPOINTS

  integer(8)      , intent(in)  :: class_KP(nKPOINTS)

  double precision, intent(in)  :: A1o
  double precision, intent(in)  :: A2o

  double precision, intent(in)  :: lattPOS(dm,nSITES)
  double precision, intent(in)  :: KPOINTS(dm,nKPOINTS)
  double precision, intent(in)  :: KP_ENRG(nKPOINTS)

  double complex  , intent(out) :: OVERLAP(nSITES,nKPOINTS)

 !-------------------------------------------------------!
 ! LOCAL VARIABLES                                       !
 !-------------------------------------------------------!

  integer(8)                  :: site
  integer(8)                  :: kpoint
  integer(8)                  :: K
  integer(8)                  :: Kp

  double precision            :: normCNST

  double complex, allocatable :: M(:,:)

!  double complex, allocatable :: pLC(:)
!  double complex, allocatable :: mLC(:)
!  logical       , allocatable :: countVEC(:)

 !-------------------------------------------------------!
 ! SUBROUTINE: GENERATION OF MATRIX M[i,j] = <R[i]|K[j]> !
 !-------------------------------------------------------!

  allocate( M(nSITES,nKPOINTS) )

  normCNST = 1.d0/dsqrt(dble(nSITES))

  do kpoint = 1,nKPOINTS,1
     do site = 1,nSITES,1
        M(site,kpoint) = normCNST * exp( ic * dot_product(lattPOS(:,site),KPOINTS(:,kpoint)) )
     enddo
  enddo

  OVERLAP(:,:) = M(:,:) 

 !-------------------------------------------------------!
 ! SUBROUTINE: ORBITAL ROTATION ACCORDING TO KP CLASS.   !
 !-------------------------------------------------------!

!  allocate( pLC(nSITES)        )
!  allocate( mLC(nSITES)        )
!  allocate( countVEC(nKPOINTS) )

!  OVERLAP(:,:) = 0.d0

!  countVEC(:)  = .FALSE.

!  do K = 1,nKPOINTS,1
!     if ( countVEC(K) .eqv. .FALSE. ) then
!        selectcase (class_KP(K))
!        case(0)
!               OVERLAP(:,K) = dble( M(:,K) )
!               countVEC(K)  = .TRUE.
!        case(1)
!              if ( KPOINTS(1,K).eq.(pi/A1o) ) then
!                 do Kp = K,nKPOINTS,1
!                    if ( (KPOINTS(2,Kp).eq.(-KPOINTS(2,K))).and.(KP_ENRG(Kp).eq.KP_ENRG(K)) ) then
!                       exit
!                    endif
!                 enddo
!              elseif (KPOINTS(2,K).eq.(pi/A2o) ) then
!                 do Kp = K,nKPOINTS,1
!                    if ( (KPOINTS(1,Kp).eq.(-KPOINTS(1,K))).and.(KP_ENRG(Kp).eq.KP_ENRG(K)) ) then
!                       exit
!                    endif
!                 enddo
!              endif
!              pLC(:) = (1/dsqrt(2.d0)) * ( M(:,K) + M(:,Kp) )
!              mLC(:) = (1/dsqrt(2.d0)) * ( M(:,K) - M(:,Kp) )
              !pLC(:) = ( M(:,K) + M(:,Kp) )
              !mLC(:) = ( M(:,K) - M(:,Kp) )
!              do site = 1,nSITES,1
!                 if ( abs(aimag(pLC(site))).lt.thrs ) then
!                    OVERLAP(site,K)  = dble(pLC(site))
!                 elseif ( abs(dble(pLC(site))).lt.thrs ) then
!                    OVERLAP(site,K)  = aimag(pLC(site))
!                 endif
!                 if ( abs(aimag(mLC(site))).lt.thrs ) then
!                    OVERLAP(site,Kp) = dble(mLC(site))
!                 elseif ( abs(dble(mLC(site))).lt.thrs ) then
!                    OVERLAP(site,Kp) = aimag(mLC(site))
!                 endif
!              enddo
!              OVERLAP(:,K ) = dble(pLC(:))
!              OVERLAP(:,Kp) = aimag(mLC(:))
!              countVEC(K)   = .TRUE.
!              countVEC(Kp)  = .TRUE.
!        case(2)
!              do Kp = K,nKPOINTS,1
!                 if ( (KPOINTS(1,Kp).eq.(-KPOINTS(1,K))).and.(KPOINTS(2,Kp).eq.(-KPOINTS(2,K)))  ) then
!                    exit
!                 endif
!              enddo
!              pLC(:) = (1/dsqrt(2.d0)) * ( M(:,K) + M(:,Kp) )
!              mLC(:) = (1/dsqrt(2.d0)) * ( M(:,K) - M(:,Kp) )
              !pLC(:) = ( M(:,K) + M(:,Kp) )
              !mLC(:) = ( M(:,K) - M(:,Kp) )
!              do site = 1,nSITES,1
!                 if ( abs(aimag(pLC(site))).lt.thrs ) then
!                    OVERLAP(site,K)  = dble(pLC(site))
!                 elseif ( abs(dble(pLC(site))).lt.thrs ) then
!                    OVERLAP(site,K)  = aimag(pLC(site))
!                 endif
!                 if ( abs(aimag(mLC(site))).lt.thrs ) then
!                    OVERLAP(site,Kp) = dble(mLC(site))
!                 elseif ( abs(dble(mLC(site))).lt.thrs ) then
!                    OVERLAP(site,Kp) = aimag(mLC(site))
!                 endif
!              enddo
!              OVERLAP(:,K ) = dble(pLC(:))
!              OVERLAP(:,Kp) = aimag(mLC(:))
!              countVEC(K)   = .TRUE.
!              countVEC(Kp)  = .TRUE.
!        endselect
!      endif
!  enddo

  deallocate( M        )
!  deallocate( pLC      )
!  deallocate( mLC      )
!  deallocate( countVEC )

  return

END SUBROUTINE get_ORB_OVERLAP_BULK

!----------------------------------------------------------------------------------------------------------

DOUBLE COMPLEX FUNCTION detM_BULK( OVERLAP,rDET,kDET,nEL,nSITES,nKPOINTS )

!-----------------------------------------------------------------------------!
! FUNCTION that calculates the PROJECTION of the TRIAL WAVE FUNCTION,         !
! the DET. in K-SPACE  onto a given spatial configuration |R> of the          !
! electrons in the lattice, for a given spin channel (up- or down-electrons). !
! It involves:                                                                !
!             |R>  = |Rup> (x) |Rdown>                                        !
!             |K>  = |Kup> (x) |Kdown>                                        !
! and:                                                                        !
!           <R|HF> = <Rup|Kup> <Rdown|Kdown>                                  !
!                  = det[M]                                                   !
!                  = det[upM] x det[downM]                                    !
! where M is matrix composed by the following elements:                       !
!           M(i,j) = <r_i|k_j> = PW (plane wave)                              !
! and is block-diagonal according to spin symmetry.                           !
!-----------------------------------------------------------------------------!

 implicit none
 
 !--------------------------------------------------!
 ! INOUT VARIABLES                                  !
 !--------------------------------------------------!

 integer(8)      , intent(in) :: nSITES
 integer(8)      , intent(in) :: nKPOINTS
 integer(8)      , intent(in) :: nEL

 integer(8)      , intent(in) :: rDET(nSITES)
 integer(8)      , intent(in) :: kDET(nKPOINTS)

 double complex  , intent(in) :: OVERLAP(nSITES,nKPOINTS)

 !--------------------------------------------------!
 ! LOCAL VARIABLES                                  !
 !--------------------------------------------------!

  double complex  , allocatable :: M(:,:)
  double complex  , allocatable :: LU(:,:)
  double complex                :: det
 
  double precision              :: normCNST
 
  integer(8), allocatable       :: IPIV(:)
  integer(8)                    :: s
  integer(8)                    :: site
  integer(8)                    :: k
  integer(8)                    :: kpoint
  
 !--------------------------------------------------!
 ! SUBROUTINE: GENERATION OF MATRIX M[sigma]        !
 !--------------------------------------------------!

  allocate( M(nEL,nEL)  )
  allocate( LU(nEL,nEL) )
  allocate( IPIV(nEL)   )

  M(:,:)   = (0.d0,0.d0)

  site = 0
  do s = 1,nSITES,1                                                                                                                      
     if ( rDET(s).eq.1 ) then
        site   = site + 1
        kpoint = 0
        do k = 1,nKPOINTS,1
           if ( kDET(k).eq.1  ) then
              kpoint         = kpoint + 1
              M(site,kpoint) = OVERLAP(s,k)
           endif
        enddo
     endif
  enddo
  
  call eval_complex_determinant( nEL,nEL,M(:,:),LU(:,:),IPIV(:),det )

  detM_BULK = det

  deallocate( M    )
  deallocate( LU   )
  deallocate( IPIV )

  return

END FUNCTION detM_BULK

!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_PROJ_TWFC_BULK( TWFC,CIindx,CI_K,OVERLAP,KupDMAT,KdownDMAT,upDMAT,downDMAT, &
                & nKPOINTS,nSITES,nELup,nELdown,nDET,nDETup,nDETdown,nTWFCDET )
        
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

  integer(8)      , intent(in)  :: nKPOINTS
  integer(8)      , intent(in)  :: nSITES
  integer(8)      , intent(in)  :: nELup , nELdown
  integer(8)      , intent(in)  :: nDET  , nDETup , nDETdown
  integer(8)      , intent(in)  :: nTWFCDET

  integer(8)      , intent(in)  :: KupDMAT(nKPOINTS,nDETup), KdownDMAT(nKPOINTS,nDETdown)
  integer(8)      , intent(in)  ::  upDMAT(nSITES,nDETup)  ,  downDMAT(nSITES,nDETdown)


  integer(8)      , intent(in)  :: CIindx(nTWFCDET)

  double precision, intent(in)  :: CI_K(nDET)

  double complex  , intent(in)  :: OVERLAP(nSITES,nKPOINTS)

  double precision, intent(out) :: TWFC(nDET)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  double complex  , allocatable :: TWFCaux(:)
  double complex                :: detM_BULK
  double complex                :: detM
  double complex                :: detMup  , detMdown

  integer(8)                    :: kDET
  integer(8)                    :: I , Iup , Idown
  integer(8)                    :: Ik, Ikup, Ikdown

 !-------------------------------------------------!
 ! ALLOCATION OF VARIABLES                         !
 !-------------------------------------------------!
 
   allocate( TWFCaux(nDET) )

 !-------------------------------------------------!
 ! SUBROUTINE: CALCULATE THE PROJECTED TRIAL WFC   !
 !-------------------------------------------------!

  TWFCaux(:) = (0.d0,0.d0)

  if ( (nELup .ne. 0) .and. (nELdown .ne. 0) ) then

     do I = 1,nDET,1
        call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
        do kDET = 1,nTWFCDET,1
           Ik          = CIindx(kDET)
           call spinDETindex( Ik,Ikup,Ikdown,nDETup,nDETdown )
           detMup      = detM_BULK( OVERLAP(:,:),upDMAT(:,Iup)    ,KupDMAT(:,Ikup)    ,nELup  ,nSITES,nKPOINTS )
           detMdown    = detM_BULK( OVERLAP(:,:),downDMAT(:,Idown),KdownDMAT(:,Ikdown),nELdown,nSITES,nKPOINTS )
           detM        = detMup * detMdown
           TWFCaux(I)  = TWFCaux(I) + CI_K(Ik) * detM
        enddo
     enddo

  elseif ( (nELup .ne. 0) .and. (nELdown .eq. 0) ) then

     do I = 1,nDET,1
        call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
        do kDET = 1,nTWFCDET,1
           Ik          = CIindx(kDET)
           call spinDETindex( Ik,Ikup,Ikdown,nDETup,nDETdown )
           detMup      = detM_BULK( OVERLAP(:,:),upDMAT(:,Iup),KupDMAT(:,Ikup),nELup,nSITES,nKPOINTS )
           detM        = detMup
           TWFCaux(I)  = TWFCaux(I) + CI_K(Ik) * detM
        enddo
     enddo

  elseif ( (nELup .eq. 0) .and. (nELdown .ne. 0) ) then

     do I = 1,nDET,1
        call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
        do kDET = 1,nTWFCDET,1
           Ik          = CIindx(kDET)
           call spinDETindex( Ik,Ikup,Ikdown,nDETup,nDETdown )
           detMdown    = detM_BULK( OVERLAP(:,:),downDMAT(:,Idown),KdownDMAT(:,Ikdown),nELdown,nSITES,nKPOINTS )
           detM        = detMdown
           TWFCaux(I)  = TWFCaux(I) + CI_K(Ik) * detM
        enddo
     enddo

  endif

  TWFC(:) = dble( TWFCaux(:) )

 !-------------------------------------------------!
 ! DEALLOCATION OF VARIABLES                       !
 !-------------------------------------------------!

  deallocate( TWFCaux )

  return

END SUBROUTINE get_PROJ_TWFC_BULK

!----------------------------------------------------------------------------------------------------------
