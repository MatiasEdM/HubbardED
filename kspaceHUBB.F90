SUBROUTINE solve_KSPACE_HUBB(BoundCond,Hk,CI_K,KupDMAT,KdownDMAT,KPOINTS,KP_ENRG,recVEC, &
                & nKPOINTS,nEL,nELup,nELdown,nDET,nDETup,nDETdown,t,U)

!----------------------------------------------------------------------!
! SUBROUTINE to construct and diagonalize the HUBBARD HAMILTONIAN in   !
! K-SPACE using HF-ORBITALS as a Single-Particle Basis Set.            !
! Returns as output the CI coefficients of the exact GS-WFC in the     !
! MULTI DET. expansion for the construction of the MULTI. DET. TWFC of !
! the FNA.                                                             !
!----------------------------------------------------------------------!

 implicit none

 !---------------------------------------------------!
 ! PARAMETERS                                        !
 !---------------------------------------------------!

  integer(8)      , parameter    :: dm = 2

 !---------------------------------------------------!
 ! INPUT VARIABLES                                   !
 !---------------------------------------------------!

  character(5)    , intent(in)  :: BoundCond

  integer(8)      , intent(in)  :: nKPOINTS
  integer(8)      , intent(in)  :: nEL   , nELup  , nELdown
  integer(8)      , intent(in)  :: nDET  , nDETup , nDETdown

  integer(8)      , intent(in)  :: KupDMAT(nKPOINTS,nDETup)
  integer(8)      , intent(in)  :: KdownDMAT(nKPOINTS,nDETdown)

  double precision, intent(in)  :: t
  double precision, intent(in)  :: U

  double precision, intent(in)  :: recVEC(dm,dm)

  double precision, intent(in)  :: KPOINTS(dm,nKPOINTS)
  double precision, intent(in)  :: KP_ENRG(nKPOINTS)

  double precision, intent(out) :: Hk(nDET,nDET)
  double precision, intent(out) :: CI_K(nDET)

 !---------------------------------------------------!
 ! LOCAL VARIABLES                                   !
 !---------------------------------------------------!

  double precision, allocatable :: CIk(:,:)
  double precision, allocatable :: Ek(:)

 !---------------------------------------------------!
 ! SUBROUTINE: GENERATION OF H(k) K-SPACE HUBB. HAM. !
 !---------------------------------------------------!

  allocate( CIk(nDET,nDET), Ek(nDET) )

  call generate_KSPACE_HAMILTONIAN( BoundCond,Hk(:,:),KupDMAT(:,:),KdownDMAT(:,:),KPOINTS(:,:),KP_ENRG(:), &
          & recVEC(:,:),nKPOINTS,nEL,nELup,nELdown,nDET,nDETup,nDETdown,t,U )
   
 !---------------------------------------------------!
 ! SUBROUTINE: DIAG. OF H(k) K-SPACE HUBB. HAM.      !
 !---------------------------------------------------!

  CIk(:,:) = Hk(:,:)

  call diagonalize_matrix( nDET,CIk(:,:),Ek(:) )

 !---------------------------------------------------!
 ! RECOVER CI(k) COEFFICIENTS OF GS-WFC              !
 !---------------------------------------------------!

  CI_K(:) = CIk(:,1)

 !---------------------------------------------------!
 ! OUTPUT WRITING: K-SPACE CALCULATION               !
 !---------------------------------------------------!

  call write_output_kspace( Hk(:,:),CIk(:,:),Ek(:),KupDMAT(:,:),KdownDMAT(:,:),KPOINTS(:,:),KP_ENRG(:), &
          & nKPOINTS,nEL,nELup,nELdown,nDET,nDETup,nDETdown,t,U )

 !---------------------------------------------------!
 ! MEMORY DEALLOCATION                               !
 !---------------------------------------------------!

  deallocate( CIk, Ek )

  return

END SUBROUTINE solve_KSPACE_HUBB

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_KSPACE_HAMILTONIAN( BoundCond,Hk,KupDMAT,KdownDMAT,KPOINTS,KP_ENRG, &
          & recVEC,nKPOINTS,nEL,nELup,nELdown,nDET,nDETup,nDETdown,t,U )

!----------------------------------------------------------------------!
! SUBROUTINE to construct the HUBBARD HAMILTONIAN in K-SPACE using     !
! HF-ORBITALS as a Single-Particle Basis Set.                          !
!----------------------------------------------------------------------!

  implicit none

 !---------------------------------------------------!
 ! PARAMETERS                                        !
 !---------------------------------------------------!

  integer(8)      , parameter    :: dm = 2

 !---------------------------------------------------!
 ! INPUT VARIABLES                                   !
 !---------------------------------------------------!

  character(5)    , intent(in)  :: BoundCond

  integer(8)      , intent(in)  :: nKPOINTS
  integer(8)      , intent(in)  :: nEL   , nELup  , nELdown
  integer(8)      , intent(in)  :: nDET  , nDETup , nDETdown

  integer(8)      , intent(in)  :: KupDMAT(nKPOINTS,nDETup)
  integer(8)      , intent(in)  :: KdownDMAT(nKPOINTS,nDETdown)

  double precision, intent(in)  :: t
  double precision, intent(in)  :: U

  double precision, intent(in)  :: recVEC(dm,dm)

  double precision, intent(in)  :: KPOINTS(dm,nKPOINTS)
  double precision, intent(in)  :: KP_ENRG(nKPOINTS)

  double precision, intent(out) :: Hk(nDET,nDET)

 !---------------------------------------------------!
 ! LOCAL VARIABLES                                   !
 !---------------------------------------------------!

  integer(8)                    :: globalDETindex
  integer(8)                    :: newDETindex
  integer(8)                    :: parity_IJ

  integer(8)      , allocatable ::     upOCC(:)
  integer(8)      , allocatable ::   upUNOCC(:)
  integer(8)      , allocatable ::   downOCC(:)
  integer(8)      , allocatable :: downUNOCC(:)

  integer(8)      , allocatable :: IupDET(:)
  integer(8)      , allocatable :: IdownDET(:)
  integer(8)      , allocatable :: JupDET(:)
  integer(8)      , allocatable :: JdownDET(:)

  integer(8)                    :: I, Iup, Idown
  integer(8)                    :: J, Jup, Jdown
  integer(8)                    :: kpoint

  integer(8)                    ::   upKocc,   upKunocc
  integer(8)                    :: downKocc, downKunocc


  integer(8)                    :: ik, jk, kk ,lk

  integer(8)                    :: diff
  integer(8)      , allocatable :: kpos(:)

  integer(8)                    ::   upPARITY
  integer(8)                    :: downPARITY

  double precision, allocatable :: kLIM1(:)
  double precision, allocatable :: kLIM2(:)

  double precision, allocatable :: q(:)
  double precision, allocatable :: kNEW(:)

  logical                       :: TR


 !---------------------------------------------------!
 ! SUBROUTINE: INTERM. STEP OF K-COMPUTATIOANL BOX   !
 !---------------------------------------------------!

  allocate( kLIM1(2) )
  allocate( kLIM2(2) )

  kLIM1(1) = minval( KPOINTS(1,:) )
  kLIM1(2) = maxval( KPOINTS(1,:) )

  kLIM2(1) = minval( KPOINTS(2,:) )
  kLIM2(2) = maxval( KPOINTS(2,:) )

 !---------------------------------------------------!
 ! SUBROUTINE: GENERATION OF H[k]                    !
 !---------------------------------------------------!

  Hk(:,:) = 0.d0

 !---------------------------------------------------!
 ! SUBROUTINE: H[k](I,I) DIAGONAL TERMS              !
 !---------------------------------------------------!

  do I = 1,nDET,1
     call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
     do kpoint = 1,nKPOINTS,1
        Hk(I,I) = Hk(I,I) + ( KupDMAT(kpoint,Iup) + KdownDMAT(kpoint,Idown) ) * KP_ENRG(kpoint)
     enddo
        Hk(I,I) = Hk(I,I) + (U/dble(nKPOINTS)) * (nELup * nELdown)
  enddo

 !---------------------------------------------------!
 ! SUBROUTINE: H[k](I,J) OFF-DIAGONAL TERMS          !
 !---------------------------------------------------!


  if ( (U.ne.0.d0).and.(nELup.ge.1).and.(nELdown.ge.1) ) then

     allocate( upOCC(nELup)                )
     allocate( upUNOCC(nKPOINTS-nELup)     )
     allocate( downOCC(nELdown)            )
     allocate( downUNOCC(nKPOINTS-nELdown) )

     allocate( IupDET(nKPOINTS)   )
     allocate( IdownDET(nKPOINTS) )
     allocate( JupDET(nKPOINTS)   )
     allocate( JdownDET(nKPOINTS) )
     allocate( kpos(2)            )

     allocate( q(dm)    )
     allocate( kNEW(dm) )

     do I = 1,nDET,1
        call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
        IupDET(:)   = KupDMAT(:,Iup)
        IdownDET(:) = KdownDMAT(:,Idown) 
        call get_OCC_UNOCC_index( upOCC(:)  ,upUNOCC(:)  ,IupDET(:)  ,nKPOINTS,nELup   )
        call get_OCC_UNOCC_index( downOCC(:),downUNOCC(:),IdownDET(:),nKPOINTS,nELdown )
        do upKocc = 1,nELup,1
           do upKunocc = 1,(nKPOINTS-nELup),1
              JupDET(:)                   = IupDET(:)
              JupDET(   upOCC(  upKocc) ) = 0
              JupDET( upUNOCC(upKunocc) ) = 1
              Jup                         = newDETindex( JupDET(:),KupDMAT(:,:),nKPOINTS,nDETup )
              call count_diff( diff,kpos(:),IupDET(:),JupDET(:),nKPOINTS )
              upPARITY                    = parity_IJ( IupDET(:),JupDET(:),kpos(:),nKPOINTS )
              q(:)                        = KPOINTS(:,upUNOCC(upKunocc)) - KPOINTS(:,upOCC(upKocc))
              do downKocc = 1,nELdown,1
                 kNEW(:)  = KPOINTS(:,downOCC(downKocc)) - q(:)
                 call transport_BZ( kNEW(:),kLIM1(:),kLIM2(:),recVEC(:,:) )
                 call check_TRANSITION( TR,downKunocc,kNEW(:),downUNOCC(:),KPOINTS(:,:),nKPOINTS,nKPOINTS-nELdown )
                 if ( TR ) then
                    JdownDET(:) = IdownDET(:)
                    JdownDET(  downOCC(  downKocc) ) = 0
                    JdownDET(downUNOCC(downKunocc) ) = 1
                    Jdown                            = newDETindex( JdownDET(:),KdownDMAT(:,:),nKPOINTS,nDETdown )
                    call count_diff( diff,kpos(:),IdownDET(:),JdownDET(:),nKPOINTS )
                    downPARITY                       = parity_IJ( IdownDET(:),JdownDET(:),kpos(:),nKPOINTS )
                    J                                = globalDETindex(Jup,Jdown,nDETup,nDETdown)
                    Hk(I,J)                          = (U/dble(nKPOINTS)) * upPARITY * downPARITY
                 endif
              enddo
           enddo
        enddo
     enddo

     deallocate( upOCC     )
     deallocate( upUNOCC   )
     deallocate( downOCC   )
     deallocate( downUNOCC )

     deallocate( IupDET   )
     deallocate( IdownDET )
     deallocate( JupDET   )
     deallocate( JdownDET )
     deallocate( kpos     )

     deallocate( q    )
     deallocate( kNEW )

  endif

  deallocate( kLIM1 )
  deallocate( kLIM2 )

  return

END SUBROUTINE generate_KSPACE_HAMILTONIAN

!----------------------------------------------------------------------------------------------------------

SUBROUTINE check_TRANSITION(TR,KPindex,KP,LIST,KPOINTS,nKPOINTS,nEL)

!----------------------------------------------------------------------!
! SUBROUTINE to check if the K-POINT KP belongs to the list LIST of    !
! K-POINTS returning the result on the logical variable TR = .T./.F.   !
! and if TR = .T. returns also the INDEX ON THE LIST to reference to   !
! the appropriate KPOINT, i.e. the INDEX KPindex references to the     !
! K-POINT as:                                                          !
!             KP(:) = KPOINTS(:, LIST(KPindex) )                       !
!----------------------------------------------------------------------!

  implicit none

 !---------------------------------------------------!
 ! PARAMETERS                                        !
 !---------------------------------------------------!

  integer(8)      , parameter    :: dm   = 2

  double precision, parameter    :: thrs = 1E-5 

 !---------------------------------------------------!
 ! INPUT VARIABLES                                   !
 !---------------------------------------------------!

  integer(8)      , intent(in)   :: nKPOINTS
  integer(8)      , intent(in)   :: nEL
  
  integer(8)      , intent(in)   :: LIST(nEL)

  integer(8)      , intent(out)  :: KPindex


  double precision, intent(in)   :: KPOINTS(dm,nKPOINTS)
  double precision, intent(in)   :: KP(dm)

  logical         , intent(out)  :: TR

 !---------------------------------------------------!
 ! LOCAL VARIABLES                                   !
 !---------------------------------------------------!

  integer(8)                    :: k
 
 !---------------------------------------------------!
 ! SUBROUTINE: CHECK THE KP INTO THE LIST OF KPs     !
 !---------------------------------------------------!

  TR = .FALSE.

  do k = 1,nEL,1
     if ( ( abs(KP(1)-KPOINTS(1,LIST(k))).lt.thrs ).and.( abs(KP(2)-KPOINTS(2,LIST(k))).lt.thrs ) ) then
        TR      = .TRUE.
        KPindex = k
        exit
     endif
  enddo
  
  return

END SUBROUTINE check_TRANSITION
