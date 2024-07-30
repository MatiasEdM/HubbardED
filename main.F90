PROGRAM main

!-----------------------------------------------------------------------------------------------------------------!
! PROGRAM for the exact diagonalization of the 2D [rectangular lattice] Fermi-Hubbard Model in Slater-Determinant !
! real-space basis.Explicit definition of up- and down-spin subspaces (each determinant is explcitly defined as i !
! the tensor product of an up-spin and down-spin determinant).                                                    !
!-----------------------------------------------------------------------------------------------------------------!

  implicit none
  

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter   :: dm = 2

 !-------------------------------------------------!
 ! VARIABLE DECLARATION                            !
 !-------------------------------------------------!

  character(6)                  :: SysType                        !System Type: LADDER (PBC[Y|X]), BULK (PBC[XY]) or OPEN (OBC)
  character(5)                  :: BoundCond                      !Boundary Conditions: Open or Periodic [X|Y|XY].

  integer(8)                    :: nSITES, nSITESx, nSITESy       !Definition of the 2D lattice. 
  integer(8)                    :: nEL, nELup, nELdown            !Number of up- and down-spin electrons.
  integer(8)                    :: nDET, nDETup, nDETdown         !Dimension of Hilbert space and spin-subspaces.

  integer(8)                    :: nTWFCDET                       !Number of DETs. in K-SPACE to be included into the MC-TWFC (short CI-EXPANSION).

  integer(8)                    :: nPDET                          !Dimension of the P-SPACE in the P/Q PARTITION (RELAXATION of FNA).

  integer(8)      , allocatable :: LMAT(:,:)                      !Connectivity MAT of the lattice points.
  integer(8)                    :: NEIGH(4)                       !Vector for neighboring lattice sites for a 2D rectangular lattice. 

  integer(8)      , allocatable :: upDMAT(:,:), downDMAT(:,:)     !Determinants MAT of the up- and down-spin electrons.  
  integer(8)      , allocatable :: hopMAPup(:,:), hopMAPdown(:,:) !Hopping MAPS for the up- and down-spin electrons.
  integer(8)      , allocatable :: dOCC(:,:)                      !The DOUBLE OCCUPANCY MAT (Ns,Ndet) of doubly occupied lattice sites.

  integer(8)                    :: row

  double precision              :: A1o, A2o                       !LATTICE PARAMETERS.

  double precision, allocatable :: lattVEC(:,:), recVEC(:,:)      !Direct Lattice and Reciprocal Lattice VEC.
  double precision, allocatable :: lattPOS(:,:)                   !List of position of lattice sites in cartesian coordinates.

  double precision              :: t ,U                           !HUBB. PARAMETERS.
  double precision              :: tk,Uk                          !HUBB. PARAMETERS for the DIAG. in K-SPACE H[k].

  double precision, allocatable :: H(:,:)                         !Hamiltonian MAT.
  double precision, allocatable :: CI(:,:)                        !MAT of CI coefficients.
  double precision, allocatable :: E(:)                           !VEC of EIGENVALUES (energies) of the many-body states.
  double precision, allocatable :: S2(:,:)                        !S2 MAT.
  double precision, allocatable :: gsVsf(:)                       ! SIGN FLIP POT. of the FNA foor the GS-WFC.

  double precision              :: K, V                           !KINETIC and INTERACTION energies.
  double precision              :: int_energy, kin_energy         !FUNCTIONS.

  double precision              :: S                              !SPIN S.
  double precision              :: tot_spin_sq                    !<S^2> FUNCTION.

  double precision              :: normDC                         !L2 NORM of the DC VEC to estimate the magnitude of the sign problem.

  character(3)                  :: stoq                           !Enables STOQUASTIZED CALCULATION.
  character(3)                  :: trunc                          !Enables TRUNCATION of the Hilbert Space (removing nodes from
                                                                  !the HF Determinant in the non-interacting limit (U=0).

  character(3)                  :: fixed_node                     !Enables FIXED NODE APPROXIMATION CALCULATION.
  character(3)                  :: fna_trunc                      !Enables FNA with TRUNCATED DIAGONALIZATION (eliminates nodes of the Trial WFC from the Hilbert Space).
  character(3)                  :: subspace_coeff                 !Enables the RECALC. of the CI-COEFF. of the TWFC in the SUPSPACE spanned by the first (highest amplitude) nTWFCDET DETs.
  character(3)                  :: fna_stoq                       !Enables STOQ. FNA CALCULATION.
  character(3)                  :: fna_relax                      !Enables FNA RELAXATION according to the P/Q PARTITION METHOD.
  
  character(3)                  :: tot_spin                       !Enables TOTAL SPIN CALCULATION.
  character(3)                  :: mat_spin                       !Enables TOTAL SPIN^2 MATRIX CALCULATION.
  character(3)                  :: tot_spin_stoq                  !Enables TOTAL SPIN CALCULATION of STOQUASTIZED system.
  character(3)                  :: tot_spin_trunc                 !Enables TOTAL SPIN CALCULATION of TRUNCATED Hilbert space.
  character(3)                  :: tot_spin_fna                   !Enables TOTAL SPIN CALCULATION of the FIXED NODE Hamiltonian.
  character(3)                  :: tot_spin_sq_fna                !Enables TOTAL SPIN CALCULATION of FNA STOQUASTIZED system.

 !-------------------------------------------------!
 ! INPUT READING AND MEMORY ALLOCATION             !
 !-------------------------------------------------!

  allocate( lattVEC(dm,dm), recVEC(dm,dm) )

  call read_input( SysType,BoundCond,nSITES,nSITESx,nSITESy,nEL,nELup,nELdown,nDET,nDETup,nDETdown,       &
          & A1o,A2o,lattVEC(:,:),t,U,stoq,trunc,fixed_node,fna_trunc,tk,Uk,nTWFCDET,nPDET,subspace_coeff, & 
          & fna_relax,fna_stoq,tot_spin,mat_spin,tot_spin_stoq,tot_spin_trunc,tot_spin_fna,tot_spin_sq_fna )

  allocate( lattPOS(dm,nSITES) )

  allocate( LMAT(nSITES,nSITES) )
  allocate( upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown) )
  allocate( hopMAPup(nDETup,nDETup), hopMAPdown(nDETdown,nDETdown) )  
  allocate( dOCC(nSITES,nDET) )
  allocate( H(nDET,nDET), CI(nDET,nDET), E(nDET) )
  allocate( S2(nDET,nDET) )
  allocate( gsVsf(nDET) )

 !-------------------------------------------------!
 ! SET UP: CONNECT. | DETS. | HOPP. MAPS | D. OCC  !
 !-------------------------------------------------!

  !call generate_lattPOS( lattPOS(:,:), lattVEC(:,:), nSITESx, nSITESy, nSITES)

  lattPOS(1,1)  = 0.d0 * lattVEC(1,1) + 0.d0 * lattVEC(1,2) ; lattPOS(2,1)  = 0.d0 * lattVEC(2,1) + 0.d0 * lattVEC(2,2)
  lattPOS(1,2)  = 0.d0 * lattVEC(1,1) + 1.d0 * lattVEC(1,2) ; lattPOS(2,2)  = 0.d0 * lattVEC(2,1) + 1.d0 * lattVEC(2,2)
  lattPOS(1,3)  = 1.d0 * lattVEC(1,1) + 1.d0 * lattVEC(1,2) ; lattPOS(2,3)  = 1.d0 * lattVEC(2,1) + 1.d0 * lattVEC(2,2)
  lattPOS(1,4)  = 2.d0 * lattVEC(1,1) + 1.d0 * lattVEC(1,2) ; lattPOS(2,4)  = 2.d0 * lattVEC(2,1) + 1.d0 * lattVEC(2,2)
  lattPOS(1,5)  = 0.d0 * lattVEC(1,1) + 2.d0 * lattVEC(1,2) ; lattPOS(2,5)  = 0.d0 * lattVEC(2,1) + 2.d0 * lattVEC(2,2)
  lattPOS(1,6)  = 1.d0 * lattVEC(1,1) + 2.d0 * lattVEC(1,2) ; lattPOS(2,6)  = 1.d0 * lattVEC(2,1) + 2.d0 * lattVEC(2,2)
  lattPOS(1,7)  = 2.d0 * lattVEC(1,1) + 2.d0 * lattVEC(1,2) ; lattPOS(2,7)  = 2.d0 * lattVEC(2,1) + 2.d0 * lattVEC(2,2)
  lattPOS(1,8)  = 0.d0 * lattVEC(1,1) + 3.d0 * lattVEC(1,2) ; lattPOS(2,8)  = 0.d0 * lattVEC(2,1) + 3.d0 * lattVEC(2,2)
  lattPOS(1,9)  = 1.d0 * lattVEC(1,1) + 3.d0 * lattVEC(1,2) ; lattPOS(2,9)  = 1.d0 * lattVEC(2,1) + 3.d0 * lattVEC(2,2)
  lattPOS(1,10) = 2.d0 * lattVEC(1,1) + 3.d0 * lattVEC(1,2) ; lattPOS(2,10) = 2.d0 * lattVEC(2,1) + 3.d0 * lattVEC(2,2)

  call generate_reciprocalVEC( recVEC(:,:), lattVEC(:,:) )

  !call generate_connectMAT( LMAT(:,:),BoundCond,nSITES,nSITESx,nSITESy  )

  open ( unit=10, file='LMAT.dat', status='old', access='sequential' )
       read(10,*)
       read(10,*)
       do row = 1,nSITES,1
          read(10,*) LMAT(row,:)
       enddo      
  close( unit=10, status='keep' ) 

  call generate_rsDETS( upDMAT(:,:)  ,nSITES,nELup  ,nDETup   )
  call generate_rsDETS( downDMAT(:,:),nSITES,nELdown,nDETdown )

  call generate_hopMAP( hopMAPup(:,:)  ,upDMAT(:,:)  ,LMAT(:,:),nSITES,nDETup )
  call generate_hopMAP( hopMAPdown(:,:),downDMAT(:,:),LMAT(:,:),nSITES,nDETdown )

  call generate_doubleOCC( dOCC(:,:),upDMAT(:,:),downDMAT(:,:),nSITES,nDET,nDETup,nDETdown )

 !-------------------------------------------------!
 ! HAMILTONIAN MATRIX AND EXACT DIAGONALIZATION    !
 !-------------------------------------------------!

  call generate_HAMILTONIAN( H(:,:),upDMAT(:,:),downDMAT(:,:),hopMAPup(:,:),hopMAPdown(:,:),dOCC(:,:),nSITES,nDET, & 
          & nDETup,nDETdown,t,U )

  CI(:,:) = H(:,:)

  call diagonalize_matrix( nDET,CI(:,:),E(:) )

 !-------------------------------------------------!
 ! GROUND STATE ENERGY DECOMPOSITION CALCULATION   !
 !-------------------------------------------------!

  V = int_energy( CI(:,1),dOCC(:,:),nSITES,nDET,U )
  K = kin_energy( E(1),V )

 !-------------------------------------------------!
 ! GROUND STATE TOTAL SPIN S CALCULATION           !
 !-------------------------------------------------!
 
  if ( tot_spin .eq. 'yes' ) then
      S = tot_spin_sq( CI(:,1),upDMAT(:,:),downDMAT(:,:),nSITES,nDET,nDETup,nDETdown )   !ATTENTION THIS IS <S^2> NOT S.
      if ( mat_spin .eq. 'yes' ) then
         call generate_S2MAT( S2(:,:),upDMAT(:,:),downDMAT(:,:),nSITES,nDET,nDETup,nDETdown )
      endif
  endif

 !-------------------------------------------------!
 ! SING-PROBLEM ESTIMATION: BY L2 NORM ||DC||2     !
 !-------------------------------------------------!

  call estimate_sign_problem( normDC,H(:,:),CI(:,1),nDET )

 !-------------------------------------------------!
 ! OUTPUT WRITING: NORMAL CALCULATION              !
 !-------------------------------------------------!

  call write_output_init( SysType,BoundCond, nSITES,nSITESx,nSITESy,nEL,nELup,nELdown,nDET,nDETup, &
          & nDETdown,t,U,A1o,A2o,lattVEC(:,:),recVEC(:,:),LMAT(:,:),upDMAT(:,:),downDMAT(:,:),     &
          & hopMAPup(:,:),hopMAPdown(:,:),dOCC(:,:),H(:,:),CI(:,:),E(:),K,V,tot_spin,mat_spin,     &
          & S,S2(:,:),normDC )

 !-------------------------------------------------!
 ! STOQ. HAMILTONIAN MATRIX EXACT DIAGONALIZATION  !
 !-------------------------------------------------!

  if ( stoq .eq. 'yes' ) then  
     call solve_stoquastized( H(:,:),E(:),nSITES,nDET,nDETup,nDETdown,upDMAT,downDMAT,tot_spin_stoq )
  endif
 
 !-------------------------------------------------!
 ! TRUNCATED HILBERT SPACE DIAGONALIZATION         !
 !-------------------------------------------------!

  if ( trunc .eq. 'yes' ) then
     call solve_truncated( upDMAT(:,:),downDMAT(:,:),H(:,:),CI(:,:),E(:),nSITES,nDET,nDETup,nDETdown,U, &
             & tot_spin_trunc )
  endif

 !-------------------------------------------------!
 ! FIXED NODE APPROXIMATION BASED ON HF-MFT DET.   !
 !-------------------------------------------------!

  call generate_Vsf_gsWFC( gsVsf(:),H(:,:),CI(:,1),nDET )

  if ( fixed_node .eq. 'yes') then
     if (SysType .eq. 'bulk') then
        call FNA_BULK( BoundCond,upDMAT(:,:),downDMAT(:,:),H(:,:),CI(:,1),lattVEC(:,:),lattPOS(:,:), &
                & recVEC(:,:),A1o,A2o,nSITES,nSITESx,nSITESy,nEL,nELup,nELdown,nDET,nDETup,nDETdown, & 
                & t,U,tk,Uk,nTWFCDET,nPDET,gsVsf(:),fna_trunc,subspace_coeff,fna_relax,fna_stoq,     &
                & tot_spin_fna,tot_spin_sq_fna )
      elseif (SysType .eq. 'ladder') then
!        call FNA_LADDER( BoundCond,upDMAT(:,:),downDMAT(:,:),H(:,:),lattVEC(:,:),lattPOS(:,:),       &
!                & recVEC(:,:),A1o,A2o,nSITES,nSITESx,nSITESy,nEL,nELup,nELdown,nDET,nDETup,nDETdown, &
!                & t,U,nTWFCDET,tk,Uk,gsVsf(:),orb_rot,fna_trunc,fna_stoq,tot_spin_fna,tot_spin_sq_fna )
         WRITE(*,*) '- FNA NOT AVAILABLE FOR LADDER SYSTEMS WITH THE NEW SCHEME (TO BE IMPLEMENTED) - '
         WRITE(*,*)
         WRITE(*,*) '  EXITING CALCULATION ...' 
      elseif( SysType .eq. 'open') then
        call FNA_ERROR()
      endif
  endif 

 !-------------------------------------------------!
 ! MEMORY DEALLOCATION                             !
 !-------------------------------------------------!

  deallocate( lattVEC, recVEC, lattPOS, LMAT ) 
  deallocate( upDMAT, downDMAT ) 
  deallocate( hopMAPup, hopMAPdown, dOCC )
  deallocate( H, CI, E, S2 )
  deallocate( gsVsf )

  stop

END PROGRAM main 
