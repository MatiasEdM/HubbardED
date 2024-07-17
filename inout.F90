!----------------------------------------------------------------------------------------------------------

SUBROUTINE read_input(SysType,BoundCond,nSITES,nSITESx,nSITESy,nEL,nELup,nELdown,nDET,nDETup,nDETdown,     &
                & A1o,A2o,lattVEC,t,U,stoq,trunc,fixed_node,fna_trunc,tk,Uk,nTWFCDET,nPDET,subspace_coeff, &
                & fna_relax,fna_stoq,tot_spin,mat_spin,tot_spin_stoq,tot_spin_trunc,tot_spin_fna,          &
                & tot_spin_sq_fna)

!----------------------------------------------------------------------------------------------------------------------!
! SUBROUTINE to read input for the 2D Fermi-HUbbard model in real space.                                               !
! EXAMPLE INPUT:                                                                                                       !
! &SYSTEM hubbard                                                                                                      !
!                                                                                                                      !
! system_type     bulk / ladder / open                                                                                 !
! boundary_cond   obc / pbcxy / pbcx / pbcy                                                                            !
! 2d_lattice      nSITESx nSITESy            - Refers to the number of sites in the A1 and A2 directions respectively. !
! up_electrons    nELup                                                                                                !
! down_electrons  nELdown                                                                                              !
! tot_spin        yes / no                                                                                             !
! mat_spin        yes / no                                                                                             !
!                                                                                                                      !
! &LATTICE PARAMETERS                                                                                                  !
!                                                                                                                      !
! A1 Ao1                                                                                                               !
! A2 Ao2                                                                                                               !
!                                                                                                                      !
! &LATTICE VECTORS                                                                                                     !
!                                                                                                                      !
! A1 a1[1] a1[2]                                                                                                       !
! A2 a2[1] a2[2]                                                                                                       !
!                                                                                                                      ! 
! &HUBBARD PARAMETERS                                                                                                  !
!                                                                                                                      !
! t Value                                                                                                              !
! U Value                                                                                                              !
!                                                                                                                      !
! &STOQUASTIZED SYSTEM                                                                                                 !
!                                                                                                                      !
! stoq            yes / no                                                                                             !
! tot_spin_stoq   yes / no                                                                                             !
!                                                                                                                      !
! &HILBERT SPACE                                                                                                       !
!                                                                                                                      !
! truncation      yes / no                                                                                             !
! tot_spin_trunc  yes / no                                                                                             !
!                                                                                                                      !
! &FIXED NODE APPROXIMATION                                                                                            !
!                                                                                                                      !
! FNA             yes / no                                                                                             !
! trunc_diag      yes / no                                                                                             !
! tot_spin_fna    yes / no                                                                                             !
!                                                                                                                      !
! &TRIAL WAVE FUNCITON                                                                                                 !
!                                                                                                                      !
! t[k]            Value                                                                                                !
! U[k]            Value                                                                                                !
! ndets           nTWFCDET                                                                                             !
! subspace_coeff  yes / no                                                                                             !
!                                                                                                                      !
! &FNA RELAXATION [P/Q PARTITION]                                                                                      !
!                                                                                                                      !
! relax           yes / no                                                                                             !
! pspace          nPDET                                                                                                !
!                                                                                                                      !
! &FNA STOQUASTIZED                                                                                                    !
!                                                                                                                      !
! stoq_fna        yes / no                                                                                             !
! tot_spin_sq_fna yes / no                                                                                             !
!                                                                                                                      !
!----------------------------------------------------------------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)     , parameter   :: dm = 2

 !-------------------------------------------------!
 ! OUTPUT VARIABLES                                !
 !-------------------------------------------------!
 
  character(6)    , intent(out) :: SysType
  character(5)    , intent(out) :: BoundCond

  integer(8)      , intent(out) :: nSITES, nSITESx, nSITESy 
  integer(8)      , intent(out) :: nEL, nELup, nELdown
  integer(8)      , intent(out) :: nDET, nDETup, nDETdown

  integer(8)      , intent(out) :: nTWFCDET

  integer(8)      , intent(out) :: nPDET
  
  double precision, intent(out) :: lattVEC(dm,dm)
  double precision, intent(out) :: A1o, A2o
  double precision, intent(out) :: t, U

  double precision, intent(out) :: tk,Uk

  character(3)    , intent(out) :: stoq
  character(3)    , intent(out) :: trunc

  character(3)    , intent(out) :: fixed_node
  character(3)    , intent(out) :: fna_trunc
  character(3)    , intent(out) :: subspace_coeff
  characteR(3)    , intent(out) :: fna_relax

  character(3)    , intent(out) :: fna_stoq

  character(3)    , intent(out) :: mat_spin
  character(3)    , intent(out) :: tot_spin
  character(3)    , intent(out) :: tot_spin_stoq
  character(3)    , intent(out) :: tot_spin_trunc
  character(3)    , intent(out) :: tot_spin_fna
  character(3)    , intent(out) :: tot_spin_sq_fna

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)       :: fact, binom

  double precision :: normCOEFF

  character(15)    :: SYSTEM
  character(2)     :: LATTICE
  character(1)     :: PARAMETERS
  character(4)     :: STOQUASTIZED
  character(10)    :: HILBERT
  character(15)    :: SPIN
  character(3)     :: FN

 !-------------------------------------------------!
 ! READING FROM STANDARD INPUT                     !
 !-------------------------------------------------!

  read(*,*)
  read(*,*)
  read(*,*) SYSTEM, SysType
            if ((SysType.ne.'bulk').and.(SysType.ne.'ladder').and.(SysType.ne.'open')) then
               write(*,*) 'ERROR WITH TYPE SYSTEM: NOT RECOGNIZED'
               stop
            endif 
  read(*,*) SYSTEM, BoundCond
            if ((BoundCond.ne.'obc').and.(BoundCond.ne.'pbcxy') &
                    & .and.(BoundCond.ne.'pbcx').and.(BoundCond.ne.'pbcy'))  then
               write(*,*) 'ERROR WITH BOUNDARY CONDITIONS: NOT RECOGNIZED'
               stop
            endif
            select case (SysType)
               case ('bulk')
                    if ((BoundCond.ne.'pbcxy')) then
                       write(*,*) 'ERROR: SYSTEM TYPE AND BOUNDARY CONDITIONS DO NOT MATCH'
                       stop
                    endif
               case ('ladder')
                    if ((BoundCond.ne.'pbcx').and.(BoundCond.ne.'pbcy')) then
                       write(*,*) 'ERROR: SYSTEM TYPE AND BOUNDARY CONDITIONS DO NOT MATCH'
                       stop
                    endif                           
               case ('open')
                    if ((BoundCond.ne.'obc')) then
                       write(*,*) 'ERROR: SYSTEM TYPE AND BOUNDARY CONDITIONS DO NOT MATCH'
                       stop
                    endif
            endselect
  read(*,*) SYSTEM, nSITESx, nSITESy
  read(*,*) SYSTEM, nELup
  read(*,*) SYSTEM, nELdown
  read(*,*) SPIN  , tot_spin
  read(*,*) SPIN  , mat_spin
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*) LATTICE , A1o
  read(*,*) LATTICE , A2o
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*) LATTICE   , lattVEC(1,1), lattVEC(2,1)
  read(*,*) LATTICE   , lattVEC(1,2), lattVEC(2,2)
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*) PARAMETERS, t
  read(*,*) PARAMETERS, U
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*) STOQUASTIZED, stoq
  read(*,*) SPIN        , tot_spin_stoq
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*) HILBERT     , trunc
  read(*,*) SPIN        , tot_spin_trunc
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*) FN          , fixed_node
  read(*,*) FN          , fna_trunc
  read(*,*) SPIN        , tot_spin_fna
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*) FN          , tk
  read(*,*) FN          , Uk
  read(*,*) FN          , nTWFCDET
  read(*,*) FN          , subspace_coeff
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*) FN          , fna_relax
  read(*,*) FN          , nPDET
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*) FN          , fna_stoq
  read(*,*) SPIN        , tot_spin_sq_fna

  lattVEC(:,1) = normCOEFF( lattVEC(:,1),dm ) * lattVEC(:,1)
  lattVEC(:,2) = normCOEFF( lattVEC(:,2),dm ) * lattVEC(:,2)

  lattVEC(:,1) = A1o * lattVEC(:,1)
  lattVEC(:,2) = A2o * lattVEC(:,2)

  nSITES       = nSITESx * nSITESy
  
  nEL          = nELup + nELdown
  
  nDETup       = binom(nSITES,nELup)
  nDETdown     = binom(nSITES,nELdown)
  nDET         = nDETup * nDETdown

  if ( fixed_node .eq. 'yes' ) then   
     if ( (nTWFCDET.lt.1).or.(nTWFCDET.gt.nDET) ) then
        write(*,*)
        write(*,*) 'ERROR WITH THE NUMBER OF DETERMINANTS OF THE CI TRIAL WFC'
        write(*,*)
        stop
     endif
     if ( fna_relax .eq. 'yes' ) then
        if ( nPDET .gt. nDET ) then
           write(*,*)
           write(*,*) 'ERROR WITH THE NUMBER OF DETERMINANTS OF THE P SUBSPACE [FNA RELAXATION | P/Q PARTITION]'
           write(*,*)
           stop
        elseif ( nPDET .eq. 0 ) then
           fna_relax = 'no'
        endif
     endif
  endif

  return

END SUBROUTINE read_input

!----------------------------------------------------------------------------------------------------------

SUBROUTINE write_output_init(SysType,BoundCond,nSITES,nSITESx,nSITESy,nEL,nELup,nELdown,nDET,nDETup,nDETdown,t,U, &
                & A1o,A2o,lattVEC,recVEC,LMAT,upDMAT,downDMAT,hopMAPup,hopMAPdown,dOCC,H,CI,E,K,V, &
                & tot_spin,mat_spin,S,S2,normDC)

!--------------------------------------------------------------------------------------------------!
! SUBROUTINE to write the output for the calculations of the 2D Fermi-Hubbard model in real space. !
!--------------------------------------------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter  :: dm   = 2
  double precision, parameter  :: thrs = 10E-8

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!
  
  character(6)    , intent(in) :: SysType
  character(5)    , intent(in) :: BoundCond

  integer(8)      , intent(in) :: nSITES, nSITESx, nSITESy 
  integer(8)      , intent(in) :: nEL, nELup, nELdown
  integer(8)      , intent(in) :: nDET, nDETup, nDETdown

  integer(8)      , intent(in) :: LMAT(nSITES,nSITES)
  integer(8)      , intent(in) :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)
  integer(8)      , intent(in) :: hopMAPup(nDETup,nDETup), hopMAPdown(nDETdown,nDETdown)
  integer(8)      , intent(in) :: dOCC(nSITES,nDET)

  double precision, intent(in) :: A1o, A2o

  double precision, intent(in) :: lattVEC(dm,dm), recVEC(dm,dm)
  
  double precision, intent(in) :: t, U

  double precision, intent(in) :: H(nDET,nDET)
  double precision, intent(in) :: CI(nDET,nDET)
  double precision, intent(in) :: E(nDET)
  double precision, intent(in) :: S2(nDET,nDET)

  double precision, intent(in) :: K, V

  double precision, intent(in) :: S

  double precision, intent(in) :: normDC 

  character(3)    , intent(in) :: tot_spin
  character(3)    , intent(in) :: mat_spin
 
  integer(8)                   :: I, Iup, Idown
  integer(8)                   :: site, ctr

 !-------------------------------------------------!
 ! WRITING STANDARD OUTPUT                         !
 !-------------------------------------------------!

  write(*,fmt=200) '---- 2-DIMENSIONAL FERMI-HUBBARD MODEL IN REAL-SPACE ----'
  write(*,*)
  write(*,fmt=200) '---- SYSTEM PROPERTIES ----------------------------------'
  write(*,*)
            if ( SysType .eq. 'bulk' ) then
               write(*,fmt=700) '- System Type                     : BULK'
            elseif( SysType .eq. 'ladder' ) then
               write(*,fmt=710) '- System Type                     : LADDER'
            elseif ( SysType .eq. 'open' ) then
               write(*,fmt=700) '- System Type                     : OPEN'
            endif
            if ( BoundCond .eq. 'obc' ) then
               write(*,fmt=500) '- Open Boundary Conditions (OBC)'
            elseif( BoundCond .eq. 'pbcxy' ) then
               write(*,fmt=510) '- XY Periodic Boundary Conditions (PBC)'
            elseif ( BoundCond .eq. 'pbcx' ) then
               write(*,fmt=520) '- X Periodic Boundary Conditions (PBC)'
            elseif (BoundCond .eq. 'pbcy' ) then
               write(*,fmt=520) '- Y Periodic Boundary Conditions (PBC)'
            endif
  write(*,*)
  write(*,fmt=210) '- Total Number of Sites (Ns)      : ', nSITES
  write(*,fmt=220) '- Lattice Geometry [x,y]          : ', nSITESX, nSITESy
  write(*,*)
  write(*,fmt=530) '- Lattice Parameters              : '
  write(*,fmt=550) '   A1o   = ', A1o
  write(*,fmt=550) '   A2o   = ', A2o
  write(*,*)
  write(*,fmt=530) '- Lattice Vectors                 :'
  write(*,fmt=540) '   a1[1] = ', lattVEC(1,1), ', a1[2] =', lattVEC(2,1)
  write(*,fmt=540) '   a2[1] = ', lattVEC(1,2), ', a2[2] =', lattVEC(2,2)
  write(*,*)
  write(*,fmt=530) '- Reciprocal Lattice Vectors      :'
  write(*,fmt=540) '   b1[1] = ', recVEC(1,1),  ', b1[2] =', recVEC(2,1)
  write(*,fmt=540) '   b2[1] = ', recVEC(1,2),  ', b2[2] =', recVEC(2,2)
  write(*,*)
  write(*,fmt=210) '- Total Number of Electrons (N)   : ', nEL
  write(*,fmt=210) '- Number of Up-Electrons   (Nu)   : ', nELup
  write(*,fmt=210) '- Number of Down-Electrons (Nd)   : ', nELdown
  write(*,*)
  write(*,fmt=210) '- Dimension of Hilbert Space      : ', nDET
  write(*,fmt=210) '- Dimension of Spin-Up Subspace   : ', nDETup
  write(*,fmt=210) '- Dimension of Spin-Down Subspace : ', nDETdown
  write(*,*)
  write(*,fmt=230) '- Hopping Parameter   (t)         : ', t
  write(*,fmt=230) '- On-Site Interaction (U)         : ', U
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- CONNECTIVITY MATRIX --------------------------------'
  write(*,*)
  call writeINTMAT( LMAT(:,:),nSITES,nSITES )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- DETERMINANTAL MATRICES -----------------------------'
  write(*,*)
  write(*,fmt=240) '----  UP-SPIN SUBSPACE  ----'
  write(*,*)
  call writeINTMAT( upDMAT(:,:),nSITES,nDETup )
  write(*,*)
  write(*,fmt=240) '---- DOWN-SPIN SUBSPACE ----'
  write(*,*)
  call writeINTMAT( downDMAT(:,:),nSITES,nDETdown )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- HOPPING MAP MATRICES -------------------------------'
  write(*,*)
  write(*,fmt=240) '----  UP-SPIN SUBSPACE  ----'
  write(*,*)
  call writeINTMAT( hopMAPup(:,:),nDETup,nDETup )
  write(*,*)
  write(*,fmt=240) '---- DOWN-SPIN SUBSPACE ----'
  write(*,*)
  call writeINTMAT( hopMAPdown(:,:),nDETdown,nDETdown )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- DOUBLE OCCUPANCY DETERMINANT MAP -------------------'
  write(*,*)
  call writeINTMAT( dOCC(:,:),nSITES,nDET )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- REAL-SPACE BASIS HAMILTONIAN MATRIX ----------------'
  write(*,*)
  call writeREALMAT( H(:,:),nDET,nDET )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  if ( mat_spin .eq. 'yes' ) then
     write(*,fmt=200) '---- REAL-SPACE BASIS SPIN SQUARED MATRIX ---------------'
     write(*,*)
     call writeREALMAT( S2(:,:),nDET,nDET )
     write(*,*)
     write(*,fmt=200) '---------------------------------------------------------'
     write(*,*)
  endif
  write(*,fmt=200) '---- DIAGONALIZATION OF THE HAMILTONIAN MATRIX ----------'
  write(*,*)
! call check_MATdiag( H(:,:),CI(:,:),nDET )
  write(*,fmt=250) '- Ground-State Energy [Eo] =', E(1)
                  ctr = 1
                  I   = 2
                  do while ( abs( E(1) - E(I) ) .lt. thrs )
                     ctr = ctr + 1
                     I   = I   + 1
                  enddo
  write(*,fmt=380) '  [Degeneracy: ', ctr, ']'
  write(*,fmt=250) '- Kinetic      Energy  [T] =', K
  write(*,fmt=250) '- Interaction  Energy  [V] =', V
  write(*,*)
  if ( tot_spin .eq. 'yes' ) then
     write(*,fmt=250) '- Total Spin   (GS)  [S^2] =', S
     write(*,*)
  endif
  write(*,fmt=260) '- Ground-State Wave Function [CI Expansion]'
  write(*,fmt=270) '  [Non-zero-amplitude Determinants]'
  write(*,fmt=280) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
            ctr = 0
            do I = 1,nDET,1
               if ( abs(CI(I,1)) .ge. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=290, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=300, advance = 'No') CI(I,1), '   | ' 
                  write(*,fmt=310, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=320, advance = 'No') ' | (x) | '
                  write(*,fmt=310, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=330                ) ' |'
               endif
            enddo
  write(*,*)
  write(*,fmt=340) '  Total number of non-zero-amplitude Determinants :', ctr
  write(*,*)
  write(*,fmt=260) '- Ground-State Wave Function [CI Expansion]'
  write(*,fmt=350) '  [Zero-amplitude Determinants]'
  write(*,fmt=280) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
             ctr = 0
             do I = 1,nDET,1
               if ( abs(CI(I,1)) .lt. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=290, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=300, advance = 'No') CI(I,1), '   | ' 
                  write(*,fmt=310, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=320, advance = 'No') ' | (x) | '
                  write(*,fmt=310, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=330                ) ' |'
               endif
            enddo
  write(*,*)
  write(*,fmt=360) '  Total number of zero-amplitude Determinants (nodes):', ctr
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- ESTIMATION OF THE SIGN-PROBLEM BY CI-VEC. NORM -----'
  write(*,*)
  write(*,fmt=370) '- ||DC||2 [Ground-State] =', normDC
  write(*,*) 

 !-------------------------------------------------!
 ! WRITING OUTPUT FILES                            !
 !-------------------------------------------------!

  open( unit=10, file='hamMAT.out', access='sequential' )
   write(10,*)'---- HAMILTONIAN MATRIX IN REAL-SPACE DETERMINANT BASIS ----'
   write(10,*)
   call matout( nDET,nDET,H(:,:),10 )
  close( unit=10, status='keep')

  open( unit=11, file='gsWFC.out', access='sequential' )
  write(11,*) '# [1] Det. INDEX | [2] Up-Spin Det. INDEX | [3] Down-Spin Det. INDEX | [4] GS-WFC CI Coeff.'
  write(11,*)
  write(11,*) 
   do I = 1,nDET,1
      call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
      write(11,fmt=400) I, Iup, Idown, CI(I,1)
   enddo
  close( unit=11, status='keep')

  open( unit=12, file='spectrum.out', access='sequential' )
  write(12,*) '---- SPRECTRUM OF THE HAMILTONIAN ----'
  write(12,*)
   do I = 1,nDET,1
      write(12,fmt=410) 'E[',I-1,'] =', E(I)
   enddo
  close( unit=12, status='keep')

  200 format(A57)
  210 format(A36, 1I4)
  220 format(A36, 2I2)
  230 format(A36, 1F5.1)
  240 format(A28)
  250 format(A28, 1F15.10)
  260 format(A43)
  270 format(A35)
  280 format(A51)
  290 format(A1,1I4,A3,1I3,A1,1I3,A2)
  300 format(1F10.5,A5)
  310 format(*(1I1))
  320 format(A9)
  330 format(A2)
  340 format(A51,1I4)
  350 format(A31)
  360 format(A54,1I4)
  370 format(A26,1F15.10)
  380 format(A15,1I2,A1)

  400 format(1I18,1I25,1I27,1F21.8)
  410 format(A2,1I3,A3,1F21.8)

  500 format(A32)
  510 format(A39)
  520 format(A38)
  530 format(A35)
  540 format(A11,1F7.3,A9,1F7.3)
  550 format(A11,1F7.3)

  700 format(A40)
  710 format(A42)

  return

END SUBROUTINE write_output_init

!----------------------------------------------------------------------------------------------------------

SUBROUTINE write_output_stoq(CI_stoq,E_stoq,DE_stoq,tot_spin_stoq,S_stoq,nSITES,nDET,nDETup, &
                & nDETdown,upDMAT,downDMAT )

!------------------------------------------------------------!
! SUBROUTINE to write the OUTPUT related to the STOQUASTIZED !
! calculations.                                              !
!------------------------------------------------------------!

 implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter  :: dm   = 2
  double precision, parameter  :: thrs = 10E-8

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!

  integer(8)      , intent(in) :: nSITES
  integer(8)      , intent(in) :: nDET, nDETup, nDETdown

  integer(8)      , intent(in) :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)

  double precision, intent(in) :: CI_stoq(nDET,nDET)
  double precision, intent(in) :: E_stoq(nDET)
  double precision, intent(in) :: DE_stoq
  double precision, intent(in) :: S_stoq

  character(3)    , intent(in) :: tot_spin_stoq

  integer(8)                   :: I, Iup, Idown
  integer(8)                   :: site, ctr

 !-------------------------------------------------!
 ! WRITING STANDARD OUTPUT                         !
 !-------------------------------------------------!

  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- DIAGONALIZATION OF THE STOQ. HAMILTONIAN MATRIX ----'
  write(*,*)
  write(*,fmt=210) '- Ground-State Stoquastized Energy [Eo] =', E_stoq(1)
  write(*,*)
  write(*,fmt=210) '- Stoquastized Gap [DEstoq.]            =', DE_stoq
  write(*,*)
  if ( tot_spin_stoq .eq. 'yes' ) then
     write(*,fmt=210) '- Total Spin (GS) Stoquastized    [S^2] =', S_stoq
     write(*,*)
  endif
  write(*,fmt=220) '- Ground-State Stoquastized Wave Function [CI Expansion]'
  write(*,fmt=230) '  [Non-zero-amplitude Determinants]'
  write(*,fmt=240) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
            ctr = 0
            do I = 1,nDET,1
               if ( abs(CI_stoq(I,1)) .ge. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=250, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=260, advance = 'No') CI_stoq(I,1), '   | ' 
                  write(*,fmt=270, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=280, advance = 'No') ' | (x) | '
                  write(*,fmt=270, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=290                ) ' |'
               endif
            enddo
  write(*,*)
  write(*,fmt=300) '  Total number of non-zero-amplitude Determinants :', ctr
  write(*,*)
  write(*,fmt=220) '- Ground-State Stoquastized Wave Function [CI Expansion]'
  write(*,fmt=310) '  [Zero-amplitude Determinants]'
  write(*,fmt=240) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
             ctr = 0
             do I = 1,nDET,1
               if ( abs(CI_stoq(I,1)) .lt. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=250, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=260, advance = 'No') CI_stoq(I,1), '   | ' 
                  write(*,fmt=270, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=280, advance = 'No') ' | (x) | '
                  write(*,fmt=270, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=290                ) ' |'
               endif
            enddo
   write(*,*)
   write(*,fmt=320) '  Total number of zero-amplitude Determinants (nodes):', ctr

  200 format(A57)
  210 format(A41, 1F15.10)
  220 format(A56)
  230 format(A35)
  240 format(A51)
  250 format(A1,1I4,A3,1I3,A1,1I3,A2)
  260 format(1F10.5,A5)
  270 format(*(1I1))
  280 format(A9)
  290 format(A2)
  300 format(A51,1I4)
  310 format(A31)
  320 format(A54,1I4)

  return

END SUBROUTINE write_output_stoq

!----------------------------------------------------------------------------------------------------------

SUBROUTINE write_output_stoq_tr(CI_stoq,E_stoq,DE_stoq,tot_spin_stoq,S_stoq,nSITES,nDET,nDET_tr,nDETup, &
                & nDETdown,upDMAT,downDMAT )

!------------------------------------------------------------!
! SUBROUTINE to write the OUTPUT related to the STOQUASTIZED !
! calculations.                                              !
!------------------------------------------------------------!

 implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter  :: dm   = 2
  double precision, parameter  :: thrs = 10E-8

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!

  integer(8)      , intent(in) :: nSITES
  integer(8)      , intent(in) :: nDET
  integer(8)      , intent(in) :: nDET_tr
  integer(8)      , intent(in) :: nDETup, nDETdown

  integer(8)      , intent(in) :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)

  double precision, intent(in) :: CI_stoq(nDET,nDET)
  double precision, intent(in) :: E_stoq(nDET_tr)
  double precision, intent(in) :: DE_stoq
  double precision, intent(in) :: S_stoq

  character(3)    , intent(in) :: tot_spin_stoq

  integer(8)                   :: I, Iup, Idown
  integer(8)                   :: site, ctr

 !-------------------------------------------------!
 ! WRITING STANDARD OUTPUT                         !
 !-------------------------------------------------!

  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- DIAGONALIZATION OF THE STOQ. HAMILTONIAN MATRIX ----'
  write(*,*)
  write(*,fmt=210) '- Ground-State Stoquastized Energy [Eo] =', E_stoq(1)
  write(*,*)
  write(*,fmt=210) '- Stoquastized Gap [DEstoq.]            =', DE_stoq
  write(*,*)
  if ( tot_spin_stoq .eq. 'yes' ) then
     write(*,fmt=210) '- Total Spin (GS) Stoquastized    [S^2] =', S_stoq
     write(*,*)
  endif
  write(*,fmt=220) '- Ground-State Stoquastized Wave Function [CI Expansion]'
  write(*,fmt=230) '  [Non-zero-amplitude Determinants]'
  write(*,fmt=240) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
            ctr = 0
            do I = 1,nDET,1
               if ( abs(CI_stoq(I,1)) .ge. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=250, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=260, advance = 'No') CI_stoq(I,1), '   | ' 
                  write(*,fmt=270, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=280, advance = 'No') ' | (x) | '
                  write(*,fmt=270, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=290                ) ' |'
               endif
            enddo
  write(*,*)
  write(*,fmt=300) '  Total number of non-zero-amplitude Determinants :', ctr
  write(*,*)
  write(*,fmt=220) '- Ground-State Stoquastized Wave Function [CI Expansion]'
  write(*,fmt=310) '  [Zero-amplitude Determinants]'
  write(*,fmt=240) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
             ctr = 0
             do I = 1,nDET,1
               if ( abs(CI_stoq(I,1)) .lt. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=250, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=260, advance = 'No') CI_stoq(I,1), '   | ' 
                  write(*,fmt=270, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=280, advance = 'No') ' | (x) | '
                  write(*,fmt=270, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=290                ) ' |'
               endif
            enddo
   write(*,*)
   write(*,fmt=320) '  Total number of zero-amplitude Determinants (nodes):', ctr

  200 format(A57)
  210 format(A41, 1F15.10)
  220 format(A56)
  230 format(A35)
  240 format(A51)
  250 format(A1,1I4,A3,1I3,A1,1I3,A2)
  260 format(1F10.5,A5)
  270 format(*(1I1))
  280 format(A9)
  290 format(A2)
  300 format(A51,1I4)
  310 format(A31)
  320 format(A54,1I4)

  return

END SUBROUTINE write_output_stoq_tr

!----------------------------------------------------------------------------------------------------------

SUBROUTINE write_output_trunc(upDMAT,downDMAT,H_tr,CI_tr,E_tr,noNODESindx, &
           & DE_error,tot_spin_trunc,S_tr,nSITES,nDET,nDETup,nDETdown,nDET_tr) 

!---------------------------------------------------------!
! SUBROUTINE to write the OUTPUT related to the TRUNCATED !
! calculations.                                           !
!---------------------------------------------------------!

 implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  double precision, parameter  :: thrs = 10E-8

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!

  integer(8)      , intent(in) :: nSITES
  integer(8)      , intent(in) :: nDET, nDETup, nDETdown
  integer(8)      , intent(in) :: nDET_tr

  integer(8)      , intent(in) :: upDMAT(nSITES,nDETup), downDMAT(nSITES,nDETdown)
  integer(8)      , intent(in) :: noNODESindx(nDET_tr)

  double precision, intent(in) :: H_tr(nDET_tr,nDET_tr)
  double precision, intent(in) :: CI_tr(nDET_tr,nDET_tr)
  double precision, intent(in) :: E_tr(nDET_tr)
  double precision, intent(in) :: DE_error
  double precision, intent(in) :: S_tr

  character(3)    , intent(in) :: tot_spin_trunc

  integer(8)                   :: I, Iup, Idown
  integer(8)                   :: site, ctr

 !-------------------------------------------------!
 ! WRITING STANDARD OUTPUT                         !
 !-------------------------------------------------!

  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- DIAGONALIZATION OF THE TRUNC. HAMILTONIAN MATRIX ---'
  write(*,*)
  write(*,fmt=210) '- Truncated Space Ground-State Energy [Eo] =', E_tr(1)
                  ctr = 1
                  I   = 2
                  do while ( abs( E_tr(1) - E_tr(I) ) .lt. thrs )
                           ctr = ctr + 1
                           I   = I   + 1
                  enddo
  write(*,fmt=400) '  [Degeneracy: ', ctr, ']'
  write(*,*)
  write(*,fmt=210) '- Energy Error [DEerr.]                    =', DE_error
  write(*,*)
  if ( tot_spin_trunc .eq. 'yes' ) then
     write(*,fmt=210) '- Total Spin (GS) Truncated Space    [S^2] =', S_tr
     write(*,*)
  endif
  write(*,fmt=220) '- Ground-State Approx. Wave Function [CI Expansion]'
  write(*,fmt=230) '  [Non-zero-amplitude Determinants]'
  write(*,fmt=240) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
            ctr = 0
            do I = 1,nDET_tr,1
               if ( abs(CI_tr(I,1)) .ge. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( noNODESindx(I),Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=250, advance = 'No') '[',noNODESindx(I),'] [',Iup,',',Idown,'] '
                  write(*,fmt=260, advance = 'No') CI_tr(I,1), '   | ' 
                  write(*,fmt=270, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=280, advance = 'No') ' | (x) | '
                  write(*,fmt=270, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=290                ) ' |'
               endif
            enddo
  write(*,*)
  write(*,fmt=300) '  Total number of non-zero-amplitude Determinants :', ctr
  write(*,*)
  write(*,fmt=220) '- Ground-State Approx. Wave Function [CI Expansion]'
  write(*,fmt=310) '  [Zero-amplitude Determinants]'
  write(*,fmt=240) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
             ctr = 0
             do I = 1,nDET_tr,1
               if ( abs(CI_tr(I,1)) .lt. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( noNODESindx(I),Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=250, advance = 'No') '[',noNODESindx(I),'] [',Iup,',',Idown,'] '
                  write(*,fmt=260, advance = 'No') CI_tr(I,1), '   | ' 
                  write(*,fmt=270, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=280, advance = 'No') ' | (x) | '
                  write(*,fmt=270, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=290                ) ' |'
               endif
            enddo
   write(*,*)
   write(*,fmt=320) '  Total number of zero-amplitude Determinants (nodes):', ctr

  200 format(A57)
  210 format(A44, 1F15.10)
  220 format(A51)
  230 format(A35)
  240 format(A51)
  250 format(A1,1I4,A3,1I3,A1,1I3,A2)
  260 format(1F10.5,A5)
  270 format(*(1I1))
  280 format(A9)
  290 format(A2)
  300 format(A51,1I4)
  310 format(A31)
  320 format(A54,1I4)

  400 format(A15,1I2,A1)

  return

END SUBROUTINE write_output_trunc

!----------------------------------------------------------------------------------------------------------

SUBROUTINE write_output_fna(upDMAT,downDMAT,KupDMAT,KdownDMAT,CI_K,CIindx,SGN,Vsf,gsVsf,H_eff,CI_eff,    &
                & E_fn,nKPOINTS,nSITES,nELup,nELdown,nDET,nDETup,nDETdown,nTWFCDET,nPDET,subspace_coeff, &
                & fna_relax,tot_spin_fna,E_T,S,normDC)

!----------------------------------------------------------!
! SUBROUTINE to write the OUTPUT related to the FIXED-NODE !
! APPROXIMATION calculations.                              !
!----------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter  :: dm   = 2

  double precision, parameter  :: thrs = 10E-8
  double precision, parameter  :: Vmax = 10E+4

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!

  integer(8)      , intent(in) :: nSITES
  integer(8)      , intent(in) :: nKPOINTS
  integer(8)      , intent(in) :: nELup, nELdown
  integer(8)      , intent(in) :: nDET , nDETup, nDETdown
  integer(8)      , intent(in) :: nTWFCDET
  integer(8)      , intent(in) :: nPDET

  integer(8)      , intent(in) :: upDMAT(nSITES,nDETup)
  integer(8)      , intent(in) :: downDMAT(nSITES,nDETdown)
  integer(8)      , intent(in) :: KupDMAT(nKPOINTS,nDETup)
  integer(8)      , intent(in) :: KdownDMAT(nKPOINTS,nDETdown)

  integer(8)      , intent(in) :: CIindx(nTWFCDET)

  double precision, intent(in) :: E_T
  double precision, intent(in) :: S
  double precision, intent(in) :: normDC
  double precision, intent(in) :: CI_K(nDET)
  double precision, intent(in) :: SGN(nDET,nDET)
  double precision, intent(in) :: Vsf(nDET)
  double precision, intent(in) :: gsVsf(nDET)
  double precision, intent(in) :: H_eff(nDET,nDET)
  double precision, intent(in) :: CI_eff(nDET,nDET)
  double precision, intent(in) :: E_fn(nDET)

  character(3)    , intent(in) :: subspace_coeff
  character(3)    , intent(in) :: fna_relax
  character(3)    , intent(in) :: tot_spin_fna

  integer(8)                   :: Ik, Ikup, Ikdown
  integer(8)                   :: I , Iup , Idown
  integer(8)                   :: kpoint, site, ctr

 !-------------------------------------------------!
 ! WRITING STANDARD OUTPUT                         !
 !-------------------------------------------------!
  
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- FIXED-NODE APPROXIMATION CALCULATIONS --------------'
   write(*,*)
            if ( fna_relax .eq. 'no' ) then
               write(*,fmt=420) '- P/Q RELAXATION [NO]'
            elseif ( fna_relax .eq. 'yes' ) then
               write(*,fmt=430) '- P/Q RELAXATION [YES]'
               write(*,fmt=440) '- Dimension of P-SPACE : ', nPDET
            endif
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- (CI EXPANSION) TRIAL WAVE FUNCTION [K-SPACE] -------'
  write(*,*)
           if ( subspace_coeff .eq. 'no' ) then
              write(*,fmt=410) '- Type of Coefficients : [EXACT]'
           elseif ( subspace_coeff .eq. 'yes' ) then
              write(*,fmt=240) '- Type of Coefficients : [SUBSPACE DIAGONALIZATION]'
           endif
  write(*,*)
            write(*,fmt=240) '  [K] [Kup,Kdown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
            do I = 1,nTWFCDET,1
               Ik = CIindx(I)
               call spinDETindex( Ik,Ikup,Ikdown,nDETup,nDETdown )
               write(*,fmt=250, advance = 'No') '[',Ik,'] [',Ikup,',',Ikdown,'] '
               write(*,fmt=260, advance = 'No') CI_K(Ik), '   | '
               write(*,fmt=270, advance = 'No') ( KupDMAT(kpoint,Ikup)    , kpoint = 1,nKPOINTS,1)
               write(*,fmt=280, advance = 'No') ' | (x) | '
               write(*,fmt=270, advance = 'No') ( KdownDMAT(kpoint,Ikdown), kpoint = 1,nKPOINTS,1)
               write(*,fmt=290                ) ' |'
            enddo
  write(*,*)           
  write(*,fmt=200) '---- REAL-SPACE SIGN-FLIP (VIOLATING) MATRIX ------------'
  write(*,*)
  call writeEREALMAT( SGN(:,:),nDET,nDET )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- FIXED-NODE APPROXIMATION EFFECTIVE HAMILTONIAN -----'
  write(*,*)
  call writeEREALMAT( H_eff(:,:),nDET,nDET )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- DIAGONALIZATION OF THE FNA EFF. HAMILTONIAN MAT. ---'
  write(*,*)
  write(*,fmt=340) '- Trial Wave Function Energy          [ET] =', E_T
  write(*,*)
  write(*,fmt=210) '- Fixed-Node App. GS Energy           [Eo] =', E_fn(1)
                  ctr = 1
                  I   = 2
                  do while ( abs( E_fn(1) - E_fn(I) ) .lt. thrs )
                           ctr = ctr + 1
                           I   = I   + 1
                  enddo
  write(*,fmt=400) '  [Degeneracy: ', ctr, ']'
  write(*,*)
  if ( tot_spin_fna .eq. 'yes' ) then
     write(*,fmt=340) '- Total Spin (GS) Fixed-Node App.    [S^2] =', S
     write(*,*)
  endif
  write(*,fmt=220) '- Ground-State   FNA   Wave Function [CI Expansion]'
  write(*,fmt=230) '  [Non-zero-amplitude Determinants]'
  write(*,fmt=240) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
            ctr = 0
            do I = 1,nDET,1
               if ( abs(CI_eff(I,1)) .ge. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=250, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=260, advance = 'No') CI_eff(I,1), '   | '
                  write(*,fmt=270, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=280, advance = 'No') ' | (x) | '
                  write(*,fmt=270, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=290                ) ' |'
               endif
            enddo
  write(*,*)
  write(*,fmt=300) '  Total number of non-zero-amplitude Determinants :', ctr
  write(*,*)
  write(*,fmt=220) '- Ground-State   FNA   Wave Function [CI Expansion]'
  write(*,fmt=310) '  [Zero-amplitude Determinants]'
  write(*,fmt=240) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
             ctr = 0
             do I = 1,nDET,1
               if ( abs(CI_eff(I,1)) .lt. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=250, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=260, advance = 'No') CI_eff(I,1), '   | '
                  write(*,fmt=270, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=280, advance = 'No') ' | (x) | '
                  write(*,fmt=270, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=290                ) ' |'
               endif
            enddo
  write(*,*)
  write(*,fmt=320) '  Total number of zero-amplitude Determinants (nodes):', ctr
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '-- ESTIMATION OF THE FNA SIGN-PROBLEM BY CI-VEC. NORM ---'
  write(*,*)
  write(*,fmt=370) '- ||DC||2 [FNA-GS] =', normDC

  !-------------------------------------------------!
  ! WRITING OUTPUT FILES                            !
  !-------------------------------------------------!

   open( unit=10, file='VsfPOT.out', access='sequential' )
   write(10,*) '# [1] Det. INDEX | [2] Up-Spin Det. INDEX | [3] Down-Spin Det. INDEX | [4] Vsf[TR]              | [5] Vsf[FCI].'

   write(10,*)
   write(10,*)
    do I = 1,nDET,1
       call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
       if ( Vsf(I) .gt. Vmax ) then
          write(10,fmt=500) I, Iup, Idown, Vmax, gsVsf(I)
       else
          write(10,fmt=500) I, Iup, Idown, Vsf(I), gsVsf(I)
       endif
    enddo
  close( unit=10, status='keep')

  200 format(A57)
  210 format(A44, 1E15.5)
  220 format(A51)
  230 format(A35)
  240 format(A51)
  250 format(A1,1I4,A3,1I3,A1,1I3,A2)
  260 format(1F10.5,A5)
  270 format(*(1I1))
  280 format(A9)
  290 format(A2)
  300 format(A51,1I4)
  310 format(A31)
  320 format(A54,1I4)
  330 format(A25)
  340 format(A44, 1F15.10)
  370 format(A20,1F15.10)
  400 format(A15,1I4,A1)
  410 format(A32)
  420 format(A21)
  430 format(A22)
  440 format(A25,1I4)

  500 format(1I18,1I25,1I27,1F21.8,1X,1F21.8)

  return

END SUBROUTINE write_output_fna


!----------------------------------------------------------------------------------------------------------

SUBROUTINE write_output_fna_tr(upDMAT,downDMAT,KupDMAT,KdownDMAT,CI_K,CIindx,SGN,Vsf,gsVsf,Heff_tr,CI_eff,       &
                & E_fn,nKPOINTS,nSITES,nELup,nELdown,nDET,nDET_tr,nDETup,nDETdown,nTWFCDET,nPDET,subspace_coeff, &
                & fna_relax,tot_spin_fna,E_T,S,normDC )

!----------------------------------------------------------!
! SUBROUTINE to write the OUTPUT related to the FIXED-NODE !
! APPROXIMATION calculations [with TRUNCATED DIAG.]        !
!----------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter  :: dm   = 2

  double precision, parameter  :: thrs = 10E-8
  double precision, parameter  :: Vmax = 10E+4

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!

  integer(8)      , intent(in) :: nSITES
  integer(8)      , intent(in) :: nKPOINTS
  integer(8)      , intent(in) :: nELup , nELdown
  integer(8)      , intent(in) :: nDET  , nDET_tr
  integer(8)      , intent(in) :: nDETup, nDETdown
  integer(8)      , intent(in) :: nTWFCDET
  integer(8)      , intent(in) :: nPDET

  integer(8)      , intent(in) :: upDMAT(nSITES,nDETup)
  integer(8)      , intent(in) :: downDMAT(nSITES,nDETdown)
  integer(8)      , intent(in) :: KupDMAT(nKPOINTS,nDETup)
  integer(8)      , intent(in) :: KdownDMAT(nKPOINTS,nDETdown)

  integer(8)      , intent(in) :: CIindx(nTWFCDET)

  double precision, intent(in) :: E_T
  double precision, intent(in) :: S
  double precision, intent(in) :: normDC
  double precision, intent(in) :: CI_K(nDET)
  double precision, intent(in) :: SGN(nDET,nDET)
  double precision, intent(in) :: Vsf(nDET)
  double precision, intent(in) :: gsVsf(nDET)
  double precision, intent(in) :: Heff_tr(nDET_tr,nDET_tr)
  double precision, intent(in) :: CI_eff(nDET)
  double precision, intent(in) :: E_fn(nDET_tr)

  character(3)    , intent(in) :: subspace_coeff
  character(3)    , intent(in) :: fna_relax
  character(3)    , intent(in) :: tot_spin_fna

  
  integer(8)                   :: Ik, Ikup, Ikdown
  integer(8)                   :: I , Iup , Idown
  integer(8)                   :: kpoint, site, ctr

 !-------------------------------------------------!
 ! WRITING STANDARD OUTPUT                         !
 !-------------------------------------------------!

  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- FIXED-NODE APPROXIMATION CALCULATIONS --------------'
  write(*,*)
            if ( fna_relax .eq. 'no' ) then
               write(*,fmt=420) '- P/Q RELAXATION [NO]'
            elseif ( fna_relax .eq. 'yes' ) then
               write(*,fmt=430) '- P/Q RELAXATION [YES]'
               write(*,fmt=440) '- Dimension of P-SPACE : ', nPDET
            endif
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- (CI EXPANSION) TRIAL WAVE FUNCTION [K-SPACE] -------'
  write(*,*)
           if ( subspace_coeff .eq. 'no' ) then
              write(*,fmt=410) '- Type of Coefficients : [EXACT]'
           elseif ( subspace_coeff .eq. 'yes' ) then
              write(*,fmt=240) '- Type of Coefficients : [SUBSPACE DIAGONALIZATION]'
           endif
  write(*,*)
            write(*,fmt=240) '  [K] [Kup,Kdown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
            do I = 1,nTWFCDET,1
               Ik = CIindx(I)
               call spinDETindex( Ik,Ikup,Ikdown,nDETup,nDETdown )
               write(*,fmt=250, advance = 'No') '[',Ik,'] [',Ikup,',',Ikdown,'] '
               write(*,fmt=260, advance = 'No') CI_K(Ik), '   | '
               write(*,fmt=270, advance = 'No') ( KupDMAT(kpoint,Ikup)    , kpoint = 1,nKPOINTS,1)
               write(*,fmt=280, advance = 'No') ' | (x) | '
               write(*,fmt=270, advance = 'No') ( KdownDMAT(kpoint,Ikdown), kpoint = 1,nKPOINTS,1)
               write(*,fmt=290                ) ' |'
            enddo
  write(*,*)
  write(*,fmt=200) '---- REAL-SPACE SIGN-FLIP (VIOLATING) MATRIX ------------'
  write(*,*)
  call writeEREALMAT( SGN(:,:),nDET,nDET )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- FIXED-NODE APP. TRUNCATED EFFECTIVE HAMILTONIAN ----'
  write(*,*)
  call writeEREALMAT( Heff_tr(:,:),nDET_tr,nDET_tr )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- DIAGONALIZATION OF THE FNA EFF. HAMILTONIAN MAT. ---'
  write(*,*)
  write(*,fmt=340) '- Trial Wave Function Energy          [ET] =', E_T
  write(*,*)
  write(*,fmt=210) '- Fixed-Node App. GS Energy           [Eo] =', E_fn(1)
                  ctr = 1
                  I   = 2
                  do while ( abs( E_fn(1) - E_fn(I) ) .lt. thrs )
                           ctr = ctr + 1
                           I   = I   + 1
                  enddo
  write(*,fmt=400) '  [Degeneracy: ', ctr, ']'
  write(*,*)
  if ( tot_spin_fna .eq. 'yes' ) then
     write(*,fmt=340) '- Total Spin (GS) Fixed-Node App.    [S^2] =', S
     write(*,*)
  endif
  write(*,fmt=220) '- Ground-State   FNA   Wave Function [CI Expansion]'
  write(*,fmt=230) '  [Non-zero-amplitude Determinants]'
  write(*,fmt=240) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
            ctr = 0
            do I = 1,nDET,1
               if ( abs(CI_eff(I)) .ge. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=250, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=260, advance = 'No') CI_eff(I), '   | '
                  write(*,fmt=270, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=280, advance = 'No') ' | (x) | '
                  write(*,fmt=270, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=290                ) ' |'
               endif
            enddo
  write(*,*)
  write(*,fmt=300) '  Total number of non-zero-amplitude Determinants :', ctr
  write(*,*)
  write(*,fmt=220) '- Ground-State   FNA   Wave Function [CI Expansion]'
  write(*,fmt=310) '  [Zero-amplitude Determinants]'
  write(*,fmt=240) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
             ctr = 0
             do I = 1,nDET,1
               if ( abs(CI_eff(I)) .lt. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=250, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=260, advance = 'No') CI_eff(I), '   | '
                  write(*,fmt=270, advance = 'No') ( upDMAT(site,Iup)    , site = 1,nSITES,1)
                  write(*,fmt=280, advance = 'No') ' | (x) | '
                  write(*,fmt=270, advance = 'No') ( downDMAT(site,Idown), site = 1,nSITES,1)
                  write(*,fmt=290                ) ' |'
               endif
            enddo
  write(*,*)
  write(*,fmt=320) '  Total number of zero-amplitude Determinants (nodes):', ctr
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '-- ESTIMATION OF THE FNA SIGN-PROBLEM BY CI-VEC. NORM ---'
  write(*,*)
  write(*,fmt=370) '- ||DC||2 [FNA-GS] =', normDC

  !-------------------------------------------------!
  ! WRITING OUTPUT FILES                            !
  !-------------------------------------------------!

   open( unit=10, file='VsfPOT.out', access='sequential' )
   write(10,*) '# [1] Det. INDEX | [2] Up-Spin Det. INDEX | [3] Down-Spin Det. INDEX | [4] Vsf[TR]              | [5] Vsf[FCI].'
   write(10,*)
   write(10,*)
    do I = 1,nDET,1
       call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
       if ( Vsf(I) .gt. Vmax ) then
          write(10,fmt=500) I, Iup, Idown, Vmax, gsVsf(I)
       else
          write(10,fmt=500) I, Iup, Idown, Vsf(I), gsVsf(I)
       endif
    enddo
  close( unit=10, status='keep')

  200 format(A57)
  210 format(A44, 1E15.5)
  220 format(A51)
  230 format(A35)
  240 format(A51)
  250 format(A1,1I4,A3,1I3,A1,1I3,A2)
  260 format(1F10.5,A5)
  270 format(*(1I1))
  280 format(A9)
  290 format(A2)
  300 format(A51,1I4)
  310 format(A31)
  320 format(A54,1I4)
  330 format(A25)
  340 format(A44, 1F15.10)
  370 format(A20,1F15.10)
  400 format(A15,1I4,A1)
  410 format(A32)
  420 format(A21)
  430 format(A22)
  440 format(A25,1I4)

  500 format(1I18,1I25,1I27,1F21.8,1X,1F21.8)

  return

END SUBROUTINE write_output_fna_tr

!----------------------------------------------------------------------------------------------------------

SUBROUTINE write_output_kspace( Hk,CIk,Ek,KupDMAT,KdownDMAT,KPOINTS,KP_ENRG, &
          & nKPOINTS,nEL,nELup,nELdown,nDET,nDETup,nDETdown,t,U )

!--------------------------------------------------------------------------------------------------!
! SUBROUTINE to write the output for the calculations of the 2D Fermi-Hubbard model in real space. !
!--------------------------------------------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter  :: dm   = 2
  double precision, parameter  :: thrs = 10E-8

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!

  integer(8)      , intent(in) :: nKPOINTS
  integer(8)      , intent(in) :: nEL, nELup, nELdown
  integer(8)      , intent(in) :: nDET, nDETup, nDETdown

  integer(8)      , intent(in) :: KupDMAT(nKPOINTS,nDETup)
  integer(8)      , intent(in) :: KdownDMAT(nKPOINTS,nDETdown)

  double precision, intent(in) :: t
  double precision, intent(in) :: U

  double precision, intent(in) :: KPOINTS(dm,nKPOINTS)
  double precision, intent(in) :: KP_ENRG(nKPOINTS)

  double precision, intent(in) :: Hk(nDET,nDET)
  double precision, intent(in) :: CIk(nDET,nDET)
  double precision, intent(in) :: Ek(nDET)

  integer(8)                   :: I, Iup, Idown
  integer(8)                   :: kpoint, ctr

 !-------------------------------------------------!
 ! WRITING STANDARD OUTPUT                         !
 !-------------------------------------------------!

  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,*)
  write(*,fmt=200) '----- 2-DIMENSIONAL FERMI-HUBBARD MODEL IN K-SPACE ------'
  write(*,*)
  write(*,fmt=200) '---- SYSTEM PROPERTIES ----------------------------------'
  write(*,*)
  write(*,fmt=210) '- Total Number of K-Points (Nk)    : ', nKPOINTS
  write(*,*)
  write(*,fmt=220) '- Hopping Parameter   (t) for H[k] : ', t
  write(*,fmt=220) '- On-Site Interaction (U) for H[k] : ', U
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- LIST OF (SORTED) K-POINTS [CARTESIAN COORD.] -------'
  write(*,*)
  call writeREALMAT( KPOINTS(:,:),dm,nKPOINTS )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- LIST OF (SORTED) K-POINTS ENERGIES -----------------'
  write(*,*)
           do kpoint = 1,nKPOINTS,1
              write(*,fmt=230) '- Ek[ ',kpoint,' ] = ', KP_ENRG(kpoint)
           enddo
  write(*,*)
  write(*,fmt=200) '---- K-SPACE DETERMINANTAL MATRICES ---------------------'
  write(*,*)
  write(*,fmt=240) '----  UP-SPIN K-SUBSPACE ---'
  write(*,*)
  call writeINTMAT( KupDMAT(:,:),nKPOINTS,nDETup )
  write(*,*)
  write(*,fmt=240) '---- DOWN-SPIN K-SUBSPACE --'
  write(*,*)
  call writeINTMAT( KdownDMAT(:,:),nKPOINTS,nDETdown )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- K-SPACE BASIS HAMILTONIAN MATRIX -------------------'
  write(*,*)
  call writeREALMAT( Hk(:,:),nDET,nDET )
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'
  write(*,*)
  write(*,fmt=200) '---- DIAGONALIZATION OF THE HAMILTONIAN [K] MATRIX ------'
  write(*,*)
  write(*,fmt=250) '- Ground-State K-Space Energy [Eo] =', Ek(1)
                  ctr = 1
                  I   = 2
                  do while ( abs( Ek(1) - Ek(I) ) .lt. thrs )
                     ctr = ctr + 1
                     I   = I   + 1
                  enddo
  write(*,fmt=380) '  [Degeneracy: ', ctr, ']'
  write(*,*)
  write(*,fmt=260) '- Ground-State K-Space Wave Function [CI Expansion]'
  write(*,fmt=270) '  [Non-zero-amplitude Determinants]'
  write(*,fmt=280) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
            ctr = 0
            do I = 1,nDET,1
               if ( abs(CIk(I,1)) .ge. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=290, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=300, advance = 'No') CIk(I,1), '   | '
                  write(*,fmt=310, advance = 'No') ( KupDMAT(kpoint,Iup)    , kpoint = 1,nKPOINTS,1)
                  write(*,fmt=320, advance = 'No') ' | (x) | '
                  write(*,fmt=310, advance = 'No') ( KdownDMAT(kpoint,Idown), kpoint = 1,nKPOINTS,1)
                  write(*,fmt=330                ) ' |'
               endif
            enddo
  write(*,*)
  write(*,fmt=340) '  Total number of non-zero-amplitude Determinants :', ctr
  write(*,*)
  write(*,fmt=260) '- Ground-State K-Space Wave Function [CI Expansion]'
  write(*,fmt=350) '  [Zero-amplitude Determinants]'
  write(*,fmt=280) '  [I] [Iup,Idown] C x { |Spin-Up| (x) |Spin-Down| }'
  write(*,*)
             ctr = 0
             do I = 1,nDET,1
               if ( abs(CIk(I,1)) .lt. thrs ) then
                  ctr = ctr + 1
                  call spinDETindex( I,Iup,Idown,nDETup,nDETdown )
                  write(*,fmt=290, advance = 'No') '[',I,'] [',Iup,',',Idown,'] '
                  write(*,fmt=300, advance = 'No') CIk(I,1), '   | '
                  write(*,fmt=310, advance = 'No') ( KupDMAT(kpoint,Iup)    , kpoint = 1,nKPOINTS,1)
                  write(*,fmt=320, advance = 'No') ' | (x) | '
                  write(*,fmt=310, advance = 'No') ( KdownDMAT(kpoint,Idown), kpoint = 1,nKPOINTS,1)
                  write(*,fmt=330                ) ' |'
               endif
            enddo
  write(*,*)
  write(*,fmt=360) '  Total number of zero-amplitude Determinants (nodes):', ctr
  write(*,*)
  write(*,fmt=200) '---------------------------------------------------------'

  200 format(A57)
  210 format(A37, 1I4)
  220 format(A37, 1F5.1)
  230 format(A6 ,1I2,A5,1F5.1)
  240 format(A28)
  250 format(A36, 1F15.10)
  260 format(A51)
  270 format(A35)
  280 format(A51)
  290 format(A1,1I4,A3,1I3,A1,1I3,A2)
  300 format(1F10.5,A5)
  310 format(*(1I1))
  320 format(A9)
  330 format(A2)
  340 format(A51,1I4)
  350 format(A31)
  360 format(A54,1I4)
  370 format(A26,1F15.10)
  380 format(A15,1I2,A1)

  return

END SUBROUTINE write_output_kspace

!----------------------------------------------------------------------------------------------------------
