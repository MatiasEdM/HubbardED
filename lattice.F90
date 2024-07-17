!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_connectMAT(LMAT,BoundCond,nSITES,nSITESx,nSITESy)

!-----------------------------------------------------------------!
! SUBROUTINE to generate the CONNECTIVITY MATRIX that defines the !
! topology of the lattice, according to the following indexing:   !
!  Nx+Ny      ...      Nx*Ny                                      !
!  ...  ...   ...  ... ...                                        !
!  Nx+1 Nx+2  Nx+3 ... 2*Nx                                       !
!  1    2     3    ... Nx                                         !
!                                                                 !
! For example a Nx=3, Ny=2 (3x2) lattice:                         !
!  4 5 6                                                          !
!  1 2 3                                                          !
!-----------------------------------------------------------------!

 implicit none

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!

  character(5), intent(in)  :: BoundCond

  integer(8)  , intent(in)  :: nSITES, nSITESx, nSITESy 
  integer(8)  , intent(out) :: LMAT(nSITES,nSITES)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)              :: i

 !-------------------------------------------------!
 ! SUBROUTINE: LMAT GENERATION                     !
 !-------------------------------------------------!

 LMAT = 0
 !------------------------------------------------------!
 ! LMAT(i,j) = LMAT(j,i) is a symmetric matrix, so it   !
 ! is only necessary to run over rows and asign at the  !
 ! same time the values of (i,j) and (j,i), and only    !
 ! up- and right neighbors or down and left neighbors.  !
 ! This is because the by symmetry, node if node (n) is !
 ! the right neighbor of (m), (m) is the left neighbor  !
 ! of (n) and the same for up and down. The latter is   !
 ! taken into account by the symmetric assignment.      !
 !------------------------------------------------------!

 do i = 1, nSITES, 1
 
    !----------------------------------------------------!
    ! Right-neighbor only for lattice points that do not !
    ! blong to the right border (index is not a multiple !
    ! of nSITESx.                                        !
    !----------------------------------------------------!

    if ( mod(i,nSITESx) /= 0 ) then
       LMAT(i, i + 1) = 1
       LMAT(i + 1 ,i) = 1
    endif
 
    !-----------------------------------------------------!
    ! Above-neighbor only for lattice points that do not  !
    ! belong to the last lattice-row, and for the lattice !
    ! point (i,j) the neighbor is (i,j+nSITESx)           !
    !-----------------------------------------------------!

    if ( i <= nSITESx*(nSITESy-1) ) then
       LMAT(i, i + nSITESx) = 1
       LMAT(i + nSITESx, i) = 1
    endif

 enddo

 !--------------------------------------------!
 ! Additional connections in the case of PBC. !
 !--------------------------------------------!

 if ( BoundCond .eq. 'pbcxy' ) then
    
    if ( nSITESx .eq. 2) then 
       do i = nSITESx,nSITES,1
          if ( mod(i,nSITESx) .eq. 0 ) then
             LMAT(i, i - nSITESx + 1) = 2
             LMAT(i - nSITESx + 1, i) = 2
          endif
       enddo
    elseif ( nSITESx .ne. 2) then
       do i = nSITESx,nSITES,1
          if ( mod(i,nSITESx) .eq. 0 ) then
             LMAT(i, i - nSITESx + 1) = 1
             LMAT(i - nSITESx + 1, i) = 1
          endif
       enddo
    endif

    if ( nSITESy .eq. 2) then
       do i = 1,nSITESx,1
          LMAT(i, i + nSITESx*(nSITESy-1)) = 2
          LMAT(i + nSITESx*(nSITESy-1), i) = 2
       enddo
    elseif ( nSITESy .ne. 2) then
       do i = 1,nSITESx,1
          LMAT(i, i + nSITESx*(nSITESy-1)) = 1
          LMAT(i + nSITESx*(nSITESy-1), i) = 1
       enddo
    endif

 elseif ( BoundCond .eq. 'pbcx' ) then

    if ( nSITESx .eq. 2) then
       do i = nSITESx,nSITES,1
          if ( mod(i,nSITESx) .eq. 0 ) then
             LMAT(i, i - nSITESx + 1) = 2
             LMAT(i - nSITESx + 1, i) = 2
          endif
       enddo
    elseif ( nSITESx .ne. 2) then
       do i = nSITESx,nSITES,1
          if ( mod(i,nSITESx) .eq. 0 ) then
             LMAT(i, i - nSITESx + 1) = 1
             LMAT(i - nSITESx + 1, i) = 1
          endif
       enddo
    endif

 elseif ( BoundCond .eq. 'pbcy' ) then

    if ( nSITESy .eq. 2) then
       do i = 1,nSITESx,1
          LMAT(i, i + nSITESx*(nSITESy-1)) = 2
          LMAT(i + nSITESx*(nSITESy-1), i) = 2
       enddo
    elseif ( nSITESy .ne. 2) then
       do i = 1,nSITESx,1
          LMAT(i, i + nSITESx*(nSITESy-1)) = 1
          LMAT(i + nSITESx*(nSITESy-1), i) = 1
       enddo
    endif

 endif

 return

END SUBROUTINE generate_connectMAT

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_connectLS(VECT,node,nSITESx,nSITESy)

!----------------------------------------------------------------!
! SUBROUTINE that given a certain node it returns a vector       !
! of length max. 4 (for a 2D rectangular lattice) which contains !
! the neighboring lattice sites or orbitals.                     !
!----------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!

  integer(8), intent(in)  :: nSITESx, nSITESy
  integer(8), intent(in)  :: node
  integer(8), intent(out) :: VECT(4)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)              :: i

 !-------------------------------------------------!
 ! SUBROUTINE: FIRST NEIGHBORS GENERATION          !
 !-------------------------------------------------!

  VECT(:) = 0

  if ( node .eq. 1) then
     VECT(1) = 2
     VECT(2) = nSITESx + 1
  elseif ( node .eq. nSITESx ) then
     VECT(1) = nSITESx - 1
     VECT(2) = 2 * nSITESx
  elseif ( node .eq. ( nSITESx*(nSITESy-1) + 1 ) ) then
     VECT(1) = nSITESx*(nSITESy-2) + 1
     VECT(2) = nSITESx*(nSITESy-1) + 2
  elseif ( node .eq. (nSITESx * nSITESy) ) then
     VECT(1) = nSITESx*(nSITESy-1)
     VECT(2) = (nSITESx * nSITESy) - 1
  elseif ( (node .gt. 1) .and. (node .lt. nSITESx) ) then
     VECT(1) = node - 1
     VECT(2) = node + 1
     VECT(3) = node + nSITESx
  elseif ( (node .gt. (nSITESx*(nSITESy-1) + 1)) .and. (node .lt. (nSITESx * nSITESy)) ) then
     VECT(1) = node - nSITESx
     VECT(2) = node - 1
     VECT(3) = node + 1
  elseif ( (mod(node,nSITESx) .eq. 0) .and. (node .ne. nSITESx) .and. (node .ne. nSITESx*nSITESy) ) then
     VECT(1) = node - nSITESx
     VECT(2) = node - 1
     VECT(3) = node + nSITESx
  elseif ( (mod(node-1,nSITESx) .eq. 0) .and. (node .ne. 1) .and. (node .ne. (nSITESx*(nSITESy-1)+1)) ) then
     VECT(1) = node - nSITESx
     VECT(2) = node + 1
     VECT(3) = node + nSITESx
  else
     VECT(1) = node - nSITESx
     VECT(2) = node - 1
     VECT(3) = node + 1
     VECT(4) = node + nSITESx
   endif

   return

END SUBROUTINE generate_connectLS

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_lattPOS(lattPOS,lattVEC,nSITESx,nSITESy,nSITES)

!----------------------------------------------------------------------------------------!
! SUBTOUINE to generate the LIST of POSITIONS of LATTICE POINTS in CARTESIAN COORDINATES !
! following the indexing:                                                                !
!  Nx+Ny      ...      Nx*Ny                                                             !
!  ...  ...   ...  ... ...                                                               !
!  Nx+1 Nx+2  Nx+3 ... 2*Nx                                                              !
!  1    2     3    ... Nx                                                                !
!                                                                                        !
! For example a Nx=3, Ny=2 (3x2) lattice:                                                !
!  4 5 6                                                                                 !
!  1 2 3                                                                                 !
!----------------------------------------------------------------------------------------!

 implicit none

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!
 
  integer(8)      , parameter :: dm = 2 

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!

  integer(8)      , intent(in) :: nSITES, nSITESx, nSITESy 
  
  double precision, intent(in) :: lattVEC(dm,dm)

  double precision, intent(out) :: lattPOS(dm,nSITES)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                    :: i , j
  integer(8)                    :: n1, n2
  integer(8)                    :: site

 !-------------------------------------------------!
 ! SUBROUTINE: LATTICE POSITION GENERATION         !
 !-------------------------------------------------!

  lattPOS(:,:) = 0.d0

  site         = 1

  !CARTESIAN COORDINATES R = Rx*Ex + Ry*Ey --------!

   do j = 1,nSITESy,1
     do i = 1,nSITESx,1
        n1 = i - 1
        n2 = j - 1
        lattPOS(1,site) = n1 * lattVEC(1,1) + n2 * lattVEC(1,2)
        lattPOS(2,site) = n1 * lattVEC(2,1) + n2 * lattVEC(2,2)
        site = site + 1
     enddo
  enddo

  return

END SUBROUTINE generate_lattPOS
 
!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_DETPOS(R,lattPOS,DET,nSITES,nDET,nEL)

!------------------------------------------------------!
! SUBROUTINE to generate the list |R> of positions in  !
! CARTESIAN COORDINATES of a given configuration (DET) !
! of up- or down-electrons. Transforms the DET |D> in  !
! OCC. NUMBER REPRESENTATION into POSITIONS.           !
!------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter   :: dm  = 2

 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!

  integer(8)      , intent(in)  :: nSITES, nDET, nEL
  integer(8)      , intent(in)  :: DET(nSITES)

  double precision, intent(in)  :: lattPOS(dm,nSITES)


  double precision, intent(out) :: R(dm,nEL)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                    :: site
  integer(8)                    :: electron

 !-------------------------------------------------!
 ! SUBROUTINE: OCC. NUMB. TO POS. TRANSFORMATION   !
 !-------------------------------------------------!

  R(:,:)   = 0.d0
  electron = 1

  do site = 1,nSITES,1
     if ( DET(site) .eq. 1 ) then
        R(1,electron)   = lattPOS(1,site) 
        R(2,electron)   = lattPOS(2,site)
        electron        = electron + 1
      endif
  enddo

  return

END SUBROUTINE get_DETPOS

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_reciprocalVEC(recVEC,lattVEC)

!------------------------------------------------------------------!
! SUBROUTINE to generate the RECIPROCAL LATTICE VECTORS {b1, b2}   !
! from the LATTICE VECTORS {a1,a2} as:                             !
! b[i] = ( 2*pi / A ) * ( a[j] x a[k] )                            !
! with A = a1 (a2 x a3) the area of the two-dimensional unit cell. !
!------------------------------------------------------------------!

 implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter   :: dm2  = 2
  integer(8)      , parameter   :: dm3  = 3
  double precision, parameter   :: pi   = dacos(-1.d0)


 !-------------------------------------------------!
 ! VARIABLES                                       !
 !-------------------------------------------------!

  double precision, intent(in)  :: lattVEC(dm2,dm2)
  double precision, intent(out) :: recVEC(dm2,dm2) 

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!
 
  double precision, allocatable :: a1(:), a2(:), a3(:)
  double precision, allocatable :: b1(:), b2(:), b3(:)
 
  double precision              :: volume_UnitCell, A

  integer(8)                    :: i

 !-------------------------------------------------!
 ! SUBROUTINE: RECIPROCAL VECTORS GENERATION       !
 !-------------------------------------------------!

  allocate( a1(dm3), a2(dm3), a3(dm3) )
  allocate( b1(dm3), b2(dm3), b3(dm3) )

  a1(1) = lattVEC(1,1) ; a1(2) = lattVEC(2,1) ; a1(3) = 0.d0
  a2(1) = lattVEC(1,2) ; a2(2) = lattVEC(2,2) ; a2(3) = 0.d0
  a3(1) = 0.d0         ; a3(2) = 0.d0         ; a3(3) = 1.d0

  b1(:) = 0.d0
  b2(:) = 0.d0
  b3(:) = 0.d0

  A     = volume_UnitCell( a1(:),a2(:),a3(:) )

  call cross_product( b1(:),a2(:),a3(:) )
  call cross_product( b2(:),a3(:),a1(:) )
  call cross_product( b3(:),a1(:),a2(:) )

  recVEC(:,:) = 0.d0
  
  recVEC(1,1) = b1(1) ; recVEC(2,1) = b1(2)
  recVEC(1,2) = b2(1) ; recVEC(2,2) = b2(2)
    
  recVEC(:,:) = ( (2*pi)/(A) ) * recVEC(:,:)

  deallocate( a1, a2, a3 )
  deallocate( b1, b2, b3 )

  return

END SUBROUTINE generate_reciprocalVEC

!----------------------------------------------------------------------------------------------------------

DOUBLE PRECISION FUNCTION volume_UnitCell(v1,v2,v3)

!-----------------------------------------------------------!
! FUNCTION to calculate the VOLUME of a UNIT CELL given the !
! lattice vectors.                                          !
!-----------------------------------------------------------!

  integer(8)      , parameter  :: dm3 = 3

  double precision, intent(in) :: v1(dm3), v2(dm3), v3(dm3)
  
  double precision             :: aux(dm3)

  call cross_product( aux(:),v2(:),v3(:) )

  volume_UnitCell = abs( dot_product( v1(:),aux(:) ))

  return

END FUNCTION volume_UnitCell 

!----------------------------------------------------------------------------------------------------------

SUBROUTINE generate_KPOINTS(KPOINTS,recVEC,nSITESx,nSITESy,mMAXx,mMAXy,nKPOINTS)

!-----------------------------------------------------!
! SUBROUTINE to generate an UNSORTED LIST of KPOINTS. !
!-----------------------------------------------------!

 implicit none
 
 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter    :: dm = 2

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)      , intent(in)  :: nSITESx, nSITESy
  integer(8)      , intent(in)  :: mMAXx  , mMAXy
  integer(8)      , intent(in)  :: nKPOINTS

  double precision, intent(in)  :: recVEC(dm,dm)
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

  !CARTESIAN COORDINATES R = Rx*Ex + Ry*Ey --------!

  if     ( ( mod(nSITESx,2).eq.0 ).and.( mod(nSITESy,2).eq.0 ) ) then
         
         do m1 = -(mMAXx) + 1,mMAXx,1
            do m2 = -(mMAXy) + 1,mMAXy,1
               KPOINTS(1,kpoint) = (dble(m1)/dble(nSITESx))*recVEC(1,1) + (dble(m2)/dble(nSITESy))*recVEC(1,2)
               KPOINTS(2,kpoint) = (dble(m1)/dble(nSITESx))*recVEC(2,1) + (dble(m2)/dble(nSITESy))*recVEC(2,2)
               kpoint            = kpoint + 1
            enddo
         enddo

  elseif ( ( mod(nSITESx,2).eq.0 ).and.( mod(nSITESy,2).ne.0 ) ) then

         do m1 = -(mMAXx) + 1,mMAXx,1
            do m2 = -(mMAXy),mMAXy,1
               KPOINTS(1,kpoint) = (dble(m1)/dble(nSITESx))*recVEC(1,1) + (dble(m2)/dble(nSITESy))*recVEC(1,2)
               KPOINTS(2,kpoint) = (dble(m1)/dble(nSITESx))*recVEC(2,1) + (dble(m2)/dble(nSITESy))*recVEC(2,2)
             kpoint            = kpoint + 1
            enddo
         enddo

  elseif ( ( mod(nSITESx,2).ne.0 ).and.( mod(nSITESy,2).eq.0 ) ) then

         do m1 = -(mMAXx),mMAXx,1
            do m2 = -(mMAXy) + 1,mMAXy,1
               KPOINTS(1,kpoint) = (dble(m1)/dble(nSITESx))*recVEC(1,1) + (dble(m2)/dble(nSITESy))*recVEC(1,2)
               KPOINTS(2,kpoint) = (dble(m1)/dble(nSITESx))*recVEC(2,1) + (dble(m2)/dble(nSITESy))*recVEC(2,2)
               kpoint            = kpoint + 1
            enddo
         enddo

  elseif ( ( mod(nSITESx,2).ne.0 ).and.( mod(nSITESy,2).ne.0 ) ) then

         do m1 = -(mMAXx),mMAXx,1
            do m2 = -(mMAXy),mMAXy,1
               KPOINTS(1,kpoint) = (dble(m1)/dble(nSITESx))*recVEC(1,1) + (dble(m2)/dble(nSITESy))*recVEC(1,2)
               KPOINTS(2,kpoint) = (dble(m1)/dble(nSITESx))*recVEC(2,1) + (dble(m2)/dble(nSITESy))*recVEC(2,2)
               kpoint            = kpoint + 1
            enddo
         enddo

  endif

  return

END SUBROUTINE generate_KPOINTS

!----------------------------------------------------------------------------------------------------------

SUBROUTINE sort_KPOINTS(KPOINTS,KP_ENRG,A1o,A2o,nKPOINTS,t)

!-------------------------------------------------------------!
! SUBROUTINE to sort the LIST OF KPOINTS based on neighboring !  
! category , GAMMA | FIRST NEIGHBORS | SECOND NEIGHBORS ... . !
!-------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter     :: dm = 2

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)      , intent(in)    :: nKPOINTS

  double precision, intent(in)    :: A1o
  double precision, intent(in)    :: A2o
  double precision, intent(in)    :: t

  double precision, intent(inout) :: KPOINTS(dm,nKPOINTS)
  double precision, intent(out)   :: KP_ENRG(nKPOINTS)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                     :: i, j, kpoint

 !-------------------------------------------------!
 ! SORTING CRITERIA: MODULE OF K POINT |K|         !
 !-------------------------------------------------!

  KP_ENRG = 0.d0

  do kpoint = 1,nKPOINTS,1
     KP_ENRG(kpoint) = -2.d0 * t *  ( dcos(KPOINTS(1,kpoint)*A1o) +  dcos(KPOINTS(2,kpoint)*A2o) )
  enddo

 !-------------------------------------------------!
 ! SUBROUTINE: BUBBLE SORTING                      !
 !-------------------------------------------------!

  do i = 1,nKPOINTS,1
     do j = 2,nKPOINTS-(i-1),1
        if ( KP_ENRG(j-1) .gt. KP_ENRG(j) ) then
           call swap( KP_ENRG(j-1)  , KP_ENRG(j)   )
           call swap( KPOINTS(1,j-1), KPOINTS(1,j) )
           call swap( KPOINTS(2,j-1), KPOINTS(2,j) )
        endif
     enddo
  enddo

  return

END SUBROUTINE sort_KPOINTS

!----------------------------------------------------------------------------------------------------------

SUBROUTINE transport_BZ( kBZ,kLIM1,kLIM2,recVEC )

!-------------------------------------------------------------!
! SUBROUTINE to TRANSPORT BACK the K-POINT kBZ back to the    !  
! COMPUTATIONAL BOX defined by:                               !
!                              kLIM1=[minK1,maxK1]            !
!                              kLIM2=[minK2,maxK2]            !
! The TRANSPORT is provided by the proper G VECTOR defined    !
! as:                                                         !
!     G = n1*b1 + n2*b2                                       !
! where the coefficients are n1,n2 = {-1,0,1} according to    !
! the K-POINT REPLICA of kNEW                                 !
!                                                             !
! - if kNEW(1) is out of the COMP. BOX. -> n1 = sgn{kNEW(1)}  !
! - if kNEW(1) is out of the COMP. BOX. -> n2 = sgn{kNEW(2)}  !
!                                                             !
! Then K-POINT kNEW is transported back into the BZ by the    !
! traslation:                                                 !
!             kBZ = kNEW - G                                  !
!-------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! PARAMETERS                                      !
 !-------------------------------------------------!

  integer(8)      , parameter     :: dm = 2

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  double precision, intent(in)    :: recVEC(dm,dm)
  double precision, intent(in)    :: kLIM1(2)
  double precision, intent(in)    :: kLIM2(2)

  double precision, intent(inout) :: kBZ(dm)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                      :: n1
  integer(8)                      :: n2

  double precision, allocatable   :: kNEW(:)
  double precision, allocatable   :: G(:)

  logical                         :: outCOMPBOX1
  logical                         :: outCOMPBOX2

 !-------------------------------------------------!
 ! SUBROUTINE: PARALLEL TRANSPORT kBZ = kNEW - G   !
 !-------------------------------------------------!

  allocate( kNEW(dm) )
  allocate(    G(dm) )

  kNEW(:) = kBZ(:)

  call check_kCOMPBOX( outCOMPBOX1,kNEW(1),kLIM1(:) )
  call check_kCOMPBOX( outCOMPBOX2,kNEW(2),kLIM2(:) )

  if     ( (outCOMPBOX1 .eqv. .FALSE.).and.(outCOMPBOX2 .eqv. .FALSE.) ) then
         n1 = 0.d0
         n2 = 0.d0
  elseif ( (outCOMPBOX1 .eqv. .FALSE.).and.(outCOMPBOX2 .eqv. .TRUE. ) ) then
         n1 = 0.d0
         n2 = sign(1.d0,kNEW(2))
  elseif ( (outCOMPBOX1 .eqv. .TRUE. ).and.(outCOMPBOX2 .eqv. .FALSE.) ) then
         n1 = sign(1.d0,kNEW(1))
         n2 = 0.d0
  elseif ( (outCOMPBOX1 .eqv. .TRUE. ).and.(outCOMPBOX2 .eqv. .TRUE. ) ) then
         n1 = sign(1.d0,kNEW(1))
         n2 = sign(1.d0,kNEW(2))
  endif

  G(:)   = dble(n1) * recVEC(:,1) + dble(n2) * recVEC(:,2)
  kBZ(:) = kNEW(:) - G(:)

  return

END SUBROUTINE transport_BZ

!---------------------------------------------------------------------------------------------------------

SUBROUTINE check_kCOMPBOX( outCOMPBOX,kNEW,kLIM )

!--------------------------------------------------------------!
! SUBROUTINE to check if the component kNEW of a given K-POINT !
! lies inside of the COMPUTATIONAL BOX of limits:              !
!                [kLIM(1),kLIM(2)]                             !
! Therefore:                                                   !
!           outCOMPBOX = .FALSE. it is NOT out of the CP if    !
!                        kLIM(1) =< kNEW =< kLIM(2)            !
!           outCOMPBOX = .TRUE.  it IS     out of the CP if    !
!                        kNEW < kLIM(1) or kLIM(2) < kNEW      !
!--------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  double precision, intent(in)    :: kNEW
  double precision, intent(in)    :: kLIM(2)

  logical         , intent(out)   :: outCOMPBOX

 !-------------------------------------------------!
 ! SUBROUTINE: CHECK                               !
 !-------------------------------------------------!

  if ( (kNEW .ge. kLIM(1)).and.(kNEW .le. kLIM(2)) ) then
     outCOMPBOX = .FALSE.
  else
     outCOMPBOX = .TRUE.
  endif

  return

END SUBROUTINE check_kCOMPBOX
