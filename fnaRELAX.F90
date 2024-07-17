SUBROUTINE generate_Hpq( Heff,Vsf,SGN,H,CI,nDET,nPDET,Vsf_MAX  )

!------------------------------------------------------------------------!
! SUBROUTINE to generate the EFFECTIVE (FNA) HAMILTONIAN Heff according  !
! to the P/Q PARTITION and FNA RELAXATION:                               !
!            Heff = PxHxP + PxHxQ + QxHxP + QxH[fn]xQ                    !
!                 = H[PP] + H[PQ] + H[QP] + H[QQ][fn]                    !
! where P and Q are PROJECTORS that project on to the P and Q SUBSPACES  !
!                            -                                           !
! - P-SUBSPACE  : of DIM. nPDET and composed by the first nPDET DETs. of !
!                 higher amplitude in the real-space representation      !
! - Q-SPACE     : the complement of P.                                   !
!                                               -                        !
! The HAM. MAT. elements <R|H[fn]|R'> in the Q-SPACE are calculated as:  !
!                             -                                          !
! - DIAGONAL    :                                                        !
!                <R|H[fn]|R>  = <R|H|R> + <R|Vsf|R>                      !
! - OFF-DIAGONAL:                                                        !
!                              -                                         !
!                             |  <R|H|R'> if <R|SGN|R'> < 0              !
!               <R|H[fn]|R'> =|                                          !
!                             |     0     if <R|SGN|R'> > 0 [OR A NODE]  !
!                              -                                         !
!------------------------------------------------------------------------!

  implicit none

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)      , intent(in)  :: nDET
  integer(8)      , intent(in)  :: nPDET

  double precision, intent(in)  :: H(nDET,nDET)
  double precision, intent(in)  :: CI(nDET)
  double precision, intent(in)  :: SGN(nDET,nDET)
  double precision, intent(in)  :: Vsf_MAX

  double precision, intent(out) :: Heff(nDET,nDET)
  double precision, intent(out) :: Vsf(nDET)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                    :: I, J
  
  double precision              :: SF_PQ

  character(1)    , allocatable :: PQ(:)

 !-------------------------------------------------!
 ! SUBROUTINE: P/Q PARTITION                       !
 !-------------------------------------------------!

 allocate( PQ(nDET) )

 call get_PQ( PQ(:),CI(:),nDET,nPDET )

 !-------------------------------------------------!
 ! SUBROUTINE: EFFECTIVE (FNA) HAMILTONIAN MAT.    !
 !-------------------------------------------------!

  Heff(:,:) = 0.d0
  Vsf(:)    = 0.d0

  do I = 1,nDET,1
     if ( PQ(I) .eq. 'P' ) then
        Heff(I,I) = H(I,I)
        do J = 1,nDET,1
           if ( J .eq. I ) cycle
           Heff(I,J) = H(I,J)
        enddo
     elseif ( PQ(I) .eq. 'Q' ) then
        Vsf(I)    = SF_PQ( I,SGN(:,:),PQ(:),H(:,:),nDET,Vsf_MAX )
        Heff(I,I) = H(I,I) + Vsf(I)
        do J = 1,nDET,1
           if ( J .eq. I ) cycle
           if ( PQ(J) .eq. 'P' ) then
              Heff(I,J)   = H(I,J)
           elseif ( PQ(J) .eq. 'Q' ) then
             if ( SGN(I,J) .lt. 0 ) then
                Heff(I,J) = H(I,J)
             endif
           endif
        enddo
     endif
  enddo

  deallocate( PQ )

  return

END SUBROUTINE generate_Hpq

!----------------------------------------------------------------------------------------------------------

SUBROUTINE get_PQ(PQ,CI,nDET,nPDET)

!------------------------------------------------------------------------!
! SUBROUTINE to generate the P/Q PARTITION of the HILBERT SPACE as       !
!                            -                                           !
! - P-SUBSPACE  : of DIM. nPDET and composed by the first nPDET DETs. of !
!                 higher amplitude in the real-space representation      !
! - Q-SPACE     : the complement of P.                                   !
!------------------------------------------------------------------------!

  implicit none

 !---------------------------------------------------!
 ! PARAMETERS                                        !
 !---------------------------------------------------!

  double precision, parameter   :: eps = 1.E-15

 !-------------------------------------------------!
 ! INPUT VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)      , intent(in)  :: nDET
  integer(8)      , intent(in)  :: nPDET

  double precision, intent(in)  :: CI(nDET)

  character(1)    , intent(out) :: PQ(nDET)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                    :: I
  integer(8)                    :: D

  double precision, allocatable :: CIaux(:)

 !-------------------------------------------------!
 ! SUBROUTINE:                      !
 !-------------------------------------------------!

  allocate( CIaux(nDET) )

  CIaux(:) = CI(:)

  do D = 1,nDET,1
     if ( CI(D) .eq. 0.d0 ) then
        CIaux(D) = eps
     endif
  enddo

  PQ(:)    = 'Q'

  do D = 1,nPDET,1
     I         = maxloc(abs(CIaux(:)),nDET)
     PQ(I)     = 'P'
     CIaux(I)  = 0.d0
  enddo
  
  deallocate( CIaux )

  return

END SUBROUTINE get_PQ

!----------------------------------------------------------------------------------------------------------

DOUBLE PRECISION FUNCTION SF_PQ( I,SGN,PQ,H,nDET,Vsf_MAX )

!--------------------------------------------------------------------------------------------!
! FUNCTION that calculates the SIGN FLIP POTENTIAL Vsf(R) contribution to the CONFIG. (DET.) !
! I (|R>) as the diagonal matrix element:                                                    !
!                                        <R|Vsf|R> = SUM_{R'} <R|H|R'> <R'|HF>/<R|HF>        !
!                                                                                            !      
! where the sum runs over all CONFIG. (DET.) for which the transition |R'> -> |R> is SIGN-   !
! VIOLATING, i.e. <R'|SGN|R> > 0, that belong to the Q-SPACE.                                !
! The projection of the TRIAL WFC (HF) corresponds to,                                       !
!                                                         <R|TWFC> = det[M]                  !
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

  character(1)    , intent(in)  :: PQ(nDET)

 !-------------------------------------------------!
 ! LOCAL VARIABLES                                 !
 !-------------------------------------------------!

  integer(8)                    :: J

 !-------------------------------------------------!
 ! SUBROUTINE: CALCULATE THE SIGN-FLIP POTENTIAL   !
 !-------------------------------------------------!

  SF_PQ = 0.d0
  do J = 1,nDET,1
     if ( PQ(J) .eq. 'Q' ) then
        if ( SGN(I,J) .eq. Vsf_MAX ) then
           SF_PQ = Vsf_MAX
           exit
        elseif ( SGN(I,J) .gt. 0.d0 ) then
           SF_PQ = SF_PQ + SGN(I,J)
        endif
     endif
  enddo

  return

END FUNCTION SF_PQ
