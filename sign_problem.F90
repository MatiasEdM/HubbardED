SUBROUTINE estimate_sign_problem(normDC,H,CI,nDET)

!--------------------------------------------------------------!
! SUBROUTINE to estimate the magnitude of the sign problem     !
! by computing the VEC containing the differences:             !
!  DC(I) = | sum_{J}[ C(J)H(J,I) ] | - sum_{J}[ |C(J)H(J,I)| ] !
! for the ground-state wave function.                          !
!--------------------------------------------------------------!

  implicit none

 !--------------------------------------------------!
 ! INOUT VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)      , intent(in)  :: nDET

  double precision, intent(in)  :: H(nDET,nDET)
  double precision, intent(in)  :: CI(nDET)
  double precision, intent(out) :: normDC

 !--------------------------------------------------!
 ! LOCAL VARIABLES                                  !
 !--------------------------------------------------!

  integer(8)                    :: I, J

  double precision              :: C, Cp
  double precision              :: L2_norm

  double precision, allocatable :: DC(:)

 !--------------------------------------------------!
 ! MEMORY ALLOCATION                                !
 !--------------------------------------------------!

  allocate( DC(nDET) )

 !--------------------------------------------------!
 ! DC VECTOR AND ||DC||2 NORM CALCULATION           !
 !--------------------------------------------------!

  DC(:)  = 0.d0
 
  do I = 1,nDET,1
     C  = 0.d0
     Cp = 0.d0
     do J = 1,nDET,1
        if ( J .eq. I ) then
           cycle
        endif
        C  = C  +      CI(J) * H(J,I)
        Cp = Cp + abs( CI(J) * H(J,I) )
     enddo
     DC(I) = Cp - abs(C)
  enddo

  normDC   = L2_norm( DC(:),nDET )

 !--------------------------------------------------!
 ! MEMORY DEALLOCATION                              !
 !--------------------------------------------------!

  deallocate( DC )

 return

END SUBROUTINE estimate_sign_problem
