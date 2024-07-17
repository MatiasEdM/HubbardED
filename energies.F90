DOUBLE PRECISION FUNCTION int_energy (CI,dOCC,nSITES,nDET,U)

!--------------------------------------------------------------!
! FUNCTION that calculates the total INTERACTION energy of the !
! many-body wave function by summing the contribution of each  !
! determinant with dubly occupied sites.                       !
!--------------------------------------------------------------!

 implicit none

 integer(8)      , intent(in) :: nSITES
 integer(8)      , intent(in) :: nDET
 integer(8)      , intent(in) :: dOCC(nSITES,nDET)

 double precision, intent(in) :: CI(nDET)

 double precision, intent(in) :: U

 integer(8)                   :: I

  int_energy = 0.d0

  do I = 1,nDET,1
     int_energy = int_energy + CI(I) * sum(dOCC(:,I))
  enddo

  int_energy = U * int_energy

  return

END FUNCTION int_energy

!----------------------------------------------------------------------------------------------------------

DOUBLE PRECISION FUNCTION kin_energy(Etot,int_energy)

!--------------------------------------------------!
! FUNCTION to calculate the KINETIC energy as the  !
! difference between the TOTAL and the INTERACTION !
! energies: T = E - U                              !
!--------------------------------------------------!

 implicit none

 double precision, intent(in) :: Etot, int_energy
 
 kin_energy = Etot - int_energy

 return

END FUNCTION kin_energy 
